#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#include "vector_filter.h"


#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* Theory on probability and scoring *ungapped* alignment
 *
 * s'(a,b) = log[P(b|a)/P(b)] = log[4P(b|a)], assuming uniform base distribution
 * s'(a,a) = log(4), s'(a,b) = log(4e/3), where e is the error rate
 *
 * Scale s'(a,b) to s(a,a) s.t. s(a,a)=x. Then s(a,b) = x*s'(a,b)/log(4), or conversely: s'(a,b)=s(a,b)*log(4)/x
 *
 * If the matching score is x and mismatch penalty is -y, we can compute error rate e:
 *   e = .75 * exp[-log(4) * y/x]
 *
 * log P(seq) = \sum_i log P(b_i|a_i) = \sum_i {s'(a,b) - log(4)}
 *   = \sum_i { s(a,b)*log(4)/x - log(4) } = log(4) * (S/x - l)
 *
 * where S=\sum_i s(a,b) is the alignment score. Converting to the phred scale:
 *   Q(seq) = -10/log(10) * log P(seq) = 10*log(4)/log(10) * (l - S/x) = 6.02 * (l - S/x)
 *
 *
 * Gap open (zero gap): q' = log[P(gap-open)], r' = log[P(gap-ext)] (see Durbin et al. (1998) Section 4.1)
 * Then q = x*log[P(gap-open)]/log(4), r = x*log[P(gap-ext)]/log(4)
 *
 * When there are gaps, l should be the length of alignment matches (i.e. the M operator in CIGAR)
 */

static const bntseq_t *global_bns = 0; // for debugging only

mem_opt_t *mem_opt_init() {
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1;
	o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	//o->zdrop = 0;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;
	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1 << 30;
	o->mapQ_coef_len = 50;
	o->mapQ_coef_fac = log(o->mapQ_coef_len);
	o->seed_type = 1;
	o->seed_intv = o->min_seed_len;
	o->dp_type = 0;
	o->opt_ext = 0;
	o->re_seed = 0;
	o->use_avx2 = 0;
	o->read_len = 0;
	o->shd_filter = 0;
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/***************************
 * Collection SA invervals *
 ***************************/

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

	typedef struct {
		bwtintv_v mem, mem1, *tmpv[2];
	} smem_aux_t;

static smem_aux_t *smem_aux_init() {
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a) {
	free(a->tmpv[0]->a);
	free(a->tmpv[0]);
	free(a->tmpv[1]->a);
	free(a->tmpv[1]);
	free(a->mem.a);
	free(a->mem1.a);
	free(a);
}

typedef struct {
	bwtint_t w;
	int bid;
} bwt_mem_width_t;

int bwt_mem_cal_width(const bwt_t *bwt, int len, const uint8_t *str, bwt_mem_width_t *width) {
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0;
	l = bwt->seq_len;
	for (i = 0; i < len; ++i) {
		uint8_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = bwt->seq_len;
			++bid;
			i--;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	//width[len].w = 0;
	//width[len].bid = ++bid;
	return bid;
}

typedef struct {
	bwtint_t w;
	int bid;
} bwt_width_t;

static void mem_collect_intv(mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a) {
	int i, k, x = 0, old_n, max, max_i;
	int start_width = 1;
	int split_len = (int) (opt->min_seed_len * opt->split_factor + .499);
	int min_miss_match = 0;
	bwt_mem_width_t *w;
	/*if (opt->seed_type == 4) {
	  w = calloc(len, sizeof(bwt_mem_width_t));
	  min_miss_match = bwt_mem_cal_width(bwt, len, seq, w);
	  }*/
	//int bowtie_seed_len = 0;
	//if (opt->seed_type == 1) bowtie_seed_len = 0;
	//else if (opt->seed_type == 2) bowtie_seed_len = opt->bowtie_seed_len;
	//else if (opt->seed_type == 3) bowtie_seed_len = 0;
	//else bowtie_seed_len = 0;
	a->mem.n = 0;
	// first pass: find all SMEMs
	//extern __thread unsigned long int bwt_extend_call_1_local, bwt_extend_call_2_local, bwt_extend_call_3_local, bwt_extend_call_4_local, bwt_extend_call_5_local, bwt_extend_call_6_local;
	//bwt_extend_call_1_local = 0, bwt_extend_call_2_local = 0, bwt_extend_call_3_local = 0, bwt_extend_call_4_local = 0, bwt_extend_call_5_local = 0, bwt_extend_call_6_local = 0;
	while (x < len) {
		if (seq[x] < 4) {
			if (opt->seed_type == 1) {
				x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			} else if (opt->seed_type == 2) {
				if (x + opt->min_seed_len <= len) {
					bwt_bowtie_seed(bwt, len, seq, x, start_width, 0, &a->mem1, opt->min_seed_len);
					x = x + opt->seed_intv;
				} else
					break;
			} else if (opt->seed_type == 3) {
				x = bwt_fwd_mem(bwt, len, seq, x, start_width, 0, &a->mem1);
			} else if (opt->seed_type == 4) {
				if (x + opt->min_seed_len <= len) {
					bwt_bowtie_seed_inexact(bwt, len, seq, x, start_width, 0, &a->mem1, 1, opt->min_seed_len);
					x = x + opt->seed_intv;
				} else
					break;
			} else
				x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				int slen = (uint32_t) p->info - (p->info >> 32); // seed length
				if (slen /*- (opt->seed_type == 4 ? (p->n_miss_match * opt->b) : 0))*/
					>= opt->min_seed_len)
						kv_push(bwtintv_t, a->mem, *p);
			}
		} else
			++x;
	}
	if (opt->seed_type == 1 && opt->re_seed == 1) {
		old_n = a->mem.n;
		for (k = 0; k < old_n; ++k) {
			bwtintv_t *p = &a->mem.a[k];
			int start = p->info >> 32, end = (int32_t) p->info;
			if (end - start < split_len || p->x[2] > opt->split_width)
				continue;
			bwt_smem1(bwt, len, seq, (start + end) >> 1, p->x[2] + 1, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i)
				if ((uint32_t) a->mem1.a[i].info - (a->mem1.a[i].info >> 32) >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
		}
		// third pass: LAST-like
		if (opt->max_mem_intv > 0) {
			x = 0;
			while (x < len) {
				if (seq[x] < 4) {
					if (1) {
						bwtintv_t m;
						x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
						if (m.x[2] > 0)
							kv_push(bwtintv_t, a->mem, m);
					} else { // for now, we never come to this block which is slower
						x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
						for (i = 0; i < a->mem1.n; ++i)
							kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
					}
				} else
					++x;
			}
		}
	}
	/*for (i = 0; i < a->mem1.n; ++i)
	  if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
	  kv_push(bwtintv_t, a->mem, a->mem1.a[i]);*/
	// second pass: find MEMs inside a long SMEM
	/*old_n = a->mem.n;
	  for (k = 0; k < old_n; ++k) {
	  bwtintv_t *p = &a->mem.a[k];
	  int start = p->info>>32, end = (int32_t)p->info;
	  if (end - start < split_len || p->x[2] > opt->split_width) continue;
	  bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
	  for (i = 0; i < a->mem1.n; ++i)
	  if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
	  kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
	  }
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
	x = 0;
	while (x < len) {
	if (seq[x] < 4) {
	if (1) {
	bwtintv_t m;
	x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
	if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
	} else { // for now, we never come to this block which is slower
	x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
	for (i = 0; i < a->mem1.n; ++i)
	kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
	}
	} else ++x;
	}
	}*/
	//sort
	ks_introsort(mem_intv, a->mem.n, a->mem.a);
}

/*static void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
  {
  int i, k, x = 0, old_n;
  int start_width = 1;
  int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
  a->mem.n = 0;
// first pass: find all SMEMs
while (x < len) {
if (seq[x] < 4) {
x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
for (i = 0; i < a->mem1.n; ++i) {
bwtintv_t *p = &a->mem1.a[i];
int slen = (uint32_t)p->info - (p->info>>32); // seed length
if (slen >= opt->min_seed_len)
kv_push(bwtintv_t, a->mem, *p);
}
} else ++x;
}
// second pass: find MEMs inside a long SMEM
old_n = a->mem.n;
for (k = 0; k < old_n; ++k) {
bwtintv_t *p = &a->mem.a[k];
int start = p->info>>32, end = (int32_t)p->info;
if (end - start < split_len || p->x[2] > opt->split_width) continue;
bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
for (i = 0; i < a->mem1.n; ++i)
if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
}
// third pass: LAST-like
if (opt->max_mem_intv > 0) {
x = 0;
while (x < len) {
if (seq[x] < 4) {
if (1) {
bwtintv_t m;
x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
} else { // for now, we never come to this block which is slower
x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
for (i = 0; i < a->mem1.n; ++i)
kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
}
} else ++x;
}
}
// sort
ks_introsort(mem_intv, a->mem.n, a->mem.a);
}*/

/************
 * Chaining *
 ************/

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int score;
} mem_seed_t; // unaligned memory

typedef struct {
	int n, m, first, rid;
	uint32_t w :29, kept :2, is_alt :1;
	float frac_rep;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain_t;

typedef struct {
	size_t n, m;
	mem_chain_t *a;
} mem_chain_v;

#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

	// return 1 if the seed is merged into the chain
	static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p, int seed_rid) {
		int64_t qend, rend, x, y;
		const mem_seed_t *last = &c->seeds[c->n - 1];
		qend = last->qbeg + last->len;
		rend = last->rbeg + last->len;
		if (seed_rid != c->rid)
			return 0; // different chr; request a new chain
		if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
			return 1; // contained seed; do nothing
		if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac)
			return 0; // don't chain if on different strand
		x = p->qbeg - last->qbeg; // always non-negtive
		y = p->rbeg - last->rbeg;
		if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
			if (c->n == c->m) {
				c->m <<= 1;
				c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
			}
			c->seeds[c->n++] = *p;
			return 1;
		}
		return 0; // request to add a new chain
	}

int mem_chain_weight(const mem_chain_t *c) {
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end)
			w += s->len;
		else if (s->qbeg + s->len > end)
			w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len ? end : s->qbeg + s->len;
	}
	tmp = w;
	w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end)
			w += s->len;
		else if (s->rbeg + s->len > end)
			w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len ? end : s->rbeg + s->len;
	}
	w = w < tmp ? w : tmp;
	return w < 1 << 30 ? w : (1 << 30) - 1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn) {
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		err_printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev)
				pos -= p->seeds[j].len - 1;
			err_printf("\t%d;%d;%d,%ld(%s:%c%ld)", p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg, (long) p->seeds[j].rbeg, bns->anns[p->rid].name,
					"+-"[is_rev], (long) (pos - bns->anns[p->rid].offset) + 1);
		}
		err_putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf) {
	int i, b, e, l_rep;
	int64_t l_pac = bns->l_pac;
	mem_chain_v chain;
	kbtree_t(chn) *tree;
	smem_aux_t *aux;

	kv_init(chain);
	if (len < opt->min_seed_len)
		return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);

	aux = buf ? (smem_aux_t*) buf : smem_aux_init();
	mem_collect_intv(opt, bwt, len, seq, aux);
	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) { // compute frac_rep
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info >> 32), se = (uint32_t) p->info;
		if (p->x[2] <= opt->max_occ)
			continue;
		if (sb > e)
			l_rep += e - b, b = sb, e = se;
		else
			e = e > se ? e : se;
	}
	l_rep += e - b;
	for (i = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int step, count, slen = (uint32_t) p->info - (p->info >> 32); // seed length
		int64_t k;
		// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
		step = p->x[2] > opt->max_occ ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			mem_chain_t tmp, *lower, *upper;
			mem_seed_t s;
			int rid, to_add = 0;
			s.rbeg = tmp.pos = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = p->info >> 32;
			s.len = slen;
			s.score = slen - ((opt->seed_type == 4) ? (p->n_miss_match * opt->b) : 0);
			rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
			if (rid < 0)
				continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
			if (kb_size(tree)) {
				kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
				if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid))
					to_add = 1;
			} else
				to_add = 1;
			if (to_add) { // add the seed as a new chain
				tmp.n = 1;
				tmp.m = 4;
				tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
				tmp.seeds[0] = s;
				tmp.rid = rid;
				tmp.is_alt = !!bns->anns[rid].is_alt;
				kb_putp(chn, tree, &tmp);
			}
		}
	}
	if (buf == 0)
		smem_aux_destroy(aux);

	kv_resize(mem_chain_t, chain, kb_size(tree));

#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
#undef traverse_func

	for (i = 0; i < chain.n; ++i)
		chain.a[i].frac_rep = (float) l_rep / len;
	if (bwa_verbose >= 4)
		printf("* fraction of repetitive seeds: %.3f\n", (float) l_rep / len);

	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)

	int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a) {
		int i, k;
		kvec_t(int)
			chains = { 0, 0, 0 }; // this keeps int indices of the non-overlapping chains
		if (n_chn == 0)
			return 0; // no need to filter
		// compute the weight of each chain and drop chains with small weight
		for (i = k = 0; i < n_chn; ++i) {
			mem_chain_t *c = &a[i];
			c->first = -1;
			c->kept = 0;
			c->w = mem_chain_weight(c);
			if (c->w < opt->min_chain_weight)
				free(c->seeds);
			else
				a[k++] = *c;
		}
		n_chn = k;
		ks_introsort(mem_flt, n_chn, a);
		// pairwise chain comparisons
		a[0].kept = 3;
		kv_push(int, chains, 0);
		for (i = 1; i < n_chn; ++i) {
			int large_ovlp = 0;
			for (k = 0; k < chains.n; ++k) {
				int j = chains.a[k];
				int b_max =
					chn_beg(a[j]) > chn_beg(a[i]) ? chn_beg(a[j]) : chn_beg(a[i]);
				int e_min =
					chn_end(a[j]) < chn_end(a[i]) ? chn_end(a[j]) : chn_end(a[i]);
				if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
					int li = chn_end(a[i]) - chn_beg(a[i]);
					int lj = chn_end(a[j]) - chn_beg(a[j]);
					int min_l = li < lj ? li : lj;
					if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
						large_ovlp = 1;
						if (a[j].first < 0)
							a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
						if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len << 1)
							break;
					}
				}
			}
			if (k == chains.n) {
				kv_push(int, chains, i);
				a[i].kept = large_ovlp ? 2 : 3;
			}
		}
		for (i = 0; i < chains.n; ++i) {
			mem_chain_t *c = &a[chains.a[i]];
			if (c->first >= 0)
				a[c->first].kept = 1;
		}
		free(chains.a);
		for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
			if (a[i].kept == 0 || a[i].kept == 3)
				continue;
			if (++k >= opt->max_chain_extend)
				break;
		}
		for (; i < n_chn; ++i)
			if (a[i].kept < 3)
				a[i].kept = 0;
		for (i = k = 0; i < n_chn; ++i) { // free discarded chains
			mem_chain_t *c = &a[i];
			if (c->kept == 0)
				free(c->seeds);
			else
				a[k++] = a[i];
		}
		return k;
	}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

	int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b,
			int *_w) {
		int w, score, q_s, r_s;
		double r;
		if (bns == 0 || pac == 0 || query == 0)
			return 0;
		assert(a->rid == b->rid && a->rb <= b->rb);
		if (a->rb < bns->l_pac && b->rb >= bns->l_pac)
			return 0; // on different strands
		if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re)
			return 0; // not colinear
		w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
		w = w > 0 ? w : -w; // l = abs(l)
		r = (double) (a->re - b->rb) / (b->re - a->rb) - (double) (a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
		r = r > 0. ? r : -r; // r = fabs(r)
		if (bwa_verbose >= 4)
			printf("* potential hit merge between [%d,%d)<=>[%ld,%ld) and [%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n", a->qb, a->qe, (long) a->rb,
					(long) a->re, b->qb, b->qe, (long) b->rb, (long) b->re, bns->anns[a->rid].name, w, r);
		if (a->re < b->rb || a->qe < b->qb) { // no overlap on query or on ref
			if (w > opt->w << 1 || r >= PATCH_MAX_R_BW)
				return 0; // the bandwidth or the relative bandwidth is too large
		} else if (w > opt->w << 2 || r >= PATCH_MAX_R_BW * 2)
			return 0; // more permissive if overlapping on both ref and query
		// global alignment
		w += a->w + b->w;
		w = w < opt->w << 2 ? w : opt->w << 2;
		if (bwa_verbose >= 4)
			printf("* test potential hit merge with global alignment; w=%d\n", w);
		bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re, &score, 0,
				0);
		q_s = (int) ((double) (b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) * (b->score + a->score) + .499); // predicted score from query
		r_s = (int) ((double) (b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) * (b->score + a->score) + .499); // predicted score from ref
		if (bwa_verbose >= 4)
			printf("* score=%d;(%d,%d)\n", score, q_s, r_s);
		if ((double) score / (q_s > r_s ? q_s : r_s) < PATCH_MIN_SC_RATIO)
			return 0;
		*_w = w;
		return score;
	}

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a) {
	int m, i, j;
	if (n <= 1)
		return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	for (i = 0; i < n; ++i)
		a[i].n_comp = 1;
	for (i = 1; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i - 1].rid || p->rb >= a[i - 1].re + opt->max_chain_gap)
			continue; // then no need to go into the loop below
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t pr, pq, mr, mq;
			int score, w;
			if (q->qe == q->qb)
				continue; // a[j] has been excluded
			pr = q->re - p->rb; // overlap length on the reference
			pq = q->qb < p->qb ? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb ? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb ? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (pr > opt->mask_level_redun * mr && pq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score) {
					p->qe = p->qb;
					break;
				} else
					q->qe = q->qb;
			} else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov ? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub ? p->sub : q->sub;
				p->csub = p->csub > q->csub ? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i)
				a[m++] = a[i];
			else
				++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i - 1].score && a[i].rb == a[i - 1].rb && a[i].qb == a[i - 1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i)
				a[m++] = a[i];
			else
				++m;
		}
	return m;
}

typedef kvec_t(int)
	int_v;

	static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z) { // similar to the loop in mem_chain_flt()
		int i, k, tmp;
		tmp = opt->a + opt->b;
		tmp = opt->o_del + opt->e_del > tmp ? opt->o_del + opt->e_del : tmp;
		tmp = opt->o_ins + opt->e_ins > tmp ? opt->o_ins + opt->e_ins : tmp;
		z->n = 0;
		kv_push(int, *z, 0);
		for (i = 1; i < n; ++i) {
			for (k = 0; k < z->n; ++k) {
				int j = z->a[k];
				int b_max = a[j].qb > a[i].qb ? a[j].qb : a[i].qb;
				int e_min = a[j].qe < a[i].qe ? a[j].qe : a[i].qe;
				if (e_min > b_max) { // have overlap
					int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb ? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
					if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
						if (a[j].sub == 0)
							a[j].sub = a[i].score;
						if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
							++a[j].sub_n;
						break;
					}
				}
			}
			if (k == z->n)
				kv_push(int, *z, i);
			else
				a[i].secondary = z->a[k];
		}
	}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id) {
	int i, n_pri;
	int_v z = { 0, 0, 0 };
	if (n == 0)
		return 0;
	for (i = n_pri = 0; i < n; ++i) {
		a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id + i);
		if (!a[i].is_alt)
			++n_pri;
	}
	ks_introsort(mem_ars_hash, n, a);
	mem_mark_primary_se_core(opt, n, a, &z);
	for (i = 0; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		p->secondary_all = i; // keep the rank in the first round
		if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt) {
			fprintf(stderr, "I am here 1");
			p->alt_sc = a[p->secondary].score;
		}
	}
	if (n_pri >= 0 && n_pri < n) {
		fprintf(stderr, "I am here 2");
		kv_resize(int, z, n);
		if (n_pri > 0)
			ks_introsort(mem_ars_hash2, n, a);
		for (i = 0; i < n; ++i)
			z.a[a[i].secondary_all] = i;
		for (i = 0; i < n; ++i) {
			if (a[i].secondary >= 0) {
				a[i].secondary_all = z.a[a[i].secondary];
				if (a[i].is_alt)
					a[i].secondary = INT_MAX;
			} else
				a[i].secondary_all = -1;
		}
		if (n_pri > 0) { // mark primary for hits to the primary assembly only
			for (i = 0; i < n_pri; ++i)
				a[i].sub = 0, a[i].secondary = -1;
			mem_mark_primary_se_core(opt, n_pri, a, &z);
		}
	} else {
		for (i = 0; i < n; ++i)
			a[i].secondary_all = a[i].secondary;
	}
	free(z.a);
	return n_pri;
}

/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_seed_t *s) {
	int qb, qe, rid;
	int64_t rb, re, mid, l_pac = bns->l_pac;
	uint8_t *rseq = 0;
	kswr_t x;

	if (s->len >= MEM_SHORT_LEN)
		return -1; // the seed is longer than the max-extend; no need to do SW
	qb = s->qbeg, qe = s->qbeg + s->len;
	rb = s->rbeg, re = s->rbeg + s->len;
	mid = (rb + re) >> 1;
	qb -= MEM_SHORT_EXT;
	qb = qb > 0 ? qb : 0;
	qe += MEM_SHORT_EXT;
	qe = qe < l_query ? qe : l_query;
	rb -= MEM_SHORT_EXT;
	rb = rb > 0 ? rb : 0;
	re += MEM_SHORT_EXT;
	re = re < l_pac << 1 ? re : l_pac << 1;
	if (rb < l_pac && l_pac < re) {
		if (mid < l_pac)
			re = l_pac;
		else
			rb = l_pac;
	}
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN)
		return -1; // the seed seems good enough; no need to do SW

	rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	x = ksw_align2(qe - qb, (uint8_t*) query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0,
			opt->use_avx2);
	free(rseq);
	return x.score;
}

void mem_shd_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a) 
{
	int i, j, k;
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			mem_seed_t *s = &c->seeds[j];
			int qb, qe, is_right_pass = 1, is_left_pass = 1, right_len = 0, left_len = 0;
			int64_t rb, re, l_pac = bns->l_pac;
			qb = s->qbeg, qe = s->qbeg + s->len;
			rb = s->rbeg, re = s->rbeg + s->len;
			uint8_t seed2bit[s->len];
			memcpy(seed2bit, query + qb, s->len);
			int l;
			char seed[s->len+1];
			for (l=0; l < s->len; l++){
				seed[l] = "ACGTN"[(int) seed2bit[l]];
			}
			seed[l] = '\0';

			if (qe < l_query) {
				uint8_t *rseq = 0;
				int rid;
				int qbeg = (qe - 5) > 0 ? (qe - 5) : 0;
				int qend = (qbeg + 128) <= l_query ? (qbeg + 128) : l_query;
				char read_display[320] __aligned__;
				{
					int q, f;
					for (q = qbeg, f = 0; f < 320; q++, f++) {
						if (q < qend)
							read_display[f] = "ACGTN"[(int) query[q]];
						else
							read_display[f] = '\0';
					}
				}
				char read_t[320] __aligned__;
				//strcpy(read_t, "TCCTCCAAGAAGATATGTAGTTGGTAAATAAACATAAGAAAAGATGCTCAAAACAATATGTTATTAGGGAACTTCAAATTAACATGATGAGATACCATTATACACCAATTAGAATGTCTAATATCTGA");
				int q, f;
				for (q = qbeg, f = 0; f < 320; q++, f++) {
					if (q < qend)
						read_t[f] = "ACGTN"[(int) query[q]];
					else
						read_t[f] = '\0';
				}
				int64_t rbeg = (re - 5) > 0 ? (re - 5) : 0;
				int64_t rend = (rbeg + 128) < (l_pac << 1) ? rbeg + 128 : (l_pac << 1);
				int64_t mid = (rbeg + rend) >> 1;
				if (rbeg < l_pac && l_pac < rend) {
					if (mid < l_pac)
						rend = l_pac;
					else
						rbeg = l_pac;
				}
				rseq = bns_fetch_seq(bns, pac, &rbeg, mid, &rend, &rid);
				char ref_display[320] __aligned__;
				{
					int r, f;
					for (r = 0, f = 0; f < 320; r++, f++) {
						if (r < (rend - rbeg)) {
							//fprintf(stderr, "%d, ", r);
							//fflush(stderr);
							ref_display[f] = "ACGTN"[(int) rseq[r]];
						} else
							ref_display[f] = '\0';
					}
				}
				char ref_t[320] __aligned__;
				//strcpy(ref_t, "TCCTCAAAGAAGATATGTAGTTGGTAAATAAACATAAGAAAAGATGCTCAAAACAATATGTTATTAGGGAACTTCAAATTAACATGATGAGATACCATTATACACCAATTAGAATGTCTAATATCTGA");
				int r;
				for (r = 0, f = 0; f < 320; r++, f++) {
					if (r < (rend - rbeg)) {
						//fprintf(stderr, "%d, ", r);
						//fflush(stderr);
						ref_t[f] = "ACGTN"[(int) rseq[r]];
					} else
						ref_t[f] = '\0';
				}
				//fprintf(stderr, "\n");
				right_len = qend - qbeg;
				//int max_err = (int) ceil(0.1 * (double) right_len);
				is_right_pass = bit_vec_filter_sse1(read_t, ref_t, (rend - rbeg), 7);
				/*if (!is_right_pass){
				  fprintf(stderr, "%d, %d, %s, %s, %s\n", max_err, rend - rbeg, seed, read_display, ref_display);
				//fprintf(stderr, "%d, %d, %s, %s, %s\n", qb, qe, seed, read_t, ref_t);
				//fprintf(stderr, "fail\n");
				}*/
				free(rseq);
			}
			if (qb > 0) {
				uint8_t *rseq = 0;
				int rid;
				int qend = (qb + 5) < l_query ? (qb + 5) : l_query;
				int qbeg = (qend - 128) > 0 ? (qend - 128) : 0;
				char read_display[320] __aligned__;
				{
					int q, f;
					for (q = qend - 1, f = 0; f < 320; q--, f++) {
						if (q >= qbeg)
							read_display[f] = "ACGTN"[(int) query[q]];
						else
							read_display[f] = '\0';
					}
				}
				char read_t[320] __aligned__;
				int q, f;
				for (q = qend - 1, f = 0; f < 320; q--, f++) {
					if (q >= qbeg)
						read_t[f] = "ACGTN"[(int) query[q]];
					else
						read_t[f] = '\0';
				}
				int64_t rend = (rb + 5) < (l_pac << 1) ? (rb + 5) : (l_pac << 1);
				int64_t rbeg = (rend - 128) > 0 ? rend - 128 : 0;
				int64_t mid = (rbeg + rend) >> 1;
				if (rbeg < l_pac && l_pac < rend) {
					if (mid < l_pac)
						rend = l_pac;
					else
						rbeg = l_pac;
				}
				rseq = bns_fetch_seq(bns, pac, &rbeg, mid, &rend, &rid);
				char ref_display[320] __aligned__;
				{
					int r, f;
					for (r = (rend - rbeg) - 1, f = 0; f < 320; r--, f++) {
						if (r >= 0)
							ref_display[f] = "ACGTN"[(int) rseq[r]];
						else
							ref_display[f] = '\0';
					}
				}
				char ref_t[320] __aligned__;
				int r;
				for (r = (rend - rbeg) - 1, f = 0; f < 320; r--, f++) {
					if (r >= 0)
						ref_t[f] = "ACGTN"[(int) rseq[r]];
					else
						ref_t[f] = '\0';
				}
				left_len = qend - qbeg;
				//int max_err = (int) ceil(0.1 * (double) left_len);
				is_left_pass = bit_vec_filter_sse1(read_t, ref_t, (rend - rbeg), 7);
				/*if (!is_left_pass) {
				  fprintf(stderr, "%d, %d, %s, %s, %s\n", max_err, rend - rbeg, seed, read_display, ref_display);
				//fprintf(stderr, "%d, %d, %s, %s, %s\n", qb, qe, seed, read_t, ref_t);
				//fprintf(stderr, "fail\n");
				}*/
				free(rseq);
			}
			//s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);

			if (is_left_pass || is_right_pass) {
				//s->score = s->score < 0 ? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
			//fprintf(stderr, "(%d,%d,%d) ", total_len, total_err,(int) ceil(0.05 * (double) total_len));
		}
		c->n = k;
	}
	//fprintf(stderr, "\n");
}
void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn,
		mem_chain_t *a) {
	double min_l = opt->min_chain_weight ?
		MEM_HSP_COEF * opt->min_chain_weight :
		MEM_MINSC_COEF * log(l_query);
	int i, j, k, min_HSP_score = (int) (opt->a * min_l + .499);
	if (min_l > MEM_SEEDSW_COEF * l_query)
		return; // don't run the following for short reads
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			mem_seed_t *s = &c->seeds[j];
			s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);
			if (s->score < 0 || s->score >= min_HSP_score) {
				s->score = s->score < 0 ? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
		}
		c->n = k;
	}
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/

static inline int cal_max_gap(const mem_opt_t *opt, int qlen) {
	int l_del = (int) ((double) (qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int) ((double) (qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins ? l_del : l_ins;
	l = l > 1 ? l : 1;
	return l < opt->w << 1 ? l : opt->w << 1;
}

#define MAX_BAND_TRY  2

typedef kvec_t(uint8_t*) seq_ptr_arr;
typedef kvec_t(int) seq_lens;
typedef kvec_t(uint8_t) seq_batch;
typedef kvec_t(int) seq_offsets;

//typedef struct {
//		gasal_gpu_storage_t *gpu_storage;
//		int batch_size;
//		int batch_start;
//		int is_active;
//		int no_extend;
//		//int32_t *max_score, *read_start, *read_end, *ref_start, *ref_end;
//		int n_query_batch, n_target_batch, n_seqs;
//	}gpu_batch;
//typedef kvec_t(int) aln_pair;

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, int *curr_read_offset, int *curr_ref_offset, gpu_batch *curr_gpu_batch)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = NULL;
	uint64_t *srt;
	if (c->n == 0) return;
	//extern uint64_t *no_of_extensions;

	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) 
	{
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) 
			max = t->len;
	}

	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) // crossing the forward-reverse boundary; then choose one side
	{
		if (c->seeds[0].rbeg < l_pac) 
			rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else 
			rmax[0] = l_pac;
	}

    /*
	    // DO NOT retrieve the reference sequence
        // ==> the sequence should be retrieved for the GPU, see below bns_fetch_seq_gpu
	    rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	    assert(c->rid == rid);
    */
   
	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) 
    {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];

        
		for (i = 0; i < av->n; ++i)  // test whether extension has been made before
        { 
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb_est || s->rbeg + s->len > p->re_est || s->qbeg < p->qb_est || s->qbeg + s->len > p->qe_est) continue; // not fully contained 
            ///ORIGINALLY: if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue;
			if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb_est; rd = s->rbeg - p->rb_est; 
            ///ORIGINALLY: qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < p->w? max_gap : p->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe_est - (s->qbeg + s->len); rd = p->re_est - (s->rbeg + s->len);
            ///ORIGINALLY: qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < p->w? max_gap : p->w;
			if (qd - rd < w && rd - qd < w) break;
		}
		if (i < av->n) 
        { 
            // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
						k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODo: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}

		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

        // building estimate values for seed filtering (estimates of the beginning and end of the alnignment on ref and read)
        
		int fwd = 0.85*(l_query - (s->qbeg + s->len));
		a->qe_est = ((s->qbeg + s->len) + fwd)  < l_query ? ((s->qbeg + s->len) + fwd) : l_query;
		a->re_est = ((s->rbeg + s->len) + fwd)  < l_pac << 1 ? ((s->rbeg + s->len) + fwd) : l_pac << 1;
		int back = 0.85*(s->qbeg + 1);
		a->qb_est =  (s->qbeg - back) > 0 ? (s->qbeg - back) : 0;
		a->rb_est =  (s->rbeg - back) > 0 ? (s->rbeg - back) : 0;
		if (a->rb_est < l_pac && l_pac < a->qe_est) // crossing the forward-reverse boundary; then choose one side
        { 
			if (s->rbeg < l_pac)
				a->re_est = l_pac; // this works because all seeds are guaranteed to be on the same strand
			else
				a->rb_est = l_pac;
		}

		// fetch data around the seed (for the reference. You always take the whole seed for read.)
        rmax[0] = l_pac<<1; rmax[1] = 0;
        rmax[0] = s->rbeg - (s->qbeg + cal_max_gap(opt, s->qbeg));
        rmax[1] = s->rbeg + s->len + ((l_query - s->qbeg - s->len) + cal_max_gap(opt, l_query - s->qbeg - s->len));
        if (s->len > max) max = s->len;
        rmax[0] = rmax[0] > 0? rmax[0] : 0;
        rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
        if (rmax[0] < l_pac && l_pac < rmax[1]) // crossing the forward-reverse boundary; then choose one side
        { 
            if (s->rbeg < l_pac)
                rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
            else 
                rmax[0] = l_pac;
        }

        /*
            // retrieve the reference sequence
            rseq = bns_fetch_seq(bns, pac, &rmax[0], s->rbeg, &rmax[1], &rid);
            assert(c->rid == rid);
        */

        // this is to check if taking "around the sequence" results in "taking the whole sequence" or not. If not... go through that if. Should be modified.
		if(s->len != l_query) {
			int64_t rmax[2];
			uint8_t *rseq = 0;
			int rid;
			rmax[0] = l_pac<<1; rmax[1] = 0;
			rmax[0] = s->rbeg - (s->qbeg + cal_max_gap(opt, s->qbeg));
			rmax[1] = s->rbeg + s->len + ((l_query - s->qbeg - s->len) + cal_max_gap(opt, l_query - s->qbeg - s->len));
			rmax[0] = rmax[0] > 0? rmax[0] : 0;
			rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
			if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
				if (s->rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
				else rmax[0] = l_pac;
			}
			// retrieve the reference sequence - from the index directly to the GPU storage.
			// ===TODO: modify this to distinguish between left and right.
            bns_fetch_seq_gpu(bns, pac, &rmax[0], s->rbeg, &rmax[1], &rid, curr_gpu_batch);
			assert(c->rid == rid);
			/*
				int rseq_beg, rseq_end;
				rseq_beg = s->rbeg - (s->qbeg + cal_max_gap(opt, s->qbeg)) - rmax[0];
				rseq_end = s->rbeg + s->len + ((l_query - s->qbeg - s->len) + cal_max_gap(opt, l_query - s->qbeg - s->len)) - rmax[0];
				rseq_beg = rseq_beg > 0 ? rseq_beg : 0;
				rseq_end = rseq_end < (rmax[1] - rmax[0]) ? rseq_end : (rmax[1] - rmax[0]);
			*/
			//uint8_t* rs = malloc(rseq_end - rseq_beg);
			int ref_l_seq_with_p = ((rmax[1] - rmax[0])%8) ? (rmax[1] - rmax[0]) + (8 - ((rmax[1] - rmax[0])%8)) : (rmax[1] - rmax[0]) ;
			int j;
			/*
				for (i = 0, j = 0; i < (rmax[1] - rmax[0]) && j < MAX_SEQ_LEN; ++i, ++j) {
					//rs[j] = rseq[i];
					//kv_push(uint8_t, *ref_seq_batch, rseq[i]);
					if (curr_gpu_batch->n_target_batch < curr_gpu_batch->gpu_storage->host_max_target_batch_bytes) curr_gpu_batch->gpu_storage->host_unpacked_target_batch[curr_gpu_batch->n_target_batch++] = rseq[i];
					else {
						fprintf(stderr, "The size of host target_batch (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_target_batch + 1, curr_gpu_batch->gpu_storage->host_max_target_batch_bytes);
							exit(EXIT_FAILURE);
					}
				}
			*/
			int ref_l_seq = rmax[1] - rmax[0];
			while (ref_l_seq < ref_l_seq_with_p) {
				//kv_push(uint8_t, *ref_seq_batch, 0);
				if (curr_gpu_batch->n_target_batch < curr_gpu_batch->gpu_storage->host_max_target_batch_bytes) 
				{
					// J.L. 2018-12-20 16:17 DONE : create some function to add a single base
					// J.L. 2019-12-20 12:35  emulating non-extensible memory host: curr_gpu_batch->gpu_storage->extensible_host_unpacked_target_batch->data[curr_gpu_batch->n_target_batch++] = 4;
                    curr_gpu_batch->n_target_batch = gasal_host_batch_addbase(curr_gpu_batch->gpu_storage, 
                                                    curr_gpu_batch->n_target_batch, 
                                                    4,
                                                    TARGET);
					//fprintf(stderr, "curr_gpu_batch->n_target_batch goes from %d to %d\n", tmp, curr_gpu_batch->n_target_batch);
				} else {
					fprintf(stderr, "The size of host target_batch (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_target_batch + 1, curr_gpu_batch->gpu_storage->host_max_target_batch_bytes);
					exit(EXIT_FAILURE);
				}
				ref_l_seq++;
			}
			a->rseq_beg = rmax[0] /*+ rseq_beg*/;
			if (bwa_verbose >= 4)
				err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long) s->len, (long) s->qbeg, (long) s->rbeg, bns->anns[c->rid].name);
           

			if (curr_gpu_batch->n_seqs < curr_gpu_batch->gpu_storage->host_max_n_alns) 
                curr_gpu_batch->gpu_storage->host_query_batch_lens[curr_gpu_batch->n_seqs] = l_query;
			else {
				fprintf(stderr, "The size of host lens1 (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_seqs + 1, curr_gpu_batch->gpu_storage->host_max_n_alns);
				exit(EXIT_FAILURE);
			}
			//kv_push(int, *read_seq_offsets, *curr_read_offset);
			if (curr_gpu_batch->n_seqs < curr_gpu_batch->gpu_storage->host_max_n_alns) 
                curr_gpu_batch->gpu_storage->host_query_batch_offsets[curr_gpu_batch->n_seqs] = *curr_read_offset;
			else {
				fprintf(stderr, "The size of host offsets1 (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_seqs + 1, curr_gpu_batch->gpu_storage->host_max_n_alns);
				exit(EXIT_FAILURE);
			}
			//kv_push(int, *ref_seq_lens, (rmax[1] - rmax[0]));
			if (curr_gpu_batch->n_seqs < curr_gpu_batch->gpu_storage->host_max_n_alns) 
                curr_gpu_batch->gpu_storage->host_target_batch_lens[curr_gpu_batch->n_seqs] = rmax[1] - rmax[0];
			else {
				fprintf(stderr, "The size of host lens2 (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_seqs + 1, curr_gpu_batch->gpu_storage->host_max_n_alns);
				exit(EXIT_FAILURE);
			}
			//kv_push(int, *ref_seq_offsets, *curr_ref_offset);
			if (curr_gpu_batch->n_seqs < curr_gpu_batch->gpu_storage->host_max_n_alns) 
                curr_gpu_batch->gpu_storage->host_target_batch_offsets[curr_gpu_batch->n_seqs] = *curr_ref_offset;
			else {
				fprintf(stderr, "The size of host offsets2 (%d) exceeds the allocation (%d)\n", curr_gpu_batch->n_seqs + 1, curr_gpu_batch->gpu_storage->host_max_n_alns);
				exit(EXIT_FAILURE);
			}
			*curr_ref_offset +=  ref_l_seq_with_p;
			curr_gpu_batch->n_seqs++;
			free(rseq);
		} else {
			a->score = a->truesc = s->score, a->qb = 0, a->rb = s->rbeg, a->qe = l_query, a->re = s->rbeg + s->len;
		}
        /*
            // Original alignment. 
            if (bwa_verbose >= 4) 
                err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
            if (s->qbeg) { // left extension
                uint8_t *rs, *qs;
                int qle, tle, gtle, gscore;

                // ===NOTE: create local copy of beginning of query, reversed (because aligning left side, so from right to left!!)
                qs = malloc(s->qbeg);
                for (i = 0; i < s->qbeg; ++i) 
                    qs[i] = query[s->qbeg - 1 - i];

                // ===NOTE: create local copy of beginning of reference, reversed (because aligning left side, so from right to left!!)
                tmp = s->rbeg - rmax[0];
                rs = malloc(tmp);
                for (i = 0; i < tmp; ++i) 
                    rs[i] = rseq[tmp - 1 - i];
                
                for (i = 0; i < MAX_BAND_TRY; ++i) {
                    int prev = a->score;
                    aw[0] = opt->w << i;
                    if (bwa_verbose >= 4) {
                        int j;
                        printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
                        printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
                    }
                    a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
                    if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
                    if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
                }
                // check whether we prefer to reach the end of the query
                if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
                    a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
                    a->truesc = a->score;
                } else { // to-end extension
                    a->qb = 0, a->rb = s->rbeg - gtle;
                    a->truesc = gscore;
                }
                free(qs); free(rs);
            } else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

            if (s->qbeg + s->len != l_query) { // right extension
                int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
                qe = s->qbeg + s->len;
                re = s->rbeg + s->len - rmax[0];
                assert(re >= 0);
                for (i = 0; i < MAX_BAND_TRY; ++i) {
                    int prev = a->score;
                    aw[1] = opt->w << i;
                    if (bwa_verbose >= 4) {
                        int j;
                        printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
                        printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
                    }
                    a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
                    if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
                    if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
                }
                // similar to the above
                if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
                    a->qe = qe + qle, a->re = rmax[0] + re + tle;
                    a->truesc += a->score - sc0;
                } else { // to-end extension
                    a->qe = l_query, a->re = rmax[0] + re + gtle;
                    a->truesc += gscore - sc0;
                }
            } else a->qe = l_query, a->re = s->rbeg + s->len;

            if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);
        */

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1] ? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(srt);
}


	/*****************************
	 * Basic hit->SAM conversion *
	 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r) {
    int w;
    if (l1 == l2 && l1 * a - score < (q + r - a) << 1)
        return 0; // to get equal alignment length, we need at least two gaps
    w = ((double) ((l1 < l2 ? l1 : l2) * a - score - q) / r + 2.);
    if (w < abs(l1 - l2))
        w = abs(l1 - l2);
    return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar) {
    int k, l;
    for (k = l = 0; k < n_cigar; ++k) {
        int op = cigar[k] & 0xf;
        if (op == 0 || op == 2)
            l += cigar[k] >> 4;
    }
    return l;
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_) {
    int i, l_name;
    mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

    if (m_)
        mtmp = *m_, m = &mtmp;
    // set flag
    p->flag |= m ? 0x1 : 0; // is paired in sequencing
    p->flag |= p->rid < 0 ? 0x4 : 0; // is mapped
    p->flag |= m && m->rid < 0 ? 0x8 : 0; // is mate mapped
    if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
        p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
    if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
        m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
    p->flag |= p->is_rev ? 0x10 : 0; // is on the reverse strand
    p->flag |= m && m->is_rev ? 0x20 : 0; // is mate on the reverse strand

    // print up to CIGAR
    l_name = strlen(s->name);
    ks_resize(str, str->l + s->l_seq + l_name + (s->qual ? s->l_seq : 0) + 20);
    kputsn(s->name, l_name, str);
    kputc('\t', str); // QNAME
    kputw((p->flag & 0xffff) | (p->flag & 0x10000 ? 0x100 : 0), str);
    kputc('\t', str); // FLAG
    if (p->rid >= 0) { // with coordinate
        kputs(bns->anns[p->rid].name, str);
        kputc('\t', str); // RNAME
        kputl(p->pos + 1, str);
        kputc('\t', str); // POS
        kputw(p->mapq, str);
        kputc('\t', str); // MAPQ
        if (p->n_cigar) { // aligned
            for (i = 0; i < p->n_cigar; ++i) {
                int c = p->cigar[i] & 0xf;
                if (!(opt->flag & MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
                    c = which ? 4 : 3; // use hard clipping for supplementary alignments
                kputw(p->cigar[i] >> 4, str);
                kputc("MIDSH"[c], str);
            }
        } else
            kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
    } else
        kputsn("*\t0\t0\t*", 7, str); // without coordinte
    kputc('\t', str);

    // print the mate position if applicable
    if (m && m->rid >= 0) {
        if (p->rid == m->rid)
            kputc('=', str);
        else
            kputs(bns->anns[m->rid].name, str);
        kputc('\t', str);
        kputl(m->pos + 1, str);
        kputc('\t', str);
        if (p->rid == m->rid) {
            int64_t p0 = p->pos + (p->is_rev ? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
            int64_t p1 = m->pos + (m->is_rev ? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
            if (m->n_cigar == 0 || p->n_cigar == 0)
                kputc('0', str);
            else
                kputl(-(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0)), str);
        } else
            kputc('0', str);
    } else
        kputsn("*\t0\t0", 5, str);
    kputc('\t', str);

    // print SEQ and QUAL
    if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
        kputsn("*\t*", 3, str);
    } else if (!p->is_rev) { // the forward strand
        int i, qb = 0, qe = s->l_seq;
        if (p->n_cigar && which && !(opt->flag & MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
            if ((p->cigar[0] & 0xf) == 4 || (p->cigar[0] & 0xf) == 3)
                qb += p->cigar[0] >> 4;
            if ((p->cigar[p->n_cigar - 1] & 0xf) == 4 || (p->cigar[p->n_cigar - 1] & 0xf) == 3)
                qe -= p->cigar[p->n_cigar - 1] >> 4;
        }
        ks_resize(str, str->l + (qe - qb) + 1);
        for (i = qb; i < qe; ++i)
            str->s[str->l++] = "ACGTN"[(int) s->seq[i]];
        kputc('\t', str);
        if (s->qual) { // printf qual
            ks_resize(str, str->l + (qe - qb) + 1);
            for (i = qb; i < qe; ++i)
                str->s[str->l++] = s->qual[i];
            str->s[str->l] = 0;
        } else
            kputc('*', str);
    } else { // the reverse strand
        int i, qb = 0, qe = s->l_seq;
        if (p->n_cigar && which && !(opt->flag & MEM_F_SOFTCLIP) && !p->is_alt) {
            if ((p->cigar[0] & 0xf) == 4 || (p->cigar[0] & 0xf) == 3)
                qe -= p->cigar[0] >> 4;
            if ((p->cigar[p->n_cigar - 1] & 0xf) == 4 || (p->cigar[p->n_cigar - 1] & 0xf) == 3)
                qb += p->cigar[p->n_cigar - 1] >> 4;
        }
        ks_resize(str, str->l + (qe - qb) + 1);
        for (i = qe - 1; i >= qb; --i)
            str->s[str->l++] = "TGCAN"[(int) s->seq[i]];
        kputc('\t', str);
        if (s->qual) { // printf qual
            ks_resize(str, str->l + (qe - qb) + 1);
            for (i = qe - 1; i >= qb; --i)
                str->s[str->l++] = s->qual[i];
            str->s[str->l] = 0;
        } else
            kputc('*', str);
    }

    // print optional tags
    if (p->n_cigar) {
        kputsn("\tNM:i:", 6, str);
        kputw(p->NM, str);
        kputsn("\tMD:Z:", 6, str);
        kputs((char*) (p->cigar + p->n_cigar), str);
    }
    if (p->score >= 0) {
        kputsn("\tAS:i:", 6, str);
        kputw(p->score, str);
    }
    if (p->sub >= 0) {
        kputsn("\tXS:i:", 6, str);
        kputw(p->sub, str);
    }
    if (bwa_rg_id[0]) {
        kputsn("\tRG:Z:", 6, str);
        kputs(bwa_rg_id, str);
    }
    if (!(p->flag & 0x100)) { // not multi-hit
        for (i = 0; i < n; ++i)
            if (i != which && !(list[i].flag & 0x100))
                break;
        if (i < n) { // there are other primary hits; output them
            kputsn("\tSA:Z:", 6, str);
            for (i = 0; i < n; ++i) {
                const mem_aln_t *r = &list[i];
                int k;
                if (i == which || (r->flag & 0x100))
                    continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
                kputs(bns->anns[r->rid].name, str);
                kputc(',', str);
                kputl(r->pos + 1, str);
                kputc(',', str);
                kputc("+-"[r->is_rev], str);
                kputc(',', str);
                for (k = 0; k < r->n_cigar; ++k) {
                    kputw(r->cigar[k] >> 4, str);
                    kputc("MIDSH"[r->cigar[k] & 0xf], str);
                }
                kputc(',', str);
                kputw(r->mapq, str);
                kputc(',', str);
                kputw(r->NM, str);
                kputc(';', str);
            }
        }
        if (p->alt_sc > 0)
            ksprintf(str, "\tpa:f:%.3f", (double) p->score / p->alt_sc);
    }
    if (p->XA) {
        kputsn("\tXA:Z:", 6, str);
        kputs(p->XA, str);
    }
    if (s->comment) {
        kputc('\t', str);
        kputs(s->comment, str);
    }
    if ((opt->flag & MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
        int tmp;
        kputsn("\tXR:Z:", 6, str);
        tmp = str->l;
        kputs(bns->anns[p->rid].anno, str);
        for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
            if (str->s[i] == '\t')
                str->s[i] = ' ';
    }
    kputc('\n', str);
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a) {
    int mapq, l, sub = a->sub ? a->sub : opt->min_seed_len * opt->a;
    double identity;
    sub = a->csub > sub ? a->csub : sub;
    if (sub >= a->score)
        return 0;
    l = a->qe - a->qb > a->re - a->rb ? a->qe - a->qb : a->re - a->rb;
    identity = 1. - (double) (l * opt->a - a->score) / (opt->a + opt->b) / l;
    if (a->score == 0) {
        mapq = 0;
    } else if (opt->mapQ_coef_len > 0) {
        double tmp;
        tmp = l < opt->mapQ_coef_len ? 1. : opt->mapQ_coef_fac / log(l);
        tmp *= identity * identity;
        mapq = (int) (6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
    } else {
        mapq = (int) (MEM_MAPQ_COEF * (1. - (double) sub / a->score) * log(a->seedcov) + .499);
        mapq = identity < 0.95 ? (int) (mapq * identity * identity + .499) : mapq;
    }
    if (a->sub_n > 0)
        mapq -= (int) (4.343 * log(a->sub_n + 1) + .499);
    if (mapq > 60)
        mapq = 60;
    if (mapq < 0)
        mapq = 0;
    mapq = (int) (mapq * (1. - a->frac_rep) + .499);
    return mapq;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
int read_no = 0;
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m) {
    // J.L. 2019-01-10 09:14 removed proto to place it in extern "C"
    //extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
    kstring_t str;
    kvec_t(mem_aln_t) aa;
    int k,
        l;
    char **XA = 0;
    //read_no++;
    //fprintf(stderr, "%d\n", read_no);
    //fflush(stderr);
    if (!(opt->flag & MEM_F_ALL))
        XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
    kv_init(aa);
    str.l = str.m = 0;
    str.s = 0;
    for (k = l = 0; k < a->n; ++k) {
        mem_alnreg_t *p = &a->a[k];
        mem_aln_t *q;
        if (p->score < opt->T)
            continue;
        if (p->secondary >= 0 && (p->is_alt || !(opt->flag & MEM_F_ALL)))
            continue;
        if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio)
            continue;
        q = kv_pushp(mem_aln_t, aa);
        *q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
        assert(q->rid >= 0); // this should not happen with the new code
        q->XA = XA ? XA[k] : 0;
        q->flag |= extra_flag; // flag secondary
        if (p->secondary >= 0)
            q->sub = -1; // don't output sub-optimal score
        if (l && p->secondary < 0) // if supplementary
            q->flag |= (opt->flag & MEM_F_NO_MULTI) ? 0x10000 : 0x800;
        if (l && !p->is_alt && q->mapq > aa.a[0].mapq)
            q->mapq = aa.a[0].mapq;
        ++l;
    }
    if (aa.n == 0) { // no alignments good enough; then write an unaligned record
        mem_aln_t t;
        t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
        t.flag |= extra_flag;
        mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);
    } else {
        for (k = 0; k < aa.n; ++k)
            mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m);
        for (k = 0; k < aa.n; ++k)
            free(aa.a[k].cigar);
        free(aa.a);
    }
    s->sam = str.s;
    if (XA) {
        for (k = 0; k < a->n; ++k)
            free(XA[k]);
        free(XA);
    }
}

void  print_seq(int length, uint8_t* seq){
    int i;
    fprintf(stderr,"seq length = %d: ", length);
    for (i = 0; i < length; ++i) {
        putc("ACGTN"[(int)seq[i]], stderr);
    }
    fprintf(stderr,"\n");
}

void mem_gasal_fill(gpu_batch *gpu_batch_arr, int gpu_batch_arr_idx, int read_l_seq, char *read_seq, int read_l_seq_with_p)
{
	 
	int i;

	for (i = 0; i < read_l_seq; ++i)
	{
		if (gpu_batch_arr[gpu_batch_arr_idx].n_query_batch < gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->host_max_query_batch_bytes) 
		{
			// J.L. 2018-12-20 16:23 DONE : add some function to add a single base
			// J.L. 2019-01-18 12:40 Emulating non-extensible host memory: gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->extensible_host_unpacked_query_batch->data[gpu_batch_arr[gpu_batch_arr_idx].n_query_batch++]=read_seq[i];                      
			gpu_batch_arr[gpu_batch_arr_idx].n_query_batch = gasal_host_batch_addbase(gpu_batch_arr[gpu_batch_arr_idx].gpu_storage, 
																gpu_batch_arr[gpu_batch_arr_idx].n_query_batch, 
																read_seq[i],
																QUERY);
		} else {
			fprintf(stderr, "The size of host query_batch (%d) exceeds the allocation (%d)\n", gpu_batch_arr[gpu_batch_arr_idx].n_query_batch + 1, gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->host_max_query_batch_bytes);
			exit(EXIT_FAILURE);
		}
		//kv_push(uint8_t, read_seq_batch, read_seq[i]);
	}
	// ===NOTE: padder for the data structure : extensible_host_unpacked_query_batch
	
	while(read_l_seq < read_l_seq_with_p) {
		//kv_push(uint8_t, read_seq_batch, 0);
		if (gpu_batch_arr[gpu_batch_arr_idx].n_query_batch < gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->host_max_query_batch_bytes)
		{
			// J.L. 2018-12-20 17:00 DONE : add some function to add a single base
			// J.L. 2019-01-18 12:40 Emulating non-extensible host memory : gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->extensible_host_unpacked_query_batch->data[gpu_batch_arr[gpu_batch_arr_idx].n_query_batch++]= 4;
			gpu_batch_arr[gpu_batch_arr_idx].n_query_batch = gasal_host_batch_addbase(gpu_batch_arr[gpu_batch_arr_idx].gpu_storage, 
																gpu_batch_arr[gpu_batch_arr_idx].n_query_batch, 
																4, 
																QUERY);
		}
		else {
			fprintf(stderr, "The size of host query_batch (%d) exceeds the allocation (%d)\n", gpu_batch_arr[gpu_batch_arr_idx].n_query_batch + 1, gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->host_max_query_batch_bytes);
			exit(EXIT_FAILURE);
		}
		read_l_seq++;
	}
}

//#define GPU_READ_BATCH_SIZE 1000
void mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *seq, void *buf, int batch_size, int batch_start_idx, mem_alnreg_v *w_regs, int tid, gasal_gpu_storage_v *gpu_storage_vec) {
    int j,  r;
    //extern double *extension_time;
    extern time_struct *extension_time;
    extern uint64_t *no_of_extensions;

	double full_mem_aln1_core;
	double full_mem_chain2aln;
	double chain_preprocess;
	double time_mem_chain, time_mem_chain_flt, time_mem_flt_chained_seeds;

	full_mem_aln1_core = realtime();

    kvec_t(mem_alnreg_v) regs_vec;
    kv_init(regs_vec);
    kv_resize(mem_alnreg_v, regs_vec, batch_size);

    int GPU_READ_BATCH_SIZE;
    if (batch_size >= 4000) GPU_READ_BATCH_SIZE = 1000;
    else {
        GPU_READ_BATCH_SIZE = (int)ceil((double)batch_size/(double)4) % 2 ? (int)ceil((double)batch_size/(double)4) + 1: (int)ceil((double)batch_size/(double)4);
    }
    int internal_batch_count = 0;
    internal_batch_count = (int)ceil(((double)batch_size)/((double)(GPU_READ_BATCH_SIZE)));
    gpu_batch gpu_batch_arr[gpu_storage_vec->n];

    for(j = 0; j < gpu_storage_vec->n; j++) {
        gpu_batch_arr[j].gpu_storage = &(gpu_storage_vec->a[j]);
    }
    
	int internal_batch_done = 0;
    int batch_processed = 0;
    //int total_internal_batches = 0;
    int internal_batch_no = 0;
    double time_extend;

    while (internal_batch_done < internal_batch_count) {
        int gpu_batch_arr_idx = 0;
        while(gpu_batch_arr_idx != gpu_storage_vec->n && gpu_batch_arr[gpu_batch_arr_idx].gpu_storage->is_free != 1) {
            gpu_batch_arr_idx++;
        }

        int internal_batch_start_idx = batch_processed;
        if (internal_batch_start_idx < batch_size && gpu_batch_arr_idx < gpu_storage_vec->n) 
		{
            gpu_batch_arr[gpu_batch_arr_idx].n_query_batch = 0;
            gpu_batch_arr[gpu_batch_arr_idx].n_target_batch = 0;
            gpu_batch_arr[gpu_batch_arr_idx].n_seqs = 0;
            int curr_read_offset = 0;
            int curr_ref_offset = 0;
            int internal_batch_size = batch_size - batch_processed >= GPU_READ_BATCH_SIZE  ? GPU_READ_BATCH_SIZE : batch_size - batch_processed;

            for (j = batch_start_idx + internal_batch_start_idx; j < (batch_start_idx + internal_batch_start_idx) + internal_batch_size; ++j) 
			{
                mem_chain_v chn;
                mem_alnreg_v regs;
                int i;
                char *read_seq = seq[j].seq;
				int read_l_seq = seq[j].l_seq;
				int read_l_seq_with_p = read_l_seq + (8 - (read_l_seq%8));
				
				// convert to 2-bit encoding if we have not done so (was before in mem_gasal_fill but must be done before. Might be a bit slower, but separates the concerns.)
				for (i = 0; i < read_l_seq; ++i)
					read_seq[i] = read_seq[i] < 4 ? read_seq[i] : nst_nt4_table[(int) read_seq[i]]; 
				
                // ===NOTE: computing chains, store them in the mem_chain_v chn
                chain_preprocess = realtime();
                time_mem_chain = realtime();
                chn = mem_chain(opt, bwt, bns, seq[j].l_seq, (uint8_t*)(read_seq), buf);
                extension_time[tid].time_mem_chain += (realtime() - time_mem_chain);
                time_mem_chain_flt = realtime();
                chn.n = mem_chain_flt(opt, chn.n, chn.a);
                extension_time[tid].time_mem_chain_flt += (realtime() - time_mem_chain_flt);

                time_mem_flt_chained_seeds = realtime();
                if (opt->shd_filter) 
					mem_shd_flt_chained_seeds(opt, bns, pac, seq[j].l_seq, (uint8_t*)(read_seq), chn.n, chn.a);
                else mem_flt_chained_seeds(opt, bns, pac, seq[j].l_seq, (uint8_t*)(read_seq), chn.n, chn.a);
                if (bwa_verbose >= 4)
                    mem_print_chain(bns, &chn);
                extension_time[tid].time_mem_flt_chained_seeds += (realtime() - time_mem_flt_chained_seeds);

                extension_time[tid].chain_preprocess += (realtime() - chain_preprocess);


                // ===NOTE: CHAINS DONE. COMPUTING ALIGNMENT

                kv_init(regs);
                for (i = 0; i < chn.n; ++i) {
                    mem_chain_t *p = &chn.a[i];
                    if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
                    //mem_chain2aln(opt, bns, pac, seq[j].l_seq, (uint8_t*)(read_seq), p, &regs, &read_seq_lens, &read_seq_offsets, &curr_read_offset, &ref_seq_batch, &ref_seq_lens, &ref_seq_offsets, &curr_ref_offset);
                    /* ===NOTE: it seems like in bwa-gasal2, mem_chain2aln has been cut down to compute fewer things. But I don't know WHAT kind of things.
                        it is probably worth it to modify that part to get back to the previous code partition:
                        - Simpler mem_align1_core, and more meaningful functions like mem_chain2aln
                        - Probably create yet another function to call GASAL2 (exactly like ksw_extend2)
				Note: that function woudln't be integrated in mem_chain2aln because they wouldn't run in the same loop. See how GASAL2 call below is out of the loop call.
				*/
					full_mem_chain2aln = realtime();
                    mem_chain2aln(opt, bns, pac, seq[j].l_seq, (uint8_t*)(read_seq), p, &regs, &curr_read_offset, &curr_ref_offset, &gpu_batch_arr[gpu_batch_arr_idx]);
                    free(chn.a[i].seeds);
					extension_time[tid].full_mem_chain2aln += (realtime() - full_mem_chain2aln);
                }
                curr_read_offset += read_l_seq_with_p;
                free(chn.a);
                kv_push(mem_alnreg_v, regs_vec, regs);

                //smem_aux_destroy((smem_aux_t*)buf);
                //buf = smem_aux_init();


				// ===NOTE: filler/padder for the data structure : extensible_host_unpacked_query_batch (we always fill the whole sequence.)
               	mem_gasal_fill(gpu_batch_arr, gpu_batch_arr_idx, read_l_seq, read_seq, read_l_seq_with_p);
            }



            // ===NOTE: Chains done, mem_chain2aln done (whatever it did)
            // ===NOTE: now, GASAL2 KERNEL LAUCNH ON BATCH

            if (/*kv_size(ref_seq_lens)*/ gpu_batch_arr[gpu_batch_arr_idx].n_seqs > 0) {
                //fprintf(stderr, "n_alns on GPU=%d\n", kv_size(ref_seq_lens));
                //no_of_extensions[tid] += kv_size(ref_seq_lens);
                no_of_extensions[tid] += gpu_batch_arr[gpu_batch_arr_idx].n_seqs;

                time_extend = realtime();
                //J.L. 2018-12-20 17:24 Added params object.
                Parameters *args;
                args = new Parameters(0, NULL);

                args->algo = LOCAL;
                args->start_pos = WITH_START;

                gasal_aln_async(gpu_batch_arr[gpu_batch_arr_idx].gpu_storage, gpu_batch_arr[gpu_batch_arr_idx].n_query_batch, gpu_batch_arr[gpu_batch_arr_idx].n_target_batch, gpu_batch_arr[gpu_batch_arr_idx].n_seqs, args);
                extension_time[tid].aln_kernel += (realtime() - time_extend);
                gpu_batch_arr[gpu_batch_arr_idx].no_extend = 0;
				gpu_batch *cur = &gpu_batch_arr[gpu_batch_arr_idx];
                //fprintf(stderr, "batch launched with batch_size=%d, n_seqs=%d alingments, with n_query_batch=%d, n_target_batch=%d\n", cur->batch_size, cur->n_seqs, cur->n_query_batch, cur->n_target_batch );
            } else {
                gpu_batch_arr[gpu_batch_arr_idx].no_extend = 1;
                //fprintf(stderr, "I am here\n");
                mem_alnreg_v regs = kv_A(regs_vec, kv_size(regs_vec) - 1);
                fprintf(stderr, "Thread no. %d is here with internal batch size %d, regs.n %d \n", tid, internal_batch_size, regs.n);
            }

            batch_processed += internal_batch_size;
            gpu_batch_arr[gpu_batch_arr_idx].batch_size = internal_batch_size;
            gpu_batch_arr[gpu_batch_arr_idx].batch_start = internal_batch_start_idx;
            gpu_batch_arr[gpu_batch_arr_idx].is_active = 1;

            assert(kv_size(regs_vec) == batch_processed);
            internal_batch_no++;
            //fprintf(stderr, "internal batch %d launched\n", internal_batch_no++);

        }

        // ===NOTE: GASAL2 GET RESULT, measure time

        //fprintf(stderr, "Current extension time of %d seeds on GPU by thread no. %d is %.3f usec\n", kv_size(ref_seq_lens), tid,  extension_time[tid]*1e6);
        int internal_batch_idx = 0;
        while (internal_batch_idx != gpu_storage_vec->n) {
            time_extend = realtime();
            int x = 0;
            if (gpu_batch_arr[internal_batch_idx].gpu_storage->is_free != 1) {
                x = (gasal_is_aln_async_done(gpu_batch_arr[internal_batch_idx].gpu_storage) == 0);
                //fprintf(stderr, "Thread no. %d stuck here with batch size %d and batch count %d. internal batch idx is %d \n", tid, batch_size, internal_batch_count, internal_batch_idx);
            }
            if (x) extension_time[tid].get_results_actual += (realtime() - time_extend);
            else if (gpu_batch_arr[internal_batch_idx].gpu_storage->is_free != 1 && x == 0) extension_time[tid].get_results_wasted += (realtime() - time_extend);
            //fprintf(stderr, "Thread no. %d stuck here with batch size %d and batch count %d. internal batch idx is %d, batches launched=%d, batches done=%d \n", tid, batch_size, internal_batch_count, internal_batch_idx, internal_batch_no, internal_batch_done);
            //fprintf(stderr, "Thread no. %d stuck here with batch size %d and batch count %d. internal batch idx is \n");
            if ((x == 1 || gpu_batch_arr[internal_batch_idx].no_extend == 1) && gpu_batch_arr[internal_batch_idx].is_active == 1)
            {
                // J.L. 2018-12-21 15:21 changed calls with best score
                int32_t *max_score = gpu_batch_arr[internal_batch_idx].gpu_storage->host_res->aln_score;
                int32_t *read_start = gpu_batch_arr[internal_batch_idx].gpu_storage->host_res->query_batch_start;
                int32_t  *read_end = gpu_batch_arr[internal_batch_idx].gpu_storage->host_res->query_batch_end;
                int32_t *ref_start = gpu_batch_arr[internal_batch_idx].gpu_storage->host_res->target_batch_start;
                int32_t   *ref_end = gpu_batch_arr[internal_batch_idx].gpu_storage->host_res->target_batch_end;

                int seq_idx=0;
                for(j = 0, r = gpu_batch_arr[internal_batch_idx].batch_start; j < gpu_batch_arr[internal_batch_idx].batch_size; ++j, ++r){
                    int i;
                    mem_alnreg_v regs = kv_A(regs_vec, r);
                    //int read_pos = kv_A(read_seq_offsets, j);
                    //int read_len = kv_A(read_seq_lens, j);
                    //uint8_t* read_seq = &(kv_A(read_seq_batch, read_pos));
                    if (gpu_batch_arr[internal_batch_idx].no_extend == 1) fprintf(stderr, "I am here too as well with regs.n %d\n", regs.n);
                    for(i = 0; i < regs.n; ++i){
                        mem_alnreg_t *a = &regs.a[i];
                        //fprintf(stderr, "I am here before\n");
                        //fprintf(stderr, "r=%d, seq[r].l_seq=%d\n", r, seq[r].l_seq);
                        //fprintf(stderr, "I am here after\n");
                        if (a->seedlen0 != seq[r].l_seq/*kv_A(read_seq_lens, seq_idx)*/) {
                            //if (gpu_batch_arr[internal_batch_idx].no_extend == 1) fprintf(stderr, "I am here too as well\n");
                            a->score = max_score[seq_idx];
                            a->qb = read_start[seq_idx];
                            a->qe = read_end[seq_idx] + 1;
                            a->rb = ref_start[seq_idx] + a->rseq_beg;
                            a->re = ref_end[seq_idx] + a->rseq_beg + 1;
                            a->truesc = max_score[seq_idx];

                            seq_idx++;
                        }
                    }
                    regs.n = mem_sort_dedup_patch(opt, bns, pac,(uint8_t*)(seq[r].seq), regs.n, regs.a);
                    if (bwa_verbose >= 4) {
                        err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
                        for (i = 0; i < regs.n; ++i) {
                            mem_alnreg_t *p = &regs.a[i];
                            printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
                        }
                    }
                    for (i = 0; i < regs.n; ++i) {
                        mem_alnreg_t *p = &regs.a[i];
                        if (p->rid >= 0 && bns->anns[p->rid].is_alt)
                            p->is_alt = 1;
                        //free(kv_A(read_seqns, i));

                    }
                    w_regs[r + batch_start_idx] = regs;
                    //free(regs.a);

                }

                gpu_batch_arr[internal_batch_idx].is_active = 0;
                internal_batch_done++;
                //fprintf(stderr, "internal batch %d done\n", internal_batch_done - 1);
            }

            internal_batch_idx++;
        }
    }
    kv_destroy(regs_vec);
    //fprintf(stderr, "--------------------------------------");
	extension_time[tid].full_mem_aln1_core += (realtime() - full_mem_aln1_core);

}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar) {
    mem_aln_t a;
    int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1 << 30), l_MD;
    score = 0;
    int64_t pos, rb, re;
    uint8_t *query;

    memset(&a, 0, sizeof(mem_aln_t));
    if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
        a.rid = -1;
        a.pos = -1;
        a.flag |= 0x4;
        return a;
    }
    qb = ar->qb, qe = ar->qe;
    rb = ar->rb, re = ar->re;
    query = malloc(l_query);
    for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
        query[i] = query_[i] < 5 ? query_[i] : nst_nt4_table[(int) query_[i]];
    a.mapq = ar->secondary < 0 ? mem_approx_mapq_se(opt, ar) : 0;
    if (ar->secondary >= 0)
        a.flag |= 0x100; // secondary alignment
    tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
    w2 = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
    w2 = w2 > tmp ? w2 : tmp;
    if (bwa_verbose >= 4)
        printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
    if (w2 > opt->w)
        w2 = w2 < ar->w ? w2 : ar->w;
    i = 0;
    a.cigar = NULL;
    do {
        free(a.cigar);
        w2 = w2 < opt->w << 2 ? w2 : opt->w << 2;
        a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*) &query[qb], rb, re,
                &score, &a.n_cigar, &NM);

        //fprintf(stderr, "in do-while: i=%d, a.cigar=%s, a.n_cigar=%d\n",i,a.cigar, a.n_cigar);

        if (bwa_verbose >= 4)
            printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
        if (score == last_sc || w2 == opt->w << 2)
            break; // it is possible that global alignment and local alignment give different scores
        last_sc = score;
        w2 <<= 1;
    } while (++i < 3 && score < ar->truesc - opt->a);
    /*
        a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*) &query[qb], rb, re,
        &score, &a.n_cigar, &NM);
        */
    if (bwa_verbose >= 4)
        fprintf(stderr, "* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
    //fflush(stderr);

    l_MD = strlen((char*) (a.cigar + a.n_cigar)) + 1;

    a.NM = NM;
    pos = bns_depos(bns, rb < bns->l_pac ? rb : re - 1, &is_rev);
    a.is_rev = is_rev;
    if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
        if ((a.cigar[0] & 0xf) == 2) {
            pos += a.cigar[0] >> 4;
            --a.n_cigar;
            memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
        } else if ((a.cigar[a.n_cigar - 1] & 0xf) == 2) {
            --a.n_cigar;
            memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
        }
    }
    if (qb != 0 || qe != l_query) { // add clipping to CIGAR
        int clip5, clip3;
        clip5 = is_rev ? l_query - qe : qb;
        clip3 = is_rev ? qb : l_query - qe;
        a.cigar = realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
        if (clip5) {
            memmove(a.cigar + 1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
            a.cigar[0] = clip5 << 4 | 3;
            ++a.n_cigar;
        }
        if (clip3) {
            memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
            a.cigar[a.n_cigar++] = clip3 << 4 | 3;
        }
    }
    a.rid = bns_pos2rid(bns, pos);
    assert(a.rid == ar->rid);
    a.pos = pos - bns->anns[a.rid].offset;
    a.score = ar->score;
    a.sub = ar->sub > ar->csub ? ar->sub : ar->csub;
    a.is_alt = ar->is_alt;
    a.alt_sc = ar->alt_sc;
    free(query);
    return a;
}

typedef struct {
    const mem_opt_t *opt;
    const bwt_t *bwt;
    const bntseq_t *bns;
    const uint8_t *pac;
    const mem_pestat_t *pes;
    smem_aux_t **aux;
    bseq1_t *seqs;
    mem_alnreg_v *regs;
    int64_t n_processed;
} worker_t;

void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec) {
    worker_t *w = (worker_t*) data;
    
    // worker truncated for gasal - see original function in the original bwamem.

    if (bwa_verbose >= 4)
        printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
    mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs, w->aux[tid], batch_size, i, w->regs, tid, gpu_storage_vec);

}

static void worker2(void *data, int i, int tid, int batch_size, int n_reads, gasal_gpu_storage_v *gpu_storage_vec) {
    extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2],
            mem_alnreg_v a[2]);
    extern void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a);
    worker_t *w = (worker_t*) data;
    if (!(w->opt->flag & MEM_F_PE)) {
        if (bwa_verbose >= 4)
            printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
        mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
        mem_reg2sam(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
        free(w->regs[i].a);
    } else {
        if (bwa_verbose >= 4)
            printf("=====> Finalizing read pair '%s' <=====\n", w->seqs[i << 1 | 0].name);
        mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed >> 1) + i, &w->seqs[i << 1], &w->regs[i << 1]);
        free(w->regs[i << 1 | 0].a);
        free(w->regs[i << 1 | 1].a);
    }
}



void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0) {
    worker_t w;
    mem_pestat_t pes[4];
    double ctime, rtime;
    int i;
    extern double *load_balance_waste_time;
    extern double total_load_balance_waste_time;
    ctime = cputime();
    rtime = realtime();
    global_bns = bns;
    w.regs = malloc(n * sizeof(mem_alnreg_v));
    w.opt = opt;
    w.bwt = bwt;
    w.bns = bns;
    w.pac = pac;
    w.seqs = seqs;
    w.n_processed = n_processed;
    w.pes = &pes[0];
    w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
    for (i = 0; i < opt->n_threads; ++i)
        w.aux[i] = smem_aux_init();

    kt_for(opt->n_threads, worker1, &w, /*(opt->flag & MEM_F_PE) ? n >> 1 : */n); // find mapping positions

    double min_exit_time = load_balance_waste_time[0];
    double max_exit_time = load_balance_waste_time[0];
    for (i = 1; i < opt->n_threads; i++) {
        if  (load_balance_waste_time[i] < min_exit_time) min_exit_time = load_balance_waste_time[i];
        if  (load_balance_waste_time[i] > max_exit_time) max_exit_time = load_balance_waste_time[i];
    }

    total_load_balance_waste_time += (max_exit_time - min_exit_time);

    for (i = 0; i < opt->n_threads; ++i)
        smem_aux_destroy(w.aux[i]);
    free(w.aux);
    if (opt->flag & MEM_F_PE) { // infer insert sizes if not provided
        if (pes0)
            memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
        else
            mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size distribution from data
    }
    kt_for(opt->n_threads, worker2, &w, (opt->flag & MEM_F_PE) ? n >> 1 : n); // generate alignment
    /*
        min_exit_time = load_balance_waste_time[0];
        max_exit_time = load_balance_waste_time[0];
        for (i = 1; i < opt->n_threads; i++) {
        if  (load_balance_waste_time[i] < min_exit_time) min_exit_time = load_balance_waste_time[i];
        if  (load_balance_waste_time[i] > max_exit_time) max_exit_time = load_balance_waste_time[i];
        }
        total_load_balance_waste_time += (max_exit_time - min_exit_time);
        */
    free(w.regs);
    if (bwa_verbose >= 3)
        fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}

