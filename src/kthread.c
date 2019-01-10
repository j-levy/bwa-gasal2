#include <pthread.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "../GASAL2/include/gasal.h"

 
// J.L. 2019-01-10 09:50 adding extern C proto
#ifdef __cplusplus
extern "C"{
	void kt_for(int n_threads, void (*func)(void*,int, int, int, int), void *data, long n);
}
#endif


/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	long i;
	long n_per_thread;
	long orig_read_batch_size;
	long curr_read_batch_size;
	long n_done;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n;
	ktf_worker_t *w;
	void (*func)(void*,int, int, int, int, gasal_gpu_storage_v*);
	void *data;
} kt_for_t;
//#define READ_BATCH_SIZE 5000
//int *is_thread_available;

//static inline long steal_work(kt_for_t *t)
//{
//   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
//        int i, min_i = -1;
//        long k, min = LONG_MAX;
//        for (i = 0; i < t->n_threads; ++i)
//                if (min > t->w[i].i) min = t->w[i].i, min_i = i;
//        if(t->func == &worker1) {
//        	int to_add =  t->n_threads * READ_BATCH_SIZE;//(__sync_add_and_fetch (&t->w[min_i].i, 0) + (t->n_threads * READ_BATCH_SIZE)) >= t->n ? READ_BATCH_SIZE/2 : (t->n_threads * READ_BATCH_SIZE);
//        	k = __sync_fetch_and_add(&t->w[min_i].i, to_add);
//        }
//        else k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads );
//
//        return k >= t->n? -1 : k;
//}

static inline long steal_work(kt_for_t *t)
{
   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
        int i, min_i = -1;
        long k, min = LONG_MAX;

        if(t->func == &worker1) {
        	 for (i = 0; i < t->n_threads; ++i)
        		 if (min > t->w[i].n_done) min = t->w[i].n_done, min_i = i;
        	//int to_add =  t->n_threads * READ_BATCH_SIZE;//(__sync_add_and_fetch (&t->w[min_i].i, 0) + (t->n_threads * READ_BATCH_SIZE)) >= t->n ? READ_BATCH_SIZE/2 : (t->n_threads * READ_BATCH_SIZE);
        	//k = __sync_fetch_and_add(&t->w[min_i].i, to_add);
        	 if ( t->w[min_i].i >= t->n ||  t->w[min_i].i >= (min_i + 1)*t->w[min_i].n_per_thread ) return -1;
        	 else return min_i;
//        	 actual_batch_size =  (k + READ_BATCH_SIZE) >=  ((min_i + 1)*t->w[min_i]->n_per_thread) ? ((min_i + 1)*t->w[min_i]->n_per_thread + 1) - k : READ_BATCH_SIZE;
//        	 actual_batch_size = (k + actual_batch_size) >= t->n ? t->n - k : actual_batch_size;
//        	 t->w[min_i]->n_done = __sync_add_and_fetch((&t->w[min_i]->n_done), actual_batch_size);
//        	 return actual_batch_size;
        }
        else  {
        	for (i = 0; i < t->n_threads; ++i)
        		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
        	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads );
        	return k >= t->n? -1 : k;
        }


}

/*static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w);
	}
	while ((i = steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w);
	pthread_exit(0);
}*/

//static void *ktf_worker(void *data)
//{
//   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
//   ktf_worker_t *w = (ktf_worker_t*)data;
//   int i;
//   extern time_struct *extension_time;
//   extern gasal_gpu_storage_v *gpu_storage_vec_arr;
//
//   //double time_extend = realtime();
//   //gasal_gpu_storage_v gpu_storage_vec = gpu_storage_vec_arr[w - w->t->w];
////   gasal_gpu_storage_v gpu_storage_vec = gasal_init_gpu_storage_v(2);
////   gasal_init_streams_new(&gpu_storage_vec, 1000*300, 40*1000*600, 250*1000, LOCAL, WITH_START);
//   //extension_time[w - w->t->w].gpu_mem_alloc += (realtime() - time_extend);
//
//   for (;;) {
//      if(w->t->func == &worker1){
//         i = __sync_fetch_and_add(&w->i, (w->t->n_threads * READ_BATCH_SIZE));
//         if (i >= w->t->n) break;
//         w->t->func(w->t->data, i, w - w->t->w, READ_BATCH_SIZE, w->t->n, &gpu_storage_vec);
//      }
//      else {
//         i = __sync_fetch_and_add(&w->i, w->t->n_threads);
//         if (i >= w->t->n) break;
//         w->t->func(w->t->data, i, w - w->t->w, READ_BATCH_SIZE, w->t->n, &gpu_storage_vec);
//      }
//   }
//   while ((i = steal_work(w->t)) >= 0){
//	  //if(w->t->func == &worker1) fprintf(stderr, "thread no. %d is stealing work\n", w - w->t->w);
//      w->t->func(w->t->data, i, w - w->t->w, READ_BATCH_SIZE, w->t->n, &gpu_storage_vec);
//   }
//   //if(w->t->func == &worker1) fprintf (stderr, "Thread %d leaving at %.3f secs\n", w - w->t->w, realtime());
//  // time_extend = realtime();
////   gasal_destroy_streams(&gpu_storage_vec);
////   gasal_destroy_gpu_storage_v(&gpu_storage_vec);
//  // extension_time[w - w->t->w].gpu_mem_free += (realtime() - time_extend);
//   if(w->t->func == &worker1) fprintf (stderr, "Thread %d leaving at %.3f secs\n", w - w->t->w, realtime());
//   pthread_exit(0);
//}

int threads_available;
static void *ktf_worker(void *data)
{
   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
   ktf_worker_t *w = (ktf_worker_t*)data;
   int i;
   extern double *load_balance_waste_time;
//   gasal_gpu_storage gpu_storage;
//   gpu_storage.max_batch1_bytes = READ_BATCH_SIZE*40*300;
//   gpu_storage.max_batch2_bytes = READ_BATCH_SIZE*40*600;
//   gpu_storage.max_n_alns = READ_BATCH_SIZE*40;

//   gasal_aln_imp_mem_alloc(&gpu_storage, LOCAL, WITH_START);

   extern gasal_gpu_storage_v *gpu_storage_vec_arr;
   gasal_gpu_storage_v gpu_storage_vec = gpu_storage_vec_arr[w - w->t->w];
   int READ_BATCH_SIZE;
   for (;;) {
      if(w->t->func == &worker1){
    	 //READ_BATCH_SIZE = (int)ceil((double)w->orig_read_batch_size/(double)(threads_available + 1)) % 2 ?  (int)ceil((double)w->orig_read_batch_size/(double)(threads_available + 1)) + 1 :  (int)ceil((double)w->orig_read_batch_size/(double)(threads_available + 1));
    	 if(w->n_done >= w->n_per_thread) break;
    	 if ((w->n_per_thread - w->n_done) <= 1600) READ_BATCH_SIZE = (w->n_per_thread - w->n_done);
    	 else READ_BATCH_SIZE = (int)ceil((double)(w->n_per_thread - w->n_done)/(double)(threads_available + 1)) % 2 ?  (int)ceil((double)(w->n_per_thread - w->n_done)/(double)(threads_available + 1)) + 1 :  (int)ceil((double)(w->n_per_thread - w->n_done)/(double)(threads_available + 1));
    	 READ_BATCH_SIZE = READ_BATCH_SIZE > w->orig_read_batch_size ? w->orig_read_batch_size : READ_BATCH_SIZE;
    	 READ_BATCH_SIZE = READ_BATCH_SIZE < 1000 ? 1000 : READ_BATCH_SIZE;
         i = __sync_fetch_and_add(&w->i, (READ_BATCH_SIZE));
         if (i >= w->t->n || i >= ((w - w->t->w) + 1)*w->n_per_thread ) break;
         __sync_add_and_fetch(&(w->n_done), (READ_BATCH_SIZE));
         int actual_batch_size =  i + READ_BATCH_SIZE >=  ((w - w->t->w) + 1)*w->n_per_thread ? ((w - w->t->w) + 1)*w->n_per_thread - i : READ_BATCH_SIZE;
         actual_batch_size = i + actual_batch_size >= w->t->n ? w->t->n - i : actual_batch_size;

         w->t->func(w->t->data, i, w - w->t->w, actual_batch_size, w->t->n, &gpu_storage_vec);
      }
      else {
         i = __sync_fetch_and_add(&w->i, w->t->n_threads);
         if (i >= w->t->n) break;
         w->t->func(w->t->data, i, w - w->t->w, READ_BATCH_SIZE, w->t->n, &gpu_storage_vec);
      }
   }
   //is_thread_available[w - w->t->w] = 1;
   if(w->t->func == &worker1){
	   __sync_add_and_fetch(&(threads_available), 1);
	   int min_i;
	   //fprintf(stderr, "Threads available = %d\n", threads_available);
	   while ((min_i = steal_work(w->t)) >= 0){
		   //is_thread_available[w - w->t->w] = 0;
		   //if(w->t->func == &worker1) fprintf(stderr, "thread no. %d is here processing i %d to %d\n", w - w->t->w, i, i + READ_BATCH_SIZE);
		   //READ_BATCH_SIZE = (int)ceil((double)w->t->w[min_i].orig_read_batch_size/(double)(threads_available + 1)) % 2 ?  (int)ceil((double)w->t->w[min_i].orig_read_batch_size/(double)(threads_available + 1)) + 1 :  (int)ceil((double)w->t->w[min_i].orig_read_batch_size/(double)(threads_available + 1));
		   if(w->t->w[min_i].n_done >= w->t->w[min_i].n_per_thread) break;
		   if ((w->t->w[min_i].n_per_thread - w->t->w[min_i].n_done) <= 1600) READ_BATCH_SIZE = w->t->w[min_i].n_per_thread - w->t->w[min_i].n_done;
		   else READ_BATCH_SIZE = (int)ceil((double)(w->t->w[min_i].n_per_thread - w->t->w[min_i].n_done)/(double)(threads_available + 1)) % 2 ?  (int)ceil((double)(w->t->w[min_i].n_per_thread - w->t->w[min_i].n_done)/(double)(threads_available + 1)) + 1 :  (int)ceil((double)(w->t->w[min_i].n_per_thread - w->t->w[min_i].n_done)/(double)(threads_available + 1));
		   READ_BATCH_SIZE = READ_BATCH_SIZE > w->t->w[min_i].orig_read_batch_size ? w->t->w[min_i].orig_read_batch_size : READ_BATCH_SIZE;
		   READ_BATCH_SIZE = READ_BATCH_SIZE < 1000 ? 1000 : READ_BATCH_SIZE;
		   i = __sync_fetch_and_add(&(w->t->w[min_i].i), (READ_BATCH_SIZE));
		   if (i >= w->t->n || i >= (min_i + 1)*w->t->w[min_i].n_per_thread) break;
		   __sync_add_and_fetch(&(w->t->w[min_i].n_done), (READ_BATCH_SIZE));
		   int actual_batch_size =  (i + READ_BATCH_SIZE) >=  ((min_i + 1)*w->t->w[min_i].n_per_thread) ? ((min_i + 1)*w->t->w[min_i].n_per_thread) - i: READ_BATCH_SIZE;
		   actual_batch_size = (i + actual_batch_size) >= w->t->n ? w->t->n - i: actual_batch_size;
		   //fprintf(stderr, "Thread %d is doing a batch of %d on the behalf of thread %d\n", w - w->t->w, actual_batch_size, min_i);
		   w->t->func(w->t->data, i, w - w->t->w, actual_batch_size, w->t->n, &gpu_storage_vec);

	   }
   }
   else {
	   while ((i = steal_work(w->t)) >= 0){
		   //is_thread_available[w - w->t->w] = 0;
		   //if(w->t->func == &worker1) fprintf(stderr, "thread no. %d is here processing i %d to %d\n", w - w->t->w, i, i + READ_BATCH_SIZE);
		   w->t->func(w->t->data, i, w - w->t->w, READ_BATCH_SIZE, w->t->n, &gpu_storage_vec);
	   }
   }

   //gasal_aln_imp_mem_free(&gpu_storage);
   //if(w->t->func == &worker1) fprintf (stderr, "Thread %d exits at %.3f secs\n", w - w->t->w, realtime());
   if(w->t->func == &worker1) load_balance_waste_time[w - w->t->w] = realtime();
   pthread_exit(0);
}

//void kt_for(int n_threads, void (*func)(void*,int, int, int, int), void *data, long n)
//{
//   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
//        int i;
//        kt_for_t t;
//        pthread_t *tid;
//        t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
//        t.w = (ktf_worker_t*)alloca(n_threads * sizeof(ktf_worker_t));
//        tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
//        for (i = 0; i < n_threads; ++i)
//                t.w[i].t = &t, t.w[i].i = (func == &worker1) ? i*READ_BATCH_SIZE : i;
//        for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
//        for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
//}


void kt_for(int n_threads, void (*func)(void*,int, int, int, int), void *data, long n)
{
   extern void worker1(void *data, int i, int tid, int batch_size, int total_reads, gasal_gpu_storage_v *gpu_storage_vec);
	int i;
	kt_for_t t;
	pthread_t *tid;
	t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
	t.w = (ktf_worker_t*)alloca(n_threads * sizeof(ktf_worker_t));
	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	threads_available = 0;
	for (i = 0; i < n_threads; ++i) {
		t.w[i].n_per_thread = ((int)ceil((double)n/(double)n_threads)) % 2 ? (int)ceil((double)n/(double)n_threads) + 1 : (int)ceil((double)n/(double)n_threads);

		t.w[i].t = &t, t.w[i].i = (func == &worker1) ? i*t.w[i].n_per_thread : i;//i*READ_BATCH_SIZE : i;
		t.w[i].n_done = 0;
		t.w[i].orig_read_batch_size = 5000;

	}
	//fprintf(stderr, "n_per_thread=%d\n", t.w[0].n_per_thread);
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}


/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

static void *ktp_worker(void *data)
{
	ktp_worker_t *w = (ktp_worker_t*)data;
	ktp_t *p = w->pl;
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		w->data = p->func(p->shared, w->step, w->step? w->data : 0); // for the first step, input is NULL

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;
		if (w->step == 0) w->index = p->index++;
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	pthread_t *tid;
	int i;

	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	pthread_mutex_init(&aux.mutex, 0);
	pthread_cond_init(&aux.cv, 0);

	aux.workers = (ktp_worker_t*)alloca(n_threads * sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

	pthread_mutex_destroy(&aux.mutex);
	pthread_cond_destroy(&aux.cv);
}
