#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "utils.h"


#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

int bwa_fa2pac(int argc, char *argv[]);
int bwa_pac2bwt(int argc, char *argv[]);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);
int bwa_index(int argc, char *argv[]);
int bwt_bwtgen_main(int argc, char *argv[]);

int bwa_aln(int argc, char *argv[]);
int bwa_sai2sam_se(int argc, char *argv[]);
int bwa_sai2sam_pe(int argc, char *argv[]);

int bwa_bwtsw2(int argc, char *argv[]);

int main_fastmap(int argc, char *argv[]);
int gase_aln(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);

int main_pemerge(int argc, char *argv[]);
int main_maxk(int argc, char *argv[]);
	
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: BWA-GASAL2 (derived from GASE: Generic Aligner for Seed-and-Extend)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "BWA-GASAL2 is an extension of BWA (version 0.7.13) which is developed by Heng Li.\n");
	fprintf(stderr, "Contact: Nauman Ahmed <n.ahmed@tudelft.nl>\n\n");
	fprintf(stderr, "Usage:   bwa-gasal2 <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         gase_aln           GASE algorithm\n");
	/*fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");*/
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: To use bwa-gasal2, you need to first index the genome with `bwa-gasal2 index'.\n\n");
	return 1;
}


time_struct *extension_time;
uint64_t *no_of_extensions;
double *load_balance_waste_time;
double total_load_balance_waste_time = 0.0;
int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	extern FILE* f_exec_time;
	f_exec_time = NULL;
	f_exec_time = fopen("time.log", "w+");
	if (f_exec_time == NULL)
	{
		fprintf(stderr, "[main] error: could not open/create the log file time.log\nAborting.\n");
		exit(EXIT_FAILURE);
	}

	kstring_t pg = {0,0,0};
	t_real = realtime();
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	extension_time = (time_struct*)calloc(100, sizeof(time_struct));
	no_of_extensions = (uint64_t*)calloc(100, sizeof(uint64_t));
	load_balance_waste_time = (double*)calloc(100, sizeof(double));
	uint64_t total_extensions = 0;

	//no_of_extensions = (uint64_t*)calloc(100, sizeof(uint64_t));
	//uint64_t total_extensions = 0;
	if (argc < 2) return usage();

	/*if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwtgen") == 0) ret = bwt_bwtgen_main(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc-1, argv+1);*/
	if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc-1, argv+1);
	/*else if (strcmp(argv[1], "aln") == 0) ret = bwa_aln(argc-1, argv+1);
	else if (strcmp(argv[1], "samse") == 0) ret = bwa_sai2sam_se(argc-1, argv+1);
	else if (strcmp(argv[1], "sampe") == 0) ret = bwa_sai2sam_pe(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtsw2") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "dbwtsw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "bwasw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc-2, argv+2, fresult);
	else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc-1, argv+1);
	else if (strcmp(argv[1], "pemerge") == 0) ret = main_pemerge(argc-1, argv+1);
	else if (strcmp(argv[1], "maxk") == 0) ret = main_maxk(argc-1, argv+1);*/
	else if (strcmp(argv[1], "gase_aln") == 0) ret = gase_aln(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(f_exec_time,"Total_time = %.3f\t", realtime() - t_real);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
		double total_time = realtime() - t_real;
		double total_extension_time = 0.0;
		int n_threads = 0;
		for (i = 0; i < 15; ++i) {
			//fprintf(stderr, "Total time spent in host_mem_alloc by thread %d = %.3f seconds\n", i, extension_time[i].host_mem_alloc);
			//fprintf(stderr, "Percentage of total time spent in extension by thread %d = %.3f\n", i, (extension_time[i]/total_time)*100);
			//fprintf(stderr, "Percentage of total time spent in extension by thread %d = %.3f\n", i, (extension_time[i]/total_time)*100);
			fprintf(stderr, "Total time spent in gpu_aln_kernel by thread %d = %.3f seconds\n", i, extension_time[i].aln_kernel);

			
			fprintf(stderr, "Total time spent in mem_aln1_core by thread %d = %.3f seconds\n", i, extension_time[i].full_mem_aln1_core);
			fprintf(stderr, "\tin chain_preprocess by thread %d = %.3f seconds\n", i, extension_time[i].chain_preprocess);
			fprintf(stderr, "\t\tin mem_chain by thread %d = %.3f seconds\n", i, extension_time[i].time_mem_chain);
			fprintf(stderr, "\t\tin mem_chain_flt by thread %d = %.3f seconds\n", i, extension_time[i].time_mem_chain_flt);
			fprintf(stderr, "\t\tin mem_flt_chained_seeds by thread %d = %.3f seconds\n", i, extension_time[i].time_mem_flt_chained_seeds);
			fprintf(stderr, "\tin mem_chain2aln by thread %d = %.3f seconds\n", i, extension_time[i].full_mem_chain2aln);
			fprintf(stderr, "\tin gpu_aln_kernel by thread %d = %.3f seconds\n", i, extension_time[i].aln_kernel);


			fprintf(stderr, "Total time spent in get_results_actual by thread %d = %.3f seconds\n", i, extension_time[i].get_results_actual);
			fprintf(stderr, "Total time spent in get_results_wasted by thread %d = %.3f seconds\n", i, extension_time[i].get_results_wasted);
			total_extensions += no_of_extensions[i];
			if (no_of_extensions[i] > 0) n_threads++;
			fprintf(stderr, "Total time spent in extension in gpu by thread %d (excluding mem_alloc and mem_free)= %.3f seconds\n", i, extension_time[i].aln_kernel + extension_time[i].get_results_actual + extension_time[i].get_results_wasted);
			total_extension_time += (extension_time[i].aln_kernel + extension_time[i].get_results_actual + extension_time[i].get_results_wasted);
		}
		fprintf(stderr, "Total time spent in gpu_mem_alloc= %.3f seconds\n", extension_time[0].gpu_mem_alloc);
		fprintf(stderr, "Total time spent in gpu_mem_free = %.3f seconds\n", extension_time[0].gpu_mem_free);
		//total_extension_time += (extension_time[0].gpu_mem_alloc + extension_time[0].gpu_mem_free);

		fprintf(f_exec_time,"Average extension time = %.3f\t", (total_extension_time/n_threads) + (extension_time[0].gpu_mem_alloc + extension_time[0].gpu_mem_free));
		fprintf(f_exec_time,"Average extension percentage = %.3f\t", ((total_extension_time/n_threads)/total_time)*100);
		fprintf(stderr,"Total no. of extensions = %llu\n", total_extensions);
		fprintf(f_exec_time, "Percentage time wasted due to imperfect load balancing = %.3f\n", (total_load_balance_waste_time/total_time)*100);


	}
	free(bwa_pg);
	free(extension_time);
	free(no_of_extensions);
	free(load_balance_waste_time);
	fclose(f_exec_time);
	return ret;
}
