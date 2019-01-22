#!/bin/bash
./bwa-gasal2 gase_aln -t 12 -l 150 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/SRR949537_1.fastq /data/work/jlevy/srr/150/SRR949537_2.fastq > /data/work/jlevy/srr/150/res_bwa_gasal2.log

./bwa-gasal2 gase_aln -t 1 -l 150 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/SRR949537_1.fastq /data/work/jlevy/srr/150/SRR949537_2.fastq > /data/work/jlevy/srr/150/res_bwa_gasal2.log

./bwa-gasal2 gase_aln -t 2 -l 150 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/SRR949537_1.fastq /data/work/jlevy/srr/150/SRR949537_2.fastq > /data/work/jlevy/srr/150/res_bwa_gasal2.log

./bwa-gasal2 gase_aln -t 4 -l 150 /data/work/jlevy/hg19.fasta /data/work/jlevy/srr/150/SRR949537_1.fastq /data/work/jlevy/srr/150/SRR949537_2.fastq > /data/work/jlevy/srr/150/res_bwa_gasal2.log
