GASE - Generic Aligner for *Seed-and-Extend*
-----------------------------------------
GASE is a DNA read aligner, developed for measuring the mapping accuracy and execution time of different combinations of seeding and extension techniques. GASE is implemented by extending BWA (version 0.7.13) which is developed by Heng Li. Currently, GASE supports 4 kinds of seeding techniques:

1. All the SMEMs (super-maximal exact matches) in a DNA read (all-SMEM). 
2. Only the non-overlapping SMEMs (nov-SMEM). 
3. Fixed length seeds with no mismatch.
4. Fixed length seeds with at most 1 mismatch. 

Three different kind of extension techniques can be used:

1. Unoptimized global alignment
2. Local alignment (unoptimized and SSE2 optimized)
3. BLAST-like seed extension (unoptimized and banded)

To compile GASE, run `make`. You may see a lot of warnings. These warning will be fixed in the future. To use GASE we first need to build the index by executing the following command:

`gase index <ref.fa>`

This commands builds the FMD index which is the same as the one used in BWA. To align DNA reads use the following command:

`gase gase_aln [options] <ref.fa> <reads.fastq>`

To find all the available options type in `gase gase_aln`. Currently GASE is only tested with single ended reads. 

Feel free to conatact Nauman Ahmed: n.ahmed@tudelft.nl for any issues or bugs.
