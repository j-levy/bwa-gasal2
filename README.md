bwa-gasal2 - an implementation of BWA-MEM using GASAL2 GPU accelerator
-----------------------------------------

**bwa-gasal2** is a modified version of BWA-MEM ([http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)) that uses a GPU-accelerated library for the extension part. The goal is to speed up the extension part and the whole program with minimum alterity to the results.

On our data sets with pair-end reads with 150 bases, we could reach a raw kernel speed-up of 5x, and with overlapped CPU-GPU execution, an effective speed-up of 16x for the extension part. A minimal overhead is although introduced. The whole program is around 1.3x faster with more than 98% of the main output results unchanged.

## Compile and run

This repository uses GASAL2 as a submodule. Clone the repository with its submodule, and switch GASAL2 to branch `bwa`. Adapt the compilation parameters in GASAL2 directory, in file `run_all.sh`. In particular, set the CUDA path, the GPU architecture. The N code should be 4 to be used with bwa-gasal2. Running file `run_all.sh` should provide a correct compilation of the library.

Then in `Makefile`, adapt the paths to the data sets you are using. We recommend adapting and using the compiling rule `srr_threads` as it provides a simple way to run bwa-gasal2 with multiple threads from a bash script.

Before running, you must provide an index for the reference genome you are using. You can use bwa-gasal2 to generate it the same way you would with BWA-MEM. In case anything goes wrong, you can simply use BWA-MEM to generate the index (I didn't test the index generation for bwa-gasal2).

## About GASAL2

Read more at [https://github.com/j-levy/GASAL2](https://github.com/j-levy/GASAL2) (which is the submodule used).

You might want to update the submodule with a more recent version developed by N. Ahmed: [https://github.com/nahmedraja/GASAL2](https://github.com/nahmedraja/GASAL2). 

## Read more

bwa-gasal2 is based on [GASAL2](https://github.com/j-levy/GASAL2) and [GASE-GASAL2](https://github.com/nahmedraja/GASE-GASAL2)

This work has been implemented as a Master thesis project, for which all the results related to it can be found in the thesis: [https://repository.tudelft.nl/islandora/object/uuid%3Abd22471f-058a-4071-95bb-5126e263124b?collection=education](https://repository.tudelft.nl/islandora/object/uuid%3Abd22471f-058a-4071-95bb-5126e263124b?collection=education)
