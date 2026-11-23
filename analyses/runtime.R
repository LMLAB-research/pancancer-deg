#### Runtime settings shared by all modules ####
#
# One BLAS/OpenMP thread per R process. The steps parallelise themselves with
# BiocParallel workers. A multi-threaded BLAS (e.g. OpenBLAS, which by default
# starts one thread per core) is much slower on the many small per-gene matrix
# products of DESeq2, because starting and synchronising the threads costs more
# than the arithmetic, and forked workers multiply the threads further. On a
# 32-core machine with OpenBLAS, step 04 did not finish the TCGA-BRCA paired
# design (114 coefficients) in more than a day; with one thread it takes minutes.
# Set here rather than through environment variables so that it applies however
# the scripts are started (run_pipeline.sh, Rscript or an interactive session).

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
