# Optional: the student tarball already includes these FASTQs.
# Run from the BDMS_sc_2026_student folder.
cd ./data/fastq_files

# Download the ~5 GB compressed test FASTQs
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar

# Extract the archive
tar -xf pbmc_1k_v3_fastqs.tar
