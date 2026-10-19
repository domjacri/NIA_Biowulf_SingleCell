# You downloaded this from: /data/CARD_datasets/acridj/BDMS_download.tar
# Set your working directory to /WHATEVER_YOU_WANT/BDMS_sc_2026_student/
# This file is in ./scripts/0_cellranger/
# The supplied output is in results/day1/cellranger_backup; this run uses cellranger_run.

module load cellranger
course_dir="$(pwd)"
cd ./results/day1

cellranger count \
  --id="cellranger_run" \
  --fastqs="${course_dir}/data/fastq_files/pbmc_1k_v3_fastqs" \
  --sample="pbmc_1k_v3" \
  --transcriptome="${CELLRANGER_REF}/refdata-gex-GRCh38-2020-A" \
  --localcores=4 \
  --localmem=16 \
  --create-bam=false \
  --disable-cell-annotation
