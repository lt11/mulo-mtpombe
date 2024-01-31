#!/bin/bash

## header  --------------------------------------------------------------------

### the one that calls small variants with strelka

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
### one chromosome per thread
### (with M threads if M > N, it takes N threads if we have N chromosomes)
n_threads=8
pll_runs=5
ref_name="${1}"
read -a popu_samp <<< "${2}"
read -a ctrl_samp <<< "${3}"

### the reference
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")

### dev off
# base_dir="/home/ltattini/prog/populombe/nc-v4-r1"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that calls small variants with strelka..."

### output folder
log_dir="${base_dir}/scr/logs"
out_dir="${base_dir}/var-calls/strelka"

### working folder
cd "${base_dir}/map-sr"

### sample loop
seq_dim=$(echo "${#popu_samp[@]}")
pll_check=$((pll_runs + 1))
for (( ind_i=0; ind_i<seq_dim; ind_i++ )); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  ### input bam files
  popu_bam="${popu_samp[ind_i]}-${ref_name}-srt-mdp-rg.bam"
  ctrl_bam="${ctrl_samp[ind_i]}-${ref_name}-srt-mdp-rg.bam"
  ### sample working folder
  samp_ctrl_pair="${popu_samp[ind_i]}-vs-${ctrl_samp[ind_i]}"
  run_dir="${out_dir}/${samp_ctrl_pair}"
  if [[ -d "${run_dir}" ]]; then rm -rf "${run_dir}"; fi
  mkdir -p "${run_dir}"
  ### configure strelka
  configureStrelkaSomaticWorkflow.py \
  --normalBam "${ctrl_bam}" \
  --tumorBam "${popu_bam}" \
  --referenceFasta "${ref_path}" \
  --runDir "${run_dir}"
  ### execution on a single local machine with "-j N" parallel jobs
  "${run_dir}/runWorkflow.py" -m local -j "${n_threads}"

  ### modify output file names
  raw_indels="${run_dir}/results/variants/somatic.indels.vcf.gz"
  bon_indels="${run_dir}/results/variants/${samp_ctrl_pair}-indels.vcf.gz"
  mv "${raw_indels}" "${bon_indels}"
  mv "${raw_indels}.tbi" "${bon_indels}.tbi"
  raw_snp="${run_dir}/results/variants/somatic.snvs.vcf.gz"
  bon_snps="${run_dir}/results/variants/${samp_ctrl_pair}-snps.vcf.gz"
  mv "${raw_snp}" "${bon_snps}"
  mv "${raw_snp}.tbi" "${bon_snps}.tbi"

  ### merge indels and SNPs, fix sample names, and index the vcf file
  sample_out_dir="${run_dir}/results/variants"
  file_merged="${sample_out_dir}/${samp_ctrl_pair}.vcf"
  bcftools concat -a -O v "${bon_snps}" "${bon_indels}" \
  -o "${file_merged}.tmp"
  grep "^##" "${file_merged}.tmp" > "${file_merged}"
  grep "^#C" "${file_merged}.tmp" | \
  sed "s|NORMAL|${ctrl_samp[ind_i]}|" | \
  sed "s|TUMOR|${popu_samp[ind_i]}|" >> "${file_merged}"
  grep -v "^#" "${file_merged}.tmp" >> "${file_merged}"
  bgzip "${file_merged}"
  tabix -f -p vcf "${file_merged}.gz"
  rm -f "${file_merged}.tmp"
  rm -f "${bon_snps}" "${bon_indels}"
  rm -f "${bon_snps}.tbi" "${bon_indels}.tbi"
  ) &
done

wait
