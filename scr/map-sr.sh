#!/bin/bash

## header  --------------------------------------------------------------------

### the one that maps short-reads and fixes the bam file (fixmate and markdup)
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=22
pll_runs=2
ref_name="${1}"
ctrl_samp="${2}"
popu_samp="${3}"

### output folder
out_dir="${base_dir}/map-sr"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that maps short-reads and \
fixes the bam file (fixmate and markdup)..."

cd "${base_dir}"
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")

all_fqs=$(echo "${ctrl_samp} ${popu_samp}" | sort | uniq)
pll_check=$((pll_runs + 1))
for ind_e in ${all_fqs}; do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  ### define read group for gatk
  str_rg="@RG\tID:${ind_e}\tSM:${ind_e}\tPL:ILLUMINA\tPU:NA\tLB:NA"
  
  ### mapping
  bwa mem -M -t "${n_threads}" "${ref_path}" \
  -R "${str_rg}" \
  "exp/${ind_e}-R1.fq.gz" \
  "exp/${ind_e}-R2.fq.gz" | \
  samtools view -b - > "${out_dir}/${ind_e}-${ref_name}.bam"
  
  ### fix mates
  samtools fixmate -O bam,level=1 -@ "${n_threads}" \
  -m "${out_dir}/${ind_e}-${ref_name}.bam" \
  "${out_dir}/${ind_e}-${ref_name}-fxm.bam"

  ### cleaning
  rm -f "${out_dir}/${ind_e}-${ref_name}.bam"
  
  ### sorting
  samtools sort -O bam,level=1 -@ "${n_threads}" \
  "${out_dir}/${ind_e}-${ref_name}-fxm.bam" \
  > "${out_dir}/${ind_e}-${ref_name}-srt.bam"
  
  ### cleaning
  rm -f "${out_dir}/${ind_e}-${ref_name}-fxm.bam"

  ### marking duplicates
  samtools markdup -O bam,level=1 -@ "${n_threads}" \
  "${out_dir}/${ind_e}-${ref_name}-srt.bam" \
  "${out_dir}/${ind_e}-${ref_name}-srt-mdp.bam"
  
  ### cleaning
  rm -f "${out_dir}/${ind_e}-${ref_name}-srt.bam"

  ### indexing
  samtools index -@ "${n_threads}" \
  "${out_dir}/${ind_e}-${ref_name}-srt-mdp.bam"
  ) &
done

wait
