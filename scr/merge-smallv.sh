#!/bin/bash

## header  --------------------------------------------------------------------

### the one that merges variant called in the evolved samples

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=4

ctrl_samp="${1}"
popu_samp="${2}"
repl_ids="${3}"

### dev off
# base_dir="/home/ltattini/prog/populombe/aka-g4-r4"
# ctrl_samp="b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86"
# popu_samp="d69 d70 d71 b79 b83 b87 d72 d73 d74 b80 b84 b88 d75 d76 d77 b81 b85 b89"
# repl_ids="V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10"

input_dir="${base_dir}/var-calls/inter-filt"
read -a arr_ctrl <<< "${ctrl_samp}"
read -a arr_popu <<< "${popu_samp}"
read -a arr_repl <<< "${repl_ids}"
unq_repl=$(echo "${repl_ids}" | tr " " "\n" | sort | uniq)

out_dir="${base_dir}/var-calls/mrg"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that merges variant called in the evolved samples..."

cd "${input_dir}"

for ind_r in ${unq_repl}; do
  echo "${ind_r}" ### tipo M3
  for ind_e in "${!arr_repl[@]}"; do
    if [[ "${arr_repl[$ind_e]}" = "${ind_r}" ]]; then
      ### remove control sample
      popu_prefix=${arr_popu["${ind_e}"]}
      popu_name=$(find . -name "${popu_prefix}-vs-*vcf.gz")
      out_popu="${popu_prefix}-${ind_r}-ptmp.vcf.gz"
      bcftools view -O z -s "${popu_prefix}" "${popu_name}" \
      > "${out_dir}/${out_popu}"
      tabix -f -p vcf "${out_dir}/${out_popu}"
      
      ### pick only the control (for the annotation in cancer mode)
      ctrl_prefix=${arr_ctrl["${ind_e}"]}
      out_ctrl="${ctrl_prefix}-${popu_prefix}-${ind_r}-ctmp.vcf.gz"
      bcftools view -O z -s "${ctrl_prefix}" "${popu_name}" \
      > "${out_dir}/${out_ctrl}"
      tabix -f -p vcf "${out_dir}/${out_ctrl}"
    fi
  done
  ### merge to one multisample vcf with bcftools
  if [[ -f "${out_dir}/multis-${ind_r}.vcf.gz" ]]; then
    rm -f "${out_dir}/multis-${ind_r}.vcf.gz"
  fi
  all_popu=( $(find "${out_dir}" -name "*${ind_r}-ptmp.vcf.gz") )
  # echo ${all_popu[@]}
  bcftools merge --output-type z --threads "${n_threads}" ${all_popu[@]} \
  > "${out_dir}/multis-${ind_r}.vcf.gz"
  tabix -f -p vcf "${out_dir}/multis-${ind_r}.vcf.gz"
  
  ### concatenate, deduplicate, sort the control variants
  all_ctrl=( $(find "${out_dir}" -name "*${ind_r}-ctmp.vcf.gz") )
  bcftools concat --remove-duplicates -a -O u ${all_ctrl[@]} |
  bcftools sort -O z -o "${out_dir}/ctrl-${ind_r}-ctmp.vcf.gz"
  tabix -f -p vcf "${out_dir}/ctrl-${ind_r}-ctmp.vcf.gz"
  ### merge "${out_dir}/multis-${ind_r}.vcf.gz" and control data 
  ### (for the annotation in cancer mode)
  bcftools merge --output-type z --threads "${n_threads}" \
  "${out_dir}/multis-${ind_r}.vcf.gz" \
  "${out_dir}/ctrl-${ind_r}-ctmp.vcf.gz" \
  > "${out_dir}/ctrl-smpl-${ind_r}.vcf.gz"
  tabix -f -p vcf "${out_dir}/ctrl-smpl-${ind_r}.vcf.gz"
done

### cleaning
cd "${out_dir}"
find . -name "*tmp.vcf.gz" | xargs rm
find . -name "*tmp.vcf.gz.tbi" | xargs rm
