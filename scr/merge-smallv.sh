#!/bin/bash

## header  --------------------------------------------------------------------

### the one that merges variant called in the evolved samples

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=4

popu_samp="${1}"
repl_ids="${2}"

### dev off
# base_dir="/Users/Lorenzo/Desktop/dev-mulo-mtpombe"
# popu_samp="a2 a3 a4 a15 a16 a21 a22 a23 a34 a35 a40 a41 a42 a53 a54 \
# a59 a60 a61 a72 a73"
# repl_ids="M0 M1 M3 D0 D2 M0 M1 M3 D0 D2 M0 M1 M3 D0 D2 \
# M0 M1 M3 D0 D2"

input_dir="${base_dir}/var-calls/inter-filt"
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
  for i in "${!arr_repl[@]}"; do
    if [[ "${arr_repl[$i]}" = "${ind_r}" ]]; then
      # echo "${i}"
      # echo ${arr_popu["${i}"]}
      file_prefix=${arr_popu["${i}"]}
      file_name=$(find . -name "${file_prefix}-vs-*vcf.gz")
      out_name="${file_prefix}-${ind_r}.vcf.gz"
      ### remove control sample
      bcftools view -O z -s "${file_prefix}" "${file_name}" \
      > "${out_dir}/${out_name}"
      tabix -f -p vcf "${out_dir}/${out_name}"
    fi
  done
  ### merge to one multisample vcf with bcftools
  if [[ -f "${out_dir}/multis-${ind_r}.vcf.gz" ]]; then rm -f "${out_dir}/multis-${ind_r}.vcf.gz"; fi
  all_files=( $(find "${out_dir}" -name "*${ind_r}.vcf.gz") )
  # echo ${all_files[@]}
  bcftools merge --output-type z --threads "${n_threads}" ${all_files[@]} \
  > "${out_dir}/multis-${ind_r}.vcf.gz"
  tabix -f -p vcf "${out_dir}/multis-${ind_r}.vcf.gz"
done

### cleaning
cd "${out_dir}"
find . -name "*vcf.gz" | grep -v "multis" | xargs rm
find . -name "*vcf.gz.tbi" | grep -v "multis" | xargs rm
