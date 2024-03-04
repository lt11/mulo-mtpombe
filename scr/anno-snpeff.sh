#!/bin/bash

## header  --------------------------------------------------------------------

### the one that annotates the variants with snpeff

## settings  ------------------------------------------------------------------

snpeff_path="${1}"
full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=12

### the reference genome MUST match the reference used to map short-reads
ref_gen_ver="Schizosaccharomyces_pombe"

in_dir="${base_dir}/var-calls/mrg"

out_dir="${base_dir}/var-calls/anno-snpeff"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that annotates the variants with snpeff..."

### parallel check
pll_check=$((pll_runs + 1))

vcf_files=$(find "${in_dir}" -name "ctrl-smpl*.vcf.gz")

for ind_f in ${vcf_files}; do
  ### get the id of the control (the last column) and the evolved samples
  ctrl_id=$(zgrep "^#C" "${ind_f}" | rev | cut -f 1 | rev)
  popu_samps=$(zgrep "^#C" "${ind_f}" | cut -f 10- | rev | cut -f 2- | rev)
  for ind_s in ${popu_samps}; do
    ### parallel samples
    ((cnt_p++))
    if (( cnt_p % pll_check == 0 )); then
      wait -n
      cnt_p=$(( pll_check - 1 ))
    fi
    
    (
    ### set output file and error file
    out_name=$(echo "${ind_s}-${ctrl_id}-anno.vcf")
    err_name=$(echo "${ind_s}-${ctrl_id}-anno.log")
    ### make the file for snpeff (-cancerSamples) with control, \t, tumour
    snpeff_file="${out_dir}/${ctrl_id}-tab-${ind_s}.txt"
    echo -e "${ctrl_id}\t${ind_s}" > "${snpeff_file}"
    ### go to sample output folder since snpeff is dumb
    ### and does not have an output folder option
    ### (the report files are created in the current folder)
    ### jeez
    samp_out_dir="${out_dir}/${ind_s}-vs-${ctrl_id}"
    if [[ -d "${samp_out_dir}" ]]; then rm -rf "${samp_out_dir}"; fi
    mkdir -p "${samp_out_dir}"
    cd "${samp_out_dir}"
    ##### run the annotator
    java -Xmx8g -jar "${snpeff_path}" -v \
    -cancer \
    -cancerSamples "${snpeff_file}" \
    "${ref_gen_ver}" \
    "${ind_f}" \
    > "${out_dir}/${out_name}" \
    2> "${samp_out_dir}/${err_name}"
    ) &
  done
done

wait

### cleaning
find "${out_dir}" -name "*-tab-*.txt" | xargs rm
