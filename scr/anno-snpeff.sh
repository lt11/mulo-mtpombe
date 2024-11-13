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

vcf_files=$(find "${in_dir}" -name "smpls-ctrl*.vcf.gz")

for ind_f in ${vcf_files}; do
  ### get the id ID the control (the last column),
  ### the evolved samples, and the replicate
  ctrl_id=$(zgrep "^#C" "${ind_f}" | rev | cut -f 1 | rev)
  popu_samps=$(zgrep "^#C" "${ind_f}" | cut -f 10- | rev | cut -f 2- | rev)
  rep_id=$(echo "${ind_f##*/}" | cut -d "-" -f 3 | cut -d "." -f 1)
  
  ### make the file for snpeff (-cancerSamples) with control, \t, tumour
  snpeff_file="${out_dir}/${ctrl_id}-tab-samples-${rep_id}.txt"
  if [[-f "${snpeff_file}" ]]; then rm -f "${snpeff_file}"; fi
  for ind_s in ${popu_samps}; do
    echo -e "${ctrl_id}\t${ind_s}" >> "${snpeff_file}"
  done
  ### go to sample output folder since snpeff is dumb
  ### and does not have an output folder option
  ### (the report files are created in the current folder)
  ### jeez
  samp_out_dir="${out_dir}/smpls-${ctrl_id}-${rep_id}"
  if [[ -d "${samp_out_dir}" ]]; then rm -rf "${samp_out_dir}"; fi
  mkdir -p "${samp_out_dir}"
  cd "${samp_out_dir}"
  ### set output file and error file
  out_name="smpls-${ctrl_id}-${rep_id}-anno.vcf"
  err_name="smpls-${ctrl_id}-${rep_id}-anno.log"

  ##### run the annotator
  java -Xmx8g -jar "${snpeff_path}" -v \
  -cancer \
  -cancerSamples "${snpeff_file}" \
  "${ref_gen_ver}" \
  "${ind_f}" \
  > "${out_dir}/${out_name}" \
  2> "${samp_out_dir}/${err_name}"

  ### change names of summary files
  mv snpEff_genes.txt "smpls-${ctrl_id}-${rep_id}-genes.txt"
  mv snpEff_summary.html "smpls-${ctrl_id}-${rep_id}-summary.html"
done

wait

### cleaning
find "${out_dir}" -name "*-tab-*.txt" | xargs rm
