#!/bin/bash

## header  --------------------------------------------------------------------

### the one that intersects the small variants 
### and keeps only the high-impact ones

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=10

v_dir="${base_dir}/var-calls/norm-vardict"
g_dir="${base_dir}/var-calls/norm-gatk"
s_dir="${base_dir}/var-calls/norm-strelka"

out_dir="${base_dir}/var-calls/inter-filt"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

### DISCLAIMER: do not confuse STATUS=Deletion with TYPE=Deletion

### the variants we keep
variants_kept="StrongLOH StrongSomatic Deletion LikelySomatic LikelyLOH"
### the TAG "STATUS=" is added in the command line below
### thus, we exclude: 
### Germline
### AFDiff

### the status is defined as in the following lines from
### https://github.com/AstraZeneca-NGS/VarDictJava (december 2022);
### when both samples have coverage for the variant:
###   Germline: detected in germline sample (pass all quality parameters)
###   StrongSomatic: detected in tumor sample only
###   LikelySomatic: the variant has at least one read support OR allele frequency < 5% (defined by -V option with default 0.05)
###   StrongLOH: detected in germline sample only, opposite of StrongSomaitc
###   LikelyLOH: detected in germline but either lost in tumor OR 20-80% in germline, but increased to 1-opt_V (95%).
###   AFDiff: detected in tumor (pass quality parameters) and present in germline but didn't pass quality parameters.
### when only one sample has coverage for the variant:
###   SampleSpecific: detected in tumor sample, but no coverage in germline sample (it's more technical than biological, as it's unlikely a tumor sample can gain a piece of sequence in reference that germline sample lacks).
###   Deletion: detected in germline sample, but no coverage in tumor sample

## clmnt  ---------------------------------------------------------------------

echo "Running the one that intersects the small variants \
and keeps only the high-impact ones..."

cd "${out_dir}"

pll_check=$((pll_runs + 1))
for ind_v in $(find "${v_dir}" -name "*vcf.gz"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  ind_g=$(echo "${ind_v}" | sed 's|vardict|gatk|')
  ind_s=$(echo "${ind_v}" | sed 's|vardict|strelka|')
  sample_id=$(basename "${ind_v}" | sed 's|-norm.vcf.gz||')
  ### strict filter:
  ### extract and write records from vardict
  ### shared by both vardict and gatk using exact allele match
  bcftools isec "${ind_v}" "${ind_g}" \
  -n =2 -w 1 -O z -o "${sample_id}-temp-1.vcf.gz"
  tabix -f -p vcf "${sample_id}-temp-1.vcf.gz"
  ### intersect the intersection with strelka
  bcftools isec "${sample_id}-temp-1.vcf.gz" "${ind_s}" \
  -n =2 -w 1 -O z -o "${sample_id}-isec.vcf.gz"
  rm -f "${sample_id}-temp-1.vcf.gz" "${sample_id}-temp-1.vcf.gz.tbi"
  
  # ### lenient filter:
  # ### shared by both vardict and strelka using exact allele match
  # bcftools isec "${ind_v}" "${ind_s}" \
  # -n =2 -w 1 -O z -o "${sample_id}-temp-2.vcf.gz"
  # tabix -f -p vcf "${sample_id}-temp-2.vcf.gz"
  # ### concatenate the shared variants
  # bcftools concat -a -D \
  # "${sample_id}-temp-1.vcf.gz" \
  # "${sample_id}-temp-2.vcf.gz" \
  # -O z -o "${sample_id}-isec.vcf.gz"
  # rm -f "${sample_id}-temp-1.vcf.gz" "${sample_id}-temp-1.vcf.gz.tbi"
  # rm -f "${sample_id}-temp-2.vcf.gz" "${sample_id}-temp-2.vcf.gz.tbi"
  
  ### make the header
  zgrep "^#" "${sample_id}-isec.vcf.gz" > "${sample_id}-isec-flt.vcf"
  ### keep only high-impact variants and keep indels
  for ind_e in ${variants_kept}; do
    zgrep -v "^#" "${sample_id}-isec.vcf.gz" | \
    grep "STATUS=${ind_e}" \
    >> "${sample_id}-isec-flt.vcf"
  done
  ### keep only high-impact variants and filter indels out
  ### for ind_e in ${variants_kept}; do
  ###   zgrep -v "^#" "${sample_id}-isec.vcf.gz" | \
  ###   grep "STATUS=${ind_e}" | \
  ###   grep -v "TYPE=Insertion" | \
  ###   grep -v "TYPE=Deletion" \
  ###   >> "${sample_id}-isec-flt.vcf"
  ### done
  rm -f "${sample_id}-isec.vcf.gz"
  
  ### sorting and indexing
  bgzip "${sample_id}-isec-flt.vcf"
  bcftools sort "${sample_id}-isec-flt.vcf.gz" \
  -O z -o "${sample_id}-isec-flt-srt.vcf.gz"
  rm -f "${sample_id}-isec-flt.vcf.gz"
  tabix -f -p vcf "${sample_id}-isec-flt-srt.vcf.gz"
  ) &
done

wait
