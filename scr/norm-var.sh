#!/bin/bash

## header  --------------------------------------------------------------------

### the one that normalises and de-duplicates the variants

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=10

### dev off
# full_dir="/home/ltattini/prog/populombe/nc-v4-r1/scr"
# base_dir="/home/ltattini/prog/populombe/nc-v4-r1"
# pll_runs=5

### output folders
log_dir="${base_dir}/scr/logs"
g_out_dir="${base_dir}/var-calls/norm-gatk"
if [[ -d "${g_out_dir}" ]]; then rm -rf "${g_out_dir}"; fi
mkdir -p "${g_out_dir}"
v_out_dir="${base_dir}/var-calls/norm-vardict"
if [[ -d "${v_out_dir}" ]]; then rm -rf "${v_out_dir}"; fi
mkdir -p "${v_out_dir}"
s_out_dir="${base_dir}/var-calls/norm-strelka"
if [[ -d "${s_out_dir}" ]]; then rm -rf "${s_out_dir}"; fi
mkdir -p "${s_out_dir}"

### input folders
g_dir="${base_dir}/var-calls/gatk"
v_dir="${base_dir}/var-calls/vardict"
s_dir="${base_dir}/var-calls/strelka"
### find the reference sequence
r_dir="${base_dir}/ref"
r_file=$(find "${r_dir}" -name "*genome.fa")

## clmnt  ---------------------------------------------------------------------

echo "Running the one that normalises and de-duplicates the variants..."

### parallel check
pll_check=$((pll_runs + 1))

### gatk (compressed) vcf files
all_gatk=$(find "${g_dir}" -name "*vcf.gz")
for ind_g in ${all_gatk}; do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  g_name=$(basename "${ind_g}" | sed 's|...$||')
  my_vcf="${g_out_dir}/${g_name}"
  gunzip -c "${ind_g}" > "${my_vcf}"
  g_out_vcf=$(echo "${my_vcf}" | sed 's|\.vcf|-norm.vcf|')
  vt normalize "${my_vcf}" -r "${r_file}" | \
  vt uniq - -o "${g_out_vcf}"
  rm -f "${my_vcf}"
  bgzip "${g_out_vcf}"
  tabix -f -p vcf "${g_out_vcf}.gz"
  ) &
done

### vardict (compressed) vcf files
all_vardict=$(find "${v_dir}" -name "*vcf.gz")
for ind_v in ${all_vardict}; do
  ### parallel samples
  ((cnt_q++))
  if (( cnt_q % pll_check == 0 )); then
    wait -n
    cnt_q=$(( pll_check - 1 ))
  fi

  (
  v_name=$(basename "${ind_v}" | sed 's|...$||')
  my_vcf="${v_out_dir}/${v_name}"
  gunzip -c "${ind_v}" > "${my_vcf}"
  v_out_vcf=$(echo "${my_vcf}" | sed 's|\.vcf|-norm.vcf|')
  vt normalize "${my_vcf}" -r "${r_file}" | \
  vt uniq - -o "${v_out_vcf}"
  rm -f "${my_vcf}"
  bgzip "${v_out_vcf}"
  tabix -f -p vcf "${v_out_vcf}.gz"
  ) &
done

### strelka (compressed) vcf files
all_strelka=$(find "${s_dir}" -name "*vcf.gz")
for ind_s in ${all_strelka}; do
  ### parallel samples
  ((cnt_r++))
  if (( cnt_r % pll_check == 0 )); then
    wait -n
    cnt_r=$(( pll_check - 1 ))
  fi

  (
  s_name=$(basename "${ind_s}" | sed 's|...$||')
  my_vcf="${s_out_dir}/${s_name}"
  gunzip -c "${ind_s}" > "${my_vcf}"
  s_out_vcf=$(echo "${my_vcf}" | sed 's|\.vcf|-norm.vcf|')
  vt normalize "${my_vcf}" -r "${r_file}" | \
  vt uniq - -o "${s_out_vcf}"
  rm -f "${my_vcf}"
  bgzip "${s_out_vcf}"
  tabix -f -p vcf "${s_out_vcf}.gz"
  ) &
done

wait
