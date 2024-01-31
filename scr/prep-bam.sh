#!/bin/bash

## header  --------------------------------------------------------------------

### the one that prepares the bam files for gatk

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=4

## clmnt  ---------------------------------------------------------------------

echo "Running the one that prepares the bam files for gatk..."

cd "${base_dir}/map-sr"
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")
ref_name=$(basename "${ref_path}" | cut -d "-" -f 1)

pll_check=$((pll_runs + 1))
for ind_s in $(find . -name "*mdp.bam"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  (
  sample_name=$(echo "${ind_s}" | sed 's|..||' | cut -d "-" -f 1)
  out_file=$(echo "${ind_s}" | sed 's|\.bam|-rg.bam|')

  gatk AddOrReplaceReadGroups \
  -I "${ind_s}" \
  -O "${out_file}" \
  --RGLB "${ref_name}-${sample_name}" \
  --RGPL "ILLUMINA" \
  --RGPU "NA" \
  --RGSM "${sample_name}" \
  --CREATE_INDEX True

  rm "${ind_s}" "${ind_s}.bai"
  ) &
done

wait
