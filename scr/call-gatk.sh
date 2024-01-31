#!/bin/bash

## header  --------------------------------------------------------------------

### the one that calls small variants with gatk

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=5
ref_name="${1}"
read -a popu_samp <<< "${2}"
read -a ctrl_samp <<< "${3}"
chr_callable="I II III"
### dev
# full_dir="/home/ltattini/ems/plato/prog/populombe/nc-v14-r1/scr"
# base_dir="/home/ltattini/ems/plato/prog/populombe/nc-v14-r1"
# pll_runs=6
# ref_name="pombe"

### output folder
log_dir="${base_dir}/scr/logs"
out_dir="${base_dir}/var-calls/gatk"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"
### gstk's intervals.list file is not a bed file,
### it can have only column 1 of a bed file
intervals_file="${out_dir}/intervals.list"
### the "dir_db" folder must not exist otherwise we get an error, 
### gatk will create it if we run it with the "panel of normals"
# dir_db="${base_dir}/var-calls/db-gatk"
# if [[ -d "${dir_db}" ]]; then rm -rf "${dir_db}"; fi

## clmnt  ---------------------------------------------------------------------

echo "Running the one that calls small variants with gatk..."

cd "${base_dir}/map-sr"
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")

### get the chromosomes to call onto from the reference genome
### (filter out chromosomes of the reference we don't want to call onto
### and check for chromosome encoding consistency, in fact an error
### is produced if the bed has contigs that are not reported in 
### the reference fasta file)
if [[ -f "${intervals_file}" ]]; then rm -f "${intervals_file}"; fi
ref_chr=$(grep "^>" "${ref_path}" | cut -d " " -f 1 | sed 's|^.||')
for ind_c in ${chr_callable}; do
  echo "${ref_chr}" | grep -w "${ind_c}" >> "${intervals_file}"
done

### run the caller
seq_dim=$(echo "${#popu_samp[@]}")
pll_check=$((pll_runs + 1))
for (( ind_i=0; ind_i<seq_dim; ind_i++ )); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi

  echo "Working on sample ${popu_samp[ind_i]} with control ${ctrl_samp[ind_i]}"
  ### call matched-sample
  gatk Mutect2 \
  -R "${ref_path}" \
  -I "${popu_samp[ind_i]}-${ref_name}-srt-mdp-rg.bam" \
  -I "${ctrl_samp[ind_i]}-${ref_name}-srt-mdp-rg.bam" \
  -normal "${ctrl_samp[ind_i]}" \
  -L "${intervals_file}" \
  --max-mnp-distance 2 \
  -O "${out_dir}/${popu_samp[ind_i]}-vs-${ctrl_samp[ind_i]}.vcf.gz" \
  &> "${log_dir}/${popu_samp[ind_i]}-vs-${ctrl_samp[ind_i]}-mutect.log" &
done

wait

### call t6 unmatched (obsolete)
# gatk Mutect2 \
# -R "${ref_path}" \
# -I t6-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/t6.vcf.gz" &> "${log_dir}/t6-mutect.log" &

### call t6 matched (start from here)
# gatk Mutect2 \
# -R "${ref_path}" \
# -I t6-rerio-srt-mdp-rg.bam \
# -I s1-rerio-srt-mdp-rg.bam \
# -normal s1 \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/t6-vs-s1.vcf.gz" &> "${log_dir}/t6-vs-s1-mutect.log" &

### make the pool from the normal sample and call t6 vs the pool
### step 1: call all the normals
# gatk Mutect2 \
# -R "${ref_path}" \
# -I s1-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s1.vcf.gz" &> "${log_dir}/s1-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s2-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s2.vcf.gz" &> "${log_dir}/s2-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s3-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s3.vcf.gz" &> "${log_dir}/s3-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s4-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s4.vcf.gz" &> "${log_dir}/s4-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s5-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s5.vcf.gz" &> "${log_dir}/s5-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s6-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s6.vcf.gz" &> "${log_dir}/s6-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s7-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s7.vcf.gz" &> "${log_dir}/s7-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s8-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s8.vcf.gz" &> "${log_dir}/s8-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s9-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s9.vcf.gz" &> "${log_dir}/s9-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s10-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s10.vcf.gz" &> "${log_dir}/s10-mutect.log" &

# gatk Mutect2 \
# -R "${ref_path}" \
# -I s11-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# -O "${out_dir}/s11.vcf.gz" &> "${log_dir}/s11-mutect.log" &

### step 2: make the pool from the skin sample
# gatk GenomicsDBImport \
# -R "${ref_path}" \
# -L "1" \
# --genomicsdb-workspace-path "${dir_db}" \
# --reader-threads 12 \
# --tmp-dir "${dir_tmp}" \
# -V "${out_dir}/s1.vcf.gz" \
# -V "${out_dir}/s2.vcf.gz" \
# -V "${out_dir}/s3.vcf.gz" \
# -V "${out_dir}/s4.vcf.gz" \
# -V "${out_dir}/s5.vcf.gz" \
# -V "${out_dir}/s6.vcf.gz" \
# -V "${out_dir}/s7.vcf.gz" \
# -V "${out_dir}/s8.vcf.gz" \
# -V "${out_dir}/s9.vcf.gz" \
# -V "${out_dir}/s10.vcf.gz" \
# -V "${out_dir}/s11.vcf.gz"

# gatk CreateSomaticPanelOfNormals \
# -R "${ref_path}" \
# -V "gendb://${dir_db}" \
# -O "${out_dir}/pooled-normals.vcf.gz"

### step 3: call t6 against the panel of normals (obsolete)
# gatk Mutect2 \
# -R "${ref_path}" \
# -I t6-rerio-srt-mdp-rg.bam \
# -L "1" \
# --max-mnp-distance 0 \
# --panel-of-normals "${out_dir}/pooled-normals.vcf.gz" \
# -O "${out_dir}/t6-vs-normals.vcf.gz" &> "${log_dir}/t6-vs-normals-mutect.log" &

### call in mitochondrial mode?
### well, with a depth of 40 it is useless

####################
# work in progress #
####################

### run FilterMutectCalls for the non-matched variants (obsolete)
# gatk FilterMutectCalls \
# -R "${ref_path}" \
# -V "${out_dir}/t6-vs-normals.vcf.gz" \
# -O "${out_dir}/t6-vs-normals-filt.vcf.gz"

### run FilterMutectCalls for the matched variants (obsolete)
# gatk FilterMutectCalls \
# -R "${ref_path}" \
# -V "${out_dir}/t6-vs-s1.vcf.gz" \
# -O "${out_dir}/t6-vs-s1-filt.vcf.gz"

# (
### call t6 matched and against the panel of normals, and make the f1r2 file
# gatk Mutect2 \
# -R "${ref_path}" \
# -I t6-rerio-srt-mdp-rg.bam \
# -I s1-rerio-srt-mdp-rg.bam \
# -normal s1 \
# --panel-of-normals "${out_dir}/pooled-normals.vcf.gz" \
# -L "1" \
# --max-mnp-distance 0 \
# --f1r2-tar-gz "${out_dir}/t6-vs-s1-vs-normals-f1r2.tar.gz" \
# -O "${out_dir}/t6-vs-s1-vs-normals.vcf.gz" \
# &> "${log_dir}/t6-vs-s1-vs-normals-mutect.log" &
### run gatk LearnReadOrientationModel (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk GetPileupSummaries (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk CalculateContamination (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk FilterMutectCalls (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
# )

# (
### call t6 unmatched but against the panel of normals, and make the f1r2 file
# gatk Mutect2 \
# -R "${ref_path}" \
# -I t6-rerio-srt-mdp-rg.bam \
# --panel-of-normals "${out_dir}/pooled-normals.vcf.gz" \
# -L "1" \
# --max-mnp-distance 0 \
# --f1r2-tar-gz "${out_dir}/t6-vs-normals-f1r2.tar.gz" \
# -O "${out_dir}/t6-vs-normals.vcf.gz" \
# &> "${log_dir}/t6-vs-normals-mutect.log" &
### run gatk LearnReadOrientationModel (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk GetPileupSummaries (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk CalculateContamination (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)

### run gatk FilterMutectCalls (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
# )
