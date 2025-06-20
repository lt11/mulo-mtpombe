#!/bin/bash

## header ---------------------------------------------------------------------

### this script is the mulo-mtpombe runner

## user's settings ------------------------------------------------------------

### internal reference name
ref_name="pombe"

### set path to snpeff.jar
snpeff_path="/home/tools/lib/snpeff/snpeff-5.1.4/snpEff.jar"

### sample IDs in the correct order, i.e. from the earlier to latest
### time-points, so that the time series plots are correctly sorted

# ### run: Akanksha's g1
# ctrl_samp="b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1"
# popu_samp="b2 b3 b4 b15 b16 c1 \
# d6 d7 d9 d21 d22 d8 \
# b21 b22 b23 b34 b35 c3 \
# d27 d28 d30 d42 d43 d29 \
# c5 b41 b42 b53 c10 c6 \
# d48 d49 d51 d63 d64 d50 \
# b59 b60 b61 b72 b73 c13"
# repl_ids="M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2 \
# M0 M1 M3 D0 D2 M2"
### run: Akanksha's g2
# ctrl_samp="b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1"
# popu_samp="b5 b6 b7 b8 b17 b18 \
# d10 d11 d12 d13 d23 d24 \
# b24 b25 b26 b27 b36 b37 \
# d31 d32 d33 d34 d44 d45 \
# b43 b44 c7 b46 c11 b56 \
# d52 d53 d54 d55 d65 d66 \
# b62 b63 b64 b65 b74 b75"
# repl_ids="M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6 \
# M4 M5 M6 M7 D4 D6"
## run: Akanksha's g3
# ctrl_samp="b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1 \
# b1 b1 b1 b1 b1 b1"
# popu_samp="b9 b10 b11 b19 b20 c2 \
# d14 d15 d16 d25 d26 d17 \
# b28 b29 b30 b38 b39 c4 \
# d35 d36 d37 d46 d47 d38 \
# b47 b48 c8 b57 c12 c9 \
# d56 d57 d58 d67 d68 d59 \
# b66 b67 b68 b76 b77 c14"
# repl_ids="M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11 \
# M8 M9 M10 D8 D9 M11"
### run: Akanksha's g4
# ctrl_samp="b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86 b78 b82 b86"
# popu_samp="d69 d70 d71 b79 b83 b87 d72 d73 d74 b80 b84 b88 d75 d76 d77 b81 b85 b89"
# repl_ids="V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10 V8 V9 V10"
### run: Daniel's g1
ctrl_samp="a1 a1 a1 a1 a1 a15 a15 a15 a15 a15"
popu_samp="a2 a3 a4 a5 a6 a16 a17 a18 a19 a20"
repl_ids="V4 V4 V4 V4 V4 V6 V6 V6 V6 V6"
### run: Daniel's g2
# ctrl_samp="a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29 a29"
# popu_samp="b100 b101 a30 b102 b103 a31 d1 d2 d3 a32 a33 c15 c16 c17 c18 c19 a34"
# repl_ids="V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14 V14"
### run: Daniel's g3
# ctrl_samp="b92 b92 b92 b92 b96 b96 b96 b96"
# popu_samp="b93 b94 b95 d4 b97 b98 b99 d5"
# repl_ids="V12 V12 V12 V12 V13 V13 V13 V13"
### run: Akanksha's g5
# ctrl_samp="b1 b1 b1 b1 b1 b1 b1"
# popu_samp="b12 d18 b31 d39 b50 d60 b69"
# repl_ids="M12 M12 M12 M12 M12 M12 M12"
### run: Akanksha's g6
# ctrl_samp="b1 b1 b1 b1 b1 b1 b1"
# popu_samp="b14 d20 b33 d41 b52 d62 b71"
# repl_ids="M15 M15 M15 M15 M15 M15 M15"
### run: Akanksha's g7
# ctrl_samp="b1 b1 b1 b1 b1 b1 b1"
# popu_samp="b13 d19 b32 d40 b51 d61 b70"
# repl_ids="M14 M14 M14 M14 M14 M14 M14"

## system's settings ----------------------------------------------------------

### check logs folder
if [[ ! -d "logs" ]]; then mkdir "logs"; fi

## clmnt ----------------------------------------------------------------------

### quality check
bash fq-check.sh \
> "logs/fq-check.out" 2> "logs/fq-check.err" &

### reference indexing
bash index-ref.sh "${ref_name}" \
> "logs/index-ref.out" 2> "logs/index-ref.err"

### mapping
bash map-sr.sh "${ref_name}" "${ctrl_samp}" "${popu_samp}" \
> "logs/map-sr.out" 2> "logs/map-sr.err"

### coverage statistics
bash depth-stats.sh \
> "logs/depth-stats.out" 2> "logs/depth-stats.err"

(
### mappability calculation
bash gem.sh "${ref_name}" "${read_len}" \
> "logs/gem.out" 2> "logs/gem.err"

### GC-content calculation
bash gc-fastas.sh "${ref_name}" \
> "logs/gc-fastas.out" 2> "logs/gc-fastas.err"

### call copy-number variants (single-sample mode)
bash cnv-freec.sh \
> "logs/cnv-freec.out" 2> "logs/cnv-freec.err" 
) &

### gatk controls which chromosomes are callable (with "chr_callable")
bash call-gatk.sh "${ref_name}" "${popu_samp}" "${ctrl_samp}" \
> "logs/call-gatk.out" 2> "logs/call-gatk.err"

### vardict sets the filters
bash call-vardict.sh "${ref_name}" "${popu_samp}" "${ctrl_samp}" \
> "logs/call-vardict.out" 2> "logs/call-vardict.err"

### calling with strelka (who controls nothing)
bash call-strelka.sh "${ref_name}" "${popu_samp}" "${ctrl_samp}" \
> "logs/call-strelka.out" 2> "logs/call-strelka.err"

### normalise and de-duplicate the variants
bash norm-var.sh \
> "logs/norm-var.out" 2> "logs/norm-var.err"

### intersect and filter
bash int-flt-var.sh \
> "logs/int-flt-var.out" 2> "logs/int-flt-var.err"

### merge the variants of evolved samples
bash merge-smallv.sh "${ctrl_samp}" "${popu_samp}" "${repl_ids}" \
> "logs/merge-smallv.out" 2> "logs/merge-smallv.err"

### parse and plot time-points (sorted as the IDs)
Rscript parse-plot-af.R \
> "logs/parse-plot-af.out" 2> "logs/parse-plot-af.err"

### annotation with snpeff
bash anno-snpeff.sh "${snpeff_path}" \
> "logs/anno-snpeff.out" 2> "logs/anno-snpeff.err"

echo "Prematura la supercazola o scherziamo?"
echo "[Conte Lello Mascetti]"
