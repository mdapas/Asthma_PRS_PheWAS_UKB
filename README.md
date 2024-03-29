This repository contains the scripts used in Dapas M, Lee YL, Wentworth-Sheilds W, Im HK, Ober C, Schoettler N. (2023) Revealing polygenic pleiotropy using genetic risk scores for asthma.


# Summary GWAS data
The polygenic risk models used in this study were derived from summary statistics from the Trans-National Asthma Genetic Consortium (TAGC) meta-analysis of asthma (23,948 cases, 118,538 controls; doi: 10.1038/s41588-017-0014-7) and from the Global Biobank Meta-analysis Initative (GBMI) meta-analysis of asthma (197,342 cases, 1,903,937 controls; doi: 10.1016/j.xgen.2022.100241). The summary statistics of the TAGC meta-analysis are available through the GWAS Catalog entry for the TAGC study on the European Bioinformatics Institute [website](https://www.ebi.ac.uk/gwas/downloads/summary-statistics). Summary statistics from TAGC were published for the entire multi-ancestry (ME) sample (23,948 cases, 118,538 controls) and for the European-ancestry (EUR) subset (19,954 cases, 107,715 controls). 

# Polygenic risk score (PRS) modeling
For TAGC, we derived polygenic risk models from the GWAS summary statistics. For GBMI, we adopted the model published by Wang and colleagues from the GBMI leave-UKB-out multi-ancestry GWAS deposited on the PGS Catalog ([PGS001787](https://www.pgscatalog.org/score/PGS001787/)). PRS models were generated using [PRS-CS](https://github.com/getian107/PRScs). We generated separate PRS models for each of the EUR and ME summary statistics reported from TAGC. Here is an example run of PRS-CS for the ME cohort using the UKB LD reference panel:

```bash
PRScs.py \
    --ref_dir=$dir/ldblk_ukbb_eur \
    --bim_prefix=$dir/ukb_bims/ukb_combined \
    --sst_file=$dir/TAGC_ME.sum_stats.txt \
    --n_gwas=142486 \
    --n_iter=10000 \
    --n_burnin=5000 \
    --out_dir=$dir/PRS.PRScs.UKB_LD.TAGC_ME/PRScs \
    --chrom=$chr
```

The posterior SNP effect sizes generated from the TAGC ME sample using the UKB LD panel are provided in [TAGC.PRScs_pst_eff.txt.gz](https://github.com/mdapas/Asthma_PRS_PheWAS_UKB/blob/main/TAGC.PRScs_pst_eff.txt.gz) file. 

Then the effect sizes were combined into aggregate scores for each chromosome using the Plink v2.0 “score” function. Effects for missing genotypes were imputed as the posterior SNP effect size multiplied by the effect allele frequency:

```bash
i=$1  # chromosome passed to script
model=PRS.PRScs.UKB_LD.TAGC_ME

input_prefix=$resource_dir/UKB__variant-qc.chr$i
score_file=$dir/$model/PRScs_pst_eff_a1_b0.5_phiauto_chr$i.txt
out_prefix=$dir/$model/UKB.PRScs_score.chr$i

plink2 \
    --bgen $input_prefix.bgen \
    --sample $input_prefix.sample \
    --score $score_file 2 4 6 ignore-dup-ids cols=+scoresums \
    --out $out_dir/$out_prefix
```

And then total scores for each chromosome were combined into an aggregate score for each individual in the UKB.

```bash
# model=PRS.PRScs.UKB_LD.TAGC_ME

# dir={root_dir}/$model

# COMPUTE SNP COUNTS & SCORES
echo "chr ukb_n prs_n scored_n" > $dir/snp_counts.txt
for i in {1..22}; do
    echo $i \
    $(grep "variants detected" $dir/UKB.PRScs_score.chr$i.log | cut -d ' ' -f 2) \
    $(wc -l $dir/PRScs_pst_eff_a1_b0.5_phiauto_chr$i.txt | cut -d ' ' -f 1) \
    $(grep processed $dir/UKB.PRScs_score.chr$i.log | cut -d ' ' -f 2); done >> $dir/snp_counts.txt

for i in {1..22}; do
    s1=$(head -n 1 $dir/UKB.PRScs_score.chr$i.sscore | awk '{print $1}')
    if [[ $s1 = "#FID" ]];
        then
            cut -f2- $dir/UKB.PRScs_score.chr$i.sscore > $dir.tmp.txt && mv $dir.tmp.txt $dir/UKB.PRScs_score.chr$i.sscore
        fi
done

# ADD ALLELE COUNTS
awk '{print $1}' $dir/UKB.PRScs_score.chr1.sscore > $dir/tmp.txt
for i in {1..22}; do awk '{print $2}' $dir/UKB.PRScs_score.chr$i.sscore | paste $dir/tmp.txt - > $dir/tmp.2.txt && mv $dir/tmp.2.txt $dir/tmp.txt; done
head -n 1 $dir/tmp.txt | awk '{print $1, "ALLELE_CT"}' > $dir/tmp.counts.txt
tail -n +2 $dir/tmp.txt | awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23}' >> $dir/tmp.counts.txt

# ADD DOSAGE SUMS
awk '{print $1}' $dir/UKB.PRScs_score.chr1.sscore > $dir/tmp.txt
for i in {1..22}; do awk '{print $3}' $dir/UKB.PRScs_score.chr$i.sscore | paste $dir/tmp.txt - > $dir/tmp.2.txt && mv $dir/tmp.2.txt $dir/tmp.txt; done
head -n 1 $dir/tmp.txt | awk '{print $1, "DOSAGE_SUM"}' > $dir/tmp.dosage_sum.txt
tail -n +2 $dir/tmp.txt | awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23}' >> $dir/tmp.dosage_sum.txt

# ADD SCORES
awk '{print $1}' $dir/UKB.PRScs_score.chr1.sscore > $dir/tmp.txt
for i in {1..22}; do awk '{print $5}' $dir/UKB.PRScs_score.chr$i.sscore | paste $dir/tmp.txt - > $dir/tmp.2.txt && mv $dir/tmp.2.txt $dir/tmp.txt; done
head -n 1 $dir/tmp.txt | awk '{print $1, "SCORE_SUM"}' > $dir/tmp.score.txt
tail -n +2 $dir/tmp.txt | awk '{print $1,$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23}' >> $dir/tmp.score.txt

# MAKE OUTPUT FILE
paste -d ' ' $dir/tmp.counts.txt $dir/tmp.dosage_sum.txt $dir/tmp.score.txt | cut -d ' ' -f 1,2,4,6 > $dir/UKB.PRScs_score.ALL.sscoreSum.txt
rm $dir/tmp.*

```

# PRS_phewas.R
This R script contains all the analyses included in the manuscript after the PRS scores were generated, including the PRS assessments and phenome-wide association (PheWAS) testing. The script is designed to be run in RStudio and is organized linearly to match the manuscript Results. Each section contains corresponding information, references, and notes.
1.  PRS Model Comparison (Figure 1)
2.  PRS Asthma Prediction (Figures 2 & 3)
3.  PheWASs by Ancestry Group (Figure 4)
4.  Comparison of PRS model PheWAS results (Figure 5)
5.  Analysis in Non-asthmatics (Figure 6)
6.  HLA Effects (Figure 7)
7.  Trait Mediation Analysis (Figure 8)

# PRS_phewas_figs.R
This R script contains all the code used to generate the figures included in the manuscript
