This repository contains the scripts used in Dapas M, Lee YL, Wentworth-Sheilds W, Im HK, Ober C, Schoettler N. (2022) Revealing polygenic pleiotropy using genetic risk scores for asthma.


# Summary GWAS data
We used summary statistics from the Trans-National Asthma Genetic Consortium (TAGC) meta-analysis of asthma (23,948 cases, 118,538 controls) as our basis for our polygenic risk model (doi: 10.1038/s41588-017-0014-7). The summary statistics of the meta-analysis are available through the GWAS Catalog entry for the TAGC study on the European Bioinformatics Institute web site (https://www.ebi.ac.uk/gwas/downloads/summary-statistics). Summary statistics from TAGC were published for the entire multi-ancestry (ME) sample (23,948 cases, 118,538 controls) and for the European-ancestry (EUR) subset (19,954 cases, 107,715 controls). 

# Polygenic risk score modeling
PRS models were generated using PRS-CS (https://github.com/getian107/PRScs). We tested two different provided LD reference panels, one generated from 1000 Genomes data and the other from the UK Biobank, both in samples of European ancestry. We generated four total PRS models, utilizing either the 1KG or UKB LD reference panels with each of the EUR and ME summary statistics reported from TAGC. Here is an example run of PRS-CS for the ME cohort using the 1000 Genomes LD reference panel:

```bash
PRScs.py \
    --ref_dir=$dir/ldblk_1kg_eur \
    --bim_prefix=$dir/ukb_bims/ukb_combined \
    --sst_file=$dir/TAGC_ME.sum_stats.txt \
    --n_gwas=142486 \
    --n_iter=10000 \
    --n_burnin=5000 \
    --out_dir=$dir/PRS.PRScs.1KG_LD.TAGC_ME/PRScs \
    --chrom=$chr
```

Posterior SNP effect sizes were combined into aggregate scores for individuals in the UKB using the Plink v2.0 “score” function. Effects for missing genotypes were imputed as the posterior SNP effect size multiplied by the effect allele frequency:

```bash
plink2 --bgen $input_prefix.bgen --sample $input_prefix.sample --score $score_file 2 4 6 ignore-dup-ids cols=+scoresums --out $out_prefix
```

# PRS_asthmaPrediction.R
This script ...

# PRS_phewas.R
This script ...

# PRS_plots.R
This script ...

# PRS_mediation.R
This script ...
