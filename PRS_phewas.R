####  SET PARAMETERS  ##############################################################################
#
#  NOTE: even if running an individual section, these parameters should be saved first
#

set.seed(0123456789)

bw <- 'White'; nbw <- 'non-British White'; afr <- 'African'
prs <- 'zscore'; covars <- paste0('+sex',paste(paste0('+pc',1:10),collapse=''))
afr_covars <- paste0('+sex',paste(paste0('+pc',1:20),collapse=''))
coa <- 'asthma_children'; aoa <- 'asthma_adults'; asthma <- 'asthma_all'
ukb_var <- 'ukb_variable'

# These paths would be user-specific
# score_dir <- [/path/to/directory/with/prs/scores]
# pheno_dir <- [/path/to/directory/with/phenotype/data]
# output_dir <- [/path/to/output/directory]


####  FUNCTIONS  ###################################################################################

rntransform <- function(x) {
  var <- ztransform(x)
  out <- rank(var) - 0.5
  out[is.na(var)] <- NA
  out <- out/(max(out,na.rm=T)+.5)
  out <- qnorm(out)
  return(out)
}

ztransform <- function(x) {
  x <- x[!is.na(x)]
  return((x-mean(x))/sd(x))
}

sfrac <- function(top, bottom, data=NULL) {
  with(data,lapply(paste0(top,"\n",bottom)))
}

get_r2 <- function(model, null_model, n, n_cases) {
  # Calculate r2 on liability scale (Lee, AJHG, 2012), assuming k=p
  LL1 <-  logLik(model)
  LL0 <-  logLik(null_model)
  r2_CS <- 1 - exp(2*(LL0[1] - LL1[1])/n)
  k <- n_cases / n
  t <- qnorm(1 - k)
  z <- dnorm(t)
  c <- k*(1-k)/(z^2)
  r2 <- c*r2_CS
  return(r2)
}


####  READ INPUT DATA  ############################################################################

# Read asthma phenotype data
# asthma_phenos <- read.csv(paste0(pheno_dir, "/asthma--data.csv"))  # asthma phenotypes
# asthma_covars <- read.csv(paste0(pheno_dir, "/asthma--covariates.csv"))  # regression covariates
# asthma_dat <- merge(asthma_phenos, asthma_covars, by='eid')

# trait_dat <- read.csv(paste0(pheno_dir, "/pwas--traits.csv"))  # ukb ethnicity and trait data

# trait_dat <- merge(asthma_dat[,!names(asthma_dat) %in% names(trait_dat)[2:ncol(trait_dat)]], trait_dat,
#                   by='eid')  # combine traits and covariates, avoid duplicate columns

# trait_dat$ancestry <- ifelse((trait_dat$ethnic_background==1001 & trait_dat$genetic_ethnic_grouping==1), bw,
#                             ifelse(trait_dat$ethnic_background %in% c(1002,1003), nbw,
#                                    ifelse((trait_dat$ethnic_background %in% c(4,4001,4002,4003) & 
#                                              trait_dat$pc1>150),afr,NA)))  # apply ancestry definitions

# write.table(trait_dat, paste0(pheno_dir, '/ukb.trait.data.csv'), quote=F, row.names=F,col.names=T, sep=",")

trait_dat <- read.csv(paste0(pheno_dir, "/ukb.trait.data.csv"))  # ukb ethnicity and trait data


###############################################################################################
### 1.  PRS Model Comparison                                                                ###
###############################################################################################
# 
# This analysis executes compares the asthma prediction results in each UKB ancestry group
#   between the PRS models
#

library(pROC)

## trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.data.csv'), header=T)

# Read PRS scores, compute asthma prediction
n_replicates <- 1000
predict_df <- data.frame(matrix(ncol=14, nrow=0))
colnames(predict_df) <- c('Study','LD', 'UKB_ancestry','Asthma_pheno','AUC','AUC_95CI_L','AUC_95CI_U',
                          'r2','r2_2.5','r2_97.5','OR','OR_95CI_L','OR_95CI_U','P')

i <- 1
for (study in c('TAGC', 'GBMI')) {
  if (study=='TAGC') {
    for (ld in c('1KG','UKB')) {
      for (pop in c('EUR','ME')) {
        print(paste0('PRScs.',ld,'_LD.TAGC_',pop))
        scores <- read.table(paste0(score_dir, '/PRScs.',ld,'_LD.TAGC_',pop,'.UKB_scores.txt'), header=T)
        scores[,prs] <- ztransform(scores$SCORE_SUM)
        scores <- merge(trait_dat, scores[,c("IID", prs)], by.x="eid", by.y="IID") 
        for (anc in c(bw, nbw, afr)) {
          print(paste0('  ',anc))
          cvs <- ifelse(anc==afr, afr_covars, covars)
          model_dat <- scores[which(scores$ancestry==anc),]
          for (case in c(asthma, coa, aoa)) {
            print(paste0('    ',case))
            case_cntl <- table(model_dat[,case])
            n <- sum(case_cntl)
            roc_auc <- roc(model_dat[,case], model_dat[,prs], plot=F, quiet=T)$auc[1]
            
            model <- glm(paste0(case,'~',prs,cvs), data=model_dat, family=binomial(logit))
            null_model <- glm(paste0(case,'~',substr(cvs,2,nchar(cvs))), data=model_dat, family=binomial(logit))
            
            r2 <- get_r2(model, null_model, n, case_cntl[2])
            
            # calculate 95% CI for r2 and AUC
            r2_values <- rep(NA, n_replicates)
            auc_values <- rep(NA, n_replicates)
            for (r in 1:n_replicates) {
              print(r)
              sample_dat <- model_dat[sample(nrow(model_dat), replace=TRUE),]
              case_cntl_sample <- table(sample_dat[,case])
              sample_model <- glm(paste0(case,'~',prs,cvs), data=sample_dat, family=binomial(logit))
              sample_null_model <- glm(paste0(case,'~',substr(cvs,2,nchar(cvs))), data=sample_dat, family=binomial(logit))
              r2_values[r] <- get_r2(sample_model, sample_null_model, n, case_cntl_sample[2])
              auc_values <- roc(sample_dat[,case], sample_dat[,prs], plot=F, quiet=T)$auc[1]
            }
            
            lower_ci_r2 <- quantile(as.numeric(r2_values), 0.025)
            upper_ci_r2 <- quantile(as.numeric(r2_values), 0.975)
            lower_ci_auc <- quantile(as.numeric(auc_values), 0.025)
            upper_ci_auc <- quantile(as.numeric(auc_values), 0.975)
            
            predict_df[i,] <- c(paste0('TAGC_',pop), ld, anc, case, roc_auc, lower_ci_auc, upper_ci_auc, 
                                r2, lower_ci_r2, upper_ci_r2,
                                exp(model$coefficients[2])[[1]], 
                                exp(model$coefficients[2]-1.96*sqrt(diag(vcov(model)))[2])[[1]],
                                exp(model$coefficients[2]+1.96*sqrt(diag(vcov(model)))[2])[[1]],
                                coef(summary(model))[2,4])
            i <- i+1
          }
        }
      }
    }
  }
}


# GBMI Model
scores <- read.table(paste0(score_dir, '/PRScs.1KG_LD.GBMI.UKB_scores.txt'), header=T)
scores[,prs] <- ztransform(scores$SCORE_SUM)

trait_dat <- merge(trait_dat, scores[,c("IID", prs)], by.x="eid", by.y="IID")
for (anc in c(bw, nbw, afr)) {
  print(paste0('  ',anc))
  cvs <- ifelse(anc==afr, afr_covars, covars)
  model_dat <- trait_dat[trait_dat$ancestry==anc,]
  for (case in c(asthma, coa, aoa)) {
    print(paste0('    ',case))
    case_cntl <- table(model_dat[,case])
    n <- sum(case_cntl)
    roc_auc <- roc(model_dat[,case], model_dat[,prs], plot=F, quiet=T)$auc[1]
    
    model <- glm(paste0(case,'~',prs,cvs), data=model_dat, family=binomial(logit))
    null_model <- glm(paste0(case,'~',substr(cvs,2,nchar(cvs))), data=model_dat, family=binomial(logit))
    
    r2 <- get_r2(model, null_model, n, case_cntl[2])
    
    # calculate 95% CI for r2 and AUC
    r2_values <- rep(NA, n_replicates)
    auc_values <- rep(NA, n_replicates)
    for (r in 1:n_replicates) {
      print(r)
      sample_dat <- model_dat[sample(nrow(model_dat), replace=TRUE),]
      case_cntl_sample <- table(sample_dat[,case])
      sample_model <- glm(paste0(case,'~',prs,cvs), data=sample_dat, family=binomial(logit))
      sample_null_model <- glm(paste0(case,'~',substr(cvs,2,nchar(cvs))), data=sample_dat, family=binomial(logit))
      r2_values[r] <- get_r2(sample_model, sample_null_model, n, case_cntl_sample[2])
      auc_values <- roc(sample_dat[,case], sample_dat[,prs], plot=F, quiet=T)$auc[1]
    }
    
    lower_ci_r2 <- quantile(as.numeric(r2_values), 0.025)
    upper_ci_r2 <- quantile(as.numeric(r2_values), 0.975)
    lower_ci_auc <- quantile(as.numeric(auc_values), 0.025)
    upper_ci_auc <- quantile(as.numeric(auc_values), 0.975)
    
    predict_df[i,] <- c('GBMI', '1KG', anc, case, roc_auc, lower_ci_auc, upper_ci_auc, 
                        r2, lower_ci_r2, upper_ci_r2,
                        exp(model$coefficients[2])[[1]], 
                        exp(model$coefficients[2]-1.96*sqrt(diag(vcov(model)))[2])[[1]],
                        exp(model$coefficients[2]+1.96*sqrt(diag(vcov(model)))[2])[[1]],
                        coef(summary(model))[2,4])
    i <- i+1
  }
}

write.table(predict_df, paste0(output_dir, '/asthma_prs_prediction_r2liability.txt'), quote=F, row.names=F,col.names=T, sep="\t")


###############################################################################################
### 2.  PRS Asthma Prediction                                                               ###
###############################################################################################
# 
# This analysis assigns individuals to PRS quantiles, computes relative risks, and generates
#   the plots for Figures 1 & 2
#

library(sjmisc)

# trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.data.csv'), header=T)

#scores <- read.table(paste0(score_dir, '/PRScs.1KG_LD.GBMI.UKB_scores.txt'), header=T)
#scores[,prs] <- ztransform(scores$SCORE_SUM)
#trait_dat <- merge(trait_dat, scores[,c("IID", prs)], by.x="eid", by.y="IID")

trait_dat <- trait_dat[!is.na(trait_dat[,prs]),]
trait_dat <- trait_dat[!is.na(trait_dat$ancestry),]

trait_dat$sex <- to_factor(ifelse(trait_dat$sex==1,'M',ifelse(trait_dat$sex==0,'F',NA)))

# Calculate PRS percentiles for each individual, calculated by ancestry group
p <- 'prs_percentile'
for (anc in c(bw, nbw, afr)) {
  trait_dat[trait_dat$ancestry==anc,p] <- 
    ecdf(trait_dat[trait_dat$ancestry==anc,prs])(trait_dat[trait_dat$ancestry==anc,prs])
}

## Define PRS quantiles
quantiles <- c('[0-1]','(1-5]','(5-10]','(10-20]','(20-40]','(40-60]',
               '(60-80]','(80-90]','(90-95]','(95-99]','(99-100]')
q <-'quantile'

trait_dat[,q] <- ifelse(trait_dat[,p]>.99,'(99-100]',
                    ifelse(trait_dat[,p]<=.99 & trait_dat[,p]>.95,'(95-99]',
                        ifelse(trait_dat[,p]<=.99 & trait_dat[,p]>.95,'(95-99]',
                            ifelse(trait_dat[,p]<=.95 & trait_dat[,p]>.9,'(90-95]',
                                ifelse(trait_dat[,p]<=.9 & trait_dat[,p]>.8,'(80-90]',
                                    ifelse(trait_dat[,p]<=.8 & trait_dat[,p]>.6,'(60-80]',
                                        ifelse(trait_dat[,p]<=.6 & trait_dat[,p]>.4,'(40-60]',
                                            ifelse(trait_dat[,p]<=.4 & trait_dat[,p]>.2,'(20-40]',
                                                ifelse(trait_dat[,p]<=.2 & trait_dat[,p]>.1,'(10-20]',
                                                    ifelse(trait_dat[,p]<=.1 & trait_dat[,p]>.05,'(5-10]',
                                                        ifelse(trait_dat[,p]<=.05 & trait_dat[,p]>.01,'(1-5]',
                                                            ifelse(trait_dat[,p]<=.01,'[0-1]',NA))))))))))))

#trait_dat$tail <- ifelse(trait_dat[,p]>=.95,'95th',ifelse(trait_dat[,p]<=.05,'5th',NA))
trait_dat$tail <- ifelse(trait_dat[,p]>.9,'90th',ifelse(trait_dat[,p]<=.1,'10th',NA))

# Compare deciles by asthma phenotype
for (phenotype in c(coa, aoa, asthma)){
  print(phenotype)
  quantile_dat <- data.frame(matrix(ncol = 4, nrow = 11))
  colnames(quantile_dat) <- c('Quantile','OR','95CI_L','95CI_U')
  quantile_dat$Quantile <- quantiles
  model <- glm(paste0(phenotype,'~',covars,'+tail'), data=trait_dat[trait_dat$ancestry==bw,], family='binomial')
  print(paste('Top 10% vs. bottom 10%:    ', round(exp(summary(model)$coefficients[13,1]),3),
              round(exp(summary(model)$coefficients[13,1]-(1.96*summary(model)$coefficients[13,2])),3),
              round(exp(summary(model)$coefficients[13,1]+(1.96*summary(model)$coefficients[13,2])),3)))
  for(quant in quantiles) {
    if (quant=='(40-60]') {
      quantile_dat[quantile_dat$Quantile==quant,c('OR','95CI_L','95CI_U')] <- 1
    } else {
      tmp_dat <- trait_dat[(trait_dat$ancestry==bw) & (trait_dat[,q] %in% c(quant, '(40-60]')),]
      tmp_dat[,q] <- factor(tmp_dat[,q])
      tmp_dat[,q] <- relevel(tmp_dat[,q], ref = '(40-60]')
      f_lm <- paste0(phenotype,'~',covars,'+',q)
      model <- glm(f_lm, data=tmp_dat, family='binomial')
      quantile_dat[quantile_dat$Quantile==quant,
                   c('OR','95CI_L','95CI_U')] <- c(round(exp(summary(model)$coefficients[13,1]),3),
                                                   round(exp(summary(model)$coefficients[13,1]-
                                                               (1.96*summary(model)$coefficients[13,2])),3),
                                                   round(exp(summary(model)$coefficients[13,1]+
                                                               (1.96*summary(model)$coefficients[13,2])),3))
    }
  }
  if (phenotype==coa) {
    quantile_dat_COA <- quantile_dat
  } else if (phenotype==aoa) {
    quantile_dat_AOA <- quantile_dat
  } else {quantile_dat_asthma <- quantile_dat}
}

# Test sex interaction effects
model <- glm(paste0(coa,'~',prs,covars,'+sex:',prs), data=trait_dat[trait_dat$ancestry==bw,], family='binomial')
summary(model)

model <- glm(paste0(aoa,'~',prs,covars,'+sex:',prs), data=trait_dat[trait_dat$ancestry==bw,], family='binomial')
summary(model)

# write.csv(trait_dat, paste0(pheno_dir, '/ukb.trait.prs.data.csv'), quote=F, row.names=F, col.names=T)
# save(quantile_dat_asthma, quantile_dat_COA, quantile_dat_AOA, file=paste0(output_dir, '/quantile_dat.RData'))


###############################################################################################
### 3.  PheWASs by ancestry group                                                           ###
###############################################################################################
# 
# Perform phenome-wide association test with the asthma PRS
#

# library(tidyverse)
# library(ggrepel)
# library(ggtext)
# library(grid)
# library(stringr)

# trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'), header=T)

pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")

quant_traits <- pwas_metadata[pwas_metadata$bin==0,ukb_var]
bin_traits <- pwas_metadata[pwas_metadata$bin==1,ukb_var]

pwas_dat <- trait_dat[,c('eid',asthma, prs,'ancestry', strsplit(afr_covars,'+',fixed=T)[[1]][-1],
                         bin_traits, quant_traits)]

for (trait in quant_traits) {
  pwas_dat[!is.na(pwas_dat[,trait]),trait] <- ztransform(pwas_dat[,trait])  # Standardize quant traits
}

# Run PheWAS in each ancestry group
pwas_results <- data.frame(matrix(ncol=13, nrow=0))
colnames(pwas_results) <- c(names(pwas_metadata),'ancestry','sex','beta','se','p','n', 'df')
i <- 1 
for (anc in c(bw, nbw, afr)) {
  print(anc)
  if (anc==afr) {
    params <- afr_covars
  } else {
    params <- covars
  }
  for (trait in c(bin_traits, quant_traits)) {
    print(c(i, trait))
    sex_check <- colSums(table(pwas_dat[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
    sexes <- names(sex_check[sex_check>0])
    test_dat <- pwas_dat[(pwas_dat$sex %in% sexes) & (pwas_dat$ancestry==anc),
                         c(trait,prs,strsplit(params,'+',fixed=T)[[1]][-1], 'ancestry')]
    if (length(sexes)<2) {
      model_covars <- substr(params, 5, nchar(params))
    } else {
      model_covars <- params
    }
    pwas_results[i,] <- tryCatch(
      {
        if (trait %in% bin_traits) {
          model <- glm(paste0(trait,'~',prs,model_covars), data=test_dat, family=binomial)
          } else {
            model <- lm(paste0(trait,'~',prs,model_covars), data=test_dat)
            }
        s <- summary(model)
        c(pwas_metadata[pwas_metadata$ukb_variable==trait,], anc, paste(sexes,collapse=' & '),
          model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals), df.residual(model))
        },
      error=function(e){
        return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], anc, rep(NA,6)))
        }
      )
    i <- i+1
  }
}


### Test replication in other ancestries
sig_trait_dat <- merge(pwas_results[(pwas_results$ancestry==bw) & 
                                      (pwas_results$p<0.05/(nrow(pwas_results)/3)), c(ukb_var,'beta','p')],
                       pwas_results[(pwas_results$ancestry==nbw), c(ukb_var,'beta','p')], 
                       by=ukb_var, all.x=T)
sig_trait_dat <- merge(sig_trait_dat,
                       pwas_results[(pwas_results$ancestry==afr), c(ukb_var,'beta','p')], 
                       by=ukb_var, all.x=T)
names(sig_trait_dat) <- c(ukb_var, 'beta.bw','p.bw','beta.nbw','p.nbw','beta.afr','p.afr')

sig_traits <- sig_trait_dat$ukb_variable

rep_traits_nbw <- na.exclude(sig_trait_dat[(sig_trait_dat$p.nbw<0.05) & 
                                            (sig_trait_dat$beta.bw*sig_trait_dat$beta.nbw>0),ukb_var])

rep_traits_afr <- na.exclude(sig_trait_dat[(sig_trait_dat$p.afr<0.05) & 
                                             (sig_trait_dat$beta.bw*sig_trait_dat$beta.afr>0),ukb_var])


### Test for African ancestry-specific effects

# Run PheWAS in each ancestry group
pwas_results_w <- data.frame(matrix(ncol=12, nrow=0))
colnames(pwas_results_w) <- c(names(pwas_metadata),'sex','beta','se','p','n', 'df')
i <- 1 

for (trait in rep_traits_afr[rep_traits_afr %in% rep_traits_nbw]) {
  print(c(i, trait))
  sex_check <- colSums(table(pwas_dat[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
  sexes <- names(sex_check[sex_check>0])
  test_dat <- pwas_dat[(pwas_dat$sex %in% sexes) & (pwas_dat$ancestry %in% c(bw, nbw)),
                       c(trait,prs,strsplit(covars,'+',fixed=T)[[1]][-1], 'ancestry')]
  if (length(sexes)<2) {
    model_covars <- substr(covars, 5, nchar(covars))
  } else {
    model_covars <- covars
  }
  pwas_results_w[i,] <- tryCatch(
    {
      if (trait %in% bin_traits) {
        model <- glm(paste0(trait,'~',prs,model_covars), data=test_dat, family=binomial)
      } else {
        model <- lm(paste0(trait,'~',prs,model_covars), data=test_dat)
      }
      s <- summary(model)
      c(pwas_metadata[pwas_metadata$ukb_variable==trait,], paste(sexes,collapse=' & '),
        model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals), df.residual(model))
    },
    error=function(e){
      return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], anc, rep(NA,6)))
    }
  )
  i <- i+1
}


# Using Wald/T statistic 
anc_comp <- merge(pwas_results_w, 
                  pwas_results[pwas_results$ancestry==afr,c('ukb_variable','beta','se','p','n','df')], 
                  by='ukb_variable')
#anc_comp <- anc_comp[anc_comp$ukb_variable %in% rep_traits_afr[rep_traits_afr %in% rep_traits_nbw],]
names(anc_comp) <- gsub(x = names(anc_comp), pattern = "\\.x", replacement = "\\.w")  
names(anc_comp) <- gsub(x = names(anc_comp), pattern = "\\.y", replacement = "\\.afr")  

anc_comp$p_diff <- ifelse(anc_comp$ukb_variable %in% bin_traits,
                                 pchisq(((anc_comp$beta.w-anc_comp$beta.afr)/
                                           sqrt(anc_comp$se.w^2 + anc_comp$se.afr^2))^2, df=1, lower.tail=F),
                                 2*pt(abs((anc_comp$beta.w-anc_comp$beta.afr)/
                                            sqrt(anc_comp$se.w^2 + anc_comp$se.afr^2)), 
                                      df=pmin(anc_comp$df.w, anc_comp$df.afr), lower.tail=F))

anc_comp$sig <- ifelse(anc_comp$p_diff<.05/length(rep_traits_afr[rep_traits_afr %in% rep_traits_nbw]),1,0)

# Using interaction model
pwas_dat$afr <- ifelse(pwas_dat$ancestry==afr,1,0)
for (trait in rep_traits_afr[rep_traits_afr %in% rep_traits_nbw]){
  sex_check <- colSums(table(pwas_dat[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
  sexes <- names(sex_check[sex_check>0])
  if (length(sexes)<2) {
    model_covars <- substr(afr_covars, 5, nchar(afr_covars))
    b <- 23
  } else {
    model_covars <- afr_covars
    b <- 24
  }
  if (trait %in% bin_traits) {
    model <- glm(paste0(trait,'~',prs,'+',model_covars,'+afr:',prs), data=pwas_dat, family='binomial')
  } else {
    model <- lm(paste0(trait,'~',prs,'+',model_covars,'+afr:',prs), data=pwas_dat)
  }
  anc_comp[anc_comp$ukb_variable==trait,'p_diff'] <- summary(model)$coefficients[b,4]
  if (summary(model)$coefficients[b,4]<0.05/length(rep_traits_afr)) {
    results <- summary(model)$coefficients[b,c(1,2,4)]
    print(c(trait, round(results[1],3), round(results[2], 3), signif(results[3],3)))
  } 
}


# Write PheWAS output tables
# write.table(pwas_dat, paste0(pheno_dir, '/ukb.trait.prs.data.csv'), quote=F, row.names=F,col.names=T, sep=",")
write.table(pwas_results, paste0(output_dir, '/phewas_results.txt'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s2 <- pwas_results[pwas_results$ancestry==bw,
#                         c('display_name',ukb_var,'ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==bw,], paste0(output_dir, '/phewas_results_bw.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s2, paste0(output_dir, '/table_s2.tsv'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s5 <- pwas_results[pwas_results$ancestry==nbw,
#                         c('display_name',ukb_var,'ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==nbw,], paste0(output_dir, '/phewas_results_nbw.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s5, paste0(output_dir, '/table_s5.tsv'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s6 <- pwas_results[pwas_results$ancestry==afr,
#                         c('display_name',ukb_var,'ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==afr,], paste0(output_dir, '/phewas_results_afr.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s6, paste0(output_dir, '/table_s6.tsv'), quote=F, row.names=F,col.names=T, sep="\t")


###############################################################################################
### 4.  Comparison of PRS model PheWAS results                                              ###
###############################################################################################
# 
# Compare PheWAS results between GBMI and TAGC models
#

# pwas_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'), header=T)
# pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))


### Run PheWAS for TAGC PRS
tagc_score <- read.table(paste0(score_dir, '/PRScs.UKB_LD.TAGC_ME.UKB_scores.txt'), header=T)
prs_tagc <- 'zscore_tagc'
tagc_score[,prs_tagc] <- ztransform(tagc_score$SCORE_SUM)

###  PheWASs in BW for TAGC ###

pwas_dat_tagc <- merge(pwas_dat, tagc_score[,c('IID', prs_tagc)], by.x='eid',by.y='IID')

pwas_results_tagc <- data.frame(matrix(ncol=13, nrow=0))
colnames(pwas_results_tagc) <- c(names(pwas_metadata),'ancestry','sex','beta','se','p','n','df')
i <- 1 

for (trait in c(bin_traits, quant_traits)) {
    print(c(i, trait))
    sex_check <- colSums(table(pwas_dat_tagc[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
    sexes <- names(sex_check[sex_check>0])
    test_dat <- pwas_dat_tagc[(pwas_dat_tagc$sex %in% sexes) & (pwas_dat_tagc$ancestry==bw),
                         c(trait,prs_tagc,strsplit(covars,'+',fixed=T)[[1]][-1])]
    if (length(sexes)<2) {
      model_covars <- substr(covars, 5, nchar(covars))
    } else {
      model_covars <- covars
    }
    pwas_results_tagc[i,] <- tryCatch(
      {
        if (trait %in% bin_traits) {
          model <- glm(paste0(trait,'~',prs_tagc,model_covars), data=test_dat, family=binomial)
        } else {
          model <- lm(paste0(trait,'~',prs_tagc,model_covars), data=test_dat)
        }
        s <- summary(model)
        c(pwas_metadata[pwas_metadata$ukb_variable==trait,], bw, paste(sexes,collapse=' & '),
          model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals), df.residual(model))
      },
      error=function(e){
        return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], anc, rep(NA,5)))
      }
    )
    i <- i+1
}

write.table(pwas_results_tagc, paste0(output_dir, '/phewas_results_tagc.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")


scatter_dat <- merge(pwas_results_tagc[,c(ukb_var,'display_name','category','invert','beta','se','p','df')],
                     pwas_results[pwas_results$ancestry==bw,c(ukb_var,'beta','se','p')], by=ukb_var)

scatter_dat$beta.x <- ifelse(scatter_dat$invert==1, scatter_dat$beta.x*-1, scatter_dat$beta.x)
scatter_dat$beta.y <- ifelse(scatter_dat$invert==1, scatter_dat$beta.y*-1, scatter_dat$beta.y)

cor(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits), c('beta.x','beta.y')]) # 0.87

scatter_dat$beta_diff <- ifelse(scatter_dat$p.x>(0.05/371) | scatter_dat$p.y>(0.05/371), NA,
                                ifelse(scatter_dat$beta.x*scatter_dat$beta.y>0,
                                       abs(scatter_dat$beta.y)-abs(scatter_dat$beta.x),
                                       abs(scatter_dat$beta.y)+abs(scatter_dat$beta.x)))

scatter_dat$p_diff <- ifelse(scatter_dat$ukb_variable %in% bin_traits,
                             pchisq(((scatter_dat$beta.y-scatter_dat$beta.x)/
                                       sqrt(scatter_dat$se.y^2 + scatter_dat$se.x^2))^2, df=1, lower.tail=F),
                             2*pt(abs((scatter_dat$beta.y-scatter_dat$beta.x)/
                                        sqrt(scatter_dat$se.y^2 + scatter_dat$se.x^2)), df=scatter_dat$df, lower.tail=F))

scatter_dat$sig <- ifelse(scatter_dat$p_diff<(0.05/371),1,0)

length(scatter_dat[(scatter_dat$p.y<(0.05/371)) & (scatter_dat$p.x<(0.05/371)), 'category'])
sum(scatter_dat$sig)


# Test for categorical biases
for (i in 1:length(unique(scatter_dat$category))) {
  b<-0; q<-0; p_bin <- 'NA'; p_quant <- 'NA'; p <- 'NA'
  cat=sub('\n',' ',unique(scatter_dat$category)[i])
  #print(cat)
  temp_bin <- scatter_dat[scatter_dat$category==cat & scatter_dat$ukb_variable %in% bin_traits &
                          scatter_dat$p.x < 0.05 & scatter_dat$p.y < 0.05,]
  temp_quant <- scatter_dat[scatter_dat$category==cat & scatter_dat$ukb_variable %in% quant_traits &
                            scatter_dat$p.x < 0.05 & scatter_dat$p.y < 0.05,]
  n_bin <- nrow(temp_bin)
  n_quant <- nrow(temp_quant)
  if (n_bin > 2) {
    b <- 1
    model <- summary(lm(I(exp(beta.y))~0+exp(beta.x), data=temp_bin))
    p_bin <- min(2*pt(-abs((coef(model)[1]-1)/coef(model)[2]), df=nrow(temp_bin)-2), 1)
  }
  if (n_quant > 2) {
    q <- 1
    model <- summary(lm(I(beta.y)~0+beta.x, data=temp_quant))
    p_quant <- min(2*pt(-abs((coef(model)[1]-1)/coef(model)[2]), df=nrow(temp_quant)-2), 1)
  }
  if (b+q==2) {
    z_b <- qnorm(1-p_bin)
    z_q <- qnorm(1-p_quant)
    z <- (n_bin*z_b + n_quant*z_q)/sqrt(n_bin^2 + n_quant^2)
    p <- 1-pnorm(z)
  } else if (b==1) {
    p <- p_bin
  } else if (q==1) {
    p <- p_quant
  }
  if (p<(0.05/13)) {print(c(cat,p))}
}

max_gbmi_bin <- scatter_dat[!is.na(scatter_dat$beta_diff) & 
                              (scatter_dat$beta_diff==max(scatter_dat[scatter_dat$ukb_variable %in% bin_traits,
                                                                      'beta_diff'], na.rm=T)),'ukb_variable']
max_gbmi_quant <- scatter_dat[!is.na(scatter_dat$beta_diff) & 
                                (scatter_dat$beta_diff==max(scatter_dat[scatter_dat$ukb_variable %in% quant_traits,
                                                                        'beta_diff'], na.rm=T)),'ukb_variable']
max_tagc_bin <- scatter_dat[!is.na(scatter_dat$beta_diff) & 
                              (scatter_dat$beta_diff==min(scatter_dat[scatter_dat$ukb_variable %in% bin_traits,
                                                                      'beta_diff'], na.rm=T)),'ukb_variable']
max_tagc_quant <- scatter_dat[!is.na(scatter_dat$beta_diff) & 
                                (scatter_dat$beta_diff==min(scatter_dat[scatter_dat$ukb_variable %in% quant_traits,
                                                                        'beta_diff'], na.rm=T)),'ukb_variable']
save(scatter_dat, max_gbmi_bin, max_gbmi_quant, max_tagc_bin, max_tagc_quant, 
     file=paste0(output_dir,'phewas_comp.Rdata'))

###############################################################################################
### 5.  Analysis of non-asthmatics                                                          ###
###############################################################################################
# 
# Perform PheWAS in non-asthmatics and evaluate differences
#

# pwas_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'), header=T)

env_dat <- read.csv(paste0(pheno_dir, '/PRS-additional.pack-years.csv'), header=T)

pwas_dat <- merge(pwas_dat, env_dat[,c('eid','townsend_dprvtn_idx','pack_years_CALC')],
                  by='eid', all.x=T)

# determine controls (make list of asthmatics and then subtract from total)
asthmatics <- unique(c(subset(pwas_dat, asthma_all==1)$eid,
                       subset(pwas_dat, asthma_children==1)$eid,
                       subset(pwas_dat, asthma_adults==1)$eid,
                       subset(pwas_dat, asthma_self_reported==1)$eid,
                       subset(pwas_dat, asthma_doctor_diagnosed0==1)$eid,
                       subset(pwas_dat, blood_clot_asthma==1)$eid,
                       subset(pwas_dat, asthma_icd10==1)$eid))

pwas_dat$asthma_allPhenos <- pwas_dat$asthma_all
pwas_dat$asthma_allPhenos[pwas_dat$eid %in% asthmatics] <- 1

# Derive probability weights for IPW
model_covars <- paste(covars, 
                      c('bmi_body_size0+townsend_dprvtn_idx+pack_years_CALC+attendance_none'),sep='+')
prop_model <- glm(paste0('asthma_allPhenos ~',prs,model_covars), data=pwas_dat, family=binomial)
prop_score <- predict(prop_model, type='response')
weights <- ifelse(pwas_dat[,asthma]==1, 1/prop_score, 1/(1-prop_score))
pwas_dat$weights <- weights

pwas_dat_na <- pwas_dat[pwas_dat$ancestry==bw,]
pwas_dat_na <- pwas_dat_na[!(pwas_dat_na$eid %in% asthmatics),]

pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")
quant_traits <- pwas_metadata[pwas_metadata$bin==0,ukb_var]
bin_traits <- pwas_metadata[pwas_metadata$bin==1,ukb_var]

# Run PheWAS 
pwas_na_results <- data.frame(matrix(ncol=12, nrow=0))
colnames(pwas_na_results) <- c(names(pwas_metadata),'sex','beta','se','p','n', 'df')
i <- 1 
for (trait in c(bin_traits, quant_traits)) {
  print(c(i, trait))
  sex_check <- colSums(table(pwas_dat_na[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
  sexes <- names(sex_check[sex_check>0])
  if (length(sexes)<2) {
    model_covars <- substr(covars, 5, nchar(covars))
  } else {
    model_covars <- covars
  }
  test_dat <- pwas_dat_na[(pwas_dat_na$sex %in% sexes),
                          c('eid',trait,prs,'weights',strsplit(model_covars,'+',fixed=T)[[1]][-1])]
  pwas_na_results[i,] <- tryCatch(
    {
      if (trait %in% bin_traits) {
        test_dat_na <- test_dat[!(test_dat$eid %in% asthmatics),]
        model <- glm(paste0(trait,'~',prs,model_covars), data=test_dat_na, family=binomial, weights=weights)
      } else {
        test_dat[!is.na(test_dat[,trait]),trait] <- ztransform(test_dat[,trait])  # standardize quant traits
        test_dat_na <- test_dat[!(test_dat$eid %in% asthmatics),]
        model <- lm(paste0(trait,'~',prs,model_covars), data=test_dat_na, weights=weights)
      }
      s <- summary(model)
      c(pwas_metadata[pwas_metadata$ukb_variable==trait,], paste(sexes,collapse=' & '),
        model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals), df.residual(model))
    },
    error=function(e){
      return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], rep(NA,6)))
    }
  )
  i <- i+1
}

pwas_na_results <- merge(pwas_na_results, pwas_results[pwas_results$ancestry==bw,
                                                       c(ukb_var,'beta','se','p','df')], by=ukb_var)

names(pwas_na_results) <- gsub(x = names(pwas_na_results), pattern = "\\.x", replacement = "\\.na")  
names(pwas_na_results) <- gsub(x = names(pwas_na_results), pattern = "\\.y", replacement = "\\.a")  

pwas_na_results$delta_beta <- (pwas_na_results$beta.na/pwas_na_results$beta.a)-1


pwas_na_results$p_diff <- ifelse(pwas_na_results$ukb_variable %in% bin_traits,
                             pchisq(((pwas_na_results$beta.a-pwas_na_results$beta.na)/
                                       sqrt(pwas_na_results$se.a^2 + pwas_na_results$se.na^2))^2, df=1, lower.tail=F),
                             2*pt(abs((pwas_na_results$beta.a-pwas_na_results$beta.na)/
                                        sqrt(pwas_na_results$se.a^2 + pwas_na_results$se.na^2)), 
                                  df=pmin(pwas_na_results$df.a, pwas_na_results$df.na), lower.tail=F))

pwas_na_results$sig <- ifelse(pwas_na_results$p_diff<(0.05/371),1,0)


#table_s4 <- pwas_na_results[,c('display_name',ukb_var,'ukb_data_field','category','beta.na',
#                                'se.na','p.na','n','delta_beta','p_diff')]

#write.table(table_s4, paste0(output_dir, '/table_s4.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
#write.table(pwas_na_results, paste0(output_dir, '/pwas_noAsthma_results.tsv'), 
#             quote=F, row.names=F,col.names=T, sep="\t")

# write.table(pwas_dat, paste0(pheno_dir, '/ukb.trait.prs.data.csv'), quote=F, row.names=F,col.names=T, sep=",")


###############################################################################################
### 6.  HLA Effects                                                                         ###
###############################################################################################
# 
# This analysis evaluates the PRS prediction and trait associations using a score generated
#   without HLA region alleles
#

# Add HLA-excluded scores
scores_noHLA <- read.table(paste0(score_dir,"/PRScs.1KG_LD.GBMI.UKB_scores.no_HLA.txt"),header=T)
names(scores_noHLA) <- c('eid','ALLELE_CT','DOSAGE_SUM','PRS_noHLA')
pwas_dat <- merge(pwas_dat, scores_noHLA, by='eid')
pwas_dat$PRS_noHLA <- ztransform(pwas_dat$PRS_noHLA)

# Correlation between PRS and HLA-removed PRS 
cor(pwas_dat[,c(prs,'PRS_noHLA')])

# Compare asthma prediction between scores
for (phenotype in c(asthma, coa, aoa)){
  model <- glm(paste0(phenotype,'~PRS_noHLA',covars), data = pwas_dat[pwas_dat$ancestry==bw,], family = 'binomial')
  roc_auc_noHLA <- roc(pwas_dat[pwas_dat$ancestry==bw,phenotype],pwas_dat[pwas_dat$ancestry==bw,'PRS_noHLA'], plot=F, quiet=T)$auc[1]
  roc_auc_all <- roc(pwas_dat[pwas_dat$ancestry==bw,phenotype],pwas_dat[pwas_dat$ancestry==bw,prs], plot=F, quiet=T)$auc[1]
  roc_diff <- roc_auc_all-roc_auc_noHLA
  print(paste(phenotype,':    ', round(exp(summary(model)$coefficients[2,1]),4),
              round(exp(summary(model)$coefficients[2,1]-(1.96*summary(model)$coefficients[2,2])),4),
              round(exp(summary(model)$coefficients[2,1]+(1.96*summary(model)$coefficients[2,2])),4),
              coef(summary(model))[2,4],
              round((summary(model)$coefficients[2,1]/
                 pwas_results[(pwas_results$ancestry==bw) & 
                                (pwas_results$ukb_variable==phenotype),'beta'])-1,4), roc_auc_noHLA, roc_diff))
}

# generate boostrap permutations
auc_diff <- data.frame(matrix(ncol=3, nrow=n_replicates))
colnames(auc_diff) <- c(coa, aoa, asthma)
model_dat <- pwas_dat[pwas_dat$ancestry==bw,]
for (i in 1:n_replicates) {
  print(i)
  sample_dat <- model_dat[sample(nrow(model_dat), replace=TRUE),]
  for (phenotype in c(coa, aoa, asthma)){
    model <- glm(paste0(phenotype,'~PRS_noHLA',covars), data = sample_dat, family = 'binomial')
    roc_auc_noHLA <- roc(sample_dat[,phenotype],sample_dat[,'PRS_noHLA'], plot=F, quiet=T)$auc[1]
    roc_auc_all <- roc(sample_dat[,phenotype],sample_dat[,prs], plot=F, quiet=T)$auc[1]
    roc_diff <- roc_auc_all-roc_auc_noHLA
    auc_diff[i,phenotype] <- roc_diff
  }
}

# for asthma_all
prs_all_asthma <- glm(paste0(phenotype,'~',prs,covars), data = pwas_dat[pwas_dat$ancestry==bw,], family = 'binomial')  
round((summary(model)$coefficients[2,1]/summary(prs_all_asthma)$coefficients[2,1])-1,4)

# Perform PheWAS using HLA-removed score 
pwas_dat_hla <- pwas_dat[pwas_dat$ancestry==bw,]

#pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")
#quant_traits <- pwas_metadata[pwas_metadata$bin==0,ukb_var]
#bin_traits <- pwas_metadata[pwas_metadata$bin==1,ukb_var]

pwas_hla_results <- data.frame(matrix(ncol=12, nrow=0))
colnames(pwas_hla_results) <- c(names(pwas_metadata),'sex','beta','se','p','n','df')
i <- 1 
for (trait in c(bin_traits, quant_traits)) {
  print(c(i, trait))
  sex_check <- colSums(table(pwas_dat_hla[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
  sexes <- names(sex_check[sex_check>0])
  test_dat <- pwas_dat_hla[(pwas_dat_hla$sex %in% sexes),
                          c(trait,'PRS_noHLA',strsplit(covars,'+',fixed=T)[[1]][-1])]
  if (length(sexes)<2) {
    model_covars <- substr(covars, 5, nchar(covars))
  } else {
    model_covars <- covars
  }
  pwas_hla_results[i,] <- tryCatch(
    {
      if (trait %in% bin_traits) {
        model <- glm(paste0(trait,'~PRS_noHLA',model_covars), data=test_dat, family=binomial)
      } else {
        test_dat[!is.na(test_dat[,trait]),trait] <- ztransform(test_dat[,trait])  # standardize quant traits
        model <- lm(paste0(trait,'~PRS_noHLA',model_covars), data=test_dat)
      }
      s <- summary(model)
      c(pwas_metadata[pwas_metadata$ukb_variable==trait,], paste(sexes,collapse=' & '),
        model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals), df.residual(model))
    },
    error=function(e){
      return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], rep(NA,6)))
    }
  )
  i <- i+1
}

pwas_hla_results <- merge(pwas_hla_results, pwas_results[pwas_results$ancestry==bw,
                                                       c(ukb_var,'beta','se','p','df')], by=ukb_var)

names(pwas_hla_results) <- gsub(x = names(pwas_hla_results), pattern = "\\.x", replacement = "\\.hla")  
names(pwas_hla_results) <- gsub(x = names(pwas_hla_results), pattern = "\\.y", replacement = "\\.a")  

pwas_hla_results$delta_beta <- (pwas_hla_results$beta.hla/pwas_hla_results$beta.a)-1

pwas_hla_results$p_diff <- ifelse(pwas_hla_results$ukb_variable %in% bin_traits,
                                 pchisq(((pwas_hla_results$beta.a-pwas_hla_results$beta.hla)/
                                           sqrt(pwas_hla_results$se.a^2 + pwas_hla_results$se.hla^2))^2, df=1, lower.tail=F),
                                 2*pt(abs((pwas_hla_results$beta.a-pwas_hla_results$beta.hla)/
                                            sqrt(pwas_hla_results$se.a^2 + pwas_hla_results$se.hla^2)), 
                                      df=pmin(pwas_hla_results$df.a, pwas_hla_results$df.hla), lower.tail=F))

pwas_hla_results$sig <- ifelse(pwas_hla_results$p_diff<.05/371,1,0)
View(pwas_hla_results)

#table_s7 <- pwas_hla_results[,c('display_name',ukb_var,'ukb_data_field','category','beta.hla',
#                                'se.hla','p.hla','n','delta_beta', 'p_diff')]
#write.table(table_s7, paste0(output_dir, '/table_s7.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
#write.table(pwas_hla_results, paste0(output_dir, '/pwas_noHLA_results.tsv'), 
#             quote=F, row.names=F,col.names=T, sep="\t")


###############################################################################################
### 7.  Trait Mediation Analysis                                                            ###
###############################################################################################

# This analysis executes mediation analyses for eosinophil percentage, FEV1/FVC ratio, 
#   age of hay fever diagnosis, and BMI, and it also produces Figure 8 from the manuscript.
#   The equations used to estimate the mediation proportions were adopted from 
#   MacKinnon 2007 (doi: 10.1177/1740774507083434) and Li 2007 (doi: 10.1002/sim.2730).
#   Mediation p-values were estimated using the ROBMED package (Alfons 2021, 
#   doi: 10.1177/1094428121999096)
#

#library(robmed)

set.seed(20220117)

# pwas_dat <- read.csv(['path/to/ukb/data/pwas_data.csv'])


# Calculate proportion mediated (we reported the a*b/(c_std) values, as from MacKinnon 2007)
traits <- c('eosinophill_count0','fev1_fvc_ratio','age_hay_fever_diagnosed')
phenos <- c('asthma_children','asthma_adults','asthma_all')

med_df <- setNames(data.frame(matrix(ncol = 3, nrow = 3)), phenos)
med_95L_df <- setNames(data.frame(matrix(ncol = 3, nrow = 3)), phenos)
med_95H_df <- setNames(data.frame(matrix(ncol = 3, nrow = 3)), phenos)
rownames(med_df) <- traits
rownames(med_95L_df) <- traits
rownames(med_95L_df) <- traits

n_replicates <- 1000

for (pheno in phenos){
  for (trait in traits){
    model_dat <- pwas_dat[pwas_dat$ancestry==bw,]
    model_covars <- paste(covars, 
                          c('bmi_body_size0+townsend_dprvtn_idx+pack_years_CALC+attendance_none'),sep='+')
    model_covars <- covars
    if (pheno=='asthma_children') {
      test_dat <- model_dat[!(is.na(model_dat$asthma_children) & model_dat$asthma_all==1),]
      test_dat <- test_dat[!is.na(model_dat$asthma_adults) & model_dat$asthma_adults!=1,]
    } else if (pheno=='asthma_adults') {
      test_dat <- model_dat[!(is.na(model_dat$asthma_adults) & model_dat$asthma_all==1),]
      test_dat <- test_dat[!is.na(model_dat$asthma_children) & model_dat$asthma_children!=1,]
    } else {
      test_dat <- model_dat
    }
    
    model.y <- glm(paste0(pheno, '~', prs, model_covars), data=test_dat, family='binomial')
    model.m <- lm(paste0(trait,'~', prs , model_covars), data=test_dat)
    model.total <- glm(paste0(pheno,'~', trait, '+', prs, model_covars), data=test_dat, family='binomial')
    
    c <- coefficients(model.y)[[prs]]
    c_prime <- coefficients(model.total)[[prs]]
    b2 <- (coefficients(model.total)[[trait]])^2
    sigma2_mx <- var(model.m$residuals)
    norm_factor <- sqrt(1+((b2*sigma2_mx)/((pi^2)/3)))
    c_std <- c*norm_factor
    
    a <- coefficients(model.m)[[prs]]
    b <- coefficients(model.total)[[trait]]
    
    
    # calculate 95% CI for r2
    prop_values <- rep(NA, n_replicates)
    for (r in 1:n_replicates) {
      print(r)
      sample_dat <- test_dat[sample(nrow(test_dat), replace=TRUE),]
      model.y <- glm(paste0(pheno, '~', prs, model_covars), data=sample_dat, family='binomial')
      model.m <- lm(paste0(trait,'~', prs , model_covars), data=sample_dat)
      model.total <- glm(paste0(pheno,'~', trait, '+', prs, model_covars), data=sample_dat, family='binomial')
      
      c <- coefficients(model.y)[[prs]]
      c_prime <- coefficients(model.total)[[prs]]
      b2 <- (coefficients(model.total)[[trait]])^2
      sigma2_mx <- var(model.m$residuals)
      norm_factor <- sqrt(1+((b2*sigma2_mx)/((pi^2)/3)))
      c_std <- c*norm_factor
      a <- coefficients(model.m)[[prs]]
      b <- coefficients(model.total)[[trait]]
      prop_values[r] <- round(a*b/(c_std),6)
    }
    
    lower_ci_r2 <- quantile(as.numeric(prop_values), 0.025)
    upper_ci_r2 <- quantile(as.numeric(prop_values), 0.975)
    
    #print(c(pheno, trait, round(a,6), round(c_std,6), round(a*b/((a*b)+c_prime),6), round(a*b/(c_std),6)))
    med_df[trait, pheno] <- round(a*b/(c_std),6)
    med_95L_df[trait, pheno] <- lower_ci_r2
    med_95H_df[trait, pheno] <- upper_ci_r2
  }
}

save(med_df, med_95L_df, med_95H_df, file=paste0(output_dir, 'mediation_results.RData'))
