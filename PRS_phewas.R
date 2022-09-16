####  LOAD LIBRARIES  ##############################################################################

library(pROC)
library(sjmisc)
library(plotrix)
library(scales)
library(emmeans)
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(grid)
#library(gridExtra)
#library(extrafont)
library(forestplot)
library(robmed)

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


####  SET PARAMETERS  ##############################################################################
#
#  NOTE: even if running an individual section, these parameters should be saved first
#

bw <- 'White'; nbw <- 'non-British White'; afr <- 'African'
prs <- 'zscore'; covars <- paste0('+sex',paste(paste0('+pc',1:10),collapse=''))
afr_covars <- paste0('+sex',paste(paste0('+pc',1:20),collapse=''))
coa <- 'asthma_children'; aoa <- 'asthma_adults'; asthma <- 'asthma_all'

# These paths are user-specific
# score_dir <- [/path/to/directory/with/prs/scores]
# pheno_dir <- [/path/to/directory/with/phenotype/data]
# output_dir <- [/path/to/output/directory]


####  READ INPUT DATA  ############################################################################

# Read asthma phenotype data
asthma_phenos <- read.csv(paste0(pheno_dir, "/asthma--data.csv"))  # asthma phenotypes
asthma_covars <- read.csv(paste0(pheno_dir, "/asthma--covariates.csv"))  # regression covariates
asthma_dat <- merge(asthma_phenos, asthma_covars, by='eid')

trait_dat <- read.csv(paste0(pheno_dir, "/pwas--traits.csv"))  # ukb ethnicity and trait data

trait_dat <- merge(asthma_dat[,!names(asthma_dat) %in% names(trait_dat)[2:ncol(trait_dat)]], trait_dat,
                  by='eid')  # combine traits and covariates, avoid duplicate columns

trait_dat$ancestry <- ifelse((trait_dat$ethnic_background==1001 & trait_dat$genetic_ethnic_grouping==1), bw,
                            ifelse(trait_dat$ethnic_background %in% c(1002,1003), nbw,
                                   ifelse((trait_dat$ethnic_background %in% c(4,4001,4002,4003) & 
                                             trait_dat$pc1>150),afr,NA)))  # apply ancestry definitions

# write.table(trait_dat, paste0(pheno_dir, '/ukb.trait.data.csv'), quote=F, row.names=F,col.names=T, sep=",")

rm(asthma_phenos, asthma_covars, asthma_dat)


###############################################################################################
### 1.  PRS Model Comparison                                                                ###
###############################################################################################
# 
# This analysis executes compares the asthma prediction results in each UKB ancestry group
#   between the PRS models
#

# library(pROC)

# trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.data.csv'), header=T)

# Read PRS scores, compute asthma prediction
predict_df <- data.frame(matrix(ncol=9, nrow=0))
colnames(predict_df) <- c('LD','TAGC_pop', 'UKB_ancestry','Asthma_pheno','AUC','OR','OR_95CI_L','OR_95CI_U','P')
i <- 1    
for (ld in c('1KG','UKB')) {
  for (pop in c('EUR','ME')) {
    print(paste0('PRScs.',ld,'_LD.TAGC_',pop))
    scores <- read.table(paste0(score_dir, '/PRScs.',ld,'_LD.TAGC_',pop,'.UKB_scores.txt'), header=T)
    scores[,prs] <- ztransform(scores$SCORE_SUM)
    scores <- merge(trait_dat, scores[,c("IID", prs)], by.x="eid", by.y="IID")
    for (anc in c(bw, nbw, afr)) {
      print(paste0('  ',anc))
      for (case in c(asthma, coa, aoa)) {
        print(paste0('    ',case))
        roc_auc <- roc(scores[scores$ancestry==anc,case],scores[scores$ancestry==anc,prs], plot=F, quiet=T)$auc[1]
        if (anc==afr) {
          model <- glm(paste0(case,'~',prs,afr_covars), data=scores[scores$ancestry==anc,], family='binomial')
        } else {
          model <- glm(paste0(case,'~',prs,covars), data=scores[scores$ancestry==anc,], family='binomial')
        }
        predict_df[i,] <- c(ld, pop, anc, case, roc_auc, exp(model$coefficients[2])[[1]], 
                            exp(model$coefficients[2]-1.96*sqrt(diag(vcov(model)))[2])[[1]],
                            exp(model$coefficients[2]+1.96*sqrt(diag(vcov(model)))[2])[[1]],
                            coef(summary(model))[2,4])
        i <- i+1
      }
    }
  }
}

write.table(predict_df, paste0(output_dir, '/asthma_prs_prediction.txt'), quote=F, row.names=F,col.names=T, sep="\t")


###############################################################################################
### 2.  PRS Asthma Prediction                                                               ###
###############################################################################################
# 
# This analysis assigns individuals to PRS quantiles, computes relative risks, and generates
#   the plots for Figures 1 & 2
#

#library(sjmisc)
#library(plotrix)
#library(scales)
#library(emmeans)
#library(ggplot2)
#library(ggnewscale)

# trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.data.csv'), header=T)

score <- read.table(paste0(score_dir, '/PRScs.UKB_LD.TAGC_ME.UKB_scores.txt'), header=T)
score[,prs] <- ztransform(score$SCORE_SUM)
trait_dat <- merge(trait_dat, score[,c("IID", prs)], by.x="eid", by.y="IID")

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
      tmp_dat <- trait_dat[trait_dat[,q] %in% c(quant, '(40-60]'),]
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

# write.table(trait_dat, paste0(pheno_dir, '/ukb.trait.prs.data.csv'), quote=F, row.names=F,col.names=T, sep=",")


#####  FIGURE 1 (PRS Asthma Prediction) #####

### Quantiles - All asthma (Figure 1A)
pdf(file=paste0(output_dir, '/Fig_1A.pdf'), width=6, height=5)
par(xpd=F, mar=c(4.5,4.5,2,2), lheight = .3)
plot(x=1:11, ylim=c(0.093,2.407), xlab='', ylab='', yaxt='n', xaxt = "n", type='n')
abline(h = seq(0,2.5,.5), lty = 2, col = 'grey80')
abline(h = 1, lty = 1, col = 'grey')
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_asthma$OR, li=quantile_dat_asthma[,c('95CI_L')], ui=quantile_dat_asthma[,c('95CI_U')], 
       xlab='\nPRS Quantile (%)', ylab='Odds Ratio (OR)', xaxt = "n", pt.bg='#9F10F8', sfrac=0.001, lwd=2.5, cex.axis=.8,
       scol='#9F10F8', cex=1, pch=21, ylim=c(0.093,2.407), yaxt='n', cex.lab=0.9)
axis(1, at=1:11, labels = FALSE)
text(1:11, par("usr")[3] - 0.11, labels = quantile_dat_COA$Quantile, srt = 35, pos = 1, xpd = TRUE, cex=0.8)

axis(2,seq(0,2.5,.5), labels = c('0.0','0.5','1.0','1.5','2.0','2.5'), xpd = TRUE, cex.axis=0.8, las=2)
dev.off()

### Quantiles -  COA & AOA (Figure 1C)
pdf(file=paste0(output_dir, '/Fig_1C.pdf'), width=6, height=5)
par(xpd=F, mar=c(4.5,4.5,2,2), lheight = .3)
plot(x=1:11, ylim=c(0.135,3.5), xlab='', ylab='', yaxt='n', xaxt = "n", type='n')
abline(h = seq(0,4,.5), lty = 2, col = 'grey80')
abline(h = 1, lty = 1, col = 'grey')
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_COA$OR, li=quantile_dat_COA[,c('95CI_L')], ui=quantile_dat_COA[,c('95CI_U')], 
       xlab='\nPRS Quantile (%)', ylab='Odds Ratio (OR)', xaxt = "n", pt.bg='#22D378', sfrac=0.001, lwd=2.5, cex.axis=.8,
       scol='#22D378', cex=0.8, pch=21, ylim=c(0.135,3.5), yaxt='n', cex.lab=0.9)
axis(1, at=1:11, labels = FALSE)
text(1:11, par("usr")[3] - 0.18, labels = quantile_dat_COA$Quantile, srt = 35, pos = 1, xpd = TRUE, cex=0.8)
axis(2,seq(0,3.5,.5), labels = c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'), xpd = TRUE, cex.axis=0.8, las=2)
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_AOA$OR, li=quantile_dat_AOA[,c('95CI_L')], ui=quantile_dat_AOA[,c('95CI_U')], 
       xlab='', ylab='', xaxt = "n", pt.bg='#FB9700', sfrac=0.001, lwd=2.5, 
       scol='#FB9700', cex=0.8, pch=21, ylim=c(0.135,3.5), yaxt='n')
points(x=6,y=1, pch=21, lwd=2.5, bg=alpha('#22D378',.4), cex=0.8)
# Key labels
legend("topleft", inset=0.01, legend=c("COA", "AOA"), lwd=2, col=c('#22D378','#FB9700'),
       pch=c(16,16), bty='n', cex=0.8, pt.cex = 1.1)
legend("topleft", inset=0.01, col='black', legend=c("", ""),  lty=F, lwd=1.5,
       pch=c(21,21), bty='n', cex=0.8, pt.cex = 1.1)
dev.off()

### DENSITY PLOT (Fig 1B)
plot_dat <- trait_dat[trait_dat$ancestry==bw,]

pdf(file=paste0(output_dir, '/Fig_1B.pdf'), width=6, height=5)
par(xpd=F, mar=c(4.1,4.5,2,2), lheight = .3)
plot(density(plot_dat[!is.na(plot_dat$asthma_adults) & (plot_dat$asthma_adults==1),'zscore']), col='orange2',lwd=1.5,
     xlim=c(-4,4), ylim=c(0.015,0.4), main='', xaxt='n', yaxt='n', xlab='PRS\n')
axis(2,seq(0,0.4,.1), labels = c('0.0','0.1','0.2','0.3','0.4'), xpd = TRUE, cex.axis=0.85, las=2)
par(new=TRUE)
plot(density(plot_dat[!is.na(plot_dat$asthma_children) & (plot_dat$asthma_children==1),'zscore']), col='seagreen3',lwd=1.5,
     xlim=c(-4,4), ylim=c(0.015,0.4), main='', xaxt='n', yaxt='n', xlab='', ylab='')
par(new=TRUE)
plot(density(plot_dat[!is.na(plot_dat$asthma_all) & (plot_dat$asthma_all==1),'zscore']), col='purple2',lwd=1.5,
     xlim=c(-4,4), ylim=c(0.015,0.4), main='', xaxt='n', yaxt='n', xlab='', ylab='')
par(new=TRUE)
plot(density(plot_dat[!is.na(plot_dat$asthma_all) & (plot_dat$asthma_all==0),'zscore']), col='black',lwd=1.5,
     xlim=c(-4,4), ylim=c(0.015,0.4), main='', xaxt='n', yaxt='n', xlab='', ylab='')
abline(v = mean(plot_dat[!is.na(plot_dat$asthma_adults) & (plot_dat$asthma_adults==1),'zscore']), lty=3, lwd=1.5, col='#FB9700')
abline(v = mean(plot_dat[!is.na(plot_dat$asthma_children) & (plot_dat$asthma_children==1),'zscore']), lty=3, lwd=1.5, col='#22D378')
abline(v = mean(plot_dat[!is.na(plot_dat$asthma_all) & (plot_dat$asthma_all==1),'zscore']), lty=3,lwd=1.5, col='#9F10F8')
abline(v = mean(plot_dat[!is.na(plot_dat$asthma_all) & (plot_dat$asthma_all==0),'zscore']), lty=3, lwd=1.5, col='black')
axis(1, at=-4:4, labels = NA)
text(-4:4, par("usr")[3] - 0.013, labels = c('-4','-3','-2','-1','0','1','2','3','4'), pos = 1, xpd = TRUE, cex=0.85)
legend("topleft", inset=0.02, legend=c("No Asthma", "Asthma","COA", "AOA"), 
       pt.bg=c('black','purple2','seagreen3','orange2'),
       pch=c(22,22,22,22), bty='n', cex=0.85, pt.cex = 1.15)
dev.off()

### AUC PLOT (Fig 1D)
pdf(file=paste0(output_dir, '/Fig_1D.pdf'), width=6, height=5)
par(mar = c(4, 4, 4, 4)+.1)
plot.roc(plot_dat$asthma_all, plot_dat$zscore, main="", xlim=c(0.963,0.037), ylim=c(0.037,0.963),
         col='purple2', asp=F, xaxt='n', yaxt='n')
legend("bottomright", title="AUC", inset=0.02, legend=c("COA: 0.640", "AOA: 0.561", "All:    0.588"), 
       col=c('seagreen3','orange2','purple2'),
       pch=c(15,15,15), bty='n', cex=0.85, pt.cex = 1.15)
plot.roc(plot_dat$asthma_children, plot_dat$zscore,col='seagreen3', add=TRUE)
plot.roc(plot_dat$asthma_adults, plot_dat$zscore ,col='orange2', add=TRUE)
axis(2, seq(0,1,.2), labels = c('0.0','0.2','0.4','0.6','0.8','1.0'), xpd = TRUE, cex.axis=0.85, las=2)
axis(1, seq(1,0,-.2), labels = NA)
text(seq(1,0,-.2), par("usr")[3] - 0.03, labels = c('1.0','0.8','0.6','0.4','0.2','0.0'), pos = 1, xpd = TRUE, cex=0.85)
dev.off()


#####  FIGURE 2 (PRS by Sex and Age of Onset) #####

### PRS VS RISK LINEAR SCALE (Fig 2A)
range <- seq(from=-6, to=6, by=0.1)

fit_coa <- glm(asthma_children ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + zscore*sex, 
               data=trait_dat[trait_dat$ancestry=='White',], family='binomial')
c_pval <- summary(fit_coa)$coefficients[nrow(summary(fit_coa)$coefficient),4]
c_emm <- emmeans(fit_coa,~zscore*sex,type="response",at= list(zscore=range, sex=c('M','F')))
c_emm2 <- as.data.frame(c_emm)
c_emm2$lower <- confint(c_emm,adjust='none',level=0.95)$asymp.LCL
c_emm2$upper <- confint(c_emm,adjust='none',level=0.95)$asymp.UCL

fit_aoa <- glm(asthma_adults ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + sex*zscore,
               data=trait_dat[trait_dat$ancestry=='White',], family='binomial')
a_pval <- summary(fit_aoa)$coefficients[nrow(summary(fit_aoa)$coefficient),4]
a_emm <- emmeans(fit_aoa,~zscore*sex,type="response",at= list(zscore=range,sex=c('M','F')))
a_emm2 <- as.data.frame(a_emm)
a_emm2$lower<-confint(a_emm,adjust='none',level=0.95)$asymp.LCL
a_emm2$upper<-confint(a_emm,adjust='none',level=0.95)$asymp.UCL

pdf(file=paste0(output_dir, '/Fig_2A.pdf'), width=6, height=5)
ggplot(a_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.9) + 
  geom_line(data=a_emm2, aes(color=as.factor(sex))) + labs(fill="Adult onset") + 
  scale_fill_manual(labels=c("Male","Female"), values = c('#BFA373','#FDC372')) + 
  scale_colour_manual(values=c('#934900','#D16E00'),guide='none') +
  coord_cartesian(ylim = c(0.014, 0.27), xlim = c(-4.5,4.5)) +
  scale_y_continuous(name ="Risk", limits=c(0,.35), breaks = c(0,0.05,0.1,0.15,0.2,0.25), 
                     labels = c('0%','5%','10%','15%','20%','25%')) +
  scale_x_continuous(name='PRS',limits = c(-5,5), breaks=c(-4:4), 
                     labels=c('-4','-3','-2','-1','0','1','2','3','4')) +
  new_scale_fill() + new_scale_colour() + theme_bw() +
  geom_ribbon(data=c_emm2, aes(ymin = lower, ymax = upper,fill=as.factor(sex)), alpha=0.75) +
  geom_line(data=c_emm2, aes(color=as.factor(sex))) + labs(fill="Childhood onset") +
  scale_fill_manual(labels=c("Male","Female"), values = c('#409348','seagreen2')) +
  scale_colour_manual(values=c('#096000','seagreen4'),guide='none') +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 12),
        strip.text = element_text(size=12), legend.text=element_text(size=10),
        legend.title=element_text(size=10, face = "bold"), plot.title = element_text(size = 18, face = "bold"),
        panel.border=element_rect(colour = "black", fill=NA, size=1.25),
        legend.justification = c(0.05, 0.94), legend.position = c(0.05, 0.94),
        legend.box.background = element_rect(fill=alpha('white',0.75), color=alpha('gray80',0.75)),
        axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=1))
dev.off()

### PRS VS COA RISK LOG SCALE (FIG 2B)
pdf(file=paste0(output_dir, '/Fig_2B.pdf'), width=5, height=4)
ggplot(c_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.75) + 
  geom_line(data=c_emm2, aes(color=as.factor(sex))) + labs(fill="Childhood onset") +
  scale_fill_manual(labels=c("Male","Female"), values = c('#409348','seagreen2')) +
  scale_colour_manual(values=c('#096000','seagreen4'),guide='none') +
  coord_cartesian(xlim = c(-3.68,3.68)) + theme_bw() +
  scale_y_continuous(name ="Risk", trans = log_trans(), limits=c(0.002,.20),
                     breaks = c(0.0025,0.005,0.01,0.025,0.05,0.15),
                     labels = c('0.25%','0.5%','1%','2.5%','5%','15%')) +
  scale_x_continuous(name='PRS', limits = c(-4.5,4.5), breaks=c(-4:4), labels=c('-4','-3','-2','-1','0','1','2','3','4')) +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 12),
        strip.text = element_text(size=12), legend.text=element_text(size=10.5),
        legend.title=element_text(size=11, face = "bold"), plot.title = element_text(size = 18, face = "bold"),
        panel.border=element_rect(colour = "black", fill=NA, size=1.25),
        legend.justification = c(0.95, 0.06), legend.position = c(0.95, 0.06),
        legend.box.background = element_rect(fill=alpha('white',0.75), color=alpha('gray80',0.75)),
        axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=1))
dev.off()

### PRS VS AOA RISK LOG SCALE (FIG 2C)
pdf(file=paste0(output_dir, '/Fig_2C.pdf'), width=5, height=4)
ggplot(a_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.9) + 
  geom_line(data=a_emm2, aes(color=as.factor(sex))) + labs(fill="Adult onset") + 
  scale_fill_manual(labels=c("Male","Female"), values = c('#BFA373','#FDC372')) + 
  scale_colour_manual(values=c('#934900','#D16E00'),guide='none') +
  coord_cartesian(xlim = c(-3.68,3.68)) +
  scale_y_continuous(name ="Risk", trans = log_trans(), limits=c(0.002,.20),
                     breaks = c(0.0025,0.005,0.01,0.025,0.05,0.15),
                     labels = c('0.25%','0.5%','1%','2.5%','5%','15%')) + theme_bw() +
  scale_x_continuous(name='PRS', limits = c(-4.5,4.5), breaks=c(-4:4), labels=c('-4','-3','-2','-1','0','1','2','3','4')) +
  theme(axis.title=element_text(size=14,face="bold"), axis.text = element_text(size = 12),
        strip.text = element_text(size=12), legend.text=element_text(size=10.5),
        legend.title=element_text(size=11, face = "bold"), plot.title = element_text(size = 18, face = "bold"),
        panel.border=element_rect(colour = "black", fill=NA, size=1.25),
        legend.justification = c(0.95, 0.06), legend.position = c(0.95, 0.06),
        legend.box.background = element_rect(fill=alpha('white',0.75), color=alpha('gray80',0.75)),
        axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=1))
dev.off()

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
# library(extrafont)

# trait_dat <- read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'), header=T)


pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")

quant_traits <- pwas_metadata[pwas_metadata$bin==0,'ukb_variable']
bin_traits <- pwas_metadata[pwas_metadata$bin==1,'ukb_variable']

pwas_dat <- trait_dat[,c('eid',asthma, prs,'ancestry', strsplit(afr_covars,'+',fixed=T)[[1]][-1],
                         bin_traits, quant_traits)]

for (trait in quant_traits) {
  pwas_dat[!is.na(pwas_dat[,trait]),trait] <- ztransform(pwas_dat[,trait])  # Standardize quant traits
}

# Run PheWAS in each ancestry group
pwas_results <- data.frame(matrix(ncol=12, nrow=0))
colnames(pwas_results) <- c(names(pwas_metadata),'ancestry','sex','beta','se','p','n')
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
          model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals))
        },
      error=function(e){
        return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], anc, rep(NA,5)))
        }
      )
    i <- i+1
  }
}

# Write PheWAS output tabeles

# write.table(pwas_dat, paste0(pheno_dir, '/ukb.pwas.data.csv'), quote=F, row.names=F,col.names=T, sep=",")
write.table(pwas_results, paste0(output_dir, '/phewas_results.txt'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s2 <- pwas_results[pwas_results$ancestry==bw,
#                         c('display_name','ukb_variable','ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==bw,], paste0(output_dir, '/phewas_results_bw.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s2, paste0(output_dir, '/table_s2.tsv'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s4 <- pwas_results[pwas_results$ancestry==nbw,
#                         c('display_name','ukb_variable','ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==nbw,], paste0(output_dir, '/phewas_results_nbw.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s4, paste0(output_dir, '/table_s4.tsv'), quote=F, row.names=F,col.names=T, sep="\t")

# table_s5 <- pwas_results[pwas_results$ancestry==afr,
#                         c('display_name','ukb_variable','ukb_data_field','category','beta','se','p','n')]
write.table(pwas_results[pwas_results$ancestry==afr,], paste0(output_dir, '/phewas_results_afr.txt'), 
            quote=F, row.names=F,col.names=T, sep="\t")
# write.table(table_s5, paste0(output_dir, '/table_s5.tsv'), quote=F, row.names=F,col.names=T, sep="\t")


#####  FIGURE 3 (PheWAS Results Plot) #####
#  Create a horizontal 'Manhattan' plot for the PheWAS.
#    NOTE: Generating a vector image of the plot may be complicated by text-encoding issues,
#          and additional manipulation was required to get all the labels readable
#

#pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'), header=T, sep='\t')

fig_data <- pwas_results[pwas_results$ancestry==bw,]

# Put traits in desired order
min_group_p <- tibble(fig_data) %>% group_by(category) %>% summarise(min_p=min(p, na.rm = T))
min_group_p <- min_group_p[order(min_group_p$min_p),]
min_group_p$order <- 1:nrow(min_group_p)
fig_data <- merge(fig_data,min_group_p[,c('category','order')], by='category')
fig_data <- fig_data[order(fig_data$order,fig_data$display_name),]

# Remove duplicate traits
n_bonf <- nrow(fig_data)
fig_data <- data.frame(tibble(fig_data) %>% group_by(display_name) %>% filter(p==min(p)))
fig_data <- fig_data[!duplicated(fig_data$display_name),]

# Iterate over traits, assigning positions
n_i <- 1
for (i in seq(1,nrow(fig_data))) {
  if (i>1){
    if (fig_data[i,'category'] != fig_data[i-1,'category'] ) {
      n_i <- n_i+9
    }
  }
  if (fig_data[i,'p'] < 0.05/n_bonf) {
      n_i <- n_i+3
  }
  fig_data[i,'Pos'] <- n_i
  n_i <- n_i+1
}

fig_data$Pos <- n_i-fig_data$Pos

# Prepare Y-axis labels
fig_data[fig_data$category=='Reproductive, Endocrine & Metabolic','category'] <- 'Reproductive,\nEndocrine & Metabolic'
fig_data[fig_data$category=='Hematologic & Blood Chemistry','category'] <- 'Hematologic &\nBlood Chemistry'

fig_tibble <- tibble(fig_data)
y_axis <- fig_tibble %>% group_by(category) %>% summarize(center=(max(Pos) + min(Pos) ) / 2 )
rects <- fig_tibble %>% group_by(category) %>% summarize(mincat=min(Pos)-5,maxcat=max(Pos)+5)
totals <- fig_tibble %>% group_by(category) %>% count()
  
# Determine whether label will be added
fig_data$label <- T
fig_data[fig_data$p>=0.05/n_bonf,'label'] <- F
fig_data$low_pvals <- NA

fig_data[fig_data$display_name=='Asthma','low_pvals'] <- 'p<4.9e-324'
fig_data[fig_data$display_name=='Childhood-onset asthma','low_pvals'] <- 'p<4.9e-324'
fig_data[fig_data$display_name=='Adult-onset asthma','low_pvals'] <- 'p=2.9e-204'

fig_data[fig_data$display_name=='Eosinophils','low_pvals'] <- 'p<4.9e-324'
fig_data[fig_data$display_name=='Hay fever, allergic rhinitis, eczema','low_pvals'] <- 'p<4.9e-324'
fig_data[fig_data$display_name=='Tx: albuterol','low_pvals'] <- 'p=1.4e-265'
fig_data[fig_data$display_name=='Wheezing','low_pvals'] <- 'p=3.0e-258'


# Flip betas if necessary
fig_data[fig_data$invert==1,'beta'] <- fig_data[fig_data$invert==1,'beta']*-1
fig_data[fig_data$display_name=='Asthma','p'] <- 5e-324
fig_data[fig_data$display_name=='Eosinophils','p'] <- 5e-324
fig_data[fig_data$display_name=='Childhood-onset asthma','p'] <- 5e-324
fig_data[fig_data$display_name=='Hay fever, allergic rhinitis, eczema','p'] <- 5e-324

fig_data[fig_data$display_name=='Adult-onset asthma','Pos'] <- n_i - 19
#fig_data[fig_data$display_name=='Age asthma diagnosed','Pos'] <- fig_data[fig_data$display_name=='Age asthma diagnosed','Pos'] + 1
fig_data[fig_data$display_name=='Asthma','Pos'] <- n_i + 3
fig_data[fig_data$display_name=='Childhood-onset asthma','Pos'] <- n_i - 30
fig_data[fig_data$display_name=='Hay fever, allergic rhinitis, eczema','Pos'] <- n_i - 41
fig_data[fig_data$display_name=='Tx: albuterol','Pos'] <- n_i - 52
fig_data[fig_data$display_name=='Wheezing','Pos'] <- fig_data[fig_data$display_name=='Wheezing','Pos'] + 1

# Plot!
group_cols <- c('Neurologic & Behavioral'='goldenrod', 'Ocular'='cornflowerblue',
                'Asthma & Allergic Disease'='indianred', 'Pulmonary'='darkseagreen',
                'Cardiovascular'='palevioletred', 'Hematologic &\nBlood Chemistry'='cadetblue',
                'Gastrointenstinal'='orange', 'Renal & Urologic'='seagreen',
                'Reproductive,\nEndocrine & Metabolic'='royalblue', 'Anthropometric'='coral',
                'Musculoskeletal & Skin'='gold', 'Neoplasms'='orchid4',
                'Other'='plum4')

y_labels <- sub('\n','<br>',
                str_c("<br>**",y_axis$category,"**<br><span style='font-size:9.5pt;color:dimgray'>",
                      paste0("^((n=",totals$n,")^)</span>")))


cairo_pdf(filename=paste0(output_dir, '/Fig_3.pdf'), width=11, height=7.5)  # requires xquartz/encodings
ggplot(fig_data, aes(x=-log10(p),y=Pos)) + 
  
  # add non-significant traits
  geom_point(data=subset(fig_data, p>0.05/n_bonf), color="grey60", size=3.1,pch=19) +
  geom_point(data=subset(fig_data, p>0.05/n_bonf), aes(bg=as.factor(category)), alpha=.3,size=3.1,pch=21) +
  scale_fill_manual(values=group_cols) +
  
  # custom axes:
  scale_y_continuous(label=y_labels, breaks= y_axis$center) +
  scale_x_continuous(expand = c(0, 0), limits=c(-1,150)) + # remove space between plot area and x axis
  geom_vline(xintercept=-log10(0.05/n_bonf),linetype='dashed',color='grey30',size=0.7) +
  coord_cartesian(ylim=c(18,n_i-18), clip="off") +
  # Add label using ggrepel to avoid overlapping
  geom_text_repel(data=subset(fig_data, label==T), aes(label=display_name), size=2.6,
                  max.overlaps=6, point.padding=unit(0.4,'lines'),nudge_x=4,
                  box.padding=unit(0.35,'lines'),segment.color='grey',min.segment.length=0.2) +
  
  # Add highlighted points
  geom_point(data=subset(fig_data, p<0.05/n_bonf & beta>0), aes(bg=as.factor(category)), alpha=0.85, size=2.8,pch=24) +
  geom_point(data=subset(fig_data, p<0.05/n_bonf & beta<0), aes(bg=as.factor(category)), alpha=0.85, size=2.8,pch=25) +
  
  # Add super-significant labels
  geom_text(data=subset(fig_data,p<1e-150), family = 'Arial Unicode MS', hjust=1,
            aes(label=paste0(display_name,', ',low_pvals,"\u2192")), x = 150, size=2.6) +

  # add color blocks
  annotate("rect", xmin=-1,xmax=-log10(0.05/n_bonf), ymin=min(rects[,'mincat'][[1]])-7,
           ymax=max(rects$maxcat)+10,
           fill='grey15', alpha=0.2) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Renal & Urologic','mincat'][[1]],
           ymax=rects[rects$category=='Renal & Urologic','maxcat'][[1]],
           fill=group_cols[['Renal & Urologic']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Neurologic & Behavioral','mincat'][[1]]-7,
           ymax=rects[rects$category=='Neurologic & Behavioral','maxcat'][[1]],
           fill=group_cols[['Neurologic & Behavioral']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Neoplasms','mincat'][[1]],
           ymax=rects[rects$category=='Neoplasms','maxcat'][[1]],
           fill=group_cols[['Neoplasms']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Ocular','mincat'][[1]],
           ymax=rects[rects$category=='Ocular','maxcat'][[1]],
           fill=group_cols[['Ocular']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Cardiovascular','mincat'][[1]],
           ymax=rects[rects$category=='Cardiovascular','maxcat'][[1]],
           fill=group_cols[['Cardiovascular']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Musculoskeletal & Skin','mincat'][[1]],
           ymax=rects[rects$category=='Musculoskeletal & Skin','maxcat'][[1]],
           fill=group_cols[['Musculoskeletal & Skin']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Other','mincat'][[1]],
           ymax=rects[rects$category=='Other','maxcat'][[1]],
           fill=group_cols[['Other']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Anthropometric','mincat'][[1]],
           ymax=rects[rects$category=='Anthropometric','maxcat'][[1]],
           fill=group_cols[['Anthropometric']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Reproductive,\nEndocrine & Metabolic','mincat'][[1]],
           ymax=rects[rects$category=='Reproductive,\nEndocrine & Metabolic','maxcat'][[1]],
           fill=group_cols[['Reproductive,\nEndocrine & Metabolic']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Gastrointenstinal','mincat'][[1]],
           ymax=rects[rects$category=='Gastrointenstinal','maxcat'][[1]],
           fill=group_cols[['Gastrointenstinal']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Pulmonary','mincat'][[1]]-3,
           ymax=rects[rects$category=='Pulmonary','maxcat'][[1]],
           fill=group_cols[['Pulmonary']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Hematologic &\nBlood Chemistry','mincat'][[1]],
           ymax=rects[rects$category=='Hematologic &\nBlood Chemistry','maxcat'][[1]]+1,
           fill=group_cols[['Hematologic &\nBlood Chemistry']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=150, ymin=rects[rects$category=='Asthma & Allergic Disease','mincat'][[1]]-2,
           ymax=max(rects$maxcat)+10,
           fill=group_cols[['Asthma & Allergic Disease']], alpha=0.15) +
  
  # Customize the theme:
  theme_bw() +
  theme( 
    text=element_text(family="Helvetica"),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_markdown(hjust=1,lineheight=1.25),
    axis.title.y = element_blank(),
    axis.title.x = element_text(vjust=-2),
    axis.line.y = element_line(),
    axis.line.x = element_line(size=0.1),
    plot.margin=unit(c(0.5,1,1,0.2),"cm")
  )
dev.off()


###############################################################################################
### 4.  Analysis of non-asthmatics                                                          ###
###############################################################################################
# 
# Perform PheWAS in non-asthmatics and evaluate differences
#

# pwas_dat <- read.csv(paste0(pheno_dir, '/ukb.pwas.data.csv'), header=T)

pwas_dat_na <- pwas_dat[pwas_dat$ancestry==bw,]
# determine controls (make list of asthmatics and then subtract from total)
asthmatics <- unique(c(subset(pwas_dat_na, asthma_all==1)$eid,
                       subset(pwas_dat_na, asthma_children==1)$eid,
                       subset(pwas_dat_na, asthma_adults==1)$eid,
                       subset(pwas_dat_na, asthma_self_reported==1)$eid,
                       subset(pwas_dat_na, asthma_doctor_diagnosed0==1)$eid,
                       subset(pwas_dat_na, blood_clot_asthma==1)$eid,
                       subset(pwas_dat_na, asthma_icd10==1)$eid))

pwas_dat_na <- pwas_dat_na[!(pwas_dat_na$eid %in% asthmatics),]

pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")
quant_traits <- pwas_metadata[pwas_metadata$bin==0,'ukb_variable']
bin_traits <- pwas_metadata[pwas_metadata$bin==1,'ukb_variable']

# Run PheWAS 
pwas_na_results <- data.frame(matrix(ncol=11, nrow=0))
colnames(pwas_na_results) <- c(names(pwas_metadata),'sex','beta','se','p','n')
i <- 1 
for (trait in c(bin_traits, quant_traits)) {
  print(c(i, trait))
  sex_check <- colSums(table(pwas_dat_na[,c(trait,'sex')])>0)>=2  # checks case/control counts for each sex
  sexes <- names(sex_check[sex_check>0])
  test_dat <- pwas_dat_na[(pwas_dat_na$sex %in% sexes),
                          c(trait,prs,strsplit(covars,'+',fixed=T)[[1]][-1])]
  if (length(sexes)<2) {
    model_covars <- substr(covars, 5, nchar(covars))
  } else {
    model_covars <- covars
  }
  pwas_na_results[i,] <- tryCatch(
    {
      if (trait %in% bin_traits) {
        model <- glm(paste0(trait,'~',prs,model_covars), data=test_dat, family=binomial)
      } else {
        test_dat[!is.na(test_dat[,trait]),trait] <- ztransform(test_dat[,trait])  # standardize quant traits
        model <- lm(paste0(trait,'~',prs,model_covars), data=test_dat)
      }
      s <- summary(model)
      c(pwas_metadata[pwas_metadata$ukb_variable==trait,], paste(sexes,collapse=' & '),
        model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals))
    },
    error=function(e){
      return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], rep(NA,5)))
    }
    )
  i <- i+1
}

pwas_na_results <- merge(pwas_na_results, pwas_results[pwas_results$ancestry==bw,
                                                       c('ukb_variable','beta','se','p')], by='ukb_variable')

pwas_na_results$delta_beta <- (pwas_na_results$beta.x/pwas_na_results$beta.y)-1

# table_s3 <- pwas_na_results[,c('display_name','ukb_variable','ukb_data_field','category','beta.x',
#                                'se.x','p.x','n','delta_beta')]

# write.table(table_s3, paste0(output_dir, '/table_s3.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
# write.table(pwas_na_results, paste0(output_dir, '/pwas_noAsthma_results.tsv'), 
#             quote=F, row.names=F,col.names=T, sep="\t")

# write.table(pwas_dat, paste0(pheno_dir, '/ukb.pwas.data.csv'), quote=F, row.names=F,col.names=T, sep=",")


#####  FIGURE 4 (PheWAS associations in non-asthmatics, binary traits) #####
#
# NOTE: category colors and text emphasis were manually added
#

# library(forestplot)

## BW PheWAS Results
fig_data <- pwas_na_results[!is.na(pwas_na_results$p.x),]
names(fig_data) <- gsub(x = names(fig_data), pattern = "\\.x", replacement = "\\.na")  
names(fig_data) <- gsub(x = names(fig_data), pattern = "\\.y", replacement = "\\.a")  

fig_data[fig_data$invert==1,'beta.a'] <- fig_data[fig_data$invert==1,'beta.a']*-1
fig_data[fig_data$bin==1,'CI_L.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] - 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'CI_U.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] + 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'beta.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'])
fig_data[fig_data$bin==0,'CI_L.a'] <- fig_data[fig_data$bin==0,'beta.a'] - 1.96*fig_data[fig_data$bin==0,'se.a']
fig_data[fig_data$bin==0,'CI_U.a'] <- fig_data[fig_data$bin==0,'beta.a'] + 1.96*fig_data[fig_data$bin==0,'se.a']

# Only plot most significant association for redundant traits (e.g. icd10 vs self-report)
fig_data <- data.frame(tibble(fig_data) %>% group_by(display_name) %>% filter(p.a==min(p.a)))
fig_data <- fig_data[!duplicated(fig_data$display_name),]

fig_data[fig_data$ukb_variable=='blood_clot_hayfever','p.a'] <- 5e-324
fig_data[fig_data$ukb_variable=='blood_clot_none','p.a'] <- 5e-324
fig_data[fig_data$ukb_variable=='eosinophill_count0','p.a'] <- 5e-324

## BW PheWAS non-Asthmatics - Results
fig_data[fig_data$bin==1,'CI_L.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'] - 1.96*fig_data[fig_data$bin==1,'se.na'])
fig_data[fig_data$bin==1,'CI_U.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'] + 1.96*fig_data[fig_data$bin==1,'se.na'])
fig_data[fig_data$bin==1,'beta.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'])
fig_data[fig_data$bin==0,'CI_L.na'] <- fig_data[fig_data$bin==0,'beta.na'] - 1.96*fig_data[fig_data$bin==0,'se.na']
fig_data[fig_data$bin==0,'CI_U.na'] <- fig_data[fig_data$bin==0,'beta.na'] + 1.96*fig_data[fig_data$bin==0,'se.na']

# Filter data for only significant findings
fig_data <- fig_data[fig_data$p.a<(0.05/n_bonf),]

# Put data in order of beta difference, filter to top results
fig_data <- merge(fig_data,min_group_p[,c('category','order')], by='category')
fig_data <- fig_data[order(fig_data$order,fig_data$p.a),]

fig_data_quant <- fig_data[fig_data$bin==0,]
fig_data_bin <- fig_data[fig_data$bin==1,]

# Text on plot
traits <- fig_data_bin$display_name
p.a <- formatC(fig_data_bin$p.a, format = "e", digits = 0)
p.na <- formatC(fig_data_bin$p.na, format = "e", digits = 0)

tabletext <- list(
  c("Trait\n",traits),
  c("P-value\n(All)    \n", p.a),
  c("P-value\n(Non-asthma)\n", p.na)
)

rows <- as.character(1:36)
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_4.pdf'), width=10, height=9)  # requires xquartz/encodings
forestplot(tabletext, 
           mean = cbind(c(NA,fig_data_bin$beta.na), c(NA,fig_data_bin$beta.a)),
           lower = cbind(c(NA,fig_data_bin$CI_L.na), c(NA,fig_data_bin$CI_L.a)), 
           upper = cbind(c(NA,fig_data_bin$CI_U.na), c(NA,fig_data_bin$CI_U.a)),
           new_page = TRUE,
           clip = c(0.8,1.8), 
           hrzl_lines=line_list,
           zero=1,
           xlog = F, xlab = "OR", 
           col = fpColors(box = c("steelblue3", "firebrick4"),
                          lines = c("steelblue", "firebrick3")),
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawDiamondCI),
           is.summary = c(TRUE,rep(FALSE,nrow(fig_data))), 
           graph.pos = 2,
           graphwidth=unit(120,'mm'),
           boxsize = 0.35, 
           line.margin = .51,
           lineheight = unit(6,"mm"),
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("WB Non-asthmatics", "WB"), 
           ci.vertices = TRUE)
dev.off()


#####  FIGURE 5 (PheWAS associations in non-asthmatics, quantitative traits) #####
#
# NOTE: category colors and text bolding were manually added
#

# Text on plot
traits <- fig_data_quant$display_name
p.a <- formatC(fig_data_quant$p.a, format = "e", digits = 0)
p.na <- formatC(fig_data_quant$p.na, format = "e", digits = 0)

p.a[p.a=="0e+00"] <- "<5e-324"

tabletext <- list(
  c("Trait\n",traits),
  c("P-value\n(All)\n", p.a),
  c("P-value\n(Non-asthma)\n", p.na)
)

rows <- as.character(1:42)
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_5.pdf'), width=11, height=10)  # requires xquartz/encodings
forestplot(tabletext, 
           mean = cbind(c(NA,fig_data_quant$beta.na), c(NA,fig_data_quant$beta.a)),
           lower = cbind (c(NA,fig_data_quant$CI_L.na), c(NA,fig_data_quant$CI_L.a)), 
           upper = cbind(c(NA,fig_data_quant$CI_U.na), c(NA,fig_data_quant$CI_U.a)),
           new_page = TRUE,
           clip = c(-0.2,.2), 
           hrzl_lines=line_list,
           xlog = F, xlab = "Beta", 
           col = fpColors(box = c("steelblue3", "firebrick4"),
                          lines = c("steelblue", "firebrick3")),
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawDiamondCI),
           is.summary = c(TRUE,rep(FALSE,nrow(fig_data))), 
           graph.pos = 2,
           graphwidth=unit(120,'mm'),
           boxsize = 0.35, 
           line.margin = .51,
           lineheight = unit(6,"mm"),
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("WB Non-asthmatics", "WB All"), 
           ci.vertices = TRUE)
dev.off()


#####  FIGURE 6 (Eosinophil effects by population) #####

trait <- 'eosinophill_percentage0'

fig_data <- c(beta.bw=pwas_results[pwas_results$ancestry==bw & pwas_results$ukb_variable==trait,'beta'],
              se.bw=pwas_results[pwas_results$ancestry==bw & pwas_results$ukb_variable==trait,'se'],
              beta.na=pwas_na_results[pwas_na_results$ukb_variable==trait,'beta.x'],
              se.na=pwas_na_results[pwas_na_results$ukb_variable==trait,'se.x'],
              beta.nbw=pwas_results[pwas_results$ancestry==nbw & pwas_results$ukb_variable==trait,'beta'],
              se.nbw=pwas_results[pwas_results$ancestry==nbw & pwas_results$ukb_variable==trait,'se'],
              beta.afr=pwas_results[pwas_results$ancestry==afr & pwas_results$ukb_variable==trait,'beta'],
              se.afr=pwas_results[pwas_results$ancestry==afr & pwas_results$ukb_variable==trait,'se'])

fig_data <- c(fig_data,
              CI_L.bw=fig_data[['beta.bw']]-(1.96*fig_data[['se.bw']]), 
              CI_U.bw=fig_data[['beta.bw']]+(1.96*fig_data[['se.bw']]),
              CI_L.na=fig_data[['beta.na']]-(1.96*fig_data[['se.na']]),
              CI_U.na=fig_data[['beta.na']]+(1.96*fig_data[['se.na']]),
              CI_L.nbw=fig_data[['beta.nbw']]-(1.96*fig_data[['se.nbw']]),
              CI_U.nbw=fig_data[['beta.nbw']]-(1.96*fig_data[['se.nbw']]),
              CI_L.afr=fig_data[['beta.afr']]-(1.96*fig_data[['se.afr']]),
              CI_U.afr=fig_data[['beta.afr']]+(1.96*fig_data[['se.afr']]))

fig_data <- round(fig_data, 3)

fig_p <- c(p.bw=pwas_results[pwas_results$ancestry==bw & pwas_results$ukb_variable==trait,'p'],
           p.na=pwas_na_results[pwas_na_results$ukb_variable==trait,'beta.x'],
           p.nbw=pwas_results[pwas_results$ancestry==nbw & pwas_results$ukb_variable==trait,'p'],
           p.afr=pwas_results[pwas_results$ancestry==afr & pwas_results$ukb_variable==trait,'p'])

## Text on plot
pops <- c('White British','White British\nnon-asthmatics','White non-British ', 'African ancestries')
betas <- c(paste0(fig_data['beta.bw'], '  [', fig_data['CI_L.bw'], '-', fig_data['CI_U.bw'],']'),
           paste0(fig_data['beta.na'], '  [', fig_data['CI_L.na'], '-', fig_data['CI_U.na'],']'),
           paste0(fig_data['beta.nbw'], '  [', fig_data['CI_L.nbw'], '-', fig_data['CI_U.nbw'],']'),
           paste0(fig_data['beta.afr'], '  [', fig_data['CI_L.afr'], '0-', fig_data['CI_U.afr'],']'))

fig_p <- formatC(fig_p, format='e', digits=2)

fig_p[fig_p=='0.00e+00'] <- '<4.94e-324'

tabletext <- list(
  c("Population",pops),
  c("\U03B2  [95% CI]", betas),
  c("P", fig_p)
)

rows <- as.character(1:5)
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
forestplot(tabletext, 
           mean = c(NA, fig_data['beta.bw'], fig_data['beta.na'], fig_data['beta.nbw'], fig_data['beta.afr']),
           lower = c(NA,fig_data['CI_L.bw'], fig_data['CI_L.na'], fig_data['CI_L.nbw'], fig_data['CI_L.afr']), 
           upper = c(NA, fig_data['CI_U.bw'], fig_data['CI_U.na'], fig_data['CI_U.nbw'], fig_data['CI_U.afr']),
           new_page = TRUE,
           clip = c(-0.2,.2), 
           hrzl_lines=line_list,
           xlog = F, xlab=expression(beta["z"]), 
           col = fpColors(box="black", lines="slategray4"),
           fn.ci_norm = fpDrawDiamondCI,
           #fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI,fpDrawDiamondCI, fpDrawDiamondCI),
           is.summary = c(TRUE,rep(FALSE,4)), 
           graph.pos = 2,
           graphwidth=unit(60,'mm'),
           boxsize = 0.2, 
           line.margin = .71,
           lineheight = unit(20,"mm"),
           lwd.zero=1.75,
           lwd.ci=1.25,
           txt_gp= fpTxtGp(label=gpar(cex=0.75), title=gpar(cex=0.85), legend=gpar(cex=0.85), 
                           xlab=gpar(cex=0.95), ticks=gpar(cex=0.7), summary=gpar(cex=0.85)),
           ci.vertices = TRUE)


###############################################################################################
### 6.  HLA Effects                                                                         ###
###############################################################################################
# 
# This analysis evaluates the PRS prediction and trait associations using a score generated
#   without HLA region alleles
#   NOTE: Figure 7 was generated externally/manually using the results generated herein
#

# Add HLA-excluded scores
scores_noHLA <- read.table(paste0(score_dir,"/PRScs.UKB_LD.TAGC_ME.UKB_scores.no_HLA.txt"),header=T)
names(scores_noHLA) <- c('eid','PRS_noHLA')
pwas_dat <- merge(pwas_dat, scores_noHLA, by='eid')
pwas_dat$PRS_noHLA <- ztransform(pwas_dat$PRS_noHLA)

# Correlation between PRS and HLA-removed PRS 
cor(pwas_dat[,c(prs,'PRS_noHLA')])

# Compare asthma prediction between scores
for (phenotype in c(coa, aoa, asthma)){
  model <- glm(paste0(phenotype,'~PRS_noHLA',covars), data = pwas_dat[pwas_dat$ancestry==bw,], family = 'binomial')
  roc_auc <- roc(pwas_dat[pwas_dat$ancestry==bw,phenotype],pwas_dat[pwas_dat$ancestry==bw,'PRS_noHLA'], plot=F, quiet=T)$auc[1]
  print(paste(phenotype,':    ', round(exp(summary(model)$coefficients[2,1]),4),
              round(exp(summary(model)$coefficients[2,1]-(1.96*summary(model)$coefficients[2,2])),4),
              round(exp(summary(model)$coefficients[2,1]+(1.96*summary(model)$coefficients[2,2])),4),
              coef(summary(model))[2,4],
              round((summary(model)$coefficients[2,1]/
                 pwas_results[(pwas_results$ancestry==bw) & 
                                (pwas_results$ukb_variable==phenotype),'beta'])-1,4), roc_auc))
}

# for asthma_all
prs_all_asthma <- glm(paste0(phenotype,'~',prs,covars), data = pwas_dat[pwas_dat$ancestry==bw,], family = 'binomial')  
round((summary(model)$coefficients[2,1]/summary(prs_all_asthma)$coefficients[2,1])-1,4)

### AUC PLOT
plot_dat <- pwas_dat[pwas_dat$ancestry==bw,]
par(mar = c(4, 4, 4, 4)+.1)
roc_auc <- roc(plot_dat[,case],scores[scores$ancestry==anc,prs], plot=F, quiet=T)$auc[1]

plot.roc(plot_dat$asthma_all, plot_dat$PRS_noHLA, main="", xlim=c(0.963,0.037), ylim=c(0.037,0.963),
         col='purple2', asp=F, xaxt='n', yaxt='n', print.auc=T)
plot.roc(plot_dat$asthma_children, plot_dat$PRS_noHLA,col='seagreen3', add=TRUE, print.auc=T)
plot.roc(plot_dat$asthma_adults, plot_dat$PRS_noHLA ,col='orange2', add=TRUE, print.auc=T)
legend("bottomright", title="AUC", inset=0.02, legend=c("COA: 0.647", "AOA: 0.556", "All:    0.587"), 
       col=c('seagreen3','orange2','purple2'),
       pch=c(15,15,15), bty='n', cex=0.85, pt.cex = 1.15)
axis(2, seq(0,1,.2), labels = c('0.0','0.2','0.4','0.6','0.8','1.0'), xpd = TRUE, cex.axis=0.85, las=2)
axis(1, seq(1,0,-.2), labels = NA)
text(seq(1,0,-.2), par("usr")[3] - 0.03, labels = c('1.0','0.8','0.6','0.4','0.2','0.0'), pos = 1, xpd = TRUE, cex=0.85)

# Perform PheWAS using HLA-removed score 
pwas_dat_hla <- pwas_dat[pwas_dat$ancestry==bw,]

pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")
quant_traits <- pwas_metadata[pwas_metadata$bin==0,'ukb_variable']
bin_traits <- pwas_metadata[pwas_metadata$bin==1,'ukb_variable']

pwas_hla_results <- data.frame(matrix(ncol=11, nrow=0))
colnames(pwas_hla_results) <- c(names(pwas_metadata),'sex','beta','se','p','n')
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
        model$coefficients[2][[1]], coef(s)[2,2], coef(s)[2,4], length(model$residuals))
    },
    error=function(e){
      return(c(pwas_metadata[pwas_metadata$ukb_variable==trait,], rep(NA,5)))
    }
  )
  i <- i+1
}

pwas_hla_results <- merge(pwas_hla_results, pwas_results[pwas_results$ancestry==bw,
                                                       c('ukb_variable','beta','se','p')], by='ukb_variable')

pwas_hla_results$delta_beta <- (pwas_hla_results$beta.x/pwas_hla_results$beta.y)-1

# table_s6 <- pwas_hla_results[,c('display_name','ukb_variable','ukb_data_field','category','beta.x',
#                                'se.x','p.x','n','delta_beta')]

# write.table(table_s6, paste0(output_dir, '/table_s6.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
# write.table(pwas_na_results, paste0(output_dir, '/pwas_noAsthma_results.tsv'), 
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

seed <- 20220117

# pwas_dat <- read.csv(['path/to/ukb/data/pwas_data.csv'])


# Calculate proportion mediated (we reported the a*b/(c_std) values, as from MacKinnon 2007)
traits <- c('eosinophill_count0','fev1_fvc_ratio','age_hay_fever_diagnosed')
phenos <- c('asthma_children','asthma_adults','asthma_all')

med_df <- setNames(data.frame(matrix(ncol = 3, nrow = 3)), phenos)
rownames(med_df) <- traits

for (pheno in phenos){
  for (trait in traits){
    model_dat <- pwas_dat[pwas_dat$ancestry==bw,]
    if (pheno=='asthma_children') {
      test_dat <- model_dat[!(is.na(model_dat$asthma_children) & model_dat$asthma_all==1),]
      test_dat <- test_dat[!is.na(model_dat$asthma_adults) & model_dat$asthma_adults!=1,]
    } else if (pheno=='asthma_adults') {
      test_dat <- model_dat[!(is.na(model_dat$asthma_adults) & model_dat$asthma_all==1),]
      test_dat <- test_dat[!is.na(model_dat$asthma_children) & model_dat$asthma_children!=1,]
    } else {
      test_dat <- model_dat
    }
    model.y <- glm(paste0(pheno, '~', prs, covars), data=test_dat, family='binomial')
    model.m <- lm(paste0(trait,'~', prs , covars), data=test_dat)
    model.total <- glm(paste0(pheno,'~', trait, '+', prs, covars), data=test_dat, family='binomial')
    
    c <- coefficients(model.y)[[prs]]
    c_prime <- coefficients(model.total)[[prs]]
    b2 <- (coefficients(model.total)[[trait]])^2
    sigma2_mx <- var(model.m$residuals)
    norm_factor <- sqrt(1+((b2*sigma2_mx)/((pi^2)/3)))
    c_std <- c*norm_factor
    
    a <- coefficients(model.m)[[prs]]
    b <- coefficients(model.total)[[trait]]
    print(c(pheno, trait, round(a,6), round(c_std,6), round(a*b/((a*b)+c_prime),6), round(a*b/(c_std),6)))
    med_df[trait, pheno] <- round(a*b/(c_std),6)
  }
}


#####  FIGURE 8 (Mediation Bar Plot) #####
names(med_df) <- c("Childhood-onset\nasthma","Adult-onset\nasthma","All\nasthma")
cols <- c('cadetblue3','darkseagreen3','coral1')

pdf(file=paste0(output_dir, '/Fig_8.pdf'), width=6, height=5)
par(xpd=F, mar=c(3,4.5,2,2))
b <- barplot(height=as.matrix(med_df), ylab="Proportion mediated", beside=TRUE, ylim=c(-0.05,0.2), col=cols,
             cex.lab=0.8, cex.names=0.8, font=2, yaxt='n')
axis(2, cex.axis=0.8)                
grid(nx = NA, ny = NULL, lty = 2, col = "gray", lwd = 1)
abline(h=0,col='black')
barplot(height=as.matrix(med_df), beside=TRUE, ylim=c(-0.05,0.2), col=cols, add=T, xaxt='n',yaxt='n', ylab='')
text(x=b, y=c(med_df[,1], med_df[,2], med_df[,3]), pos=c(1,1,1,1,1,3,1,1,1), 
     label=round(c(med_df[,1], med_df[,2], med_df[,3]),2), cex=0.6, col='black')
legend('topright', legend=c('Eosinophil count', expression('FEV'[1]*'/FVC'), 'Age hay fever diagnosed'), cex=0.7,
       fill=cols, box.col='gray75', bty=1, inset=c(0.02,0.01), x.intersp=0.75, y.intersp=0.9, text.width=3.5)
dev.off()

# Estimate asthma mediation of BMI association, proportion mediated from Li 2007
bmi_med_df <- setNames(data.frame(matrix(ncol = 2, nrow = 3)), c('Proportion','P-value'))
rownames(bmi_med_df) <- phenos

for (pheno in phenos){
  model_dat <- pwas_dat[pwas_dat$ancestry==bw,]
  if (pheno=='asthma_children') {
    test_dat <- model_dat[!(is.na(model_dat$asthma_children) & model_dat$asthma_all==1),]
    test_dat <- test_dat[!is.na(model_dat$asthma_adults) & model_dat$asthma_adults!=1,]
  } else if (pheno=='asthma_adults') {
    test_dat <- model_dat[!(is.na(model_dat$asthma_adults) & model_dat$asthma_all==1),]
    test_dat <- test_dat[!is.na(model_dat$asthma_children) & model_dat$asthma_children!=1,]
  } else {
    test_dat <- model_dat
  }
  tm_sobel <- test_mediation(test_dat, prs, 'bmi_body_size0',pheno, 
                             test='sobel', covariates=strsplit(covars,'\\+')[[1]][-1], robust=F)
  
  model.y <- lm(paste0('bmi_body_size0 ~', prs, covars), data=test_dat)
  model.m <- glm(paste0(pheno, '~', prs, covars), data=test_dat, family='binomial')
  model.total <- lm(paste0('bmi_body_size0 ~', pheno, '+', prs, covars), data=test_dat)
  
  a_hat <- coefficients(model.m)[[prs]]
  b_hat <- coefficients(model.total)[[pheno]]
  b1 <- coefficients(model.y)[[prs]]
  delt_AL <- a_hat*b_hat*(exp(sum(coefficients(model.m)))/((1 + exp(sum(coefficients(model.m))))^2))
  
  a <- coefficients(model.m)[[prs]]
  b <- coefficients(model.total)[[pheno]]
  c_prime <- coefficients(model.total)[[prs]]
  c <- coefficients(model.y)[[prs]]
  
  prop_AL <- delt_AL/(delt_AL + c_prime)
  
  print(c(pheno, prop_AL, tm_sobel$p_value))
  bmi_med_df[pheno, ] <- c(prop_AL, tm_sobel$p_value)
}
