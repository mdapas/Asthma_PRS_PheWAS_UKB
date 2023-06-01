####  SET PARAMETERS  ##############################################################################
#
#  NOTE: even if running an individual section, these parameters should be saved first
#

set.seed(0123456789)

bw <- 'White'; prs <- 'zscore'; ukb_var <- 'ukb_variable'
coa <- 'asthma_children'; aoa <- 'asthma_adults'; asthma <- 'asthma_all'

# These paths would be user-specific
# pheno_dir <- [/path/to/directory/with/phenotype/data]
# output_dir <- [/path/to/output/directory]

pheno_dir <- '/Users/mdn578/Documents/Research/Ober/Asthma_PRS/Phenotypes/'
output_dir <- '/Users/mdn578/Documents/Research/Ober/Asthma_PRS/GBMI_Results/'

group_cols <- c('Neurologic & Behavioral'='royalblue', 'Ocular'='cornflowerblue',
                'Asthma & Allergic Disease'='indianred', 'Pulmonary'='orange',
                'Cardiovascular'='palevioletred', 'Hematologic &\nBlood Chemistry'='cadetblue',
                'Gastrointenstinal'='coral', 'Renal & Urologic'='seagreen',
                'Reproductive,\nEndocrine & Metabolic'='orchid4', 'Anthropometric'='darkseagreen',
                'Musculoskeletal & Skin'='goldenrod', 'Neoplasms'='plum4',
                'Other'='gold')

#####  FIGURE 1 (PRS Asthma R2)  #############################################################################

library(ggplot2)

predict_df <- read.delim(paste0(output_dir, '/asthma_prs_prediction_r2liability.txt'), header=T)

# compare r2 best GBMI and TAGC models (using r2 covariance from 10.1016/j.ajhg.2023.01.004)
for (i in c(5:11)) {predict_df[,i] <- as.numeric(predict_df[,i])}
predict_df$model <- paste(predict_df$Study, predict_df$LD,sep='_')
model_comp <- predict_df[predict_df$model %in% c('GBMI_1KG','TAGC_ME_UKB'),]
model_comp$model <- substr(model_comp$model,1,4)
model_comp$UKB_ancestry <- factor(gsub("White", "white", model_comp$UKB_ancestry), 
                                  levels=c('white','non-British white','African'))
model_comp$Asthma_pheno <- factor(model_comp$Asthma_pheno, levels=c(asthma,coa,aoa))


pdf(file=paste0(output_dir, '/Fig_1.pdf'), width=8, height=4)
ggplot(model_comp, 
       aes(factor(model), y=r2, ymin=r2_2.5, ymax=r2_97.5)) + 
  geom_linerange(aes(colour=factor(Asthma_pheno)), position=position_dodge(width=0.5), linewidth=1.2) + 
  geom_point(aes(fill=factor(Asthma_pheno)), shape=21, position=position_dodge(width=0.5), size=1.8, stroke=1.1) +
  theme_classic() + ylab(bquote({italic(R)^2}[liability])) +
  scale_y_continuous(limits=c(0.0,0.08)) + facet_wrap(~UKB_ancestry) + 
  scale_color_manual(values=c('#9F10F8','#22D378','#FB9700'), labels=c('Asthma (All)', 'COA','AOA')) + 
  scale_fill_manual(values=c('#9F10F8','#22D378','#FB9700'), labels=c('Asthma (All)', 'COA','AOA')) + 
  theme(axis.title=element_text(size=11, face='bold'), axis.title.x=element_blank(),
        legend.title=element_blank(), legend.text=element_text(size=8))
dev.off()


#####  FIGURE 2 (PRS Asthma Prediction)  #####################################################################

library(plotrix)
library(scales)
library(pROC)

trait_dat <-  read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'))
load(paste0(output_dir, '/quantile_dat.RData'))
plot_dat <- trait_dat[trait_dat$ancestry==bw,]

### Quantiles - All asthma (Figure 2A)
y_range <- c(0.142,3.7)
pdf(file=paste0(output_dir, '/Fig_2A.pdf'), width=6, height=5)
par(xpd=F, mar=c(4.5,4.5,2,2), lheight = .3)
plot(x=1:11, ylim=y_range, xlab='', ylab='', yaxt='n', xaxt = "n", type='n')
abline(h = seq(0,3.5,.5), lty = 2, col = 'grey80')
abline(h = 1, lty = 1, col = 'grey')
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_asthma$OR, li=quantile_dat_asthma[,c('95CI_L')], ui=quantile_dat_asthma[,c('95CI_U')], 
       xlab='\nPRS Quantile (%)', ylab='Odds Ratio (OR)', xaxt = "n", pt.bg='#9F10F8', sfrac=0.001, lwd=2.5, cex.axis=.8,
       scol='#9F10F8', cex=1, pch=21, ylim=y_range, yaxt='n', cex.lab=0.9)
axis(1, at=1:11, labels = FALSE)
text(1:11, par("usr")[3] - 0.11, labels = quantile_dat_COA$Quantile, srt = 35, pos = 1, xpd = TRUE, cex=0.8)

axis(2,seq(0,3.5,.5), labels = c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'), xpd = TRUE, cex.axis=0.8, las=2)
dev.off()


### DENSITY PLOT (Fig 2B)

pdf(file=paste0(output_dir, '/Fig_2B.pdf'), width=6, height=5)
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


### Quantiles -  COA & AOA (Figure 2C)
y_range <- c(0.21,5.5)
pdf(file=paste0(output_dir, '/Fig_2C.pdf'), width=6, height=5)
par(xpd=F, mar=c(4.5,4.5,2,2), lheight = .3)
plot(x=1:11, ylim=y_range, xlab='', ylab='', yaxt='n', xaxt = "n", type='n')
abline(h = seq(0,6,.5), lty = 2, col = 'grey80')
abline(h = 1, lty = 1, col = 'grey')
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_COA$OR, li=quantile_dat_COA[,c('95CI_L')], ui=quantile_dat_COA[,c('95CI_U')], 
       xlab='\nPRS Quantile (%)', ylab='Odds Ratio (OR)', xaxt = "n", pt.bg='#22D378', sfrac=0.001, lwd=2.5, cex.axis=.8,
       scol='#22D378', cex=0.8, pch=21, ylim=y_range, yaxt='n', cex.lab=0.9)
axis(1, at=1:11, labels = FALSE)
text(1:11, par("usr")[3] - 0.18, labels = quantile_dat_COA$Quantile, srt = 35, pos = 1, xpd = TRUE, cex=0.8)
axis(2,seq(0,5.5,.5), labels = c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5'), xpd = TRUE, cex.axis=0.8, las=2)
par(new=TRUE)
plotCI(x=1:11, y=quantile_dat_AOA$OR, li=quantile_dat_AOA[,c('95CI_L')], ui=quantile_dat_AOA[,c('95CI_U')], 
       xlab='', ylab='', xaxt = "n", pt.bg='#FB9700', sfrac=0.001, lwd=2.5, 
       scol='#FB9700', cex=0.8, pch=21, ylim=y_range, yaxt='n')
points(x=6,y=1, pch=21, lwd=2.5, bg=alpha('#22D378',.4), cex=0.8)
# Key labels
legend("topleft", inset=0.01, legend=c("COA", "AOA"), lwd=2, col=c('#22D378','#FB9700'),
       pch=c(16,16), bty='n', cex=0.8, pt.cex = 1.1)
legend("topleft", inset=0.01, col='black', legend=c("", ""),  lty=F, lwd=1.5,
       pch=c(21,21), bty='n', cex=0.8, pt.cex = 1.1)
dev.off()


### AUC PLOT (Fig 2D)
pdf(file=paste0(output_dir, '/Fig_2D.pdf'), width=6, height=5)
par(mar = c(4, 4, 4, 4)+.1)
plot.roc(plot_dat$asthma_all, plot_dat$zscore, main="", xlim=c(0.963,0.037), ylim=c(0.037,0.963),
         col='purple2', asp=F, xaxt='n', yaxt='n')
legend("bottomright", title="AUC", inset=0.02, legend=c("COA: 0.665", "AOA: 0.600", "All:    0.623"), 
       col=c('seagreen3','orange2','purple2'),
       pch=c(15,15,15), bty='n', cex=0.85, pt.cex = 1.15)
plot.roc(plot_dat$asthma_children, plot_dat$zscore,col='seagreen3', add=TRUE)
plot.roc(plot_dat$asthma_adults, plot_dat$zscore ,col='orange2', add=TRUE)
axis(2, seq(0,1,.2), labels = c('0.0','0.2','0.4','0.6','0.8','1.0'), xpd = TRUE, cex.axis=0.85, las=2)
axis(1, seq(1,0,-.2), labels = NA)
text(seq(1,0,-.2), par("usr")[3] - 0.03, labels = c('1.0','0.8','0.6','0.4','0.2','0.0'), pos = 1, xpd = TRUE, cex=0.85)
dev.off()


#####  FIGURE 3 (PRS by Sex and Age of Onset)  ###############################################################

library(emmeans)
library(ggplot2)
library(ggnewscale)
library(scales)

# trait_dat <-  read.csv(paste0(pheno_dir, '/ukb.trait.prs.data.csv'))
# plot_dat <- trait_dat[trait_dat$ancestry==bw,]


### PRS VS RISK LINEAR SCALE (Fig 3A)
range <- seq(from=-7, to=7, by=0.05)

fit_coa <- glm(asthma_children ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + zscore*sex, 
               data=plot_dat, family='binomial')
c_pval <- summary(fit_coa)$coefficients[nrow(summary(fit_coa)$coefficient),4]
c_emm <- emmeans(fit_coa,~zscore*sex,type="response",at= list(zscore=range, sex=c('M','F')))
c_emm2 <- as.data.frame(c_emm)
c_emm2$lower <- confint(c_emm,adjust='none',level=0.95)$asymp.LCL
c_emm2$upper <- confint(c_emm,adjust='none',level=0.95)$asymp.UCL

fit_aoa <- glm(asthma_adults ~ pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + sex*zscore,
               data=plot_dat, family='binomial')
a_pval <- summary(fit_aoa)$coefficients[nrow(summary(fit_aoa)$coefficient),4]
a_emm <- emmeans(fit_aoa,~zscore*sex,type="response",at= list(zscore=range,sex=c('M','F')))
a_emm2 <- as.data.frame(a_emm)
a_emm2$lower<-confint(a_emm,adjust='none',level=0.95)$asymp.LCL
a_emm2$upper<-confint(a_emm,adjust='none',level=0.95)$asymp.UCL

pdf(file=paste0(output_dir, '/Fig_3A.pdf'), width=6, height=5)
ggplot(a_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.9) + 
  geom_line(data=a_emm2, aes(color=as.factor(sex))) + labs(fill="Adult onset") + 
  scale_fill_manual(labels=c("Male","Female"), values = c('#BFA373','#FDC372')) + 
  scale_colour_manual(values=c('#934900','#D16E00'),guide='none') +
  coord_cartesian(ylim = c(0.017, 0.37), xlim = c(-4.475,4.475)) +
  scale_y_continuous(name ="Risk", limits=c(0,.4), breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35), 
                     labels = c('0%','5%','10%','15%','20%','25%','30%','35%')) +
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

### PRS VS COA RISK LOG SCALE (FIG 3B)
pdf(file=paste0(output_dir, '/Fig_3B.pdf'), width=5, height=4)
ggplot(c_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.75) + 
  geom_line(data=c_emm2, aes(color=as.factor(sex))) + labs(fill="Childhood onset") +
  scale_fill_manual(labels=c("Male","Female"), values = c('#409348','seagreen2')) +
  scale_colour_manual(values=c('#096000','seagreen4'),guide='none') +
  coord_cartesian(xlim = c(-3,3)) + theme_bw() +
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

### PRS VS AOA RISK LOG SCALE (FIG 3C)
pdf(file=paste0(output_dir, '/Fig_3C.pdf'), width=5, height=4)
ggplot(a_emm2, aes(x=zscore, y=prob,fill=as.factor(sex),)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.9) + 
  geom_line(data=a_emm2, aes(color=as.factor(sex))) + labs(fill="Adult onset") + 
  scale_fill_manual(labels=c("Male","Female"), values = c('#BFA373','#FDC372')) + 
  scale_colour_manual(values=c('#934900','#D16E00'),guide='none') +
  coord_cartesian(xlim = c(-3,3)) +
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


#####  FIGURE 4 (PheWAS Results Plot)  #######################################################################
#
#  Create a horizontal 'Manhattan' plot for the PheWAS.
#    NOTE: Generating a vector image of the plot may be complicated by text-encoding issues,
#          and additional manipulation was required to get all the labels readable
#

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggtext)

# pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))

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
      n_i <- n_i+14
      print(fig_data[i,'category'])
    }
    else if (fig_data[i,'p'] < 0.05/n_bonf) {
      n_i <- n_i+5
    }
  }
  fig_data[i,'Pos'] <- n_i
  n_i <- n_i+2
}

fig_data$Pos <- n_i-fig_data$Pos

# Prepare Y-axis labels
fig_data[fig_data$category=='Reproductive, Endocrine & Metabolic','category'] <- 'Reproductive,\nEndocrine & Metabolic'
fig_data[fig_data$category=='Hematologic & Blood Chemistry','category'] <- 'Hematologic &\nBlood Chemistry'

fig_tibble <- tibble(fig_data)
y_axis <- fig_tibble %>% group_by(category) %>% summarize(center=(max(Pos) + min(Pos) ) / 2 )
rects <- fig_tibble %>% group_by(category) %>% summarize(mincat=min(Pos)-8,maxcat=max(Pos)+8)
totals <- fig_tibble %>% group_by(category) %>% count()

# Determine whether label will be added
fig_data$label <- T
#fig_data[fig_data$p>=0.05/n_bonf,'label'] <- F
#fig_data[fig_data$bin==0 & abs(fig_data$beta)<0.01, 'label'] <- F
#fig_data[fig_data$bin==1 & abs(fig_data$beta)<log(1.05), 'label'] <- F
fig_data[fig_data$p>1e-10, 'label'] <- F

fig_data$low_pvals <- NA

low_pval_traits <- fig_data[fig_data$p==0, 'display_name']

fig_data[fig_data$display_name %in% low_pval_traits, 'low_pvals'] <- 'p<4.9e-324'
fig_data[fig_data$display_name %in% low_pval_traits,'p'] <- 5e-324

#fig_data[fig_data$display_name=='Tx: fluticasone/salmeterol','low_pvals'] <- 'p=8.9e-243'
#fig_data[fig_data$display_name=='FEV1','low_pvals'] <- 'p=9.5e-234'

# Flip betas if necessary
fig_data[fig_data$invert==1,'beta'] <- fig_data[fig_data$invert==1,'beta']*-1

fig_data[fig_data$display_name=='Adult-onset asthma','Pos'] <- n_i - 19
#fig_data[fig_data$display_name=='Age asthma diagnosed','Pos'] <- fig_data[fig_data$display_name=='Age asthma diagnosed','Pos'] + 1
fig_data[fig_data$display_name=='Asthma','Pos'] <- n_i + 3
fig_data[fig_data$display_name=='Childhood-onset asthma','Pos'] <- n_i - 30
fig_data[fig_data$display_name=='Hay fever, allergic rhinitis, eczema','Pos'] <- n_i - 41
fig_data[fig_data$display_name=='Tx: albuterol','Pos'] <- n_i - 52
fig_data[fig_data$display_name=='Wheezing','Pos'] <- fig_data[fig_data$display_name=='Wheezing','Pos'] + 1


y_labels <- sub('\n','<br>',
                str_c("<br>**",y_axis$category,"**<br><span style='font-size:9.5pt;color:dimgray'>",
                      paste0("^((n=",totals$n,")^)</span>")))

xmax=250

cairo_pdf(filename=paste0(output_dir, '/Fig_3.pdf'), width=11, height=7.5)  # requires xquartz/encodings
ggplot(fig_data, aes(x=-log10(p),y=Pos)) + 
  
  # add non-significant traits
  geom_point(data=subset(fig_data, p>0.05/n_bonf), color="grey60", size=3.1,pch=19) +
  geom_point(data=subset(fig_data, p>0.05/n_bonf), aes(bg=as.factor(category)), alpha=.3,size=3.1,pch=21) +
  scale_fill_manual(values=group_cols) +
  
  # custom axes:
  scale_y_continuous(label=y_labels, breaks= y_axis$center) +
  scale_x_continuous(expand = c(0, 0), limits=c(-1,xmax)) + # remove space between plot area and x axis
  geom_vline(xintercept=-log10(0.05/n_bonf),linetype='dashed',color='grey30',size=0.7) +
  coord_cartesian(ylim=c(28,n_i-28), clip="off") +
  # Add label using ggrepel to avoid overlapping
  geom_text_repel(data=subset(fig_data, label==T), aes(label=display_name), size=2.6,
                  max.overlaps=6, point.padding=unit(0.4,'lines'),nudge_x=4,
                  box.padding=unit(0.35,'lines'),segment.color='grey',min.segment.length=0.2) +
  
  # Add highlighted points
  geom_point(data=subset(fig_data, p<0.05/n_bonf & beta>0), aes(bg=as.factor(category)), alpha=0.85, size=2.8,pch=24) +
  geom_point(data=subset(fig_data, p<0.05/n_bonf & beta<0), aes(bg=as.factor(category)), alpha=0.85, size=2.8,pch=25) +
  geom_point(data=subset(fig_data, label==T), color='black', size=1,pch=16) +
  
  # Add super-significant labels
  geom_text(data=subset(fig_data,p<1e-250), family = 'Arial Unicode MS', hjust=1,
            aes(label=paste0(display_name,', ',low_pvals,"\u2192")), x = xmax, size=2.6) +
  
  # add color blocks
  annotate("rect", xmin=-1,xmax=-log10(0.05/n_bonf), ymin=min(rects[,'mincat'][[1]])-32,
           ymax=max(rects$maxcat)+32,
           fill='grey15', alpha=0.2) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Renal & Urologic','mincat'][[1]]-32,
           ymax=rects[rects$category=='Renal & Urologic','maxcat'][[1]],
           fill=group_cols[['Renal & Urologic']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Neurologic & Behavioral','mincat'][[1]],
           ymax=rects[rects$category=='Neurologic & Behavioral','maxcat'][[1]],
           fill=group_cols[['Neurologic & Behavioral']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Neoplasms','mincat'][[1]],
           ymax=rects[rects$category=='Neoplasms','maxcat'][[1]],
           fill=group_cols[['Neoplasms']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Ocular','mincat'][[1]],
           ymax=rects[rects$category=='Ocular','maxcat'][[1]],
           fill=group_cols[['Ocular']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Cardiovascular','mincat'][[1]],
           ymax=rects[rects$category=='Cardiovascular','maxcat'][[1]],
           fill=group_cols[['Cardiovascular']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Musculoskeletal & Skin','mincat'][[1]],
           ymax=rects[rects$category=='Musculoskeletal & Skin','maxcat'][[1]],
           fill=group_cols[['Musculoskeletal & Skin']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Other','mincat'][[1]],
           ymax=rects[rects$category=='Other','maxcat'][[1]],
           fill=group_cols[['Other']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Anthropometric','mincat'][[1]],
           ymax=rects[rects$category=='Anthropometric','maxcat'][[1]],
           fill=group_cols[['Anthropometric']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Reproductive,\nEndocrine & Metabolic','mincat'][[1]],
           ymax=rects[rects$category=='Reproductive,\nEndocrine & Metabolic','maxcat'][[1]],
           fill=group_cols[['Reproductive,\nEndocrine & Metabolic']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Gastrointenstinal','mincat'][[1]],
           ymax=rects[rects$category=='Gastrointenstinal','maxcat'][[1]],
           fill=group_cols[['Gastrointenstinal']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Pulmonary','mincat'][[1]],
           ymax=rects[rects$category=='Pulmonary','maxcat'][[1]],
           fill=group_cols[['Pulmonary']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Hematologic &\nBlood Chemistry','mincat'][[1]],
           ymax=rects[rects$category=='Hematologic &\nBlood Chemistry','maxcat'][[1]],
           fill=group_cols[['Hematologic &\nBlood Chemistry']], alpha=0.15) +
  annotate("rect", xmin=-1,xmax=xmax, ymin=rects[rects$category=='Asthma & Allergic Disease','mincat'][[1]],
           ymax=max(rects$maxcat)+32,
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


#####  FIGURE 5 (PheWAS association scatter plot GBMI vs TAGC) ###############################################

# pwas_results_tagc <- read.delim(paste0(output_dir, '/phewas_results_tagc.txt'))
# pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))
pwas_metadata <- read.table(paste0(pheno_dir,'/pwas_metadata.tsv'), header=T, sep='\t', quote="")
quant_traits <- pwas_metadata[pwas_metadata$bin==0,ukb_var]
bin_traits <- pwas_metadata[pwas_metadata$bin==1,ukb_var]

load(paste0(output_dir,'phewas_comp.Rdata'))


# BINARY TRAITS (5a)
r2 <- round(cor(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits), c('beta.x','beta.y')]))[2]^2,2)
sig_bin_traits <- scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & ((scatter_dat$p.x<(0.05/371))|(scatter_dat$p.y<(0.05/371))), ukb_var]

range=c(-0.3, 1)
pdf(file=paste0(output_dir, '/Fig_5a.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1, type='n', xlim=range, ylim=range, xlab='TAGC PRS (OR)', 
     ylab='GBMI PRS (OR)', cex.axis=0.8, cex.lab=0.9, xaxt='n', yaxt='n')
axis(1, c(0,log2(1.5),log2(2)), labels = c('1.0','1.5','2'), xpd = TRUE, cex.axis=0.85)
axis(2, c(0,log2(1.5),log2(2)), labels = c('1.0','1.5','2'), xpd = TRUE, cex.axis=0.85, las=2)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3); 
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) &
                            (!scatter_dat$ukb_variable %in% sig_bin_traits), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
#mtext('r=0.87',side=3,adj=0.05,padj=2)
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)
# significant in both datasets
#arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)
arrows(log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y']-
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.y'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y']+
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.y'])), length=0.01, angle=90, code=3, col='grey50')
arrows(log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x']-
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x']+
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y'])), length=0.01, angle=90, code=3, col='grey50')

# significant in both datasets
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & (scatter_dat$p.x<(0.05/371)) & 
                            (scatter_dat$p.y<(0.05/371)), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
     xlim=range, ylim=range,xlab='', ylab='',pch=21, bg='forestgreen', cex=0.8)
#significant in GBMI but not TAGC
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & (scatter_dat$p.y<(0.05/371)) & 
                            (scatter_dat$p.x>(0.05/371)), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=21, bg='deepskyblue', cex=0.8)
#significant in TAGC but not GBMI
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & (scatter_dat$p.y>(0.05/371)) & 
                            (scatter_dat$p.x<(0.05/371)), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=21, bg='goldenrod2', cex=0.8)
#par(new=T)
#plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & (scatter_dat$sig==1), c('beta.x','beta.y')])), 
#     xaxt='n',yaxt='n', 
#     xlim=range, ylim=range, xlab='', ylab='', pch=1, col='red', cex=0.8, lwd=1)
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% c(max_tagc_bin, max_gbmi_bin, 'asthma_icd10')), c('beta.x','beta.y')])), 
     xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='', pch=1, col='red', cex=0.8, lwd=1.5)
text(x=log2(exp(scatter_dat[scatter_dat$ukb_variable== max_gbmi_bin, 'beta.x']))-.03, 
     y=log2(exp(scatter_dat[scatter_dat$ukb_variable==max_gbmi_bin, 'beta.y']))+.06,
     label='Breathing\nproblems', cex=0.4, col='black', adj=1)
text(x=log2(exp(scatter_dat[scatter_dat$ukb_variable=='asthma_icd10', 'beta.x']))-.03, 
     y=log2(exp(scatter_dat[scatter_dat$ukb_variable=='asthma_icd10', 'beta.y']))+.06,
     label='Asthma\n(ICD10)', cex=0.4, col='black', adj=1)
text(x=log2(exp(scatter_dat[scatter_dat$ukb_variable== max_tagc_bin, 'beta.x']))+.04, 
     y=log2(exp(scatter_dat[scatter_dat$ukb_variable==max_tagc_bin, 'beta.y']))-.06,
     label='Celiac\ndisease', cex=0.4, col='black', adj=0)

legend('bottomright',legend=c('TAGC & GBMI','GBMI only','TAGC only'), bty='n', cex=0.7,
       title='Significant traits', title.cex=0.7, title.font=2, title.adj=-0.02,
       pch=21,pt.bg=c('forestgreen','deepskyblue','goldenrod2'))
dev.off()


# QUANT TRAITS (5b)
cor(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & (scatter_dat$p.x<(0.05/371)) & 
                  (scatter_dat$p.y<(0.05/371)), c('beta.x','beta.y')]) # 0.92

r2 <- round(cor(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits), c('beta.x','beta.y')])[2]^2,2)
sig_quant_traits <- scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & ((scatter_dat$p.x<(0.05/371))|(scatter_dat$p.y<(0.05/371))), ukb_var]

range=c(-0.14,0.14)
pdf(file=paste0(output_dir, '/Fig_5b.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1,type='n',xlim=range, ylim=range, xlab=expression(paste('TAGC PRS (',beta,')')), 
     ylab=expression(paste('GBMI PRS (',beta,')')), cex.axis=0.7, cex.lab=0.8)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3)
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) &
                   (!scatter_dat$ukb_variable %in% sig_quant_traits), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
#mtext('r=0.87',side=3,adj=0.05,padj=2)
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)
# significant in both datasets
#arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)
arrows(scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y']-
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.y'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y']+
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.y'], length=0.01, angle=90, code=3, col='grey50')
arrows(scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x']-
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x']+
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y'], length=0.01, angle=90, code=3, col='grey50')
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & (scatter_dat$p.x<(0.05/371)) & 
                   (scatter_dat$p.y<(0.05/371)), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=21, bg='forestgreen', cex=0.8)
#significant in GBMI but not TAGC
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & (scatter_dat$p.y<(0.05/371)) & 
                   (scatter_dat$p.x>(0.05/371)), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=21, bg='deepskyblue', cex=0.8)
#significant in TAGC but not GBMI
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & (scatter_dat$p.y>(0.05/371)) & 
                   (scatter_dat$p.x<(0.05/371)), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=21, bg='goldenrod2', cex=0.8)
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% c(max_tagc_quant, max_gbmi_quant)), c('beta.x','beta.y')], 
     xaxt='n',yaxt='n', xlim=range, ylim=range, xlab='', ylab='', pch=1, col='red', cex=0.8, lwd=1.5)
text(x=scatter_dat[scatter_dat$ukb_variable== max_tagc_quant, 'beta.x']-.004, 
     y=scatter_dat[scatter_dat$ukb_variable==max_tagc_quant, 'beta.y']+.011,
     label='Age asthma\ndiagnosed', cex=0.4, col='black', adj=1)
text(x=scatter_dat[scatter_dat$ukb_variable== max_gbmi_quant, 'beta.x']+.006, 
     y=scatter_dat[scatter_dat$ukb_variable==max_gbmi_quant, 'beta.y']-.011,
     label='FEV1\n(% predicted)', cex=0.4, col='black', adj=0)

legend('bottomright',legend=c('TAGC & GBMI','GBMI only','TAGC only'), bty='n', cex=0.7,
       title='Significant traits', title.cex=0.7, title.font=2, title.adj=-0.02, 
       pch=21,pt.bg=c('forestgreen','deepskyblue','goldenrod2'))
dev.off()


# BINARY TRAITS, with categorical differences (5c)
r2 <- round(cor(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits), c('beta.x','beta.y')]))[2]^2,2)
sig_bin_traits <- scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & ((scatter_dat$p.x<(0.05/371))|(scatter_dat$p.y<(0.05/371))), ukb_var]

range=c(-0.3, 1)
pdf(file=paste0(output_dir, '/Fig_5c.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1, type='n', xlim=range, ylim=range, xlab='TAGC PRS (OR)', 
     ylab='GBMI PRS (OR)', cex.axis=0.8, cex.lab=0.9, xaxt='n', yaxt='n')
axis(1, c(0,log2(1.5),log2(2)), labels = c('1.0','1.5','2'), xpd = TRUE, cex.axis=0.85)
axis(2, c(0,log2(1.5),log2(2)), labels = c('1.0','1.5','2'), xpd = TRUE, cex.axis=0.85, las=2)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3); 
par(new=T)
plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) &
                            (!scatter_dat$ukb_variable %in% sig_bin_traits), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
#mtext('r=0.87',side=3,adj=0.05,padj=2)
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)
# significant in both datasets
#arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)
arrows(log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y']-
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.y'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y']+
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.y'])), length=0.01, angle=90, code=3, col='grey50')
arrows(log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x']-
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.x']+
                  1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'se.x'])),
       log2(exp(scatter_dat[scatter_dat$ukb_variable %in% sig_bin_traits, 'beta.y'])), length=0.01, angle=90, code=3, col='grey50')

i_cats <- c()
for (i in 1:length(group_cols)) {
  par(new=T)
  cat=sub('\n',' ',names(group_cols)[i])
  print(cat)
  temp_dat <- scatter_dat[scatter_dat$category==cat & scatter_dat$ukb_variable %in% bin_traits,]
  if (nrow(temp_dat) > 2) {
    model <- summary(lm(I(exp(beta.y))~0+exp(beta.x), data=temp_dat))
    p <- 2*pt((1-coef(model)[1])/coef(model)[2], df=nrow(temp_dat)-2)
  } else {
    p <- 1
  }
  if (p<(0.05/13)) {
    print(p)
    i_cats <- c(i_cats, i)
    abline(0,log2(exp(coef(model)[1])), col=group_cols[i], lty=5, lwd=2)
  } else {
    plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & 
                                (scatter_dat$category==cat), c('beta.x','beta.y')])), xaxt='n',yaxt='n', 
         xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
  }
  
}
for (i in i_cats) {
  par(new=T)
  plot(log2(exp(scatter_dat[(scatter_dat$ukb_variable %in% bin_traits) & 
                              (scatter_dat$category==sub('\n',' ',names(group_cols)[i])), c('beta.x','beta.y')])), 
       xaxt='n',yaxt='n', xlim=range, ylim=range, xlab='', ylab='',pch=21, bg=group_cols[i], cex=0.8)
}

legend('bottomright',legend=names(group_cols[i_cats]), bty='n', cex=0.7,
       title='Significant categories', title.cex=0.7, title.font=2, title.adj=-0.02,
       pch=21,pt.bg=group_cols[i_cats])

dev.off()



# QUANT TRAITS, color by category (5d)
cor(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & (scatter_dat$p.x<(0.05/371)) & 
                  (scatter_dat$p.y<(0.05/371)), c('beta.x','beta.y')]) # 0.92

r2 <- round(cor(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits), c('beta.x','beta.y')])[2]^2,2)
sig_quant_traits <- scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & ((scatter_dat$p.x<(0.05/371))|(scatter_dat$p.y<(0.05/371))), ukb_var]

range=c(-0.14,0.14)
pdf(file=paste0(output_dir, '/Fig_5d.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1,type='n',xlim=range, ylim=range, xlab=expression(paste('TAGC PRS (',beta,')')), 
     ylab=expression(paste('GBMI PRS (',beta,')')), cex.axis=0.7, cex.lab=0.8)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3)
par(new=T)
plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) &
                   (!scatter_dat$ukb_variable %in% sig_quant_traits), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
     xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
#mtext('r=0.87',side=3,adj=0.05,padj=2)
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)
# significant in both datasets
#arrows(x, avg-sdev, x, avg+sdev, length=0.02, angle=90, code=3)
arrows(scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y']-
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.y'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y']+
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.y'], length=0.01, angle=90, code=3, col='grey50')
arrows(scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x']-
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.x']+
         1.96*scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'se.x'],
       scatter_dat[scatter_dat$ukb_variable %in% sig_quant_traits, 'beta.y'], length=0.01, angle=90, code=3, col='grey50')
i_cats <- c()
for (i in 1:length(group_cols)) {
  par(new=T)
  cat=sub('\n',' ',names(group_cols)[i])
  print(cat)
  temp_dat <- scatter_dat[scatter_dat$category==cat & scatter_dat$ukb_variable %in% quant_traits,]
  if (nrow(temp_dat) > 2) {
    model <- summary(lm(I(beta.y)~0+beta.x, data=temp_dat))
    p <- 2*pt((1-coef(model)[1])/coef(model)[2], df=nrow(temp_dat)-2)
  } else {
    p <- 1
  }
  if (p<(0.05/13)) {
    print(p)
    i_cats <- c(i_cats, i)
    abline(0,coef(model)[1], col=group_cols[i], lty=5, lwd=2)
  } else {
    plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & 
                       (scatter_dat$category==cat), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
         xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.8)
  }
  
}
for (i in i_cats) {
  par(new=T)
  plot(scatter_dat[(scatter_dat$ukb_variable %in% quant_traits) & 
                     (scatter_dat$category==sub('\n',' ',names(group_cols)[i])), c('beta.x','beta.y')], xaxt='n',yaxt='n', 
       xlim=range, ylim=range, xlab='', ylab='',pch=21, bg=group_cols[i], cex=0.8)
}
legend('bottomright',legend=names(group_cols[i_cats]), bty='n', cex=0.7,
       title='Significant categories', title.cex=0.7, title.font=2, title.adj=-0.02,
       pch=21,pt.bg=group_cols[i_cats])
dev.off()


#####  FIGURE 6 (Scatter plot of asthma vs. non-asthmatic associations) #####################################

library(tidyverse)
library(forestplot)
library(grid)

#pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))
#pwas_na_results <- read.delim(paste0(output_dir, '/pwas_noAsthma_results.tsv'))

scatter_dat <- pwas_na_results[!is.na(pwas_na_results$beta.na) & pwas_na_results$p.a<(0.05/371),]
scatter_dat$beta.a <- ifelse(scatter_dat$invert==1, scatter_dat$beta.a*-1, scatter_dat$beta.a)
scatter_dat$beta.na <- ifelse(scatter_dat$invert==1, scatter_dat$beta.na*-1, scatter_dat$beta.na)

scatter_dat$beta_diff <- ifelse(scatter_dat$p.a>(0.05/371) | scatter_dat$p.na>(0.05/371), NA,
                                ifelse(scatter_dat$beta.a*scatter_dat$beta.na>0,
                                       abs(scatter_dat$beta.na)-abs(scatter_dat$beta.a),
                                       abs(scatter_dat$beta.na)+abs(scatter_dat$beta.a)))

scatter_quant <- scatter_dat[scatter_dat$ukb_variable %in% quant_traits,]
scatter_bin <- scatter_dat[scatter_dat$ukb_variable %in% bin_traits,]


# BINARY TRAITS, with categorical differences (6a)
scatter_bin <- scatter_bin[scatter_bin$se.na <1, ]  # removing traits with exceptionally wide CIs
r2 <- round(cor(exp(scatter_bin[,c('beta.a','beta.na')]))[2]^2,2)
sig_bin_traits <- scatter_bin[scatter_bin$sig==1, ukb_var]

xrange <- c(-0.26,1)
yrange <- c(-0.34,1.4)
pdf(file=paste0(output_dir, '/Fig_6a.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1, type='n', xlim=xrange, ylim=yrange, xlab='All (OR)', 
     ylab='No Asthma (OR)', cex.axis=0.8, cex.lab=0.9, xaxt='n', yaxt='n')
axis(1, c(0,log2(1.5),log2(2)), labels = c('1.0','1.5','2'), xpd = TRUE, cex.axis=0.85)
axis(2, c(0,log2(1.5),log2(2), log2(2.5)), labels = c('1.0','1.5','2','2.5'), xpd = TRUE, cex.axis=0.85, las=2)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3); 
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)
# significant in both datasets
arrows(log2(exp(scatter_bin[, 'beta.a'])), log2(exp(scatter_bin[,'beta.na']-1.96*scatter_bin[,'se.na'])),
       log2(exp(scatter_bin[,'beta.a'])), log2(exp(scatter_bin[,'beta.na']+1.96*scatter_bin[,'se.na'])),
       length=0.01, angle=90, code=3, col='grey50')
arrows(log2(exp(scatter_bin[,'beta.a']-1.96*scatter_bin[,'se.a'])), log2(exp(scatter_bin[,'beta.na'])),
       log2(exp(scatter_bin[,'beta.a']+1.96*scatter_bin[,'se.a'])), log2(exp(scatter_bin[,'beta.na'])),
       length=0.01, angle=90, code=3, col='grey50')
par(new=T)
plot(log2(exp(scatter_bin[!scatter_bin$ukb_variable %in% sig_bin_traits, c('beta.a','beta.na')])), 
     xaxt='n',yaxt='n', xlim=xrange, ylim=yrange, xlab='', ylab='',pch=16, col='grey70', cex=0.6)

cat_indx <- c()
for (cat in names(table(scatter_bin[scatter_bin$ukb_variable %in% sig_bin_traits,'category']))) {
  cat_match <- substr(cat,1,5)
  par(new=T)
  plot(log2(exp(scatter_bin[(scatter_bin$ukb_variable %in% sig_bin_traits) & (scatter_bin$category==cat), 
                            c('beta.a','beta.na')])), xaxt='n',yaxt='n', xlim=xrange, ylim=yrange, xlab='', ylab='',
       pch=21, bg=group_cols[substr(names(group_cols),1,5)==cat_match], cex=0.8)
  cat_indx <- c(cat_indx, cat_match)
}

names <- c()
for (trait in sig_bin_traits) {
  row <- scatter_dat[scatter_dat$ukb_variable==trait,]
  if (trait %in% c('blood_clot_hayfever','blood_clot_none','hayfever_doctor_diagnosed0')) {
    row['display_name'] <- 'Hay fever'
  }
  if (!row['display_name'] %in% names) {
    if (row$beta.a>0) {
      text(x=log2(exp(row['beta.a']))+.1, y=log2(exp(row['beta.na']))-.03, label=row['display_name'], cex=0.5, col='black', adj=1)
    } else {
      text(x=log2(exp(row['beta.a']))-.04, y=log2(exp(row['beta.na']))+.067, label=row['display_name'], cex=0.5, col='black', adj=1)
    }
    names <- c(names, row['display_name'])
  }
}

legend('bottomright',legend=names(group_cols[substr(names(group_cols),1,5) %in% cat_indx]), bty='n', cex=0.7,
       title.cex=0.7, title.font=2, title.adj=-0.02,
       pch=21,pt.bg=group_cols[substr(names(group_cols),1,5) %in% cat_indx])
dev.off()


# QUANT TRAITS, color by category (6b)
cor(scatter_quant[,  c('beta.a','beta.na')]) # 0.92

r2 <- round(cor(scatter_quant[, c('beta.a','beta.na')])[2]^2,2)
sig_quant_traits <- scatter_quant[(scatter_quant$sig==1), ukb_var]

range=c(-0.14,0.14)
pdf(file=paste0(output_dir, '/Fig_6b.pdf'), width=5, height=5)
par(xpd=F, mgp=c(2.5,0.75,0), mar=c(4.1,4.5,2,2), lheight = .3)
plot(1,type='n',xlim=range, ylim=range, xlab=expression(paste('All (',beta,')')), 
     ylab=expression(paste('No asthma (',beta,')')), cex.axis=0.7, cex.lab=0.8)
abline(0,1, col='pink', lty=2); abline(h=0, lty=3); abline(v=0, lty=3)
mtext(parse(text=paste("r^2==",r2)), side=3, adj=0.05, padj=2, cex=0.8)

arrows(scatter_quant$beta.a, scatter_quant$beta.na - 1.96*scatter_quant$se.na,
       scatter_quant$beta.a, scatter_quant$beta.na + 1.96*scatter_quant$se.na, 
       length=0.01, angle=90, code=3, col='grey50')
arrows(scatter_quant$beta.a - 1.96*scatter_quant$se.a, scatter_quant$beta.na, 
       scatter_quant$beta.a + 1.96*scatter_quant$se.a, scatter_quant$beta.na, 
       length=0.01, angle=90, code=3, col='grey50')
par(new=T)
plot(scatter_quant[!scatter_quant$ukb_variable %in% sig_quant_traits, c('beta.a','beta.na')], 
     xaxt='n',yaxt='n', xlim=range, ylim=range, xlab='', ylab='',pch=16, col='grey70', cex=0.6)

cat_indx <- c()
for (cat in names(table(scatter_quant[scatter_quant$ukb_variable %in% sig_quant_traits,'category']))) {
  cat_match <- substr(cat,1,5)
  par(new=T)
  plot(scatter_quant[(scatter_quant$ukb_variable %in% sig_quant_traits) & scatter_quant$category==cat, 
                     c('beta.a','beta.na')], pch=21, cex=0.8, xaxt='n',yaxt='n', 
       xlim=range, ylim=range, xlab='', ylab='', bg=group_cols[substr(names(group_cols),1,5)==cat_match])
  cat_indx <- c(cat_indx, cat_match)
}

names <- c()
for (trait in sig_quant_traits) {
  row <- scatter_dat[scatter_dat$ukb_variable==trait,]
  if (!row$display_name %in% names) {
    if (row$beta.a>0) {
      text(x=row$beta.a+.032, y=row$beta.na-.008, label=row$display_name, cex=0.5, col='black', adj=1)
    } else {
      text(x=row$beta.a-.004, y=row$beta.na+.009, label=row$display_name, cex=0.5, col='black', adj=1)
    }
    names <- c(names, row$display_name)
  }
}

legend('bottomright',legend=names(group_cols[substr(names(group_cols),1,5) %in% cat_indx]), bty='n', cex=0.7,
       title.cex=0.7, title.font=2, title.adj=-0.02,
       pch=21,pt.bg=group_cols[substr(names(group_cols),1,5) %in% cat_indx])
dev.off()


# FIGURE 6c (PheWAS associations in non-asthmatics, forest plot binary traits)
#
# NOTE: category colors and text emphasis were manually added
#

pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))
fig_data <- pwas_results[pwas_results$ancestry==bw,]

# Put traits in desired order
min_group_p <- tibble(fig_data) %>% group_by(category) %>% summarise(min_p=min(p, na.rm = T))
min_group_p <- min_group_p[order(min_group_p$min_p),]
min_group_p$order <- 1:nrow(min_group_p)

## BW PheWAS Results
fig_data <- pwas_na_results[!is.na(pwas_na_results$p.na),]

fig_data[fig_data$ukb_variable %in% c('blood_clot_hayfever','blood_clot_none','hayfever_doctor_diagnosed0'),
         'display_name'] <- 'Hay fever'

# Filter data for only significant findings
fig_data <- fig_data[fig_data$sig==1,]

fig_data[fig_data$invert==1,'beta.a'] <- fig_data[fig_data$invert==1,'beta.a']*-1
fig_data[fig_data$bin==1,'CI_L.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] - 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'CI_U.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] + 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'beta.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'])
fig_data[fig_data$bin==0,'CI_L.a'] <- fig_data[fig_data$bin==0,'beta.a'] - 1.96*fig_data[fig_data$bin==0,'se.a']
fig_data[fig_data$bin==0,'CI_U.a'] <- fig_data[fig_data$bin==0,'beta.a'] + 1.96*fig_data[fig_data$bin==0,'se.a']

# Only plot most significant association for redundant traits (e.g. icd10 vs self-report)
fig_data <- data.frame(tibble(fig_data) %>% group_by(display_name) %>% filter(p.a==min(p.a)))
fig_data <- fig_data[!duplicated(fig_data$display_name),]

fig_data[fig_data$p.a==0,'p.a'] <- 5e-324
fig_data[fig_data$p.na==0,'p.na'] <- 5e-324

## BW PheWAS non-Asthmatics - Results
fig_data[fig_data$invert==1,'beta.na'] <- fig_data[fig_data$invert==1,'beta.na']*-1
fig_data[fig_data$bin==1,'CI_L.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'] - 1.96*fig_data[fig_data$bin==1,'se.na'])
fig_data[fig_data$bin==1,'CI_U.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'] + 1.96*fig_data[fig_data$bin==1,'se.na'])
fig_data[fig_data$bin==1,'beta.na'] <- exp(fig_data[fig_data$bin==1,'beta.na'])
fig_data[fig_data$bin==0,'CI_L.na'] <- fig_data[fig_data$bin==0,'beta.na'] - 1.96*fig_data[fig_data$bin==0,'se.na']
fig_data[fig_data$bin==0,'CI_U.na'] <- fig_data[fig_data$bin==0,'beta.na'] + 1.96*fig_data[fig_data$bin==0,'se.na']

# Put data in order of beta difference, filter to top results
fig_data <- merge(fig_data, min_group_p[,c('category','order')], by='category')
fig_data <- fig_data[order(fig_data$order,fig_data$p.a, 1/fig_data$beta.a),]

fig_data_quant <- fig_data[fig_data$bin==0,]
fig_data_bin <- fig_data[fig_data$bin==1,]

# Text on plot
traits <- fig_data_bin$display_name
p.a <- formatC(fig_data_bin$p.a, format = "e", digits = 0)
p.na <- formatC(fig_data_bin$p.na, format = "e", digits = 0)
p.diff <- formatC(fig_data_bin$p_diff, format = "e", digits = 0)

tabletext <- list(
  c("Trait\n",traits),
  c("Association P\n(All)    \n", p.a),
  c("Association P\n(Non-asthma)\n", p.na),
  c('Effect difference\nP\n', p.diff)
)

rows <- as.character(1:(nrow(fig_data_bin)+2))
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_6c.pdf'), width=12, height=5)  # requires xquartz/encodings
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
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("WB Non-asthmatics", "WB"), 
           ci.vertices = TRUE)
dev.off()


# FIGURE 6d (PheWAS associations in non-asthmatics, forest plot quantitative traits)
#
# NOTE: category colors and text bolding were manually added. Extra NA rows were created to match
#       dimensions in 6c and then subsequently removed manually.

# Text on plot
traits <- fig_data_quant$display_name
p.a <- formatC(fig_data_quant$p.a, format = "e", digits = 0)
p.na <- formatC(fig_data_quant$p.na, format = "e", digits = 0)
p.diff <- formatC(fig_data_quant$p_diff, format = "e", digits = 0)

#p.a[p.a=="0e+00"] <- "<5e-324"

tabletext <- list(
  c("Trait\n",traits,NA,NA),
  c("P-value\n(All)\n", p.a,NA,NA),
  c("P-value\n(Non-asthma)\n", p.na,NA,NA),
  c('Effect difference\nP\n', p.diff,NA,NA)
)

rows <- as.character(1:(nrow(fig_data_quant)+4))
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_6d.pdf'), width=12, height=5)  # requires xquartz/encodings
forestplot(tabletext, 
           mean = cbind(c(NA,fig_data_quant$beta.na,NA,NA), c(NA,fig_data_quant$beta.a,NA,NA)),
           lower = cbind (c(NA,fig_data_quant$CI_L.na,NA,NA), c(NA,fig_data_quant$CI_L.a,NA,NA)), 
           upper = cbind(c(NA,fig_data_quant$CI_U.na,NA,NA), c(NA,fig_data_quant$CI_U.a,NA,NA)),
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
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("WB Non-asthmatics", "WB All"), 
           ci.vertices = TRUE)
dev.off()


#####  FIGURE 7 (Signif diff traits w/o HLA, forest plots) ###################################################
#
# NOTE: Binary and quantitative plots were manually aligned and special characters were subsequently added

library(tidyverse)
library(forestplot)
library(grid)
library(scales)

# pwas_results <- read.delim(paste0(output_dir, '/phewas_results.txt'))
# pwas_hla_results <- read.delim(paste0(output_dir, '/pwas_noHLA_results.tsv'))

fig_data <- pwas_results[pwas_results$ancestry==bw,]

# Put traits in desired order
min_group_p <- tibble(fig_data) %>% group_by(category) %>% summarise(min_p=min(p, na.rm = T))
min_group_p <- min_group_p[order(min_group_p$min_p),]
min_group_p$order <- 1:nrow(min_group_p)

## BW PheWAS Results
fig_data <- pwas_hla_results[(pwas_hla_results$sig==1) | (pwas_hla_results$ukb_variable %in% c(coa,aoa)),]

fig_data[fig_data$invert==1,'beta.a'] <- fig_data[fig_data$invert==1,'beta.a']*-1
fig_data[fig_data$bin==1,'CI_L.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] - 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'CI_U.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'] + 1.96*fig_data[fig_data$bin==1,'se.a'])
fig_data[fig_data$bin==1,'beta.a'] <- exp(fig_data[fig_data$bin==1,'beta.a'])
fig_data[fig_data$bin==0,'CI_L.a'] <- fig_data[fig_data$bin==0,'beta.a'] - 1.96*fig_data[fig_data$bin==0,'se.a']
fig_data[fig_data$bin==0,'CI_U.a'] <- fig_data[fig_data$bin==0,'beta.a'] + 1.96*fig_data[fig_data$bin==0,'se.a']

# Only plot most significant association for redundant traits (e.g. icd10 vs self-report)
fig_data <- data.frame(tibble(fig_data) %>% group_by(display_name) %>% filter(p_diff==min(p_diff)))
fig_data <- fig_data[!duplicated(fig_data$display_name),]

fig_data[fig_data$p.a==0,'p.a'] <- 5e-324
fig_data[fig_data$p.hla==0,'p.hla'] <- 5e-324

## BW PheWAS non-Asthmatics - Results
fig_data[fig_data$invert==1,'beta.hla'] <- fig_data[fig_data$invert==1,'beta.hla']*-1
fig_data[fig_data$bin==1,'CI_L.hla'] <- exp(fig_data[fig_data$bin==1,'beta.hla'] - 1.96*fig_data[fig_data$bin==1,'se.hla'])
fig_data[fig_data$bin==1,'CI_U.hla'] <- exp(fig_data[fig_data$bin==1,'beta.hla'] + 1.96*fig_data[fig_data$bin==1,'se.hla'])
fig_data[fig_data$bin==1,'beta.hla'] <- exp(fig_data[fig_data$bin==1,'beta.hla'])
fig_data[fig_data$bin==0,'CI_L.hla'] <- fig_data[fig_data$bin==0,'beta.hla'] - 1.96*fig_data[fig_data$bin==0,'se.hla']
fig_data[fig_data$bin==0,'CI_U.hla'] <- fig_data[fig_data$bin==0,'beta.hla'] + 1.96*fig_data[fig_data$bin==0,'se.hla']

# Put data in order of beta difference, filter to top results
fig_data <- merge(fig_data, min_group_p[,c('category','order')], by='category')
fig_data <- fig_data[order(fig_data$delta_beta,fig_data$p.a, 1/fig_data$beta.a),]

fig_data_quant <- fig_data[fig_data$bin==0,]
fig_data_bin <- fig_data[fig_data$bin==1,]

# Text on plot
traits <- fig_data_bin$display_name
delta.beta <- percent(fig_data_bin$delta_beta)
p.diff <- formatC(fig_data_bin$p_diff, format = "e", digits = 0)

tabletext <- list(
  c("Trait\n",traits),
  c('Effect difference\ndBeta\n', delta.beta),
  c('Effect difference\nP\n', p.diff)
)

rows <- as.character(1:(nrow(fig_data_bin)+2))
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_7a.pdf'), width=12, height=5)  # requires xquartz/encodings
forestplot(tabletext, 
           mean = cbind(c(NA,fig_data_bin$beta.hla), c(NA,fig_data_bin$beta.a)),
           lower = cbind(c(NA,fig_data_bin$CI_L.hla), c(NA,fig_data_bin$CI_L.a)), 
           upper = cbind(c(NA,fig_data_bin$CI_U.hla), c(NA,fig_data_bin$CI_U.a)),
           new_page = TRUE,
           clip = c(0.8,2), 
           hrzl_lines=line_list,
           zero=1,
           xlog = F, xlab = "OR", 
           col = fpColors(box = c("slateblue3", "firebrick4"),
                          lines = c("slateblue", "firebrick3")),
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawDiamondCI),
           is.summary = c(TRUE,rep(FALSE,nrow(fig_data))), 
           graph.pos = 2,
           graphwidth=unit(120,'mm'),
           boxsize = 0.35, 
           line.margin = .51,
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("HLA-excluded PRS", "GBMI PRS"), 
           ci.vertices = TRUE)
dev.off()


# FIGURE 7b (Signif diff traits w/o HLA, forest plot quant traits)

# Text on plot
traits <- fig_data_quant$display_name
delta.beta <- percent(fig_data_quant$delta_beta)
p.diff <- formatC(fig_data_quant$p_diff, format = "e", digits = 0)

tabletext <- list(
  c("Trait\n",traits),
  c('Effect difference\ndBeta\n', delta.beta),
  c('Effect difference\nP\n', p.diff)
)

rows <- as.character(1:(nrow(fig_data_quant)+2))
line_list <- list("1"=gpar(col=rgb(1,1,1,0),lty='dotted'), "2"=gpar(col='gray20'))
for (i in c(3:length(rows))) {
  line_list[[rows[i]]] <- gpar(col='gray25',lty=3)
}

# Plot
cairo_pdf(filename=paste0(output_dir, '/Fig_7b.pdf'), width=12, height=2)  # requires xquartz/encodings
forestplot(tabletext, 
           mean = cbind(c(NA,fig_data_quant$beta.hla), c(NA,fig_data_quant$beta.a)),
           lower = cbind (c(NA,fig_data_quant$CI_L.hla), c(NA,fig_data_quant$CI_L.a)), 
           upper = cbind(c(NA,fig_data_quant$CI_U.hla), c(NA,fig_data_quant$CI_U.a)),
           new_page = TRUE,
           clip = c(-0.2,.2), 
           hrzl_lines=line_list,
           xlog = F, xlab = "Beta", 
           col = fpColors(box = c("slateblue3", "firebrick4"),
                          lines = c("slateblue", "firebrick3")),
           fn.ci_norm = c(fpDrawDiamondCI, fpDrawDiamondCI),
           is.summary = c(TRUE,rep(FALSE,nrow(fig_data))), 
           graph.pos = 2,
           graphwidth=unit(120,'mm'),
           boxsize = 0.35, 
           line.margin = .51,
           lwd.zero=1.75,
           txt_gp= fpTxtGp(label=gpar(cex=0.8), title=gpar(cex=0.8), legend=gpar(cex=0.8), 
                           xlab=gpar(cex=0.8), ticks=gpar(cex=0.8),summary=gpar(cex=0.75)),
           legend = c("HLA-excluded PRS", "GBMI PRS"), 
           ci.vertices = TRUE)
dev.off()


#####  FIGURE 8 (Mediation Bar Plot) ########################################################################

load(paste0(output_dir, 'mediation_results.RData'))

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
arrows(x0=b,x1=b,y0=c(med_95L_df[,1], med_95L_df[,2], med_95L_df[,3]),
       y1=c(med_95H_df[,1], med_95H_df[,2], med_95H_df[,3]), angle=90, length=0.02, code=3)
dev.off()



