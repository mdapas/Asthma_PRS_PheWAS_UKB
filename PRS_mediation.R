###########################################################################################
###  TRAIT MEDIATION ANALYSIS                                                           ###
###########################################################################################

# This script executes mediation analyses for eosinophil percentage, FEV1/FVC ratio, 
#   age of hay fever diagnosis, and BMI, and it also produces Figure 8 from the manuscript.
#   The equations used to estimate the mediation proportions were adopted from 
#   MacKinnon 2007 (doi: 10.1177/1740774507083434) and Li 2007 (doi: 10.1002/sim.2730).
#   Mediation p-values were estimated using the ROBMED package (Alfons 2021, 
#   doi: 10.1177/1094428121999096)
#


library(robmed)

seed <- 20220117
bw <- 'White'; nbw <- 'non-British White'; afr <- 'African'
prs <- 'zscore'; covars <- '+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10'

# trait_dat <- read.csv(['path/to/ukb/data/trait_data.csv'])
trait_dat <- trait_dat[!is.na(trait_dat[,prs]),]
trait_dat <- trait_dat[!is.na(trait_dat$ancestry),]

# Normalize trait values within each ancestry group
for (g in c(bw, nbw, afr)) {
  trait_dat[trait_dat$ancestry==g,'eosinophil_z'] <- rntransform(trait_dat[trait_dat$ancestry==g,'eosinophill_percentage0'])
  trait_dat[trait_dat$ancestry==g,'fev1_fvc_z'] <- ztransform(trait_dat[trait_dat$ancestry==g,'fev1_fvc_ratio'])
  trait_dat[trait_dat$ancestry==g,'age_hay_z'] <- ztransform(trait_dat[trait_dat$ancestry==g,'age_hay_fever_diagnosed'])
  trait_dat[trait_dat$ancestry==g,'bmi_z'] <- ztransform(trait_dat[trait_dat$ancestry==g,'bmi_body_size0'])
}


# Calculate proportion mediated (we reported the a*b/(c_std) values, as from MacKinnon 2007)
traits <- c('eosinophil_z','fev1_fvc_z','age_hay_z')
phenos <- c('asthma_children','asthma_adults','asthma_all')

med_df <- setNames(data.frame(matrix(ncol = 3, nrow = 3)), phenos)
rownames(med_df) <- traits

for (pheno in phenos){
  for (trait in traits){
    model_dat <- trait_dat[trait_dat$ancestry==bw,]
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
    # print(c(pheno, trait, round(a,6), round(c_std,6), round(a*b/((a*b)+c_prime),6), round(a*b/(c_std),6)))
    med_df[trait, pheno] <- round(a*b/(c_std),6)
  }
}


#####  Bar Plot (Figure 8) #####
names(med_df) <- c("Childhood-onset\nasthma","Adult-onset\nasthma","All\nasthma")
cols <- c('cadetblue3','darkseagreen3','coral1')

par(xpd=F, mar=c(3,4.5,2,2), lheight = .3)
b <- barplot(height=as.matrix(med_df), ylab="Proportion mediated", beside=TRUE, ylim=c(-0.05,0.2), col=cols,
             cex.lab=0.8, cex.names=0.8, font=2, yaxt='n')
axis(2, cex.axis=0.8)                
grid(nx = NA, ny = NULL, lty = 2, col = "gray", lwd = 1)
abline(h=0,col='black')
barplot(height=as.matrix(med_df), beside=TRUE, ylim=c(-0.05,0.2), col=cols, add=T, xaxt='n',yaxt='n', ylab='')
text(x=b, y=c(coa,aoa,all), pos=c(1,1,1,1,1,3,1,1,1), label=round(c(coa,aoa,all),2), cex=0.6, col='black')
legend('topright', legend=c('Eosinophil count', expression('FEV'[1]*'/FVC'), 'Age hay fever diagnosed'), cex=0.7,
       fill=cols, box.col='gray75', bty=1, inset=c(0.02,0.01), x.intersp=0.75, y.intersp=0.9, text.width=3.5)


# Estimate asthma mediation of BMI association, proportion mediated from Li 2007
bmi_med_df <- setNames(data.frame(matrix(ncol = 2, nrow = 3)), c('Proportion','P-value'))
rownames(bmi_med_df) <- phenos

for (pheno in phenos){
  model_dat <- trait_dat[trait_dat$ancestry==bw,]
  if (pheno=='asthma_children') {
    test_dat <- model_dat[!(is.na(model_dat$asthma_children) & model_dat$asthma_all==1),]
    test_dat <- test_dat[!is.na(model_dat$asthma_adults) & model_dat$asthma_adults!=1,]
  } else if (pheno=='asthma_adults') {
    test_dat <- model_dat[!(is.na(model_dat$asthma_adults) & model_dat$asthma_all==1),]
    test_dat <- test_dat[!is.na(model_dat$asthma_children) & model_dat$asthma_children!=1,]
  } else {
    test_dat <- model_dat
  }
  tm_sobel <- test_mediation(test_dat, prs, 'bmi_z',pheno, 
                             test='sobel', covariates=strsplit(covars,'\\+')[[1]][-1], robust=F)

  model.y <- lm(paste0('bmi_z ~', prs, covars), data=test_dat)
  model.m <- glm(paste0(pheno, '~', prs, covars), data=test_dat, family='binomial')
  model.total <- lm(paste0('bmi_z ~', pheno, '+', prs, covars), data=test_dat)
  
  a_hat <- coefficients(model.m)[[prs]]
  b_hat <- coefficients(model.total)[[pheno]]
  b1 <- coefficients(model.y)[[prs]]
  delt_AL <- a_hat*b_hat*(exp(sum(coefficients(model.m)))/((1 + exp(sum(coefficients(model.m))))^2))
  
  a <- coefficients(model.m)[[prs]]
  b <- coefficients(model.total)[[pheno]]
  c_prime <- coefficients(model.total)[[prs]]
  c <- coefficients(model.y)[[prs]]
  
  prop_AL <- delt_AL/(delt_AL + c_prime)
  
  # print(c(pheno, prop_AL, tm_sobel$p_value))
  bmi_med_df[pheno, ] <- c(prop_AL, tm_sobel$p_value)
}
