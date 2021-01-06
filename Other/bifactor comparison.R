library(lavaan)
library(simstandard)
library(tidyverse)

bimod <- "G =~ .55*x1 + .55*x2 + .55*x3 + .55*x4 + .55*x5 + .55*x6 + .55*x7 + .55*x8 + .55*x9 + .55*x10 + .55*x11 + .55*x12 
          F1 =~ .55*x1 + .55*x2 + .55*x3 + .55*x4
          F2 =~ .55*x5 + .55*x6 + .55*x7 + .55*x8
          F3 =~ .55*x9 + .55*x10 + .55*x11 + .55*x12
          G ~~ 0*F1
          G ~~ 0*F2
          G ~~ 0*F3
          F1 ~~ 0*F2
          F1 ~~ 0*F3
          F2 ~~ 0*F3"

cfmod <- "F1 =~ .55*x1 + .55*x2 + .55*x3 + .55*x4
          F2 =~ .55*x5 + .55*x6 + .55*x7 + .55*x8
          F3 =~ .55*x9 + .55*x10 + .55*x11 + .55*x12
          F1 ~~ .55*F2
          F1 ~~ .55*F3
          F2 ~~ .55*F3"

dfbi <- sim_standardized(bimod, n=1000, latent = FALSE, errors = FALSE)

cor(dfbi)

dfcf <- sim_standardized(cfmod, n=1000, latent = FALSE, errors = FALSE)

cor(dfcf)

bimod2 <- "G =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12
           F1 =~ x1 + x2 + x3 + x4
           F2 =~ x5 + x6 + x7 + x8
           F3 =~ x9 + x10 + x11 + x12
           G ~~ 0*F1
           G ~~ 0*F2
           G ~~ 0*F3
           F1 ~~ 0*F2
           F1 ~~ 0*F3
           F2 ~~ 0*F3"

cfmod2 <- "F1 =~ x1 + x2 + x3 + x4
           F2 =~ x5 + x6 + x7 + x8
           F3 =~ x9 + x10 + x11 + x12"




bifit <- sem(bimod2, dfbi, std.lv=TRUE)
bt <- fitMeasures(bifit, fit.measures = c("SRMR", "RMSEA", "CFI"))
#summary(bifit)
bifit_cf <- sem(cfmod2, dfbi, std.lv=TRUE)
bm <- fitMeasures(bifit_cf, fit.measures = c("SRMR", "RMSEA", "CFI"))
#summary(bifit_cf)

cffit <- sem(cfmod2, dfcf, std.lv=TRUE)
ct <- fitMeasures(cffit, fit.measures = c("SRMR", "RMSEA", "CFI"))
#summary(cffit)
cffit_bi <- sem(bimod2, dfcf, std.lv=TRUE)
cm <- fitMeasures(cffit_bi, fit.measures = c("SRMR", "RMSEA", "CFI"))
#summary(cffit_bi)

true <- round(c(bt, ct),3)
misspec <- round(c(bm, cm),3)

rbind(true,misspec)

###########

bt <- as.data.frame(matrix(nrow=200,ncol=6))
for(i in 1:200){
        dfbi <- sim_standardized(bimod, n=1000, latent = FALSE, errors = FALSE)
        bifit <- sem(bimod2, dfbi, std.lv=TRUE)
        bt[i,] <- fitMeasures(bifit, fit.measures = c("chisq", "df", "pvalue", "SRMR", "RMSEA", "CFI"))
}

bm <- as.data.frame(matrix(nrow=200,ncol=6))
for(i in 1:200){
        dfbi <- sim_standardized(bimod, n=1000, latent = FALSE, errors = FALSE)
        bifit_cf <- sem(cfmod2, dfbi, std.lv=TRUE)
        bm[i,] <- fitMeasures(bifit_cf, fit.measures = c("chisq", "df", "pvalue", "SRMR", "RMSEA", "CFI"))
}

ct <- as.data.frame(matrix(nrow=200,ncol=6))
for(i in 1:200){
        dfcf <- sim_standardized(cfmod, n=1000, latent = FALSE, errors = FALSE)
        cffit <- sem(cfmod2, dfcf, std.lv=TRUE)
        ct[i,] <- fitMeasures(cffit, fit.measures = c("chisq", "df", "pvalue", "SRMR", "RMSEA", "CFI"))
}

cm <- as.data.frame(matrix(nrow=200,ncol=6))
for(i in 1:200){
        dfcf <- sim_standardized(cfmod, n=1000, latent = FALSE, errors = FALSE)
        cffit_bi <- sem(bimod2, dfcf, std.lv=TRUE)
        cm[i,] <- fitMeasures(cffit_bi, fit.measures = c("chisq", "df", "pvalue", "SRMR", "RMSEA", "CFI"))
}

btm <- round(colMeans(bt),3)
bmm <- round(colMeans(bm),3)
ctm <- round(colMeans(ct),3)
cmm <- round(colMeans(cm),3)

rbind(btm,bmm)
rbind(ctm,cmm)


###########


dfbila <- simulateData(model=bimod,sample.nobs = 1000, model.type="cfa",
                       orthogonal = TRUE)
bifitla <- sem(bimod2, data=dfbila, std.lv=TRUE)
fitMeasures(bifitla, fit.measures = c("SRMR", "RMSEA", "CFI"))



###########

bimisc <- "G =~ .44*x1 + .44*x2 + .44*x3 + .44*x4 + .44*x5 + .44*x6 + .44*x7 + .44*x8 + .44*x9 + .44*x10 + .44*x11 + .44*x12 
           F1 =~ .44*x1 + .44*x2 + .44*x3 + .44*x4 
           F2 =~ .44*x5 + .44*x6 + .44*x7 + .44*x8
           F3 =~ .44*x9 + .44*x10 + .44*x11 + .44*x12
           x1 ~~ .44*x6
           G ~~ 0*F1
           G ~~ 0*F2
           G ~~ 0*F3
           F1 ~~ 0*F2
           F1 ~~ 0*F3
           F2 ~~ 0*F3"

bi_mis <- as.data.frame(matrix(nrow=250,ncol=3))
for(i in 1:250){
        bdat <- sim_standardized(bimisc,500,latent = FALSE,errors = FALSE)
        bmfit <- sem(bimod2,bdat,std.lv=TRUE)
        bi_mis[i,]<- fitMeasures(bmfit, fit.measures = c("SRMR","RMSEA","CFI"))  
}

bi_mis %>% 
        summarise(SRMR=quantile(V1,.05),
                  RMSEA=quantile(V2,.05),
                  CFI=quantile(V3,.95))

cfmisc <- "F1 =~ .44*x1 + .44*x2 + .44*x3 + .44*x4 
           F2 =~ .44*x5 + .44*x6 + .44*x7 + .44*x8
           F3 =~ .44*x9 + .44*x10 + .44*x11 + .44*x12
           F1 =~ .44*x6"

cf_mis <- as.data.frame(matrix(nrow=250,ncol=3))
for(i in 1:250){
        cdat <- sim_standardized(cfmisc,500,latent = FALSE,errors = FALSE)
        cmfit <- cfa(cfmod2,cdat,std.lv=TRUE)
        cf_mis[i,]<- fitMeasures(cmfit, fit.measures = c("SRMR","RMSEA","CFI"))  
}

cf_mis %>% 
        summarise(SRMR=quantile(V1,.05),
                  RMSEA=quantile(V2,.05),
                  CFI=quantile(V3,.95))
