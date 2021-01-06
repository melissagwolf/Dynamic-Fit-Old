#https://community.rstudio.com/t/background-images-in-shiny/12261
#https://support.rstudio.com/hc/en-us/articles/360001157793-Static-Content-on-RStudio-Connect
#https://support.rstudio.com/hc/en-us/articles/360025989313-RStudio-Connect-Custom-Landing-Page


library(psych)

model <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6"

loads <- c(.602,.805,.857,.631,.345)

d <- combn(loads, 2, FUN = function(x) x[1] * x[2])
mean(d)

clean <- cleanmodel(model)

repid <- rep(1:500,500)

#######

mod642 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
          Y3 ~~ .642*Y2"

.91*.705

mdat642 <- sim_standardized(m=mod642,n=500*500,latent=FALSE,errors=FALSE)

mdat642rep <- cbind(mdat642,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa642 <- lapply(mdat642rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit642 <- map_dfr(cfa642, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit642,quant = c(.05,.10,.90,.95))

######

mod249 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
          Y6 ~~ .249*Y5"

.456*.545

mdat249 <- sim_standardized(m=mod249,n=500*500,latent=FALSE,errors=FALSE)

mdat249rep <- cbind(mdat249,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa249 <- lapply(mdat249rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit249 <- map_dfr(cfa249, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit249,quant = c(.05,.10,.90,.95))

#####

mod249642 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
          Y6 ~~ .249*Y5
          Y3 ~~ .642*Y2"

mdat249642 <- sim_standardized(m=mod249642,n=500*500,latent=FALSE,errors=FALSE)

mdat249642rep <- cbind(mdat249642,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa249642 <- lapply(mdat249642rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit249642 <- map_dfr(cfa249642, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit249642,quant = c(.05,.10,.90,.95))

######

.631*.601

mod379 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
          Y4 ~~ .379*Y1"

mdat379 <- sim_standardized(m=mod379,n=500*500,latent=FALSE,errors=FALSE)

mdat379rep <- cbind(mdat379,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa379 <- lapply(mdat379rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit379 <- map_dfr(cfa379, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit379,quant = c(.05,.10,.90,.95))

#####

mod249379 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
           Y6 ~~ .249*Y5          
           Y4 ~~ .379*Y1"


mdat249379 <- sim_standardized(m=mod249379,n=500*500,latent=FALSE,errors=FALSE)

mdat249379rep <- cbind(mdat249379,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa249379 <- lapply(mdat249379rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit249379 <- map_dfr(cfa249379, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit249379,quant = c(.05,.10,.90,.95))

####

mod249379642 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
           Y6 ~~ .249*Y5          
           Y4 ~~ .379*Y1
           Y3 ~~ .642*Y2"


mdat249379642 <- sim_standardized(m=mod249379642,n=500*500,latent=FALSE,errors=FALSE)

mdat249379642rep <- cbind(mdat249379642,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa249379642 <- lapply(mdat249379642rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit249379642 <- map_dfr(cfa249379642, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit249379642,quant = c(.05,.10,.90,.95))

####

mod415 <- "F1 =~ .602*Y1 + .705*Y2 + .910*Y3 + .631*Y4 + .545*Y5 + .456*Y6
           Y3 ~~ .415*Y6"


mdat415 <- sim_standardized(m=mod415,n=500*500,latent=FALSE,errors=FALSE)

mdat415rep <- cbind(mdat415,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa415 <- lapply(mdat415rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

fit415 <- map_dfr(cfa415, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

describe(fit415,quant = c(.05,.10,.90,.95))