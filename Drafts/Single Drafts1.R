

co <- function(x,y){

  c <- x*y
  
  return(c)
}

co(.8,.8)


model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5"
.345*.602
.602*.631
.631*.805
.345*.631
.805*.857
.857*.631

model <- "F1 =~ .85*Y1 + .85*Y2 + .85*Y3 + .85*Y4 + 85*Y5"

loads <- c(.602,.805,.857,.631,.345)

vec <- numeric(10)

for (i in loads){
  vec[i] <- loads[1]*loads[2]
}

d <- combn(loads, 2, FUN = function(x) x[1] * x[2])
mean(d)

(.85*.85)*.95

mod_85 <- "F1 =~ .85*Y1 + .85*Y2 + .85*Y3 + .85*Y4 + .85*Y5
           Y1 ~~ .72*Y2"

mod_8585 <- "F1 =~ .85*Y1 + .85*Y2 + .85*Y3 + .85*Y4 + .85*Y5
             Y1 ~~ .72*Y2
             Y3 ~~ .72*Y4"

mod_8585d <- "F1 =~ .85*Y1 + .85*Y2 + .85*Y3 + .85*Y4 + .85*Y5
             Y1 ~~ .72*Y2
             Y2 ~~ .72*Y4"



mod_low12 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
              Y5 ~~ .208*Y1" 

mod_low23 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
              Y1 ~~ .380*Y4"

mod_low1223 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
              Y5 ~~ .208*Y1
              Y1 ~~ .380*Y4"

mod_low45 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
              Y3 ~~ .690*Y2"
mod_low4535 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y3 ~~ .690*Y2
                Y3 ~~ .541*Y4"

mod_low1234 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .208*Y1
                Y4 ~~ .508*Y2" 

mod_low1245 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .208*Y1
                Y3 ~~ .690*Y2"

mod_low1212 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .208*Y1
                Y4 ~~ .208*Y2" 

mod_low1213 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .208*Y1
                Y4 ~~ .218*Y2" 

mod_mean <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .412*Y1" 

mod_mean2 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5
                Y5 ~~ .412*Y1
                Y4 ~~ .412*Y2" 

clean <- cleanmodel(model)

mdat85 <- sim_standardized(m=mod_85,n=500*500,latent=FALSE,errors=FALSE)
mdat8585 <- sim_standardized(m=mod_8585,n=500*500,latent=FALSE,errors=FALSE)
mdat8585d <- sim_standardized(m=mod_8585d,n=500*500,latent=FALSE,errors=FALSE)

mdat12 <- sim_standardized(m=mod_low12,n=500*500,latent=FALSE,errors=FALSE)
mdat23 <- sim_standardized(m=mod_low23,n=500*500,latent=FALSE,errors=FALSE)
mdat45 <- sim_standardized(m=mod_low45,n=500*500,latent=FALSE,errors=FALSE)
mdat4535 <- sim_standardized(m=mod_low4535,n=500*500,latent=FALSE,errors=FALSE)
mdat1234 <- sim_standardized(m=mod_low1234,n=500*500,latent=FALSE,errors=FALSE)
mdat1245 <- sim_standardized(m=mod_low1245,n=500*500,latent=FALSE,errors=FALSE)
mdat1212 <- sim_standardized(m=mod_low1212,n=500*500,latent=FALSE,errors=FALSE)
mdat1213 <- sim_standardized(m=mod_low1213,n=500*500,latent=FALSE,errors=FALSE)
mdat1223 <- sim_standardized(m=mod_low1223,n=500*500,latent=FALSE,errors=FALSE)

mdatmean <- sim_standardized(m=mod_mean,n=500*500,latent=FALSE,errors=FALSE)
mdatmean2 <- sim_standardized(m=mod_mean2,n=500*500,latent=FALSE,errors=FALSE)



repid <- rep(1:500,500)

mdat85rep <- cbind(mdat85,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat8585rep <- cbind(mdat8585,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat12rep <- cbind(mdat12,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat23rep <- cbind(mdat23,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat45rep <- cbind(mdat45,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat4535rep <- cbind(mdat4535,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat1234rep <- cbind(mdat1234,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat1245rep <- cbind(mdat1245,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat1212rep <- cbind(mdat1212,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat1213rep <- cbind(mdat1213,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdat1223rep <- cbind(mdat1223,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdatmeanrep <- cbind(mdatmean,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

mdatmean2rep <- cbind(mdatmean2,repid) %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

cfa85 <- lapply(mdat85rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))
cfa8585 <- lapply(mdat8585rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))


cfa12 <- lapply(mdat12rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa23 <- lapply(mdat23rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa45 <- lapply(mdat45rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa4535 <- lapply(mdat4535rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa1234 <- lapply(mdat1234rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa1245 <- lapply(mdat1245rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa1212 <- lapply(mdat1212rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa1213 <- lapply(mdat1213rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfa1223 <- lapply(mdat1223rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfamean <- lapply(mdatmeanrep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))

cfamean2 <- lapply(mdatmean2rep[[2]], function(x) cfa(model=clean,data=x,std.lv=TRUE))


fit85 <- map_dfr(cfa85, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit8585 <- map_dfr(cfa8585, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

fit12 <- map_dfr(cfa12, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit23 <- map_dfr(cfa23, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit45 <- map_dfr(cfa45, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit4535 <- map_dfr(cfa45, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit1234 <- map_dfr(cfa1234, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit1245 <- map_dfr(cfa1245, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))


fit1212 <- map_dfr(cfa1212, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit1213 <- map_dfr(cfa1213, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fit1223 <- map_dfr(cfa1223, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))


fitmean <- map_dfr(cfamean, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
fitmean2 <- map_dfr(cfamean2, ~fitMeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

library(psych)

describe(fit85,quant = c(.05,.10,.90,.95))
describe(fit8585,quant = c(.05,.10,.90,.95))

describe(fit12,quant = c(.05,.10,.90,.95))
describe(fit23,quant = c(.05,.10,.90,.95))
describe(fit45,quant = c(.05,.10,.90,.95))
describe(fit4535,quant = c(.05,.10,.90,.95))
describe(fit1234,quant = c(.05,.10,.90,.95))
describe(fit1245,quant = c(.05,.10,.90,.95))


describe(fit1212,quant = c(.05,.10,.90,.95))
describe(fit1213,quant = c(.05,.10,.90,.95))
describe(fit1223,quant = c(.05,.10,.90,.95))

describe(fitmean,quant = c(.05,.10,.90,.95))
describe(fitmean2,quant = c(.05,.10,.90,.95))





