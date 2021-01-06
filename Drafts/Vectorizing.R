# This is a bad loop with 'growing' data
set.seed(42)
m=10
n=10

# Create matrix of normal random numbers
mymat <- replicate(m, rnorm(n))

# Transform into data frame
mydframe <- data.frame(mymat)

for (i in 1:m) {
  for (j in 1:n) {
    mydframe[i,j]<-mydframe[i,j] + 10*sin(0.75*pi)
    print(mydframe)
  }
}

#Vectorized

set.seed(42)
m=10
n=10
mymat <- replicate(m, rnorm(n)) 
mydframe <- data.frame(mymat)
mydframe <- mydframe + 10*sin(0.75*pi)
mydframe

##

misspecified_model_fit.v <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM(model)
  
  #Use max sample size of 10000
  n <- min(n,10000)
  
  #Simulate data and return fit stats
  set.seed(649364)
  
  misspec_data <- replicate(n=500,sim_standardized(m=misspec_dgm,n = n,
                                        latent = FALSE,
                                        errors = FALSE),simplify = FALSE)
  
  misspec_cfa <- map(misspec_data,~cfa(model = mod, data=., std.lv=TRUE))
  
  misspec_fit_sum <- map_dfr(misspec_cfa,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
    `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
    dplyr::mutate(Type_M="Misspecified")
  
  set.seed(NULL)
  
  #return(misspec_data)
  return(misspec_fit_sum)
  
}

system.time({misspecified_model_fit.v(mod.test,500)})

system.time({misspecified_model_fit(mod.test,500)})

library(purrr)



mod.test <- "f1 =~ .5*x1 + .44*x2 + .643*x3 + .454*x4"
mod.test2 <- "f1 =~ x1 + x2 + x3 + x4"

dat <- replicate(n=5,sim_standardized(m=mod.test,n = 200,
                           latent = FALSE,
                           errors = FALSE),simplify = FALSE)

dat2 <- map(dat,~cfa(model = mod.test2, data=., std.lv=TRUE))

dat3 <- map_dfr(dat2,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
  dplyr::mutate(Type_M="Misspecified")

## lapply

datl2 <- lapply(dat,cfa(model = mod.test2, data=dat[[1]], std.lv=TRUE))

names(dat)

##Loop (data only)

misspecified_model_fit <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM(model)
  
  #Create empty df to put fit stats in
  #misspec_fit <- data.frame(matrix(nrow=n*500,ncol=4))
  misspec_data2 <- data.frame(matrix(nrow=n*500,ncol=4))
  
  #Use max sample size of 10000
  n <- min(n,10000)
  
  #Simulate data through loop 500 times
  set.seed(649364)
  for (i in 1:500){
    misspec_data <- simstandard::sim_standardized(m=misspec_dgm,n = n,
                                                  latent = FALSE,
                                                  errors = FALSE)
    misspec_data2[i,] <- misspec_data
    #misspec_cfa <- base::withCallingHandlers(
    #lavaan::cfa(model = mod, data = misspec_data, std.lv=TRUE), warning = hide_ov)
    #misspec_fits <- lavaan::fitMeasures(misspec_cfa, c("srmr","rmsea","cfi"))
    #misspec_fit[i,] <- misspec_fits
  }
  set.seed(NULL)
  
  #Clean up data
  #misspec_fit_sum <- misspec_fit %>%
  #`colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
  #dplyr::mutate(Type_M="Misspecified")
  
  return(misspec_data2)
  #return(misspec_fit_sum)
  
}

### n = 250000

misspec_data <- sim_standardized(m=mod.test,n = 500*500,
                                 latent = FALSE,
                                 errors = FALSE)

rep_id <- rep(1:500,500)

dat_rep <- cbind(misspec_data,rep_id)

nested <- dat_rep %>% 
  group_by(rep_id) %>% 
  nest() %>% 
  as.list()

system.time({misspec_cfa <- map(nested$data,~cfa(model = mod.test2, data=., std.lv=TRUE))})

misspec_fit_sum <- map_dfr(misspec_cfa,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
  dplyr::mutate(Type_M="Misspecified")

### by

misspec_data <- sim_standardized(m=mod.test,n = 500*500,
                                 latent = FALSE,
                                 errors = FALSE)

rep_id <- rep(1:500,500)

dat_rep <- cbind(misspec_data,rep_id)

system.time({by(dat_rep,rep_id,function(x) cfa(model = mod.test2, data=x, std.lv=TRUE))})




### try in function form

misspecified_model_fit.d <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM(model)
  
  #Use max sample size of 10000
  n <- min(n,10000)
  
  #Set seed
  set.seed(649364)

  #Simulate data  
  all_data <- sim_standardized(m=misspec_dgm,n = n*500,
                                   latent = FALSE,
                                   errors = FALSE)
  
  #Create indicator
  rep_id <- rep(1:500,n)
  
  #Combine
  dat_rep <- cbind(all_data,rep_id)
  
  #Group
  misspec_data <- dat_rep %>% 
    group_by(rep_id) %>% 
    nest() %>% 
    as.list()
  
  #Run 500 cfa
  misspec_cfa <- map(misspec_data$data,~cfa(model = mod, data=., std.lv=TRUE))
  
  #Extract fit stats from each rep
  misspec_fit_sum <- map_dfr(misspec_cfa,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
    `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
    dplyr::mutate(Type_M="Misspecified")
  
  set.seed(NULL)
  
  #return(misspec_data)
  return(misspec_fit_sum)
  
}

misspecified_model_fit.d(mod.test,200)

system.time({misspecified_model_fit.d(mod.test,500)})

## Attempting to run on mult cores

library(parallel)
library(snow)

cl <- makeCluster(rep("localhost", 4), type = "SOCK")
clusterCall(cl,misspec_cfa2 <- map(nested$data,~cfa(model = mod.test2, data=., std.lv=TRUE)))
stopCluster(cl)

## furrr (multicore of purr)

library(furrr)

system.time({misspec_cfa3 <- future_map(nested$data,~cfa(model = mod.test2, data=., std.lv=TRUE))})
system.time({misspec_cfa4 <- map(nested$data,~cfa(model = mod.test2, data=., std.lv=TRUE))})