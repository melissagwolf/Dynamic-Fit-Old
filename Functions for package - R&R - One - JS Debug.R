################
##### Test #####
################

library(lavaan)
library(tidyverse)
library(readr)

scs_data <- read_csv("C:/Users/missg/Dropbox (Personal)/Dynamic Fit/scs_data.csv")

scs.mod <- '
sc =~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10
'

scs.clean <- scs_data[,1:10]

test <- complete.cases(scs.clean)

scs.fit <- cfa(scs.mod, data = scs.clean,
               estimator = "mlr", #robust ml estimation
               std.lv = TRUE, #using fixed-factor id/ss (latent var. constrained to 1)
               missing = "ml")#FIML for missingness



########## Run #############

summary(scs.fit)

cfaOne(scs.fit)

########## Test each function #############
n <- cfa_n(scs.fit)
model <- cfa_lavmod(scs.fit)

cleanmodel(model)

#Extract standardized solution from lavaan object
lav <- lavaan::standardizedSolution(scs.fit)

#Create model statement
ss_mod <- suppressMessages(lav %>% 
                             dplyr::filter(lhs != rhs) %>% 
                             dplyr::group_by(lhs,op) %>% 
                             dplyr::select(lhs,op,rhs,est.std) %>% 
                             dplyr::mutate(est.std=round(est.std,digits=4)) %>% 
                             dplyr::summarise(rhs=paste(est.std,"*",rhs,collapse=" + ")) %>% 
                             dplyr::arrange(desc(op)) %>% 
                             tidyr::unite("mod",lhs,op,rhs,sep="") %>% 
                             dplyr::pull(mod))

ss_mod

#Collapse into one string because my other functions expect that
mod <- base::paste(ss_mod, sep="", collapse="\n") 
mod


########################################
##### Functions Below (Don't edit) #####
########################################

################
### Packages ###
################

library(lavaan)
library(tidyverse)
library(simstandard)
library(tools)
library(patchwork)
library(purrr)
library(stringr)

#######################################
########## HELPER FUNCTIONS ###########
#######################################

#### Function to create model statement without numbers from user model (for input) ####
#Copy from OG
cleanmodel <- function(model){
  
  suppressMessages(model %>%
    lavaan::lavaanify(fixed.x = FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs) %>%
    dplyr::group_by(.data$lhs, .data$op) %>%
    dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    dplyr::arrange(dplyr::desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    dplyr::pull(.data$l))
  
}

#### Function for Number of Factors ####
#Copy from OG

number_factor <- function(model){
  
  #prep the model
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  #isolate factors
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()
  
  #identify number of factors in model
  num_factors <- base::nrow(factors)
  
  return(num_factors)
}

#Did they enter unstandardized loadings?  Aka, do they have any loadings = 1?
#Copy from OG
#Used for error message

unstandardized <- function(model){
  
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1) %>%
    base::nrow()
  
  return(one_plus)
}

#### Function to calculate degrees of freedom ####
# Used for error message

defre <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Rename 
  true_dgm <- model
  
  #Run one simulation
  dat <- simstandard::sim_standardized(true_dgm,n=n,latent=FALSE,errors=FALSE)
  fit <- lavaan::cfa(model=mod,data=dat,std.lv=TRUE)
  
  #Number of freely estimated paths
  paths <- base::max(lavaan::parTable(fit)$free)
  
  #Number of unique values in input matrix
  parms <- base::nrow(lavaan::lavInspect(fit,"std.lv")$theta)
  tot.parms <- (parms*(1+parms))/2
  
  #Subtract
  return(tot.parms-paths)
}

### One-factor: Function to see which items are available ###
## This name is new!!!

one_num <- function(model){
  
  #Rename (just to be consistent with shiny app)
  Mod_C <- model 
  
  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()
  
  #Identify any items that already have an error covariance
  items_covariance <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(-type,type) %>%
    dplyr::select(lhs,op,rhs,type) %>%
    dplyr::filter(op=="=~" | is.na(type)) %>%
    dplyr::filter(is.na(type)) %>%
    dplyr::select(-type) %>%
    tidyr::pivot_longer(-op,names_to = "test", values_to = "rhs") %>%
    dplyr::select(-op,-test) %>%
    dplyr::mutate(lhs=NA,op=NA,ustart=NA)
  
  #Isolate the items that do not already have an error covariance
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup() %>% 
    arrange(abs(ustart))
  
  return(solo_items)
}



#### One-Factor: Function to create misspecification statement ####
#This function name is new!!!!

one_add <- function(model){
  
  #Read in available items
  itemoptions <- one_num(model)
  
  #Count number of available items
  num_i <- base::nrow(itemoptions)
  
  #Select items for misspecification depending on number of available items
  if(num_i==4){
    num_m <- itemoptions %>% 
      dplyr::slice(1:2)
  }else if(num_i==5){
    num_m <- itemoptions %>% 
      dplyr::slice(1:4)
  }else{
    num_m <- itemoptions %>% 
      dplyr::slice(1:(floor(num_i/2)*2))
  }
  
  #Identifiers to separate odds and even rows
  evenindex <- base::seq(2,base::nrow(num_m),2)
  oddindex <- base::seq(1,base::nrow(num_m),2)
  
  #Separate
  left <- num_m[evenindex,]
  right <- num_m[oddindex,] %>% 
    `colnames<-`(c("lhs_1","op_1","rhs_1","ustart_1","n_1"))
  
  #Create misspecification statements
  Residual_Correlation <- base::cbind(left,right) %>% 
    dplyr::mutate(cor=.3,
                  opp="~~",
                  star="*") %>% 
    tidyr::unite(V1,c("rhs","opp","cor","star","rhs_1"),sep=" ") %>% 
    dplyr::select(V1)
  
  return(Residual_Correlation)
}

#### One-factor: Function to create Misspecified DGM ####

DGM_one <- function(model){
  
  #Count number of available items for number of misspecifications
  num_m<- base::nrow(one_num(model))
  
  #Figure out number of levels given number of available items
  if(num_m==4){
    L1 <- 1
    levels <- L1
  }else if(num_m==5){
    L1 <- 1
    L2 <- 2
    levels <- base::rbind(L1,L2)
  }else{
    L3 <- base::floor(num_m/2)
    L2 <- base::floor((2*L3)/3)
    L1 <- base::floor(L3/3)
    levels <- base::rbind(L1,L2,L3)
  }
  
  #Read in misspecifications
  mod <- one_add(model)
  
  #Get parameters for true dgm
  Mods <- model  
  #Mod_C <- base::as.character(Mods$V1) 
  
  #single_mod <- base::lapply(levels, function(x) base::rbind(Mod_C,mod[base::seq(x), ,drop = FALSE]) %>%
  #                       base::data.frame() %>% 
  #                       dplyr::pull(V1))
  #This made you miserable. Shiny was struggling with \n at the end of strings here, for some reason.
  
  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  single_mod <- base::lapply(levels, function(x) base::rbind(Mods,mod[base::seq(x), ,drop = FALSE]) %>%
                               base::data.frame() %>% 
                               dplyr::pull(V1))
  
  return(single_mod)
  
}

#From OG
#Catch regular warning "some estimated ov variances are negative"
#Use in multi_fit_HB function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings
#SPECIFIC TO R PACKAGE

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
}

### One-factor: Simulate fit indices for misspecified model for all levels ###
## This name is new!!!

one_fit <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)
  
  #Use max sample size of 10000
  n <- base::min(n,2000)
  
  #Set seed
  set.seed(649364)
  
  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*50,
                                                                            latent=FALSE,errors=FALSE))
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:50,n)
  
  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))
  
  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>% 
                               tidyr::nest())
  
  #Grab data level of the list
  data <- purrr::map(misspec_data,2)
  
  #Run 500 cfa
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) base::withCallingHandlers(lavaan::cfa(model = mod, data=y, std.lv=TRUE),
                                                                                                  warning=hide_ov)))
  
  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, c("srmr","rmsea","cfi"))) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))
  
  set.seed(NULL)
  
  return(misspec_fit_sum)
  
}

#### One_Factor: Function to create True DGM (aka, just the model the user read in) ####
## This name is new!!

true_fit_one <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model
  
  #Use max sample size of 10000
  n <- base::min(n,2000)
  
  #Set Seed
  set.seed(326267)
  
  #Simulate one large dataset  
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*50,
                                                 latent = FALSE,
                                                 errors = FALSE)
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:50,n)
  
  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)
  
  #Group and list
  true_data <- dat_rep_true %>% 
    dplyr::group_by(rep_id_true) %>% 
    tidyr::nest() %>% 
    base::as.list()
  
  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,~base::withCallingHandlers(lavaan::cfa(model = mod, data=., std.lv=TRUE),
                                                                   warning=hide_ov))
  
  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")
  
  set.seed(NULL)
  
  return(true_fit_sum)
  
}

#### One-Factor: Function to combine both model fit stats for all levels into one dataframe ####
## New name!!

one_df <- function(model,n){
  
  #Use max sample size of 2000
  n <- min(n,2000)
  
  #Get fit stats for misspecified model
  misspec_fit <- one_fit(model,n)
  
  #Get fit stats for correctly specified model
  true_fit <- true_fit_one(model,n)
  
  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))
  
  #Final table
  return(Table)
}


##### NEW: Extract model statement from lavaan object #####
#SPECIFIC TO R PACKAGE

cfa_lavmod <- function(model){
  
  #Extract standardized solution from lavaan object
  lav <- lavaan::standardizedSolution(model)
  
  #Create model statement
  ss_mod <- suppressMessages(lav %>% 
                               dplyr::filter(lhs != rhs) %>% 
                               dplyr::group_by(lhs,op) %>% 
                               dplyr::select(lhs,op,rhs,est.std) %>% 
                               dplyr::mutate(est.std=round(est.std,digits=4)) %>% 
                               dplyr::summarise(rhs=paste(est.std,"*",rhs,collapse=" + ")) %>% 
                               dplyr::arrange(desc(op)) %>% 
                               tidyr::unite("mod",lhs,op,rhs,sep="") %>% 
                               dplyr::pull(mod))
  
  #Collapse into one string because my other functions expect that
  mod <- base::paste(ss_mod, sep="", collapse="\n") 
  
  return(mod)
  
}

##### NEW: Extract n from lavaan object #####
#SPECIFIC TO R PACKAGE

cfa_n <- function(model){
  
  #Extract n from lavaan object
  n <- base::unlist(model@SampleStats@nobs)
  return(n)
}

#############################################
############ cfaOne.R FUNCTION ###############
#############################################

cfaOne <- function(model,n=NULL,plot=FALSE,string=FALSE){
  
  #If string, expect string (a la Shiny app)
  if(string){
    model=model
    n=n
  }else{
  #Use these functions to convert to string (input is a lavaan object)
  #Probably what we should expect for people using R
  #need 'n' first because otherwise model will overwrite  
    n <- cfa_n(model)
    model <- cfa_lavmod(model)
  }
  
  if (unstandardized(model)>0){
    stop("dynamic Error: Your model has loadings greater than or equal to 1 (an impossible value). Please use standardized loadings.")
  }
  
  if (number_factor(model)>1){
    stop("dynamic Error: You entered a multi-factor model.  Use cfaHB instead.")
  }
  
  if (defre(model,n)==0){
    stop("dynamic Error: It is impossible to add misspecifications to a just identified model.")
  }
  
  if ( nrow(one_num(model)) < (number_factor(model)-1)){
    stop("dynamic Error: There are not enough free items to produce all misspecification levels.")
  }
  
  #Create list to store outputs (table and plot)
  res <- list(input=as.list(environment),
              output=list())
  
  #Run simulation  
  results <- one_df(model,n)
  
  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(.05,.1)),
                                                      RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
                                                      CFI_M=stats::quantile(CFI_M, c(.95,.9))))
  
  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
                                                   RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
                                                   CFI_T=stats::quantile(CFI_T, c(.05,.1))))
  
  #Bind each of the misspecified cutoffs to the true cutoffs, listwise
  Table <- purrr::map(misspec_sum,~base::cbind(.,true_sum[[1]]) %>% 
                        dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                                      RMSEA_R=base::round(RMSEA_M,3),
                                      CFI_R=base::round(CFI_M,3),
                                      SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                                      RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                                      CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>% 
                        dplyr::select(SRMR,RMSEA,CFI)) 
  
  #This is to clean up the table for presentation
  #list is a function within mutate to apply function lead across each element
  Row2 <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR_1=SRMR,
                                              RMSEA_1=RMSEA,
                                              CFI_1=CFI) %>%
                           dplyr::mutate_at(c("SRMR_1","RMSEA_1","CFI_1"),base::list(lead)) %>% 
                           dplyr::slice(1) %>% 
                           dplyr::mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR_1,"--"),
                                         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA_1,"--"),
                                         CFI=base::ifelse(base::is.character(CFI),CFI_1,"--"),
                                         SRMR=stringr::str_replace_all(base::as.character(SRMR),"0\\.","."),
                                         RMSEA=stringr::str_replace_all(base::as.character(RMSEA),"0\\.","."),
                                         CFI=stringr::str_replace_all(base::as.character(CFI),"0\\.",".")) %>% 
                           dplyr::select(SRMR,RMSEA,CFI)) 
  
  #Still cleaning
  #Unlist Table
  Table_C <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR=stringr::str_replace_all(base::as.character(SRMR),"0\\.","."),
                                                 RMSEA=stringr::str_replace_all(base::as.character(RMSEA),"0\\.","."),
                                                 CFI=stringr::str_replace_all(base::as.character(CFI),"0\\.",".")))

  #Cleaning  
  Table_C[base::seq(2,nrow(Table_C),by=2),] <- Row2 
  
  #Create row names for level
  Table_C$levelnum <- base::paste("Level", base::rep(1:(base::nrow(Table_C)/2),each=2))
  
  #Create row names for proportions
  Table_C$cut <- base::rep(c("95/5","90/10"))
  
  #Add rownames to final table
  Final_Table <- Table_C %>% 
    tidyr::unite(Cut,levelnum,cut,sep=": ") %>% 
    column_to_rownames(var='Cut')
  
  #Put into list
  res$output$Cutoffs <- Final_Table
  
  #If user selects plot = T
  if(plot){
    
    #For each list element (misspecification) compute the cutoffs
    misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(.05,.1)),
                                                        RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
                                                        CFI_M=stats::quantile(CFI_M, c(.95,.9))))
    
    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
                                                     RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
                                                     CFI_T=stats::quantile(CFI_T, c(.05,.1))))
    
    #Select just those variables and rename columns to be the same
    Misspec_dat <- purrr::map(results,~dplyr::select(.,SRMR_M:Type_M) %>% 
                                `colnames<-`(c("SRMR","RMSEA","CFI","Model")))
    
    #Select just those variables and rename columns to be the same
    True_dat <- purrr::map(results,~dplyr::select(.,SRMR_T:Type_T) %>% 
                             `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #For each element in the list, bind the misspecified cutoffs to the true cutoffs
    #rbind doesn't work well with lists (needs do.call statement)    
    plot <- base::lapply(base::seq(base::length(Misspec_dat)),function(x) dplyr::bind_rows(Misspec_dat[x],True_dat[x]))
    
    #Plot SRMR. Need map2 and data=.x (can't remember why).
    SRMR_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=SRMR,fill=Model))+
                               geom_histogram(position="identity",
                                              alpha=.5, bins=30)+
                               scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                               geom_vline(aes(xintercept=.y$SRMR_M[1],
                                              linetype="misspec_sum$SRMR_M[1]",color="misspec_sum$SRMR_M[1]"),
                                          size=.6)+
                               geom_vline(aes(xintercept=.08,
                                              linetype=".08",color=".08"),
                                          size=.75)+
                               scale_color_manual(name="Cutoff Values",
                                                  labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                  values=c("misspec_sum$SRMR_M[1]"="black",
                                                           ".08"="black"))+
                               scale_linetype_manual(name="Cutoff Values",
                                                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                     values=c("misspec_sum$SRMR_M[1]"="longdash",
                                                              ".08"="dotted"))+
                               theme(axis.title.y = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     panel.background = element_blank(),
                                     axis.line = element_line(color="black"),
                                     legend.position = "none",
                                     legend.title = element_blank(),
                                     legend.box = "vertical"))
    
    #Plot RMSEA.  Need map2 and data=.x (can't remember why).
    RMSEA_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
                                geom_histogram(position="identity",
                                               alpha=.5, bins=30)+
                                scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                                geom_vline(aes(xintercept=.y$RMSEA_M[1],
                                               linetype="misspec_sum$RMSEA_M[1]",color="misspec_sum$RMSEA_M[1]"),
                                           size=.6)+
                                geom_vline(aes(xintercept=.06,
                                               linetype=".06",color=".06"),
                                           size=.75)+
                                scale_color_manual(name="Cutoff Values",
                                                   labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                   values=c("misspec_sum$RMSEA_M[1]"="black",
                                                            ".06"="black"))+
                                scale_linetype_manual(name="Cutoff Values",
                                                      labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                      values=c("misspec_sum$RMSEA_M[1]"="longdash",
                                                               ".06"="dotted"))+
                                theme(axis.title.y = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      panel.background = element_blank(),
                                      axis.line = element_line(color="black"),
                                      legend.position = "none",
                                      legend.title = element_blank(),
                                      legend.box = "vertical"))
    
    #Plot CFI. Need map2 and data=.x (can't remember why).
    CFI_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=CFI,fill=Model))+
                              geom_histogram(position="identity",
                                             alpha=.5, bins=30)+
                              scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                              geom_vline(aes(xintercept=.y$CFI_M[1],
                                             linetype="misspec_sum$CFI_M[1]",color="misspec_sum$CFI_M[1]"),
                                         size=.6)+
                              geom_vline(aes(xintercept=.95,
                                             linetype=".95",color=".95"),
                                         size=.75)+
                              scale_color_manual(name="Cutoff Values",
                                                 labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                 values=c("misspec_sum$CFI_M[1]"="black",
                                                          ".95"="black"))+
                              scale_linetype_manual(name="Cutoff Values",
                                                    labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                    values=c("misspec_sum$CFI_M[1]"="longdash",
                                                             ".95"="dotted"))+
                              theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    panel.background = element_blank(),
                                    axis.line = element_line(color="black"),
                                    legend.position = "none",
                                    legend.title = element_blank(),
                                    legend.box = "vertical"))
    
    #Create a list with the plots combined for each severity level
    plots_combo <- base::lapply(base::seq(base::length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))
    
    #Add a collective legend and title with the level indicator
    plots <- base::lapply(base::seq(base::length(plots_combo)), function(x) patchwork::wrap_plots(plots_combo[[x]])+
                            plot_layout(guides = "collect")+
                            plot_annotation(title=paste("Level", x))
                          & theme(legend.position = 'bottom'))
  
    #Put into list
    res$output$Plots <- plots
    
  }
  
  #Create object (necessary for subsequent print statement)
  class(res) <- 'cfaOne'
  
  return(res)
  
}

#Print suppression/organization statement for list
#Needs same name as class, not function name
print.cfaOne <- function(res){
  
  base::cat("Your DFI cutoffs: \n")
  base::print(res$output$Cutoffs)
  
  if(!is.null(res$output$Plots)){
    
    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(res$output$Plots)
  }
  
  #Hides this function
  base::invisible()
}