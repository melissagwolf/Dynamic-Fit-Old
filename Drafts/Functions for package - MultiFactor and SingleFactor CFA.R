library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(simstandard)
library(tools)
library(lemon)
library(tibble)
library(magrittr)
library(purrr)

#### Function for Number of Factors ####

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
unstandardized <- function(model){

  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1) %>%
    base::nrow()

  return(one_plus)
}

#### Function to create model statement without numbers from user model (for input) ####

cleanmodel <- function(model){
  
  model %>%
    lavaan::lavaanify(fixed.x = FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs) %>%
    dplyr::group_by(.data$lhs, .data$op) %>%
    dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    dplyr::arrange(dplyr::desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    dplyr::pull(.data$l)
  
}

#### Function to compute DF to know how many misspecifications we can add ####

defre <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for true dgm
  true_dgm <- model
  
  #Run one simulation
  dat <- simstandard::sim_standardized(true_dgm,n=n,latent=FALSE,errors=FALSE)
  fit <- lavaan::cfa(model=mod,data=dat,std.lv=TRUE)
  
  #Number of freely estimated paths
  paths <- max(lavaan::parTable(fit)$free)
  
  #Number of unique values in input matrix
  parms <- nrow(lavaan::lavInspect(fit,"std.lv")$theta)
  tot.parms <- (parms*(1+parms))/2
  
  #Subtract
  return(tot.parms-paths)
}

#### Function for Single Factor Misspecification (Correlation) ####

single_factor <- function(model){

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  #identify the factor name
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
    dplyr::ungroup()

  #Select the item with the lowest loading
  lowest1 <- solo_items %>%
    dplyr::arrange(abs(ustart)) %>%
    dplyr::slice(1)

  #Select the item with the second lowest loading
  lowest2 <- solo_items %>%
    dplyr::arrange(abs(ustart)) %>%
    dplyr::slice(2) %>%
    `colnames<-`(c("lhs2","op2","rhs2","ustart2","n2"))

  #Compute the model implied residual correlation (simstandard)
  U1 <- lowest1$ustart
  U2 <- lowest2$ustart2
  Cor <- U1*U2
  #U1E <- (1-(U1^2))              #for dg by lavaan
  #U2E <- (1-(U2^2))
  #Cov <- Cor*sqrt(U1E*U2E)

  #Create model DF
  Residual_Correlation <- Cor %>%
    base::as.data.frame() %>%
    `colnames<-`("Cor") %>%
    dplyr::mutate(Cor=round(Cor,4)) %>%
    base::cbind(lowest1,lowest2) %>%
    dplyr::mutate(operator="~~") %>%
    dplyr::mutate(times="*") %>%
    dplyr::select(rhs,operator,Cor,times,rhs2) %>%
    tidyr::unite("V1",sep=" ")

  return(Residual_Correlation)
}

#### Function for multi-factor misspecification (Cross-loading) ####

multi_factor <- function(model){

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Compute Coefficient H for each factor
  Coef_H <- lavaan::lavaanify(model, fixed.x = FALSE) %>%
    dplyr::filter(lhs != rhs) %>%
    dplyr::filter(op == "=~") %>%
    dplyr::mutate(L_Sq=ustart^2) %>%
    dplyr::mutate(E_Var=1-L_Sq) %>%
    dplyr::mutate(Div=L_Sq/E_Var) %>%
    dplyr::group_by(lhs) %>%
    dplyr::summarise(Sum=sum(Div)) %>%
    dplyr::mutate(H=((1+(Sum^-1))^-1)) %>%
    dplyr::select(-Sum) %>%
    dplyr::arrange(-H) %>%
    `colnames<-`(c("rhs","H"))

  #Identify number of items per factor
  num_items <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::group_by(lhs) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Original"))

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

  #Isolate the items that do not already have an error covariance or cross-loading
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()

  #Count number of items remaining per factor
  remaining <- solo_items %>%
    dplyr::group_by(lhs) %>%
    dplyr::select(-n) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::full_join(num_items,by="lhs") %>% 
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Remaining","Original"))


  #Add in factor loadings, group by number of items per factor (>2 or 2)
  #And sort factor loadings magnitude within group 
  itemoptions <- solo_items %>%  
    dplyr::full_join(remaining,by="lhs") %>% 
    dplyr::mutate(priority=ifelse(Original>2 & Remaining !="NA","Three","Two")) %>% 
    dplyr::group_by(priority) %>% 
    dplyr::arrange(abs(ustart), .by_group=TRUE) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(lhs,rhs,ustart,priority) %>% 
    base::as.data.frame() %>%
    dplyr::as_tibble() %>%
    `colnames<-`(c("lhs","Item","Loading","Priority"))
  
  #Identify lowest loadings for severity of misspecification
  if(defre(model,n)>2){
    crosses <- itemoptions %>% 
      slice(1:3)
  }else if(defre(model,n)>1){
    crosses <- itemoptions %>% 
      slice(1:2)
  }else if(defre(model,n)>0){
    crosses <- itemoptions %>% 
      slice(1)
  }

  #isolate factors and factor correlations
  factcor1 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::mutate(type=recode(type, .missing ="Error Correlation")) %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor")

  #flip in reverse so we get a list of all factors in one column
  factcor2 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor") %>%
    `colnames<-`(c("rhs","op","lhs","ustart","type")) %>%
    dplyr::select(lhs,op,rhs,ustart,type)

  #combine and clean
  modinfo <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
    dplyr::full_join(crosses,by="lhs") %>% 
    dplyr::full_join(Coef_H,by="rhs") %>% 
    dplyr::filter(Item != "NA") %>% 
    dplyr::group_by(Item) %>% 
    dplyr::arrange(-H, .by_group=TRUE) %>% 
    dplyr::slice(1, .preserve = TRUE) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(operator="=~") %>% 
    dplyr::arrange(Priority,Loading)

  #Compute maximum allowable cross loading value
  F1 <- modinfo$ustart
  F1_Sq <- F1^2                       #.95 = monte carlo correction
  L1 <- modinfo$Loading               #to ensure that the dgm matrix
  L1_Sq <- L1^2                       #is positive definite
  E <- 1-L1_Sq
  MaxAllow <- ((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95)

  #extract value of loading
  Final_Loading <- round(apply(cbind(L1,MaxAllow),1,FUN=min),4)

  #Create model DF
  Cross_Loading <- modinfo %>%
    dplyr::select(rhs,Item,operator) %>%
    base::cbind(Final_Loading) %>%
    dplyr::mutate(times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")

  #return value to append to model statement
  return(Cross_Loading)
}

#### Function to create Misspecified DGM given the number of factors ####

### Need to add in different levels**

Misspecified_DGM_Multi <- function(model){

  if(defre(model,n)>2){
    L1 <- multi_factor(model)[1,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1")) 
    L2 <- multi_factor(model)[2,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1")) 
    L3 <- multi_factor(model)[3,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1")) 
    
    L1_misspec <- rbind(model,L1)
    L2_misspec <- rbind(model,L1,L2)
    L3_misspec <- rbind(model,L1,L2,L3)
    
    L1_misspec_c <- L1_misspec$V1
    L2_misspec_c <- L2_misspec$V1
    L3_misspec_c <- L3_misspec$V1
    
    multi_mod <- list(L1_misspec_c,L2_misspec_c,L3_misspec_c)
  }else if(defre(model,n)>1){
    L1 <- multi_factor(model)[1,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1")) 
    L2 <- multi_factor(model)[2,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1")) 
    
    L1_misspec <- rbind(model,L1)
    L2_misspec <- rbind(model,L1,L2)
    
    L1_misspec_c <- L1_misspec$V1
    L2_misspec_c <- L2_misspec$V1
    
    multi_mod <- list(L1_misspec_c,L2_misspec_c)
  }else if(defre(model,n)>0){
    L1 <- multi_factor(model)[1,] %>% 
      data.frame() %>% 
      `colnames<-`(c("V1"))
    
    L1_misspec <- rbind(model,L1)
    
    L1_misspec_c <- L1_misspec$V1
    
    multi_mod <- list(L1_misspec_c)
  }
  
  return(multi_mod)

}

#Catch regular warning "some estimated ov variances are negative"
#Use in misspecified_model_fit function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
}

#Simulate misspecified data fit stats

#Gonna need if/else statements here, too

misspecified_model_fit <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM_Multi(model)

  #Use max sample size of 10000
  n <- min(n,10000)
  
  #Set seed
  set.seed(649364)
  
  #Simulate one large dataset for each misspecification
  all_data_misspec <- map(misspec_dgm,~sim_standardized(m=.,n=n*500,
                                                        latent=FALSE,errors=FALSE))
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:500,n)
  
  #Combine indicator with dataset
  dat_rep_misspec <- map(all_data_misspec,~cbind(.,rep_id_misspec))
  
  #Group and list
  misspec_data <- map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>% 
                        nest())
  
  #Grab data level of the list
  data <- map(misspec_data,2)
  
  #Run 500 cfa
  misspec_cfa <- base::withCallingHandlers(map(data, function(x) map(x, function(y) 
    cfa(model = mod, data=y, std.lv=TRUE))), warning = hide_ov)
  
  #Extract fit stats from each rep (list) into a data frame and clean
  misspec_fit_sum <- map(misspec_cfa, function(x) map_dfr(x, function(y) fitMeasures(y, c("srmr","rmsea","cfi"))) %>% 
                           `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                           dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#Simulate true data fit stats

#Gonna need if/else statements here, too

true_model_fit <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true dgm
  true_dgm <- model
  
  #Use max sample size of 10000
  n <- min(n,10000)
  
  #Set Seed
  set.seed(326267)
  
  #Simulate one large dataset  
  all_data_true <- sim_standardized(m=true_dgm,n = n*500,
                                       latent = FALSE,
                                       errors = FALSE)
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- rep(1:500,n)
  
  #Combine indicator with dataset
  dat_rep_true <- cbind(all_data_true,rep_id_true)
  
  #Group and list
  true_data <- dat_rep_true %>% 
    group_by(rep_id_true) %>% 
    nest() %>% 
    as.list()
  
  #Run 500 cfa
  true_cfa <- map(true_data$data,~cfa(model = mod, data=., std.lv=TRUE))
  
  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- map_dfr(true_cfa,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
    `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
    dplyr::mutate(Type_M="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#Combine into one dataframe
dynamic_fit <- function(model,n){

  #Probably need some sort of grouping statement
  
  #Use max sample size of 10000
  n <- min(n,10000)

  #Get fit stats for misspecified model
  misspec_fit <- misspecified_model_fit(model,n)

  #Get fit stats for correctly specified model
  true_fit <- true_model_fit(model,n)

  #Create groups
  levels <- rep(1:3,each=500)
  
  #Produce final table by level
  Table <- map_dfr(misspec_fit,~cbind(.,true_fit)) %>% 
    cbind(levels) %>% 
    `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M","Type_M",
                    "SRMR_T","RMSEA_T","CFI_T","Type_T",
                    "levels"))

  #Final table
  return(Table)
}

#Generate dynamic model fit index cutoffs and table
#Need warning for text file
#Need warning for non-standardized loadings

multiCFA <- function(model,n){
  
  #Probably need some sort of grouping statement
  #Gonna need if/else statements here, too

  if (unstandardized(model)>0){
    stop("dynamic Error: Your model has loadings greater than or equal to 1 (an impossible value). Please use standardized loadings.")
  }
  
  if (number_factor(model)<2){
    stop("dynamic Error: You entered a single factor model. Use singleCFA.")
  }

  system.time(results <- dynamic_fit(model,n))

  misspec_sum <- results %>%
    dplyr::group_by(levels) %>% 
    dplyr::summarise(SRMR_M=quantile(SRMR_M, c(.05,.1)),
                     RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                     CFI_M=quantile(CFI_M, c(.95,.9))) %>% 
    ungroup()

  true_sum <- results %>%
    dplyr::group_by(levels) %>% 
    dplyr::summarise(SRMR_T=quantile(SRMR_T, c(.95,.9)),
                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                     CFI_T=quantile(CFI_T, c(.05,.1))) %>% 
    ungroup() %>% 
    select(-levels)

  Table <- cbind(misspec_sum,true_sum) %>%
    dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                  RMSEA_R=base::round(RMSEA_M,3),
                  CFI_R=base::round(CFI_M,3),
                  SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                  RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                  CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
    dplyr::select(SRMR,RMSEA,CFI)


  T_R1 <- Table %>%
    dplyr::slice(1) %>%
    dplyr::rename(SRMR_1=SRMR,
                  RMSEA_1=RMSEA,
                  CFI_1=CFI)

  T_R2 <- Table%>%
    dplyr::slice(2) %>%
    dplyr::rename(SRMR_2=SRMR,
                  RMSEA_2=RMSEA,
                  CFI_2=CFI)

  T_C <- cbind(T_R1,T_R2) %>%
    dplyr::mutate(SRMR=base::ifelse(base::is.character(SRMR_1),SRMR_2,"--"),
                  RMSEA=base::ifelse(base::is.character(RMSEA_1),RMSEA_2,"--"),
                  CFI=base::ifelse(base::is.character(CFI_1),CFI_2,"--")) %>%
    dplyr::select(SRMR,RMSEA,CFI)

  Table_Final <- Table %>%
    dplyr::slice(1) %>%
    base::rbind(T_C) %>%
    dplyr::mutate(Cut=c("95/5","90/10")) %>%
    dplyr::select(Cut,SRMR,RMSEA,CFI) %>%
    tibble::column_to_rownames(var="Cut")

  Misspec_dat <- results %>%
    dplyr::select(SRMR_M:Type_M) %>%
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))

  True_dat <- results %>%
    dplyr::select(SRMR_T:Type_T) %>%
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))

  misspec_sum <- results %>%
    dplyr::summarise(SRMR_M=quantile(SRMR_M, c(.05,.1)),
                     RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                     CFI_M=quantile(CFI_M, c(.95,.9)))

  true_sum <- results %>%
    dplyr::summarise(SRMR_T=quantile(SRMR_T, c(.95,.9)),
                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                     CFI_T=quantile(CFI_T, c(.05,.1)))

  plot <- base::rbind(Misspec_dat,True_dat)

  SRMR_plot <- plot %>%
    ggplot(aes(x=SRMR,fill=Model))+
    geom_histogram(position="identity",
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$SRMR_M[1],
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
          legend.box = "vertical")

  RMSEA_plot <- plot %>%
    ggplot(aes(x=RMSEA,fill=Model))+
    geom_histogram(position="identity",
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$RMSEA_M[1],
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
          legend.box = "vertical")

  CFI_plot <- plot %>%
    ggplot(aes(x=CFI,fill=Model))+
    geom_histogram(position="identity",
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$CFI_M[1],
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
          legend.box = "vertical")

  Plot_Final <- base::suppressMessages(lemon::grid_arrange_shared_legend(SRMR_plot, RMSEA_plot, CFI_plot,
                                    ncol=3,nrow=1))

  return(list(Table_Final, Plot_Final))
}
