model <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .497*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .479*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3

x1 ~~ .55*x12"

n <- 200

cfaHB(model,n,plot=T)

### Packages ###
library(lavaan)
library(tidyverse)
library(simstandard)
library(tools)
library(patchwork)
library(purrr)
library(stringr)

########## HELPER FUNCTIONS ###########

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

### Multi-factor: Function to see which items are available ###
## This name is new!!!

multi_num_HB <- function(model){
  
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
  
  return(itemoptions)
}

#### Multi-Factor: Function to identify available items and loading magnitude ####
#This function name is new!!!!

multi_add_HB <- function(model){
  
  #read in the model
  Mod_C <- model
  
  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  #read in number of factors
  num_fact <- number_factor(model)
  
  #read in viable items from each factor
  itemoptions <- multi_num_HB(model)
  
  #select lowest loading from each factor, in order of magnitude
  crosses <- itemoptions %>% 
    dplyr::group_by(lhs) %>% 
    dplyr::slice_min(base::abs(Loading)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(Priority,base::abs(Loading)) %>% 
    dplyr::slice(1:(num_fact-1))
  
  #identify all factor names (again)
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()
  
  #Compute Coefficient H for each factor
  suppressMessages(Coef_H <- lavaan::lavaanify(Mod_C, fixed.x = FALSE) %>%
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
    `colnames<-`(c("rhs","H")))
  
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
  
  #Isolate items
  dup1 <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
    dplyr::full_join(crosses,by="lhs") %>% 
    dplyr::full_join(Coef_H,by="rhs") %>% 
    dplyr::filter(Item != "NA") %>% 
    dplyr::arrange(abs(Loading))
  
  #Run twice for cleaning later
  dup2 <- dup1
  
  #Manipulate to create model statement
  #Need to add factor correlation statement where lowest comes first
  #So that we can remove from contention once a factor correlation is added
  setup <- base::rbind(dup1,dup2) %>%  
    dplyr::mutate(lhs_1=lhs,
                  rhs_1=rhs,
                  f_min=base::pmin(lhs_1,rhs_1),
                  f_max=base::pmax(lhs_1,rhs_1)) %>% 
    tidyr::unite(facts,c("f_min","f_max")) %>% 
    dplyr::select(-lhs_1,-rhs_1) %>% 
    dplyr::distinct(lhs,op,rhs,ustart,type,Item,Loading,Priority,H,.keep_all = TRUE)
  
  
  #Rename for iteration
  setup_copy <- setup
  
  #Create empty dataframe
  cleaned <- base::data.frame(base::matrix(nrow=0,ncol=10)) %>% 
    `colnames<-`(names(setup)) %>% 
    dplyr::mutate_if(is.logical, as.character)
  
  #Cycle through to grab F-1 misspecifications
  #Select the highest H for first item I (for crossloading)
  #Use anti_join to remove that factor correlation from the list for the next item
  for (i in unique(setup_copy$Item)){
    cleaned[i,] <- setup_copy %>% 
      dplyr::filter(Item==i) %>% 
      dplyr::slice_max(H)
    setup_copy <- dplyr::anti_join(setup_copy,cleaned,by="facts")
  }
  
  #Prep dtaframe for model statement
  modinfo <- cleaned %>% 
    dplyr::mutate(operator="=~",
                  H=base::as.numeric(H),
                  Loading=base::as.numeric(Loading),
                  ustart=base::as.numeric(ustart)) %>% 
    dplyr::arrange(Priority,Loading,-H) 
  
  #Compute maximum allowable cross loading value
  Cross_Loading <- modinfo %>% 
    dplyr::mutate(F1=ustart,
                  F1_Sq=F1^2,
                  L1=Loading,
                  L1_Sq=L1^2,
                  E=1-L1_Sq) %>% 
    dplyr::mutate(MaxAllow=((base::sqrt(((L1_Sq*F1_Sq)+E))-(L1*F1))*.95),
                  MaxAllow2=base::round(MaxAllow,digits=4),
                  Final_Loading=base::pmin(Loading,MaxAllow2),
                  times="*") %>% 
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>% 
    tidyr::unite("V1",sep=" ")
  
  #return value to append to model statement
  return(Cross_Loading)
}


#### Multi-factor: Function to create Misspecified DGM given the number of factors ####
## This name is new!!!

DGM_Multi_HB <- function(model){
  
  mod <- multi_add_HB(model)
  
  #Get parameters for true dgm
  Mods <- model   
  #Mod_C <- base::as.character(Mods$V1) 
  
  #multi_mod <- lapply(mod, function(x) rbind(Mods,mod[seq(x), ,drop = FALSE]) %>%
  #                      data.frame() %>% 
  #                      pull(V1))
  
  #This made you miserable. Shiny was struggling with \n at the end of strings here, for some reason.
  
  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  multi_mod <- lapply(seq(nrow(mod)), function(x) rbind(Mods,mod[seq(x), ,drop = FALSE]) %>%
                        base::data.frame() %>% 
                        dplyr::pull(V1))
  
  return(multi_mod)
  
}

#From OG
#Catch regular warning "some estimated ov variances are negative"
#Use in misspecified_model_fit function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
}

### Multi-factor: Simulate fit indices for misspecified model for all levels ###
## This name is new!!!

multi_fit_HB <- function(model,n){
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for misspecified dgm (this is a list)
  misspec_dgm <- DGM_Multi_HB(model)
  
  #Use max sample size of 2000
  n <- min(n,2000)
  
  #Set seed
  set.seed(269854)
  
  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*50,
                                                                            latent=FALSE,errors=FALSE))
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:50,n)
  
  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))
  
  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>% 
                               tidyr::nest())
  
  #Grab data level of the list
  data <- purrr::map(misspec_data,2)
  
  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) base::withCallingHandlers(lavaan::cfa(model = mod, data=y, std.lv=TRUE),
                                                                                                  warning = hide_ov)))
  
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

#### Multi_Factor: Function to create True DGM (aka, just the model the user read in) ####
## This name is new!!

true_fit_HB <- function(model,n){
  
  #Can make this faster by only doing it once
  #Would need to change table. Not sure what would happen to plot.
  #Already did this
  
  #Get clean model equation
  mod <- cleanmodel(model)
  
  #Get parameters for true DGM
  true_dgm <- model
  
  #Use max sample size of 10000
  n <- base::min(n,2000)
  
  #Set Seed
  set.seed(267326)
  
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

#### Multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####
## New name!!

multi_df_HB <- function(model,n){
  
  #Use max sample size of 2000
  n <- min(n,2000)
  
  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_HB(model,n)
  
  #Get fit stats for correctly specified model
  true_fit <- true_fit_HB(model,n)
  
  #Produce final table of fit indices for each level (as a list)
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))
  
  #Final table
  return(Table)
}


##### NEW: Extract model statement from lavaan object #####

cfa_lavmod <- function(model){
  
  lav <- lavaan::standardizedSolution(model)
  ss_mod <- suppressMessages(lav %>% 
                               dplyr::filter(lhs != rhs) %>% 
                               dplyr::group_by(lhs,op) %>% 
                               dplyr::select(lhs,op,rhs,est.std) %>% 
                               dplyr::mutate(est.std=round(est.std,digits=4)) %>% 
                               dplyr::summarise(rhs=paste(est.std,"*",rhs,collapse=" + ")) %>% 
                               dplyr::arrange(desc(op)) %>% 
                               tidyr::unite("mod",lhs,op,rhs,sep="") %>% 
                               dplyr::pull(mod))
  
  mod <- paste(ss_mod, sep="", collapse="\n") 
  
  return(mod)
  
}

##### NEW: Extract n from lavaan object #####

cfa_n <- function(model){
  n <- unlist(model@SampleStats@nobs)
  return(n)
}

#############################################
############ cfaHB.R FUNCTION ###############
#############################################

cfaHB <- function(model,n=NULL,plot=FALSE,string=F){
  
  res <- list(input=as.list(environment),
              output=list())
  
  if(string){
    model=model
    n=n
  }else{
    n <- cfa_n(model)
    model <- cfa_lavmod(model)
  }
  
  results <- multi_df_HB(model,n)
  
  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=quantile(SRMR_M, c(.05,.1)),
                                                      RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                                                      CFI_M=quantile(CFI_M, c(.95,.9))))
  
  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=quantile(SRMR_T, c(.95,.9)),
                                                   RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                                                   CFI_T=quantile(CFI_T, c(.05,.1))))
  
  #Bind each of the misspecified cutoffs to the true cutoffs, listwise
  Table <- purrr::map(misspec_sum,~cbind(.,true_sum[[1]]) %>% 
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
                           dplyr::mutate_at(c("SRMR_1","RMSEA_1","CFI_1"),list(lead)) %>% 
                           dplyr::slice(1) %>% 
                           dplyr::mutate(SRMR=ifelse(is.character(SRMR),SRMR_1,"--"),
                                         RMSEA=ifelse(is.character(RMSEA),RMSEA_1,"--"),
                                         CFI=ifelse(is.character(CFI),CFI_1,"--"),
                                         SRMR=str_replace_all(as.character(SRMR),"0\\.","."),
                                         RMSEA=str_replace_all(as.character(RMSEA),"0\\.","."),
                                         CFI=str_replace_all(as.character(CFI),"0\\.",".")) %>%
                           dplyr::select(SRMR,RMSEA,CFI)) 
  
  #Still cleaning
  #Unlist Table
  Table_C <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR=stringr::str_replace_all(as.character(SRMR),"0\\.","."),
                                                 RMSEA=stringr::str_replace_all(as.character(RMSEA),"0\\.","."),
                                                 CFI=stringr::str_replace_all(as.character(CFI),"0\\.",".")))
  
  #Cleaning
  Table_C[seq(2,nrow(Table_C),by=2),] <- Row2 
  
  #For row names
  num_fact <- (number_factor(model)-1)
  
  #Create row names for level
  Table_C$levelnum <- paste("Level", rep(1:num_fact,each=2))
  
  #Create row names for proportions
  Table_C$cut <- rep(c("95/5","90/10"))
  
  #Add cross-loading magnitude
  suppressMessages(mag <- multi_add_HB(model) %>% 
    tidyr::separate(V1,into=c("a","b","Magnitude","d","e"),sep=" ") %>% 
    select(Magnitude) %>% 
    mutate(Magnitude=as.numeric(Magnitude),
           Magnitude=round(Magnitude,digits=3)) %>% 
    slice(rep(1:n(),each=2)))
  
  even <- seq_len(nrow(mag))%%2
  mag2 <- cbind(mag,even) %>% 
    mutate(Magnitude=ifelse(even==0," ",Magnitude)) %>%
    mutate(Magnitude=str_replace_all(as.character(Magnitude),"0\\.",".")) %>% 
    select(Magnitude)
  
  #Add to table
  Table_C <- cbind(Table_C,mag2)
      
  #Add rownames to final table
  Final_Table <- Table_C %>% 
    tidyr::unite(Cut,levelnum,cut,sep=": ") %>% 
    tibble::column_to_rownames(var='Cut')
  
  res$output$Cutoffs <- Final_Table
  
  if(plot){
    #For each list element (misspecification) compute the cutoffs
    misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=quantile(SRMR_M, c(.05,.1)),
                                                        RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                                                        CFI_M=quantile(CFI_M, c(.95,.9))))
    
    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=quantile(SRMR_T, c(.95,.9)),
                                                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                                                     CFI_T=quantile(CFI_T, c(.05,.1))))
    
    #Select just those variables and rename columns to be the same
    Misspec_dat <- purrr::map(results,~dplyr::select(.,SRMR_M:Type_M) %>% 
                                `colnames<-`(c("SRMR","RMSEA","CFI","Model")))
    
    #Select just those variables and rename columns to be the same
    True_dat <- purrr::map(results,~dplyr::select(.,SRMR_T:Type_T) %>% 
                             `colnames<-`(c("SRMR","RMSEA","CFI","Model")))
    
    #For each element in the list, bind the misspecified cutoffs to the true cutoffs
    #rbind doesn't work well with lists (needs do.call statement)
    plot <- lapply(seq(length(Misspec_dat)),function(x) dplyr::bind_rows(Misspec_dat[x],True_dat[x]))
    
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
    plots_combo <- lapply(seq(length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))
    
    #Add a collective legend and title with the level indicator
    plots <- lapply(seq(length(plots_combo)), function(x) wrap_plots(plots_combo[[x]])+
                      plot_layout(guides = "collect")+
                      plot_annotation(title=paste("Level", x))
                    & theme(legend.position = 'bottom'))
    
    res$output$Plots <- plots
    
  }
  
  class(res) <- 'cfaHB'
  
  return(res)
  
}

print.cfaHB <- function(res){
  
  cat("Your DFI cutoffs: \n")
  print(res$output$Cutoffs)
  
  if(!is.null(res$output$Plots)){
    
    cat("\n The distributions for each level are in the Plots tab \n")
    print(res$output$Plots)
  }
  
  invisible()
}

cfaHB(model,n,string=T)
cfaHB(obj)

sim_standardized(cfa_lavmod(obj),200)
#rerun function in purrr looks useful

ml <- cfa_lavmod(obj)
ml
nl <- cfa_n(obj)
nl

cfaHB(ml,nl)

class(cleanmodel(model))
class(model_l)
model_l

cleanmodel(ml)
number_factor(ml)
defre(ml,nl)
multi_num_HB(ml)
multi_add_HB(ml)
DGM_Multi_HB(m2)
multi_df_HB(ml,nl)

as.character(ml)

ml
model

m2 <- paste(ml, sep="", collapse="\n") 

cleanmodel(m2)