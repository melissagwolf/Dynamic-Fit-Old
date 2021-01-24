#Test Model

model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5 + .646*Y6 + .732*Y7 + .44*Y8 + .38*Y9 + .66*Y10 + .57*Y11 + .77*Y12"

model <- "F1 =~ .602 * Y1 + .805 * Y2 + .857 * Y3 + .631 * Y4"

n <- 500

single_factor_num(model)
single_factor(model)
Misspecified_DGM_Single(model)
singleCFA(model,n)

######

library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(simstandard)
library(tools)
library(patchwork)
library(tibble)
library(magrittr)
library(purrr)
library(stringr)

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
  
  clean <- model %>%
    lavaan::lavaanify(fixed.x = FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs) %>%
    dplyr::group_by(.data$lhs, .data$op) %>%
    dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    dplyr::arrange(dplyr::desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    dplyr::pull(.data$l)
  
  return(clean)
  
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

#### Function for number of items in the model ####

num_items <- function(model){
  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  #Count number of items
  nums <- lav_file %>% 
    filter(op=="=~")
  
  num <- nrow(nums)
  
  return(num)
}

#### Function for single-factor misspecification cross-loading) ####

single_factor_num <- function(model){
  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
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
    arrange(-abs(ustart))
  
  return(solo_items)
}

single_factor <- function(model){
  
  itemoptions <- single_factor_num(model)
  
  num_i <- nrow(itemoptions)
  
  if(num_i==4){
    num_m <- itemoptions %>% 
      slice(1:2)
  }else if(num_i==5){
    num_m <- itemoptions %>% 
      slice(1:4)
  }else{
    num_m <- itemoptions %>% 
      slice(1:(floor(num_i/2)*2))
  }
  
  evenindex <- seq(2,nrow(num_m),2)
  oddindex <- seq(1,nrow(num_m),2)
  
  left <- num_m[evenindex,]
  right <- num_m[oddindex,] %>% 
    `colnames<-`(c("lhs_1","op_1","rhs_1","ustart_1","n_1"))
  
  Residual_Correlation <- cbind(left,right) %>% 
    mutate(cor=.3,
           opp="~~",
           star="*") %>% 
    unite(V1,c("rhs","opp","cor","star","rhs_1"),sep=" ") %>% 
    select(V1)
  
  return(Residual_Correlation)
}

#### Function to create Misspecified DGM given the number of factors ####

### Need to add in different levels**

Misspecified_DGM_Single <- function(model){

  num_m<- nrow(single_factor_num(model))
  
  if(num_m==4){
    L1 <- 1
    levels <- L1
  }else if(num_m==5){
    L1 <- 1
    L2 <- 2
    levels <- rbind(L1,L2)
  }else{
    L3 <- floor(num_m/2)
    L2 <- floor((2*L3)/3)
    L1 <- floor(L3/3)
    levels <- rbind(L1,L2,L3)
  }
  
  mod <- single_factor(model)
  
  single_mod <- lapply(levels, function(x) rbind(model,mod[seq(x), ,drop = FALSE]) %>%
                         data.frame() %>% 
                         pull(V1))
  
  return(single_mod)

}

#Catch regular warning "some estimated ov variances are negative"
#Use in misspecified_model_fit function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
}

#Simulate misspecified data fit stats

misspec_fit_single <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM_Single(model)

  #Use max sample size of 10000
  n <- min(n,2000)
  
  #Set seed
  set.seed(649364)
  
  #Simulate one large dataset for each misspecification
  all_data_misspec <- map(misspec_dgm,~sim_standardized(m=.,n=n*50,
                                                        latent=FALSE,errors=FALSE))
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:50,n)
  
  #Combine indicator with dataset
  dat_rep_misspec <- map(all_data_misspec,~cbind(.,rep_id_misspec))
  
  #Group and list
  misspec_data <- map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>% 
                        nest())
  
  #Grab data level of the list
  data <- map(misspec_data,2)
  
  #Run 500 cfa
  misspec_cfa <- map(data, function(x) map(x, function(y) 
    base::withCallingHandlers(cfa(model = mod, data=y, std.lv=TRUE), warning = hide_ov)))
  
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
  n <- min(n,2000)
  
  #Set Seed
  set.seed(326267)
  
  #Simulate one large dataset  
  all_data_true <- sim_standardized(m=true_dgm,n = n*50,
                                       latent = FALSE,
                                       errors = FALSE)
  
  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- rep(1:50,n)
  
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
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#Combine into one dataframe
single_df<- function(model,n){

  #Probably need some sort of grouping statement
  
  #Use max sample size of 10000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- misspec_fit_single(model,n)

  #Get fit stats for correctly specified model
  true_fit <- true_model_fit(model,n)

  #Produce final table by level
  Table <- map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

#Generate dynamic model fit index cutoffs and table
#Need warning for text file
#Need warning for non-standardized loadings

singleCFA <- function(model,n){
  
  #Probably need some sort of grouping statement
  #Gonna need if/else statements here, too

  if (unstandardized(model)>0){
    stop("dynamic Error: Your model has loadings greater than or equal to 1 (an impossible value). Please use standardized loadings.")
  }
  
  if (number_factor(model)>1){
    stop("dynamic Error: You entered a multi-factor model. Use multiCFA.")
  }
  
  if (defre(model,n)==0){
    stop("dynamic Error: It is impossible to add misspecifications to a just identified model.")
  }

  if (nrow(single_factor_num(model))<4){
    stop("dynamic Error: There are not enough free items to produce misspecification levels.")
  }

  results <- single_df(model,n)

  misspec_sum <- map(results,~summarise(.,SRMR_M=quantile(SRMR_M, c(.05,.1)),
                                        RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                                        CFI_M=quantile(CFI_M, c(.95,.9))))
  true_sum <- map(results,~dplyr::summarise(.,SRMR_T=quantile(SRMR_T, c(.95,.9)),
                                            RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                                            CFI_T=quantile(CFI_T, c(.05,.1))))
  
  Table <- map(misspec_sum,~cbind(.,true_sum[[1]]) %>% 
                 mutate(SRMR_R=base::round(SRMR_M,3),
                        RMSEA_R=base::round(RMSEA_M,3),
                        CFI_R=base::round(CFI_M,3),
                        SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                        RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                        CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>% 
                 dplyr::select(SRMR,RMSEA,CFI)) 
  
  Row2 <- map_dfr(Table,~mutate(.,SRMR_1=SRMR,
                                RMSEA_1=RMSEA,
                                CFI_1=CFI) %>%
                    mutate_at(c("SRMR_1","RMSEA_1","CFI_1"),list(lead)) %>% 
                    slice(1) %>% 
                    mutate(SRMR=ifelse(is.character(SRMR),SRMR_1,"--"),
                           RMSEA=ifelse(is.character(RMSEA),RMSEA_1,"--"),
                           CFI=ifelse(is.character(CFI),CFI_1,"--"),
                           SRMR=stringr::str_replace_all(as.character(SRMR),"0\\.","."),
                           RMSEA=stringr::str_replace_all(as.character(RMSEA),"0\\.","."),
                           CFI=stringr::str_replace_all(as.character(CFI),"0\\.",".")) %>% 
                    select(SRMR,RMSEA,CFI)) 
  
  Table_C <- map_dfr(Table,~mutate(.,SRMR=stringr::str_replace_all(as.character(SRMR),"0\\.","."),
                                   RMSEA=stringr::str_replace_all(as.character(RMSEA),"0\\.","."),
                                   CFI=stringr::str_replace_all(as.character(CFI),"0\\.",".")))
  
  Table_C[seq(2,nrow(Table_C),by=2),] <- Row2 
  
  Table_C$levelnum <- paste("Level", rep(1:(nrow(Table_C)/2),each=2))
  
  Table_C$cut <- rep(c("95/5","90/10"))
  
  Final_Table <- Table_C %>% 
    unite(Cut,levelnum,cut,sep=": ") %>% 
    column_to_rownames(var='Cut') 
  
  ##Plots
  
  repid1 <- rep(1:length(results),each=50)
  repid2 <- rep(1:length(results),each=50)
  repid <- c(repid1,repid2)
  
  df <- do.call(rbind.data.frame,results)
  
  df_1 <- df[,c(1:4)] %>% 
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))
  
  df_2 <- df[,c(5:8)] %>% 
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))
  
  df2 <- rbind(df_1,df_2)
  
  df3 <- cbind(df2,repid)
  
  rm(SRMR_plot)
  
  SRMR_plot <- df3 %>%
    group_by(repid) %>% 
    ggplot(aes(x=SRMR,fill=Model))+
    geom_histogram(position="identity",
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=SRMR[1],
                   linetype="SRMR[1]",color="SRMR[1]"),
               size=.6)+
    geom_vline(aes(xintercept=.08,
                   linetype=".08",color=".08"),
               size=.75)+
    scale_color_manual(name="Cutoff Values",
                       labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                       values=c("SRMR[1]"="black",
                                ".08"="black"))+
    scale_linetype_manual(name="Cutoff Values",
                          labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                          values=c("SRMR[1]"="longdash",
                                   ".08"="dotted"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color="black"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.box = "vertical")
  
  SRMR_Mean <- df3 %>% 
    group_by(repid,Model) %>% 
    summarise(mean=mean(SRMR)) %>% 
    filter(Model=="Misspecified")
  
  SRMR_plot
  
  
  RMSEA_plot <- map2(plot,misspec_sum,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
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
  
  CFI_plot <- map2(plot,misspec_sum,~ggplot(data=.x,aes(x=CFI,fill=Model))+
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
  
  
  plots_combo <- lapply(seq(length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))
  
  Final_Plots <- suppressMessages(lapply(seq(length(plots_combo)), function(x) wrap_plots(plots_combo[[x]])+
           plot_layout(guides = "collect")+
           plot_annotation(title=paste("Level", x))
         & theme(legend.position = 'bottom')))
  
  return(list(Final_Table, Final_Plots))
}

