################
##### Test #####
################

####### String (Normal app input) ########

p <- 36
T_ml <- 31923.829
df <- 639
T_mli <- 229809.310
n <- 12000

####### Lavaan ########

df <- sim_standardized(model,n,latent=F,errors=F)
modclean <- cleanmodel(model)
obj <- cfa(modclean,df)

########## Run #############

cfaHB(model,n,string=T,plot=T)
cfaHB(obj,plot=T)

########################################
##### Functions Below (Don't edit) #####
########################################

################
### Packages ###
################

#######################################
########## HELPER FUNCTIONS ###########
#######################################

##### NCP Chi-Sq #####

ncp_chi2 <- function(alpha,T_ml,df){
  
  z=qnorm(1-alpha)
  z2=z*z 
  z3=z2*z 
  z4=z3*z 
  z5=z4*z
  sig2=2*(2*T_ml-df+2)
  sig=sqrt(sig2)
  sig3=sig*sig2
  sig4=sig2*sig2
  sig5=sig4*sig
  sig6=sig2*sig4
  
  delta=T_ml-df+2+sig*
    (
      z+(z2-1)/sig-z/sig2 + 2*(df-1)*(z2-1)/(3*sig3)
      +( -(df-1)*(4*z3-z)/6+(df-2)*z/2 )/sig4
      +4*(df-1)*(3*z4+2*z2-11)/(15*sig5)
      +(
        -(df-1)*(96*z5+164*z3-767*z)/90-4*(df-1)*(df-2)*(2*z3-5*z)/9
        +(df-2)*z/2
      )/sig6
    )
  
  delta=max(delta,0)
  
  return(delta)
}

######## Compute Values #########

results <- function(p,T_ml,df,T_mli,n){
  
  #Set parms
  df_i <- p*(p-1)/2
  alpha <- .05
  
  #T-size RMSEA#;
  delta_t_r <- ncp_chi2(alpha,T_ml,df)
  RMSEA_t <- sqrt(delta_t_r/(df*n))
  
  #T-size CFI
  
  delta_t_c <- ncp_chi2(alpha/2, T_ml,df)
  delta_it <- ncp_chi2(1-alpha/2, T_mli,df_i)
  CFI_t <- 1-max(delta_t_c,0)/max(delta_t_c,delta_it,0)
  
  #Recalculate Bins based on Model Characteristics - RMSEA#
  
  RMSEA_e01=exp(
    1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n))
  
  
  RMSEA_e05=exp(2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
                +.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n))
  
  
  RMSEA_e08=exp(2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
                +.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n))
  
  
  RMSEA_e10=exp(2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
                +.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n))
  
  ## Recalculate - CFI
  
  CFI_e99=1-exp(
    4.67603-.50827*log(df)+.87087*(df^(1/5))-.59613*((df_i)^(1/5))-1.89602*log(n)
    + .10190*((log(n))^2)+ .03729*log(df)*log(n)
  );
  
  #corresponding to R-square=.9836;
  
  CFI_e95=1-exp(
    4.12132-.46285*log(df)+.52478*(df^(1/5))-.31832*((df_i)^(1/5))-1.74422*log(n)
    +.13042*((log(n))^2)-.02360*(n^(1/2))+.04215*log(df)*log(n)
  );
  
  #corresponding to R-square=.9748;
  
  CFI_e92=1-exp(
    6.31234-.41762*log(df)+.01554*((log(df))^2)-.00563*((log(df_i))^2)-1.30229*log(n)
    +.19999*((log(n))^2)-2.17429*(n^(1/5))+.05342*log(df)*log(n)-.01520*log(df_i)*log(n)
  );
  
  #corresponding to R-square=.9724
  
  CFI_e90=1-exp(
    5.96633-.40425*log(df)+.01384*((log(df))^2)-.00411*((log(df_i))^2)-1.20242*log(n)
    +.18763*((log(n))^2)-2.06704*(n^(1/5))+.05245*log(df)*log(n)-.01533*log(df_i)*log(n)
  );
  
  ## Create bins
  
  cutoff_rmsea <- cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10, RMSEA_t)
  cutoff_cfi <- cbind(CFI_e90, CFI_e92, CFI_e95, CFI_e99, CFI_t)
  cutoff_combo <- rbind(cutoff_rmsea,cutoff_cfi)
  cutoff_3 <- round(cutoff_combo,3)
  colnames(cutoff_3) <- c("Cut_1","Cut_2","Cut_3","Cut_4","T")
  
  return(cutoff_3)
}

equivTest <- function(p,T_ml,df,T_mli,n){
  
  #Create list to store outputs (table and plot)
  res <- list(input=as.list(environment),
              output=list())
  
  #Get data in
  cut <- results(p,T_ml,df,T_mli,n)
  
  ## Extract T-Size and save to list
  rmsea <- cut[1,5]
  res$output$rmsea <- rmsea
  
  cfi <- cut[2,5]
  res$output$cfi <- cfi
  
  ## Label cutoffs
  good <- c("Excellent:","Close:","Fair:","Mediocre:","Poor:")
  
  ## Create RMSEA bins
  one_r <- paste(cut[1,1],"or below")
  two_r <- paste(cut[1,1],"to",cut[1,2])
  three_r <- paste(cut[1,2],"to",cut[1,3])
  four_r <- paste(cut[1,3],"to",cut[1,4])
  five_r <- paste(cut[1,4],"or above")
  
  vals_r <- rbind(one_r,two_r,three_r,four_r,five_r)
  
  bins_r <- paste(good,vals_r, sep=" ")
  
  colnames(bins_r) <- NULL
  rownames(bins_r) <- NULL
  
  #Save to list
  res$output$bins_r <- bins_r
  
  #Create CFI Bins
  one_c <- paste(cut[2,1],"or below")
  two_c <- paste(cut[2,1],"to",cut[2,2])
  three_c <- paste(cut[2,2],"to",cut[2,3])
  four_c <- paste(cut[2,3],"to",cut[2,4])
  five_c <- paste(cut[2,4],"or above")
  
  vals_c <- rbind(five_c,four_c,three_c,two_c,one_c)
  
  bins_c <- paste(good,vals_c, sep=" ")
  
  colnames(bins_c) <- NULL
  rownames(bins_c) <- NULL
  
  #Save to list
  res$output$bins_c <- bins_c
  
  class(res) <- 'equivTest'
  
  return(res)
}

####### Print #########

print.equivTest <- function(res){
  
  ufs::cat0("Your T-size RMSEA: ", res$output$rmsea,
            "\n",
            "Adjusted RMSEA cutoff values: \n",
            res$output$bins_r,
            "\n","\n",
            "Your T-size CFI: ", res$output$cfi,
            "\n",
            "Adjusted CFI cutoff values: \n",
            res$output$bins_c)
  
  #Hides this function
  base::invisible()
}

equivTest(p,T_ml,df,T_mli,n)



library(ufs)

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
############ cfaHB.R FUNCTION ###############
#############################################

cfaHB <- function(model,n=NULL,plot=FALSE,string=FALSE){
  
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
  
  if (number_factor(model)<2){
    stop("dynamic Error: You entered a one-factor model.  Use cfaOne instead.")
  }
  
  if (defre(model,n)==0){
    stop("dynamic Error: It is impossible to add misspecifications to a just identified model.")
  }
  
  if ( nrow(multi_num_HB(model)) < (number_factor(model)-1)){
    stop("dynamic Error: There are not enough free items to produce all misspecification levels.")
  }
  
  #Create list to store outputs (table and plot)
  res <- list(input=as.list(environment),
              output=list())
  
  #Run simulation  
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
  
  #Clean cross-loading magnitude
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
  
  #Put into list
  res$output$Cutoffs <- Final_Table
  
  #If user selects plot = T
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
    plots <- lapply(seq(length(plots_combo)), function(x) patchwork::wrap_plots(plots_combo[[x]])+
                      plot_layout(guides = "collect")+
                      plot_annotation(title=paste("Level", x))
                    & theme(legend.position = 'bottom'))
    
    #Put into list
    res$output$Plots <- plots
    
  }
  
  #Create object (necessary for subsequent print statement)
  class(res) <- 'cfaHB'
  
  return(res)
  
}

#Print suppression/organization statement for list
#Needs same name as class, not function name
print.cfaHB <- function(res){
  
  base::cat("Your DFI cutoffs: \n")
  base::print(res$output$Cutoffs)
  
  if(!is.null(res$output$Plots)){
    
    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(res$output$Plots)
  }
  
  #Hides this function
  base::invisible()
}