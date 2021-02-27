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

library(lavaan)
library(simstandard)
model <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .497*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .479*x13 + .570*x14 + .628*x15
F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3
x1 ~~ .55*x12"
ns <- 200

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

dat <- sim_standardized(model,ns,latent=F,errors=F)
modclean <- cleanmodel(model)
obj <- cfa(modclean,dat)

########## Run #############

equivTest(obj,plot=T)
equivTest(p,T_ml,df,T_mli,n,string=T,plot=T)

########################################
##### Functions Below (Don't edit) #####
########################################

################
### Packages ###
################

library(stringr)
library(ggplot2)

#######################################
########## HELPER FUNCTIONS ###########
#######################################

##### NCP Chi-Sq #####

equiv_ncp_chi2 <- function(alpha,T_ml,df){
  
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

equiv_cutoffs <- function(p,T_ml,df,T_mli,n){
  
  #Set parms
  df_i <- p*(p-1)/2
  alpha <- .05
  
  #T-size RMSEA#;
  delta_t_r <- equiv_ncp_chi2(alpha,T_ml,df)
  RMSEA_t <- sqrt(delta_t_r/(df*n))
  
  #T-size CFI
  
  delta_t_c <- equiv_ncp_chi2(alpha/2, T_ml,df)
  delta_it <- equiv_ncp_chi2(1-alpha/2, T_mli,df_i)
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

#### Lavaan extraction ####

equiv_n <- function(obj){
  n <- base::unlist(obj@SampleStats@nobs)
  return(n)
}

equiv_p <- function(obj){
  p <- base::unlist(obj@Model@nvar)
  return(p)
}

equiv_T_ml <- function(obj){
  T_ml <- base::unlist(obj@test[["standard"]][["stat"]])
  return(T_ml)
}

equiv_T_mli <- function(obj){
  T_mli <- base::unlist(obj@baseline[["test"]][["standard"]][["stat"]])
  return(T_mli)
}

equiv_df <- function(obj){
  df <- base::unlist(obj@test[["standard"]][["df"]])
  return(df)
}


#############################################
########## equivTest.R FUNCTION #############
#############################################

equivTest <- function(p,T_ml,df,T_mli,n,string=FALSE,plot=FALSE){
  
  #if string, expect string (a la shiny app)
  if(string){
    p=p
    T_ml=T_ml
    df=df
    T_mli=T_mli
    n=n
  }else{
    #Use these functions to convert to string (input is a lavaan object)
    #Probably what we should expect for people using R
    #need 'p' last because otherwise model will overwrite  
    T_ml <- equiv_T_ml(p)
    df <- equiv_df(p)
    T_mli <- equiv_T_mli(p)
    n <- equiv_n(p)
    p <- equiv_p(p)
  }
  
  #Create list to store outputs (table and plot)
  res <- list(input=as.list(environment),
              output=list())
  
  #Get data in
  dat <- equiv_cutoffs(p,T_ml,df,T_mli,n)
  
  #Remove 0's
  clean <- lapply(dat, function(x) stringr::str_replace_all(x,"0\\.","."))
  ul <- unlist(clean)
  cut <- matrix(ul,nrow=2,ncol=5)
  
  ## Extract T-Size and save to list
  rmsea <- cut[1,5]
  res$output$rmsea <- rmsea
  
  cfi <- cut[2,5]
  res$output$cfi <- cfi
  
  ## Label cutoffs
  excellent <- "Excellent: "
  close <- "Close: "
  fair <- "Fair: "
  mediocre <- "Mediocre: "
  poor <- "Poor: "
  
  ## Create RMSEA bins
  one_r <- paste(cut[1,1],"or below")
  two_r <- paste(cut[1,1],"to",cut[1,2])
  three_r <- paste(cut[1,2],"to",cut[1,3])
  four_r <- paste(cut[1,3],"to",cut[1,4])
  five_r <- paste(cut[1,4],"or above")
  
  #Combine
  eo_r <- paste(excellent,one_r,sep="")
  ct_r <- paste(close,two_r,sep="")
  ft_r <- paste(fair,three_r,sep="")
  mf_r <- paste(mediocre,four_r,sep="")
  pf_r <- paste(poor,five_r,sep="")
  
  #Save to list
  res$output$eo_r <- eo_r
  res$output$ct_r <- ct_r
  res$output$ft_r <- ft_r
  res$output$mf_r <- mf_r
  res$output$pf_r <- pf_r
  
  #Create CFI Bins
  one_c <- paste(cut[2,1],"or below")
  two_c <- paste(cut[2,1],"to",cut[2,2])
  three_c <- paste(cut[2,2],"to",cut[2,3])
  four_c <- paste(cut[2,3],"to",cut[2,4])
  five_c <- paste(cut[2,4],"or above")
  
  #Combine
  eo_c <- paste(poor,one_c,sep="")
  ct_c <- paste(mediocre,two_c,sep="")
  ft_c <- paste(fair,three_c,sep="")
  mf_c <- paste(close,four_c,sep="")
  pf_c <- paste(excellent,five_c,sep="")
  
  #Save to list
  res$output$eo_c <- eo_c
  res$output$ct_c <- ct_c
  res$output$ft_c <- ft_c
  res$output$mf_c <- mf_c
  res$output$pf_c <- pf_c
  
  if (plot){
    
    #RMSEA
    e <- max(dat[1,4],dat[1,5])    
    
    x <- dat[1,1:4]
    
    m <- e+(dat[1,1]-0)
    
    ex <- mean(c(.00001,x[1]))
    cl <- mean(c(x[1],x[2]))
    fa <- mean(c(x[2],x[3]))
    me <- mean(c(x[3],x[4]))
    po <- mean(c(x[4],m))
    
    rmsea_plot <- ggplot(data.frame(x), aes(x=x, y=0)) +
      geom_point(alpha=0)  +
      annotate("segment",x=0,xend=m, y=0, yend=0, size=1,col="grey50") +
      annotate("segment",x=0,xend=0, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=dat[1,5],xend=dat[1,5],y=-0.1,yend=.25, size=1, col="tomato4")+
      annotate("text",x=dat[1,5],y=.6,label=paste("T-size \n RMSEA \n",cut[1,5]),
               col="tomato4",size=3.5)+
      annotate("text",x=ex,y=-.5,label="Excellent",size=3.5)+
      annotate("text",x=cl,y=-.5,label="Close",size=3.5)+
      annotate("text",x=fa,y=-.5,label="Fair",size=3.5)+
      annotate("text",x=me,y=-.5,label="Mediocre",size=3.5)+
      annotate("text",x=po,y=-.5,label="Poor",size=3.5)+
      geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=3.5) +
      scale_x_continuous(limits = c(0,m)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) + 
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    
    #CFI
    e <- min(dat[2,1],dat[2,5])
    
    x <- dat[2,1:4]
    
    m <- e-(1-dat[2,4])
    
    ex <- mean(c(1,x[4]))
    cl <- mean(c(x[4],x[3]))
    fa <- mean(c(x[3],x[2]))
    me <- mean(c(x[2],x[1]))
    po <- mean(c(x[1],m))
    
    cfi_plot <- ggplot(data.frame(x), aes(x=x, y=0)) +
      geom_point(alpha=0)  +
      annotate("segment",x=m,xend=1, y=0, yend=0, size=1,col="grey50") +
      annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=dat[2,5],xend=dat[2,5],y=-0.1,yend=.25, size=1, col="tomato4")+
      annotate("text",x=dat[2,5],y=.6,label=paste("T-size \n CFI \n",cut[2,5]),
               col="tomato4", size=3.5)+
      annotate("text",x=ex,y=-.5,label="Excellent",size=3.5)+
      annotate("text",x=cl,y=-.5,label="Close",size=3.5)+
      annotate("text",x=fa,y=-.5,label="Fair",size=3.5)+
      annotate("text",x=me,y=-.5,label="Mediocre",size=3.5)+
      annotate("text",x=po,y=-.5,label="Poor",size=3.5)+
      geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=3.5) +
      scale_x_continuous(limits = c(m,1)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) + 
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())
    
    res$output$R_Plot <- rmsea_plot
    res$output$C_Plot <- cfi_plot
  }

  class(res) <- 'equivTest'
  
  return(res)
}

#Print suppression/organization statement for list
#Needs same name as class, not function name

print.equivTest <- function(res){
  
  cat("\n",
      "Your T-size RMSEA:", 
      res$output$rmsea,
      "\n","\n",
      "Adjusted RMSEA cutoff values: \n",
      res$output$eo_r,"\n",
      res$output$ct_r,"\n",
      res$output$ft_r,"\n",
      res$output$mf_r,"\n",
      res$output$pf_r,
      "\n","\n",
      "Your T-size CFI:", res$output$cfi,
      "\n","\n",
      "Adjusted CFI cutoff values: \n",
      res$output$pf_c,"\n",
      res$output$mf_c,"\n",
      res$output$ft_c,"\n",
      res$output$ct_c,"\n",
      res$output$eo_c,"\n")
  
  if(!is.null(res$output$R_Plot)){
    
    cat("\n",
        "The plots for each fit index are in the Plots tab")
    print(res$output$R_Plot)
    print(res$output$C_Plot)
  }
  
  #Hides this function
  base::invisible()
}
