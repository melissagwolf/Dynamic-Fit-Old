T_ml <- 5412.017 
df <- 50
n <- 12752
p <- 12
T_mli <- 67019.332

# everything after this is does not change #

#The formula is from Venables 1975 for obtaining the noncentrality #of a non-central chi-square distribution;

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

##### Combined Function

results <- function(T_ml,df,n,p,T_mli){
  
  n <- n-1
  alpha <- .05
  df_i <- p*(p-1)/2
  
  #T-size RMSEA#;
  delta_t_r <- ncp_chi2(.05,T_ml,df)
  RMSEA_t <- sqrt(delta_t_r/(df*n))
  
  #T-size CFI
  
  delta_t_c <- ncp_chi2(alpha/2, T_ml,df)
  delta_it <- ncp_chi2(1-alpha/2, T_mli,df_i)
  CFI_t=1-max(delta_t_c,0)/max(delta_t_c,delta_it,0)
  
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
  #corresponding to R-square=.9724;
  
  
  CFI_e90=1-exp(
    5.96633-.40425*log(df)+.01384*((log(df))^2)-.00411*((log(df_i))^2)-1.20242*log(n)
    +.18763*((log(n))^2)-2.06704*(n^(1/5))+.05245*log(df)*log(n)-.01533*log(df_i)*log(n)
  );
  
  ##
  
  
  cutoff_rmsea <- cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10, RMSEA_t)
  cutoff_cfi <- cbind(CFI_e90, CFI_e92, CFI_e95, CFI_e99, CFI_t)
  cutoff_combo <- rbind(cutoff_rmsea,cutoff_cfi)
  cutoff_3 <- round(cutoff_combo,3)
  colnames(cutoff_3) <- c("Cut_1","Cut_2","Cut_3","Cut_4","T")
  
  return(cutoff_3)
}

cut <- results(T_ml,df,n,p,T_mli)
cut
cut[1,5]
cut[2,5]
cut

##RMSEA Table

good <- c("Excellent:","Close:","Fair:","Mediocre:","Poor:")

one_r <- paste(cut[1,1],"or below")
two_r <- paste(cut[1,1],"to",cut[1,2])
three_r <- paste(cut[1,2],"to",cut[1,3])
four_r <- paste(cut[1,3],"to",cut[1,4])
five_r <- paste(cut[1,4],"or above")

vals_r <- rbind(one_r,two_r,three_r,four_r,five_r)

as.data.frame(cbind(good,vals_r))

##CFI Table

good <- c("Excellent:","Close:","Fair:","Mediocre:","Poor:")

one_c <- paste(cut[2,1],"or below")
two_c <- paste(cut[2,1],"to",cut[2,2])
three_c <- paste(cut[2,2],"to",cut[2,3])
four_c <- paste(cut[2,3],"to",cut[2,4])
five_c <- paste(cut[2,4],"or above")

vals_c <- rbind(five_c,four_c,three_c,two_c,one_c)

as.data.frame(cbind(good,vals_c))

cut
test <- cut[1,1:4]
test


### RMSEA PLOT
e <- max(cut[1,4],cut[1,5])
e

x <- cut[1,1:4]
x

cut

m <- e+(cut[1,1]-0)
m

ex <- mean(c(.00001,x[1]))
cl <- mean(c(x[1],x[2]))
fa <- mean(c(x[2],x[3]))
me <- mean(c(x[3],x[4]))
po <- mean(c(x[4],m))
po

ggplot(data.frame(x), aes(x=x, y=0)) +
  geom_point(alpha=0)  +
  annotate("segment",x=0,xend=m, y=0, yend=0, size=1,col="grey50") +
  annotate("segment",x=0,xend=0, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=cut[5],xend=cut[5],y=-0.1,yend=.3, size=1, col="tomato4")+
  annotate("text",x=cut[5],y=.6,label=paste("T-size \n RMSEA \n",cut[5]),
           col="tomato4", size=4.5)+
  annotate("text",x=ex,y=-.5,label="Excellent",size=4.5)+
  annotate("text",x=cl,y=-.5,label="Close",size=4.5)+
  annotate("text",x=fa,y=-.5,label="Fair",size=4.5)+
  annotate("text",x=me,y=-.5,label="Mediocre",size=4.5)+
  annotate("text",x=po,y=-.5,label="Poor",size=4.5)+
  geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=4.5) +
  scale_x_continuous(limits = c(0,m)) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_color_manual(values = unname(colours)) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


##CFI

cut

e <- min(cut[2,1],cut[2,5])
e

x <- cut[2,1:4]
x

m <- e-(1-cut[2,4])
m

ex <- mean(c(1,x[4]))
cl <- mean(c(x[4],x[3]))
fa <- mean(c(x[3],x[2]))
me <- mean(c(x[2],x[1]))
po <- mean(c(x[1],m))
ex
cl
fa
me
po
m

x

cut

ggplot(data.frame(x), aes(x=x, y=0)) +
  geom_point(alpha=0)  +
  annotate("segment",x=m,xend=1, y=0, yend=0, size=1,col="grey50") +
  annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=cut[2,5],xend=cut[2,5],y=-0.1,yend=.3, size=1, col="tomato4")+
  annotate("text",x=cut[2,5],y=.6,label=paste("T-size \n CFI \n",cut[2,5]),
           col="tomato4", size=4.5)+
  annotate("text",x=ex,y=-.5,label="Excellent",size=4.5)+
  annotate("text",x=cl,y=-.5,label="Close",size=4.5)+
  annotate("text",x=fa,y=-.5,label="Fair",size=4.5)+
  annotate("text",x=me,y=-.5,label="Mediocre",size=4.5)+
  annotate("text",x=po,y=-.5,label="Poor",size=4.5)+
  geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=4.5) +
  scale_x_continuous(limits = c(m,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_color_manual(values = unname(colours)) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())