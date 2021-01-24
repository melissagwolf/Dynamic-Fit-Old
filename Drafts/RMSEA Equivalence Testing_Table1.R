#Input Values#
#this are the only thing that change#
N=487;#sample size
T_ml= 469.29; #model chi-square
df=149;#model degrees of freedom

# everything after this is does not change #

#The formula is from Venables 1975 for obtaining the noncentrality #of a non-central chi-square distribution;

ncp_chi2 <- function(T_ml,df){
  
  alpha <- .05    
  
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

ncp_chi2(T_ml,df)


rmsea_results <- function(T_ml,df,n){
  
  n=N-1
  alpha=.05
  
  #T-size RMSEA#;
  delta_t=ncp_chi2(T_ml,df)
  RMSEA_t=sqrt(delta_t/(df*n))
  
  #Recalculate Bins based on Model Characteristics#
  
  RMSEA_e01=exp(
    1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n))
  
  
  RMSEA_e05=exp(2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
                +.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n))
  
  
  RMSEA_e08=exp(2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
                +.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n))
  
  
  RMSEA_e10=exp(2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
                +.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n))
  
  
  cutoff=cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10, RMSEA_t)
  cutoff_3=round(cutoff,3)
  
  return(cutoff_3)
}



cutoff_3 <- rmsea_results(T_ml,df,n)

good <- c("Excellent:","Close:","Fair:","Mediocre:","Poor:")

one <- paste(cutoff_3[1],"or below")
two <- paste(cutoff_3[1],"to",cutoff_3[2])
three <- paste(cutoff_3[2],"to",cutoff_3[3])
four <- paste(cutoff_3[3],"to",cutoff_3[4])
five <- paste(cutoff_3[4],"or above")

vals <- rbind(one,two,three,four,five)

tt <- as.data.frame(cbind(good,vals))


#####

cutoff_3

e <- max(cutoff_3[4],cutoff_3[5])
e

x <- c(cutoff_3[1:3],e)
x

m <- e+(cutoff_3[1]-0)
m

ex <- mean(c(.00001,x[1]))
cl <- mean(c(x[1],x[2]))
fa <- mean(c(x[2],x[3]))
me <- mean(c(x[3],x[4]))
po <- mean(c(x[4],m))

ggplot(data.frame(x), aes(x=x, y=0)) +
  geom_point(alpha=0)  +
  annotate("segment",x=0,xend=m, y=0, yend=0, size=1,col="grey50") +
  annotate("segment",x=0,xend=0, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
  annotate("segment",x=cutoff_3[5],xend=cutoff_3[5],y=-0.1,yend=.3, size=1, col="tomato4")+
  annotate("text",x=cutoff_3[5],y=.6,label=paste("T-size \n RMSEA \n",cutoff_3[5]),col="tomato4", size=3.5)+
  annotate("text",x=ex,y=-.5,label="Excellent")+
  annotate("text",x=cl,y=-.5,label="Close")+
  annotate("text",x=fa,y=-.5,label="Fair")+
  annotate("text",x=me,y=-.5,label="Mediocre")+
  annotate("text",x=po,y=-.5,label="Poor")+
  geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=3.5) +
  scale_x_continuous(limits = c(0,m)) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_color_manual(values = unname(colours)) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


######

#The code is to evaluate the cutoff values CFI_e in equivalence testing
#corresponding to the conventional cutoff values of CFI=.99, .95, .92, and .90, respectively;
#as described in the article
#"Confirm Structural Equation Models by Equivalence Testing with Adjusted Fit Indices"
#by Yuan, Chan, Marcoulides and Bentler;

#needed inputs are degrees of freedom (df), sample size (N), and number of observed variables (p);

df=23; N=145; p=9; 




n=N-1; df_i=p*(p-1)/2;

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
#corresponding to R-square=.9713;

cutoff=cbind(CFI_e90, CFI_e92, CFI_e95, CFI_e99);






####

cutoff_3

e <- max(cutoff_3[4],cutoff_3[5])
e

x <- c(0,cutoff_3[1:4],e+(cutoff_3[1]-0))
x


plot(0,xlim=c(x[1],x[6]),axes=FALSE, type = "n", xlab = "", ylab = "")
axis(1, at = x, labels = x)
text(cutoff_3[5],.5,cutoff_3[5],col="blue")
clip((cutoff_3[5]-.01),(cutoff_3[5]+.01),-5,0)
abline(v=cutoff_3[5],col="blue")



