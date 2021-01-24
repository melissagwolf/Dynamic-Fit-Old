#Input Values#
#this are the only thing that change#
N=487;#sample size
p=5;#number of observed variables
T_ml= 469.29; #model chi-square
df=149;#model degrees of freedom

# everything after this is does not change #

#The formula is from Venables 1975 for obtaining the noncentrality #of a non-central chi-square distribution;
n=N-1;
alpha=.05;
ncp_chi2=function(alpha, T_ml,df){
z=qnorm(1-alpha);
z2=z*z; z3=z2*z; z4=z3*z; z5=z4*z;
sig2=2*(2*T_ml-df+2);
sig=sqrt(sig2); sig3=sig*sig2; sig4=sig2*sig2;sig5=sig4*sig;
sig6=sig2*sig4;

delta=T_ml-df+2+sig*
(
  z+(z2-1)/sig-z/sig2 + 2*(df-1)*(z2-1)/(3*sig3)
  +( -(df-1)*(4*z3-z)/6+(df-2)*z/2 )/sig4
  +4*(df-1)*(3*z4+2*z2-11)/(15*sig5)
  +(
  -(df-1)*(96*z5+164*z3-767*z)/90-4*(df-1)*(df-2)*(2*z3-5*z)/9
      +(df-2)*z/2
  )/sig6
);
delta=max(delta,0);
return(delta)
}


#Raw RMSEA#
delta_c=max(0,T_ml-df);
RMSEA_c=sqrt(delta_c/((N-1)*df));

#T-size RMSEA#;
delta_t=ncp_chi2(alpha, T_ml,df);
RMSEA_t=sqrt(delta_t/(df*(N-1)));

#Recalculate Bins based on Model Characteristics#

RMSEA_e01=exp(
1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n));


RMSEA_e05=exp(2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
+.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n));


RMSEA_e08=exp(2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
+.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n));


RMSEA_e10=exp(2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
+.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n));


cutoff=cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10);
cutoff_3=round(cutoff,3);

class(cutoff_3)

good <- c("Excellent:","Close:","Fair:","Mediocre:","Poor:")

one <- paste(cutoff_3[1],"or below")
two <- paste(cutoff_3[1],"to",cutoff_3[2])
three <- paste(cutoff_3[2],"to",cutoff_3[3])
four <- paste(cutoff_3[3],"to",cutoff_3[4])
five <- paste(cutoff_3[4],"or above")

vals <- rbind(one,two,three,four,five)

as.data.frame(cbind(good,vals))

cat("\n",'Model Raw RMSEA:', RMSEA_c,"\n",'Model T-Size RMSEA:', RMSEA_t,"\n")

obj <- cat('T-Size RMSEA Bins',"\n",'Excellent:',cutoff_3[1], 'or below',"\n",
'Close:', cutoff_3[1], 'to', cutoff_3[2],"\n",
'Fair:', cutoff_3[2], 'to', cutoff_3[3], "\n",
'Mediocre:', cutoff_3[3], 'to', cutoff_3[4], "\n",
'Poor: above', cutoff_3[4], "\n")



