############ Website Model

model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
        F2 =~ .413*Y5 + -.516*Y6
        F1 ~~ .443*F2
        Y4 ~~ .301*Y5"

cleanmod <- cleanmodel(model)

multi_factor(model,4)

Misspecified_DGM_Multi(model,4)

misspec1 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
             F2 =~ .413*Y5 + -.516*Y6
             F1 ~~ .443*F2        
             Y4 ~~ .301*Y5
             F2 =~ 0.5464 * Y1"

misspec3 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
             F2 =~ .413*Y5 + -.516*Y6
             F1 ~~ .443*F2
             Y4 ~~ .301*Y5
             F2 =~ 0.5464 * Y1
             F2 =~ 0.3188 * Y2                                                                                                          
             F2 =~ 0.2474 * Y3"  

misspec4 <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
             F2 =~ .413*Y5 + -.516*Y6
             F1 ~~ .443*F2
             Y4 ~~ .301*Y5
             F2 =~ 0.5464 * Y1
             F2 =~ 0.3188 * Y2                                                                                                          
             F2 =~ 0.2474 * Y3
             F1 =~ -0.516 * Y6"

defre(model,500)

############# H&B Model

model <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .479*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .497*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3"

coef

cleanmod <- cleanmodel(model)



Misspecified_DGM_Multi(model,4)

defre(model,500)

misspec1 <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .479*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .497*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3

F2 =~ 0.445 * x2"

misspec2 <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .479*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .497*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3

F2 =~ 0.445 * x2
F3 =~ .479 * x5"

misspec4 <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .497*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .479*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3

F2 =~ 0.445 * x2
F2 =~ 0.479 * x13
F2 =~ 0.497 * x5
F3 =~ 0.507 * x7"

#############

m1_all <- sim_standardized(misspec1,500*500,latent=FALSE,errors=FALSE)
m2_all <- sim_standardized(misspec2,500*500,latent=FALSE,errors=FALSE)
m4_all <- sim_standardized(misspec4,500*500,latent=FALSE,errors=FALSE)

repid <- rep(1:500,500)

m1_rep <- cbind(m1_all,repid)
m2_rep <- cbind(m2_all,repid)
m4_rep <- cbind(m4_all,repid)

m1_dat <- m1_rep %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

m2_dat <- m2_rep %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

m4_dat <- m4_rep %>% 
  group_by(repid) %>% 
  nest() %>% 
  as.list()

m1_dat2 <- m1_dat[[2]]
m2_dat2 <- m2_dat[[2]]
m4_dat2 <- m4_dat[[2]]

m1_cfa <- map(m1_dat2,~cfa(cleanmod,data=.,std.lv=TRUE))
m2_cfa <- map(m2_dat2,~cfa(cleanmod,data=.,std.lv=TRUE))
m4_cfa <- map(m4_dat2,~cfa(cleanmod,data=.,std.lv=TRUE))

m1_fit <- map_dfr(m1_cfa,~fitmeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
m2_fit <- map_dfr(m2_cfa,~fitmeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))
m4_fit <- map_dfr(m4_cfa,~fitmeasures(.,c("srmr","rmsea","cfi","fmin","chisq","df","pvalue")))

library(psych)

d1 <- describe(m1_fit,quant=c(.05,.10,.90,.95))
d2 <- describe(m2_fit,quant=c(.05,.10,.90,.95))

print(d1,digits=3)
print(d2,digits=3)


describe(m4_fit)

m1_fit %>% 
  ggplot(aes(x=srmr))+
  geom_histogram()

m4_fit %>% 
  ggplot(aes(x=srmr))+
  geom_histogram()

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

Coef_H

#.733 and .712 H going into F1
#.719 and .676 H into F





#########

mismod <- Misspecified_DGM_Multi(model)

mismod

sim_standardized(modmod,500*500,latent=FALSE,errors=FALSE)

modmod <- "F1 =~ .602*Y1 +.805*Y2 + .857*Y3 + .631*Y4
F2 =~ .413*Y5 + -.516*Y6 + .602*Y1
F1 ~~ .443*F2
Y4 ~~ .301*Y5"


webmod <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
F2 =~ .413*Y5 + -.516*Y6
F1 ~~ .443*F2
Y4 ~~ .301*Y5
F2 =~ .602*Y1"

((sqrt((.602^2*.443^2)+(1-(.602^2)))) - (.602*.443))*.95

modinfo

F1 <- modinfo$ustart
F1_Sq <- F1^2                       #.95 = monte carlo correction
L1 <- modinfo$Loading               #to ensure that the dgm matrix
L1_Sq <- L1^2                       #is positive definite
E <- 1-L1_Sq
MaxAllow <- ((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95)


sqrt((L1_Sq*F1_Sq)+E-L1*F1)*.95