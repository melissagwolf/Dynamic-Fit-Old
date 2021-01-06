#OG

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

multi_mod


#Lower/Upper

lower <- 2
upper <- 3

mod2 <- multi_factor(model,lower,upper)

multi_mod2 <- lapply(seq(nrow(mod2)), function(x) rbind(model,mod2[seq_len(x), ,drop = FALSE]) %>%
                      pull(V1))

multi_mod_final<- multi_mod2[lower:upper]

Misspecified_DGM_Multi(model,2,3)

multi_factor(model)

Misspecified_DGM_Multi(model,2,3)

system.time(misspecified_model_fit(model,500,1,3))
#user  system elapsed 
#102.04    2.03  105.42

system.time(misspecified_model_fit(model,500))
#user  system elapsed 
#102.20    2.68  106.17 

system.time(misspecified_model_fit(model,500,2,3))
#user  system elapsed 
#67.56    1.33   69.92 

multiCFA(model,n,lower,upper)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` ungrouping output (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#SRMR RMSEA   CFI
#Level 2: 95/5   NONE  NONE  NONE
#Level 2: 90/10  NONE  NONE  NONE
#Level 3: 95/5  0.036 0.038 0.986
#Level 3: 90/10 0.039 0.042 0.982

system.time(multiCFA(model,n,lower,upper))
#user  system elapsed 
#64.23    1.49   66.54 

multiCFA(model,n)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` ungrouping output (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` ungrouping output (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` ungrouping output (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#SRMR RMSEA   CFI
#Level 1: 95/5   NONE  NONE  NONE
#Level 1: 90/10  NONE  NONE  NONE
#Level 2: 95/5   NONE  NONE  NONE
#Level 2: 90/10  NONE  NONE  NONE
#Level 3: 95/5  0.037 0.039 0.984
#Level 3: 90/10    --    --    --

system.time(multiCFA(model,n))
#user  system elapsed 
#93.98    2.08   97.38 

##Seed fix##

multiCFA(model,n,lower,upper)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#`summarise()` ungrouping output (override with `.groups` argument)
#`summarise()` regrouping output by 'lhs' (override with `.groups` argument)
#SRMR RMSEA   CFI
#Level 2: 95/5   NONE  NONE  NONE
#Level 2: 90/10  NONE  NONE  NONE
#Level 3: 95/5  0.037 0.039 0.984
#Level 3: 90/10 0.039 0.043 0.981

system.time(multiCFA(model,n,lower,upper))
#user  system elapsed 
#66.37    1.63   69.00 