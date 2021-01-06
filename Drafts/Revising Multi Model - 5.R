

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
                     SRMR=as.character(SRMR),
                     RMSEA=as.character(RMSEA),
                     CFI=as.character(CFI)) %>% 
              select(SRMR,RMSEA,CFI)) 

Table_C <- map_dfr(Table,~mutate(.,SRMR=as.character(SRMR),
                                 RMSEA=as.character(RMSEA),
                                 CFI=as.character(CFI)))

Table_C[seq(2,nrow(Table_C),by=2),] <- Row2 

Table_C$level <- paste("Level",Table_C$num)

Table_C$num <- rep(1:nrow(Row2),each=2)

Table_C$cut <- rep(c("95/5","90/10"),nrow(Row2))

Final_Table <- Table_C %>% 
  unite(a,level,num,sep="") %>% 
  unite(Cut,a,cut,sep=": ") %>% 
  column_to_rownames(var='Cut')

Final_Table

##########






model <- "F1 =~ .9*X1 + .9*X2 + .9*X3
F2 =~ .9*X4 + .9*X5 + .9*X6
F1 ~~ .1*F2"

model <- "F1 =~ .8*X1 + .5*X2
F2 =~ .8*X4 + .7*X5
F1 ~~ .3*F2"

defre(model,2000)


defre(model,250)

multi_factor(model)

mod <- cleanmodel(model)

dat <- sim_standardized(m=model,n=250,latent=FALSE,errors=FALSE)

cfa(mod,data=dat)


multiCFA(model,2000)
