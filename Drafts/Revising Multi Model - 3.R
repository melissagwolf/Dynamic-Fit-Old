lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
  dplyr::filter(.data$lhs != .data$rhs)

lav_file

factors <- lav_file %>%
  dplyr::filter(op=="=~") %>%
  dplyr::select(lhs) %>%
  base::unique()

factors

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

num_items <- lav_file %>%
  dplyr::filter(op=="=~") %>%
  dplyr::group_by(lhs) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  base::as.data.frame() %>%
  `colnames<-`(c("lhs","Original"))

num_items

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

items_covariance

solo_items <- lav_file %>%
  dplyr::select(lhs,op,rhs,ustart) %>%
  base::rbind(items_covariance) %>%
  dplyr::filter(op=="=~"|is.na(op)) %>%
  dplyr::group_by(rhs) %>%
  dplyr::add_tally() %>%
  dplyr::filter(n==1) %>%
  dplyr::ungroup()

solo_items

remaining <- solo_items %>%
  dplyr::group_by(lhs) %>%
  dplyr::select(-n) %>%
  dplyr::count() %>%
  dplyr::ungroup() %>%
  dplyr::full_join(num_items,by="lhs") %>% 
  base::as.data.frame() %>%
  `colnames<-`(c("lhs","Remaining","Original"))

remaining

## Probably do something here where we grab the number of items eligible per factor
## Could also remove the lowest loading and then redo the if/then statement without that
## Like, slice (2)

## minor_eligible
## moderate_eligible
## severe_eligible

## num items eligible overall on eligible factors.  Need to check both options.
## But also need to start with df

## Begin by separating into two datasets (one with available loadings that were on factors > 2 and one without)

###########

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

itemoptions

#Identify lowest loadings for severity of misspecification
if(defre(model,n)>2){
  crosses <- itemoptions %>% 
    slice(1:3)
}else if(defre(model,n)>1){
  crosses <- itemoptions %>% 
    slice(1:2)
}else{
  crosses <- itemoptions %>% 
    slice(1)
} 

crosses

#isolate factors and factor correlations
factcor1 <- factors %>%
  dplyr::mutate(type="Factor") %>%
  dplyr::full_join(lav_file, by = "lhs") %>%
  dplyr::mutate(type=recode(type, .missing ="Error Correlation")) %>%
  dplyr::select(lhs,op,rhs,ustart,type) %>%
  dplyr::filter(op=="~~" & type=="Factor")

factcor1

#flip in reverse so we get a list of all factors in one column
factcor2 <- factors %>%
  dplyr::mutate(type="Factor") %>%
  dplyr::full_join(lav_file, by = "lhs") %>%
  dplyr::select(lhs,op,rhs,ustart,type) %>%
  dplyr::filter(op=="~~" & type=="Factor") %>%
  `colnames<-`(c("rhs","op","lhs","ustart","type")) %>%
  dplyr::select(lhs,op,rhs,ustart,type)

factcor2

Coef_H

F1 ~~ .443*F2
F1 ~~ .555*F3
F2 ~~ .333*F3

#Combine and clean?
modinfo <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::group_by(Item) %>% 
  dplyr::arrange(-H, .by_group=TRUE) %>% 
  dplyr::slice(1, .preserve = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(operator="=~") %>% 
  dplyr::arrange(Priority,Loading)
  
modinfo

#Compute maximum allowable cross loading value
F1 <- modinfo$ustart
F1_Sq <- F1^2                       #.95 = monte carlo correction
L1 <- modinfo$Loading               #to ensure that the dgm matrix
L1_Sq <- L1^2                       #is positive definite
E <- 1-L1_Sq
MaxAllow <- ((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95)
MaxAllow

Final_Loading <- round(apply(cbind(L1,MaxAllow),1,FUN=min),3)

Final_Loading

#Create model DF
Cross_Loading <- modinfo %>%
  dplyr::select(rhs,Item,operator) %>%
  base::cbind(Final_Loading) %>%
  dplyr::mutate(times="*") %>%
  dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
  tidyr::unite("V1",sep=" ")

Cross_Loading


multi_factor(model)


###################


## If/else for defre

if(defre(model,n)>2){
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
}else if(defre(model,n)>1){
  L1 <- multi_factor(model)[1,] %>% 
    data.frame() %>% 
    `colnames<-`(c("V1")) 
  L2 <- multi_factor(model)[2,] %>% 
    data.frame() %>% 
    `colnames<-`(c("V1")) 
  
  L1_misspec <- rbind(model,L1)
  L2_misspec <- rbind(model,L1,L2)
  
  L1_misspec_c <- L1_misspec$V1
  L2_misspec_c <- L2_misspec$V1
  
  multi_mod <- list(L1_misspec_c,L2_misspec_c)
  }else{
  L1 <- multi_factor(model)[1,] %>% 
    data.frame() %>% 
    `colnames<-`(c("V1"))
  
  L1_misspec <- rbind(model,L1)
  
  L1_misspec_c <- L1_misspec$V1
  
  multi_mod <- list(L1_misspec_c)
  }

###################

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


####################

mod <- cleanmodel(model)

#Get parameters for misspecified dgm
misspec_dgm <- Misspecified_DGM_Multi(model)

#Use max sample size of 10000
n <- min(n,10000)

#Set seed
set.seed(649364)

#Simulate one large dataset for each misspecification

all_data_misspec <- map(misspec_dgm,~sim_standardized(m=.,n=n*10,
                                  latent=FALSE,errors=FALSE))

rep_id_misspec <- rep(1:10,n)

dat_rep_misspec <- map(all_data_misspec,~cbind(.,rep_id_misspec))

rm(misspec_data)

misspec_data <- map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>% 
      nest())

data <- map(misspec_data,2)


misspec_cfa <- base::withCallingHandlers(map(data, function(x) map(x, function(y) cfa(model = mod, data=y, std.lv=TRUE))), warning = hide_ov)

misspec_fit_sum <- map(misspec_cfa, function(x) map_dfr(x, function(y) fitMeasures(y, c("srmr","rmsea","cfi"))) %>% 
                                                          `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                                          dplyr::mutate(Type_M="Misspecified"))
                       
                       
                       
                       

  



#########

mat1 <- replicate(n=10,data.frame(matrix(data=rnorm(20,0,1),nrow=5,ncol=5)),simplify=FALSE)
mat2 <- replicate(n=10,data.frame(matrix(data=rnorm(20,0,1),nrow=5,ncol=5)),simplify=FALSE)
mat3 <- replicate(n=10,data.frame(matrix(data=rnorm(20,0,1),nrow=5,ncol=5)),simplify=FALSE)
combined <- list(mat1,mat2,mat3)

map(combined[[i]],~length(.))

test <- map(combined[[i]],2)

map(combined,length)

#########

test <- rapply(combined, length, how = 'list')

tes <- lapply(combined, function(x) lapply(x, function(y) length(y)))



tes <- map(data, function(x) map(x, function(y) cfa(model = mod, data=y, std.lv=TRUE)))


#Run 500 cfa
misspec_cfa <- base::withCallingHandlers(map(misspec_data$data,~cfa(model = mod, data=., std.lv=TRUE)), warning = hide_ov)



#Extract fit stats from each rep (list) into a data frame and clean
misspec_fit_sum <- map_dfr(misspec_cfa,~fitMeasures(., c("srmr","rmsea","cfi"))) %>% 
  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
  dplyr::mutate(Type_M="Misspecified")

set.seed(NULL)

system.time(results <- misspecified_model_fit(model,500))

system.time(results2 <- true_model_fit(model,500))

########

n <- 500

n <- min(n,10000)

#Get fit stats for misspecified model
misspec_fit <- misspecified_model_fit(model,n)

#Get fit stats for correctly specified model
true_fit <- true_model_fit(model,n)

#Create groups
levels <- rep(1:3,each=500)

#ifelse statements to produce final table
Table <- map_dfr(misspec_fit,~cbind(.,true_fit)) %>% 
                   cbind(levels)

Table

########

cbind(misspec_sum,true_sum)


#########

Table <- cbind(misspec_sum,true_sum) %>%
  dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                RMSEA_R=base::round(RMSEA_M,3),
                CFI_R=base::round(CFI_M,3),
                SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI)

Table

Table <- cbind(misspec_sum,true_sum) %>%
  dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                RMSEA_R=base::round(RMSEA_M,3),
                CFI_R=base::round(CFI_M,3))

L1 <- Table %>% 
  slice(1,2)

L2 <- Table %>% 
  slice(3,4)

L3 <- Table %>% 
  slice(5,6)

L1
L2
L3


L190 <- L1 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  mutate(SRMR2=SRMR,
         RMSEA2=RMSEA,
         CFI2=CFI) %>% 
  mutate_at(c("SRMR2","RMSEA2","CFI2"),list(lead)) %>% 
  slice(1) %>% 
  mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR,"--"),
         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA,"--"),
         CFI=base::ifelse(base::is.character(CFI),CFI,"--")) %>% 
  select(SRMR,RMSEA,CFI)

L290 <- L2 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  mutate(SRMR2=SRMR,
         RMSEA2=RMSEA,
         CFI2=CFI) %>% 
  mutate_at(c("SRMR2","RMSEA2","CFI2"),list(lead)) %>% 
  slice(1) %>% 
  mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR,"--"),
         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA,"--"),
         CFI=base::ifelse(base::is.character(CFI),CFI,"--")) %>% 
  select(SRMR,RMSEA,CFI)

L390 <- L3 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  mutate(SRMR2=SRMR,
         RMSEA2=RMSEA,
         CFI2=CFI) %>% 
  mutate_at(c("SRMR2","RMSEA2","CFI2"),list(lead)) %>% 
  slice(1) %>% 
  mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR,"--"),
         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA,"--"),
         CFI=base::ifelse(base::is.character(CFI),CFI,"--")) %>% 
  select(SRMR,RMSEA,CFI)

L190
L290
L390

L195 <- L1 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  slice(1)

L295 <- L2 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  slice(1)

L395 <- L3 %>% 
  mutate(SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
         RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
         CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  slice(1)

L195
L295
L395

Table_Final <- rbind(L195,L190,L295,L290,L395,L390) %>% 
  dplyr::mutate(Cut=c("Level 1: 95/5","Level 1: 90/10","Level 2: 95/5","Level 2: 90/10","Level 3: 95/5","Level 3: 90/10")) %>%
  dplyr::select(Cut,SRMR,RMSEA,CFI) %>%
  tibble::column_to_rownames(var="Cut")

Table_Final

##Gonna need if/else statements since I'm slicing by length and defre will change length



















#### New (doesn't follow from above table - within multiCFA fuction)
Table2 <- Table %>% 
  mutate(SRMR2=SRMR,
         RMSEA2=RMSEA,
         CFI2=CFI) %>% 
  mutate_at(c("SRMR2","RMSEA2","CFI2"),list(lead)) %>% 
  slice(5) %>% 
  mutate(id=seq(1:3)) %>% 
  rowwise(id) %>% 
  mutate(SRMR1=base::ifelse(base::is.character(SRMR),SRMR2,"--"),
         RMSEA1=base::ifelse(base::is.character(RMSEA),RMSEA2,"--"),
         CFI1=base::ifelse(base::is.character(CFI),CFI2,"--")) %>%
  ungroup() %>% 
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  data.frame()

class(Table2$CFI)



Table %>% 
  mutate(SRMR2=SRMR,
         RMSEA2=RMSEA,
         CFI2=CFI) %>% 
  mutate_at(c("SRMR2","RMSEA2","CFI2"),list(lead)) %>% 
  mutate(id=rep(1:3,each=2)) %>% 
  rowwise(id) %>% 
  mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR2,"--"),
         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA2,"--"),
         CFI=base::ifelse(base::is.character(CFI),CFI2,"--")) %>%
  ungroup() %>% 
  dplyr::select(SRMR,RMSEA,CFI) %>% 
  data.frame()


  
#### New (doesn't follow from above table - within multiCFA fuction)










########

model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
       F2 =~ .413*Y5 + -.516*Y6
       F1 ~~ .443*F2
       Y4 ~~ .301*Y5"

model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
       F2 =~ .413*Y5 + -.516*Y6 + .333*Y11
       F3 =~ .584*Y7 + .683*Y8 + .382*Y9 + .777*Y10
       F1 ~~ .443*F2
       F1 ~~ .555*F3
       F2 ~~ .333*F3
       Y4 ~~ .301*Y5"

n <- 500

lower <- 2

upper <- 3

model <- "F1 =~ .602*Y1 + .805*Y2 + .631*Y4
       F2 =~ .413*Y5 + -.516*Y6
       F1 ~~ .443*F2
       Y4 ~~ .301*Y5"

model <- "F1 =~ .602*Y1 + .805*Y2 + .631*Y4
          F2 =~ .413*Y5 + -.516*Y6 + .666*Y3"

single_factor(model)

mod <- cleanmodel(model)

#Compute DF

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




### Not using

threeplus <- solo_items %>% 
  full_join(factitems,by="lhs") %>% 
  filter(Original>2 & Remaining != "NA") %>% 
  arrange(ustart) %>% 
  select(-n) 

threeplus

two <- solo_items %>% 
  full_join(factitems,by="lhs") %>% 
  filter(Original<3) %>% 
  arrange(ustart) %>% 
  select(-n)

two

#Combine info about number of items available per factor
factitems <- remaining %>%
  dplyr::full_join(num_items,by="lhs") %>%
  dplyr::select(lhs,Remaining,Original)

factitems



