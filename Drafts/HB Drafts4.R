itemoptions



crosses <- itemoptions %>% 
  distinct_at(vars(lhs,Loading), .keep_all = T) %>% 
  group_by(lhs) %>% 
  slice_min(Loading) %>%
  ungroup() %>% 
  arrange(Loading) %>% 
  slice(1:(num_fact-1))

crosses <- itemoptions %>% 
  slice(1:(num_fact-1))

crosses

factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) %>% 
  distinct_at(vars(lhs,rhs), .keep_all=T) %>%
  group_by(Item) %>% 
  slice_max(H) %>% 
  ungroup() %>% 
  mutate(operator="=~") %>% 
  arrange(Priority,Loading,-H)

####


model <- "F1 =~ .43*X1 + .75*X2 + .76*X3
          F2 =~ .44*X4 + .85*X5 + .86*X6 
          F3 =~ .45*X7 + .55*X8 + .56*X9 
          F4 =~ .46*X10 + .65*X11 + .66*X12

          F1 ~~ .5*F2
          F1 ~~ .5*F3
          F1 ~~ .5*F4
          F2 ~~ .5*F3
          F2 ~~ .5*F4
          F3 ~~ .5*F4"  

num_fact <- number_factor(model)

itemoptions <- multi_factor_num(model)

crosses <- itemoptions %>% 
  distinct_at(vars(lhs,Loading), .keep_all = T) %>% 
  group_by(lhs) %>% 
  slice_min(Loading) %>%
  ungroup() %>% 
  arrange(Loading) %>% 
  slice(1:(num_fact-1))

crosses

Coef_H

#F2 to F1
#F1 to F2 < shouldn't happen
#F4 to F2
#F2 to F3

factcor1
factcor2


factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) %>% 
  distinct_at(vars(lhs,rhs), .keep_all=T) %>%
  group_by(Item) %>% 
  slice_max(H) %>% 
  ungroup() %>% 
  mutate(operator="=~") %>% 
  arrange(Priority,Loading,-H)


working <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) 
working 

working %>% 
  separate(lhs, into = c("lhs_T","lhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  separate(rhs, into = c("rhs_T","rhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  select(-lhs_T,-rhs_T) %>% 
  mutate(f_min=pmin(lhs_N,rhs_N),
         f_max=pmax(lhs_N,rhs_N)) %>% 
  unite(facts,c("f_min","f_max")) %>% 
  distinct(facts,.keep_all=T) %>% 
  group_by(Item) %>% 
  slice_max(H)

crosses

Coef_H

#F2 to F1
#F1 to F2 < shouldn't happen
#F4 to F2
#F2 to F3

working

####### Checking again

model <- "F1 =~ .43*X1 + .75*X2 + .76*X3
          F2 =~ .44*X4 + .85*X5 + .86*X6 
          F3 =~ .45*X7 + .95*X8 + .88*X9 
          F4 =~ .46*X10 + .65*X11 + .66*X12

          F1 ~~ .5*F2
          F1 ~~ .5*F3
          F1 ~~ .5*F4
          F2 ~~ .5*F3
          F2 ~~ .5*F4
          F3 ~~ .5*F4"  

crosses
Coef_H

#F3 to F1
#F3 to F2
#F4 to F3

####### Seems to work - NOPE

working %>% 
  separate(lhs, into = c("lhs_T","lhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  separate(rhs, into = c("rhs_T","rhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  select(-lhs_T,-rhs_T) %>% 
  mutate(f_min=pmin(lhs_N,rhs_N),
         f_max=pmax(lhs_N,rhs_N)) %>% 
  unite(facts,c("f_min","f_max")) %>%
  select(-lhs_N,-rhs_N) %>% 
  group_by(facts) %>% 
  slice_max(H) %>% 
  ungroup() %>% 
  group_by(Item) %>% 
  slice_max(H) %>% 
  ungroup() %>% 
  arrange(Loading)

#######


test <- working %>% 
  separate(lhs, into = c("lhs_T","lhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  separate(rhs, into = c("rhs_T","rhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  select(-lhs_T,-rhs_T) %>% 
  mutate(f_min=pmin(lhs_N,rhs_N),
         f_max=pmax(lhs_N,rhs_N)) %>% 
  unite(facts,c("f_min","f_max")) %>%
  select(-lhs_N,-rhs_N)

test

vec <- unique(test$Item)
vec
result <- numeric(length(vec))

result <- matrix(nrow=length(vec),ncol=2)

for (i in seq_along(vec)){
  tmp <- cbind(test$H[test$Item==vec[i]],test$facts[test$Item==vec[i]])
  result[i,] <- cbind(max(tmp[,1][!tmp[,1] %in% result], na.rm=T),tmp[,2][!tmp[,2]%in% result])
}


###

result <- matrix(nrow=3,ncol=2)

for (i in seq_along(vec)){
  tmp[,1] <- test$H[test$Item==vec[i]]
  tmp[,2] <- test$facts[test$Item==vec[i]]
  result[i,] <- c(max(tmp[,1][!tmp[,1] %in% result[,1]], na.rm=T),tmp[,2][!tmp[,2]%in% result[,2]])
}

###

test

result <- matrix(nrow=3,ncol=2)

for (i in seq_along(vec)){
  tmp[,1] <- test$H[test$Item==vec[i]]
  tmp[,2] <- test$facts[test$Item==vec[i]]
  result[i,] <- which.max(tmp[,1][!tmp[,2] %in% result[,2]])
}


#F3 > F1 (.92)
#F3 > F2 (.92)
#F4 > F3 (.64)


#### dplyr attempt

names(test)

results <- data.frame(matrix(nrow=0,ncol=10)) %>% 
  `colnames<-`(names(test)) %>% 
  mutate_if(is.logical, as.character)



for (i in unique(test$Item)){
  results[i,] <- test %>% 
    filter(Item==i) %>% 
    slice_max(H)
  test <- anti_join(test,results,by="facts")
}

########



setup <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) %>% 
  separate(lhs, into = c("lhs_T","lhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  separate(rhs, into = c("rhs_T","rhs_N"),sep = "(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>% 
  select(-lhs_T,-rhs_T) %>% 
  mutate(f_min=pmin(lhs_N,rhs_N),
         f_max=pmax(lhs_N,rhs_N)) %>% 
  unite(facts,c("f_min","f_max")) %>%
  select(-lhs_N,-rhs_N)

setup_copy <- setup

cleaned <- data.frame(matrix(nrow=0,ncol=10)) %>% 
  `colnames<-`(names(setup)) %>% 
  mutate_if(is.logical, as.character)

for (i in unique(setup_copy$Item)){
  cleaned[i,] <- setup_copy %>% 
    filter(Item==i) %>% 
    slice_max(H)
  setup_copy <- anti_join(setup_copy,cleaned,by="facts")
}

modinfo <- cleaned %>% 
  mutate(operator="=~",
         H=as.numeric(H),
         Loading=as.numeric(Loading),
         ustart=as.numeric(ustart)) %>% 
  arrange(Priority,Loading,-H) 

Cross_Loading <- modinfo %>% 
  mutate(F1=ustart,
         F1_Sq=F1^2,
         L1=Loading,
         L1_Sq=L1^2,
         E=1-L1_Sq,
         MaxAllow=((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95),
         Final_Loading=pmin(Loading,MaxAllow),
         times="*") %>% 
  select(rhs,operator,Final_Loading,times,Item) %>% 
  unite("V1",sep=" ")

######



multi_factor(model)


#####





test %>% 
  filter(Item=="X1") %>% 
  slice_max(H) %>% 
  anti_join(test,results,by="facts")


class(results$facts)








tmp

tmp

tmp <- matrix(nrow=3,ncol=2)

tmp[,1] <- test$H[test$Item==vec[1]]
tmp[,2] <- test$facts[test$Item==vec[2]]


tp <- test$H[test$Item==vec[1]]

tc <- cbind(max(tmp[,1][!tmp[,1] %in% result], na.rm=T),tmp[,2][!tmp[,2]%in% result])

tj <- tmp[,2][!tmp[,2]%in% result]
tj

tt <- max(tmp[,1][!tmp[,1] %in% result], na.rm=T)
length(tt)

class(tj)
length(tj)

class(tc)
dim(tc)

class(tp)
dim(tp)

class(tmp)
dim(tmp)

tmp[,1]
###









tmp

tmp$V2

tmp[,2]

for (i in seq_along(vec)){
  tmp <- test$H[test$Item==vec[i]]
  test$load[i] <- max(tmp[!tmp %in% test], na.rm=T)
}


