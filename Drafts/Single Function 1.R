j <- seq(4,30)
j


#Need to change for 4 and 5
l3 <- floor(j/2)
l3

l2 <- floor((2*l3)/3)
l2

l1 <- round(l3/3)
l1

cbind(j,l1,l2,l3)

###

itemoptions <- single_factor_num(model)

num_i <- nrow(itemoptions)

if(num_i==4){
  L1 <- 1
  
  num_m <- itemoptions %>% 
    slice(1:2)
  
}else if(num_i==5){
  L1 <- 1
  L2 <- 2
  
  num_m <- itemoptions %>% 
    slice(1:4)
  
}else{
  L3 <- floor(num_i/2)
  L2 <- floor((2*L3)/3)
  L1 <- floor(L3/3)
  
  num_m <- itemoptions %>% 
    slice(1:(L3*2))
}

###

if(num_i==4){
  num_m <- itemoptions %>% 
    slice(1:2)
}else if(num_i==5){
  num_m <- itemoptions %>% 
    slice(1:4)
}else{
  num_m <- itemoptions %>% 
    slice(1:(floor(num_i/2)*2))
}

evenindex <- seq(2,nrow(num_m),2)
oddindex <- seq(1,nrow(num_m),2)

left <- num_m[evenindex,]
right <- num_m[oddindex,] %>% 
  `colnames<-`(c("lhs_1","op_1","rhs_1","ustart_1","n_1"))

Residual_Correlation <- cbind(left,right) %>% 
  mutate(cor=.35,
         opp="~~",
         star="*") %>% 
  unite(V1,c("rhs","opp","cor","star","rhs_1"),sep=" ") %>% 
  select(V1)

Residual_Correlation

single_factor(model)

#######

num_m<- nrow(single_factor_num(model))

if(num_m==4){
  L1 <- 1
  levels <- L1
}else if(num_m==5){
  L1 <- 1
  L2 <- 2
  levels <- rbind(L1,L2)
}else{
  L3 <- floor(num_m/2)
  L2 <- floor((2*L3)/3)
  L1 <- floor(L3/3)
  levels <- rbind(L1,L2,L3)
}

mod <- single_factor(model)

single_mod <- lapply(levels, function(x) rbind(model,mod[seq(x), ,drop = FALSE]) %>%
                      data.frame() %>% 
                      pull(V1))

single_mod

Misspecified_DGM_Single(model)

misspec_fit_single(model,500)