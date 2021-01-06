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
  base::as.data.frame() %>%
  `colnames<-`(c("lhs","Remaining"))

remaining

## Probably do something here where we grab the number of items eligible per factor
## Could also remove the lowest loading and then redo the if/then statement without that
## Like, slice (2)

## minor_eligible
## moderate_eligible
## severe_eligible

## num items eligible overall on eligible factors.  Need to check both options.
## But also need to start with df

eligible <- remaining %>%
  dplyr::full_join(num_items,by="lhs") %>%
  dplyr::filter(Original>2 & Remaining != "NA") %>%
  dplyr::mutate(Eligible=1) %>%
  dplyr::select(lhs,Eligible)

eligible

Num_Eligible <- base::nrow(eligible)

Num_Eligible

## If/else for defre

n <- 100

if(defre(model,n)>2){
  if(n>5){
    print(6)
  }else{(print(4))}
}else if(defre(model,n)>1){
  print(2)
}else{
  print(3)
}

class(solo_items)

if(defre(model,n)>2){
  for(i in 1:3){
    
    #Create place to put misspecification
    m <- data.frame(matrix(nrow=3,ncol=5))
    
    if(Num_Eligible>0){
      m[i,] <- solo_items %>%
        dplyr::filter(op=="=~") %>%
        dplyr::select(lhs,op,rhs,ustart) %>%
        dplyr::group_by(rhs) %>%                   #this is where we remove
        dplyr::add_tally() %>%                     #any items that load on
        dplyr::filter(n==1) %>%                    #more than one factor
        dplyr::full_join(eligible,by="lhs") %>%    #add factor eligibility status
        dplyr::filter(Eligible == 1) %>%           #select only factors that are eligible
        dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
        base::as.data.frame() %>%                  #from a factor with more than 2 items
        dplyr::select(-Eligible) %>% 
        dplyr::slice(1) 
    } else {
      crosses <- solo_items %>%
        dplyr::filter(op=="=~") %>%
        dplyr::select(lhs,op,rhs,ustart) %>%
        dplyr::group_by(rhs) %>%                   #this is where we remove
        dplyr::add_tally() %>%                     #any items that load on
        dplyr::filter(n==1) %>%                    #more than one factor
        dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
        base::as.data.frame() %>%                  #from any factor
        dplyr::slice(1) %>%
        dplyr::select(-n,-op) %>%
        dplyr::as_tibble() %>%
        `colnames<-`(c("Factor","Item","Loading"))
    }
  }
  
}



solo_items[-1,]



for(i in 1:3){
  
  #Create place to put misspecification
  m <- data.frame(matrix(nrow=3,ncol=5))
  
  m[i,] <- solo_items %>%
      dplyr::filter(op=="=~") %>%
      dplyr::select(lhs,op,rhs,ustart) %>%
      dplyr::group_by(rhs) %>%                   #this is where we remove
      dplyr::add_tally() %>%                     #any items that load on
      dplyr::filter(n==1) %>%                    #more than one factor
      dplyr::full_join(eligible,by="lhs") %>%    #add factor eligibility status
      dplyr::filter(Eligible == 1) %>%           #select only factors that are eligible
      dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
      base::as.data.frame() %>%                  #from a factor with more than 2 items
      dplyr::select(-Eligible) %>% 
      dplyr::slice(1) 
  }


lowest <- function(model,n){
  solo_items %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    dplyr::group_by(rhs) %>%                   #this is where we remove
    dplyr::add_tally() %>%                     #any items that load on
    dplyr::filter(n==1) %>%                    #more than one factor
    dplyr::full_join(eligible,by="lhs") %>%    #add factor eligibility status
    dplyr::filter(Eligible == 1) %>%           #select only factors that are eligible
    dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
    base::as.data.frame() %>%                  #from a factor with more than 2 items
    dplyr::select(-Eligible) %>% 
    dplyr::slice(1)
}

for(i in 1:3){
  
  #Create place to put misspecification
  m <- data.frame(matrix(nrow=3,ncol=5))
  
  m[i,] <- lowest(model,n) 
  
  n <- m %>% 
    `colnames<-`(c("lhs","op","rhs","ustart","n"))
  
  solo_items <- anti_join(solo_items,n,by=c('lhs','rhs','ustart'))
  
  return(n)
}



m[i,]<- lowest(model,n) 

n <- m %>% 
  `colnames<-`(c("lhs","op","rhs","ustart","n"))

solo_items2 <- anti_join(solo_items,n,by=c('lhs','rhs','ustart'))
solo_items2

m



if(Num_Eligible>0){
  crosses <- solo_items %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    dplyr::group_by(rhs) %>%                   #this is where we remove
    dplyr::add_tally() %>%                     #any items that load on
    dplyr::filter(n==1) %>%                    #more than one factor
    dplyr::full_join(eligible,by="lhs") %>%    #add factor eligibility status
    dplyr::filter(Eligible == 1) %>%           #select only factors that are eligible
    dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
    base::as.data.frame() %>%                  #from a factor with more than 2 items
    dplyr::slice(1) %>%
    dplyr::select(-n,-op, -Eligible) %>%
    dplyr::as_tibble() %>%
    `colnames<-`(c("Factor","Item","Loading"))
} else {
  crosses <- solo_items %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    dplyr::group_by(rhs) %>%                   #this is where we remove
    dplyr::add_tally() %>%                     #any items that load on
    dplyr::filter(n==1) %>%                    #more than one factor
    dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
    base::as.data.frame() %>%                  #from any factor
    dplyr::slice(1) %>%
    dplyr::select(-n,-op) %>%
    dplyr::as_tibble() %>%
    `colnames<-`(c("Factor","Item","Loading"))
}

crosses

factcor1 <- factors %>%
  dplyr::mutate(type="Factor") %>%
  dplyr::full_join(lav_file, by = "lhs") %>%
  dplyr::mutate(type=recode(type, .missing ="Error Correlation")) %>%
  dplyr::select(lhs,op,rhs,ustart,type) %>%
  dplyr::filter(op=="~~" & type=="Factor")

factcor1

factcor2 <- factors %>%
  dplyr::mutate(type="Factor") %>%
  dplyr::full_join(lav_file, by = "lhs") %>%
  dplyr::select(lhs,op,rhs,ustart,type) %>%
  dplyr::filter(op=="~~" & type=="Factor") %>%
  `colnames<-`(c("rhs","op","lhs","ustart","type")) %>%
  dplyr::select(lhs,op,rhs,ustart,type)

factcor2

modinfo <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>%
  base::cbind(crosses) %>%
  dplyr::full_join(Coef_H,by="rhs") %>%
  dplyr::filter(lhs == Factor) %>%
  dplyr::arrange(-H) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(operator="=~")

modinfo


######











V2
df




for(i in 1:3){
  
  V1 <- c(5,6,7,8,9,10)
  df <- data.frame(V1)
  
  V2 <- as.data.frame(matrix(nrow=3,ncol=1))
  maximum <- function(x){
    max(x)
  }
  V2[i,]<- maximum(df)
  df <- anti_join(df,V2,by='V1')
  
  return(V2)
}

V2

for(i in 1:3){
  
  #Create place to put misspecification
  m <- data.frame(matrix(nrow=3,ncol=5))
  
  m[i,] <- lowest(model,n) 
  
  n <- m %>% 
    `colnames<-`(c("lhs","op","rhs","ustart","n"))
  
  solo_items <- anti_join(solo_items,n,by=c('lhs','rhs','ustart'))
  
  return(n)
}


########

model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
       F2 =~ .413*Y5 + -.516*Y6
       F1 ~~ .443*F2
       Y4 ~~ .301*Y5"

mod <- "F1 =~ Y1 + Y2 + Y3 + Y4
        F2 =~ Y5 + Y6
        F1 ~~ F2
        Y4 ~~ Y5"

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





### shouldn't need this
hide_se <- function(z){
  if(any(grepl("Could not compute standard errors!",z)))
    invokeRestart("muffleWarning")
}




