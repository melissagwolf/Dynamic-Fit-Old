

mod <- multi_factor(model,lower,upper)

multi_mod <- lapply(seq(nrow(mod)), function(x) rbind(mod[seq_len(x), ,drop = FALSE],model) %>% 
                pull(V1))

huh


Misspecified_DGM_Multi(model,2,4)


#############

testing <- function(lower=1){
  paste(lower)
}

testing(1:10)


########


trying <- function(model,lower=1,upper=lower){
  
  dif <- (upper+1)-lower
  one <- (upper+1)-upper
  mod <- multi_factor(model,lower,upper)
  
  m <- data.frame(matrix(ncol=1,nrow=dif))
  
  for(i in one:dif){
    
    m[i,] <- mod[i,]
  }
  return(m)
}

trying(model,3,4)

#########

mlist <- list()

trying <- function(model,lower=1,upper=lower){
  
  dif <- (upper+1)-lower
  one <- (upper+1)-upper
  
  for(i in one:dif){
    
    og <- model
    
    mod <- multi_factor(model,lower,upper)
    
    mlist[[i]] <- mod[i,]
  }
  return(mlist)
}

n <- trying(model,3,4)

n


#######

blah <- list()

for(i in 1:3){
  blah[[i]]<- i
}

blah



##########

Misspecified_DGM_Multi <- function(model,lower,upper){
  
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

  return(multi_mod)
  
}

Misspecified_DGM_Multi(model)

class(multi_factor(model,5,7))

defre(model,500)

#########

df <- c("one","two","three","four")
df <- as.data.frame(df)
df

# Coerce df vector of data.frame to character, store as new data.frame: str_df => data.frame 
str_df <- transform(df, df = as.character(df))

# Allocate some memory in order to split data into a list:  df_list => empty list
df_list <- vector("list", nrow(str_df))

# Split the string version of the data.frame into a list as required: 
# df_list => list of character vectors
df_list <- lapply(seq_len(nrow(str_df)), function(i){
  str_df[if(i == 1){1}else{1:i}, grep("df", names(str_df))]
}
)

names(str_df)

df_list

#########

trying2 <- function(model,lower=1,upper=lower){
  
  mod <- multi_factor(model,lower,upper)
  
  str_mod <- transform(mod, mod =as.character(mod$V1))
  
  multi_levels <- vector("list",nrow(str_mod))
  
  multi_levels <- lapply(seq_len(nrow(str_mod)), function(i){
    str_mod[if(i==1){1}else{1:i}, grep("mod", names(str_mod))]
  })
  
  return(multi_levels)
  
}

hm <- trying2(model,1,5)

cleanmodel(hm[[2]])


######

mod <- multi_factor(model,lower,upper)

mod_c <- mod$V1

class(mod_c)

class(mod)

######

huh <- lapply(seq(nrow(mod)), function(x) rbind(mod[seq_len(x), , drop = FALSE],model) %>% 
                pull(V1))

huh



cleanmodel(huh[[3]])

cleanmodel(huh[[3]]$V1)

huh2 <- huh[1]

str(huh)


multi_factor(model,1,2)



######

View(huh)

maybe <- map(huh,~rbind(.,model))

maybe



cleanmodel(maybe[[2]]$V1)

model

mod

seq_len(mod)

function(x) mod[seq_len(x), , drop = FALSE]

lapply(seq(nrow(df)), function(x) df[seq_len(x), , drop = FALSE])

#######


trying <- function(model,lower=1,upper=lower){
  
  dif <- (upper+1)-lower
  one <- (upper+1)-upper
  mod <- multi_factor(model,lower,upper)
  
  m <- data.frame(matrix(ncol=1,nrow=dif))
  
  for(i in one:dif){
    
    m[i,] <- mod[i,]
  }
  return(m)
}

trying(model,3,4)


















