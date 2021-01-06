model <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4
       F2 =~ .413*Y5 + -.516*Y6
       F1 ~~ .443*F2
       Y4 ~~ .301*Y5"

multi_factor(model)

multi_factor <- function(model){
  
  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)
  
  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()
  
  #Compute Coefficient H for each factor
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
  
  #Identify number of items per factor
  num_items <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::group_by(lhs) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Original"))
  
  #Identify any items that already have an error covariance
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
  
  #Isolate the items that do not already have an error covariance or cross-loading
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()
  
  #Count number of items remaining per factor
  remaining <- solo_items %>%
    dplyr::group_by(lhs) %>%
    dplyr::select(-n) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Remaining"))
  
  #Figure out which items are ideally eligible for cross-loading misspecification
  eligible <- remaining %>%
    dplyr::full_join(num_items,by="lhs") %>%
    dplyr::filter(Original>2 & Remaining != "NA") %>%
    dplyr::mutate(Eligible=1) %>%
    dplyr::select(lhs,Eligible)
  
  #Compute number that are ideally eligible
  Num_Eligible <- base::nrow(eligible)
  
  #Identify item and factor with lowest loading (that doesn't have an existing error cov)
  #We prefer to select a factor with 3+ items for identification purposes
  #If there is an eligible factor with at least 3 items, pick the lowest loading from one of those
  #If there isn't, then you can still proceed
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
      base::as.data.frame() %>%                 #from any factor
      dplyr::slice(1) %>%
      dplyr::select(-n,-op) %>%
      dplyr::as_tibble() %>%
      `colnames<-`(c("Factor","Item","Loading"))
  }
  
  #isolate factors and factor correlations
  factcor1 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::mutate(type=recode(type, .missing ="Error Correlation")) %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor")
  
  #flip in reverse so we get a list of all factors in one column
  factcor2 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor") %>%
    `colnames<-`(c("rhs","op","lhs","ustart","type")) %>%
    dplyr::select(lhs,op,rhs,ustart,type)
  
  #combine and clean
  modinfo <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>%
    base::cbind(crosses) %>%
    dplyr::full_join(Coef_H,by="rhs") %>%
    dplyr::filter(lhs == Factor) %>%
    dplyr::arrange(-H) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(operator="=~")
  
  #Compute maximum allowable cross loading value
  F1 <- modinfo$ustart
  F1_Sq <- F1^2                       #.95 = monte carlo correction
  L1 <- modinfo$Loading               #to ensure that the dgm matrix
  L1_Sq <- L1^2                       #is positive definite
  E <- 1-L1_Sq
  MaxAllow <- ((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95)
  
  #extract value of loading
  Final_Loading <- round(min(L1,MaxAllow),4)
  
  #Create model DF
  Cross_Loading <- modinfo %>%
    dplyr::select(rhs,Item,operator) %>%
    base::cbind(Final_Loading) %>%
    dplyr::mutate(times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")
  
  #return value to append to model statement
  return(Cross_Loading)
}