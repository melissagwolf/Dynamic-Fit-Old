num_fact <- number_factor(model)

itemoptions <- multi_factor_num(model)

facts <- unique(itemoptions$lhs)

results_f <- numeric(length(facts))

for (i in seq_along(facts)){
  
  
}

test

items <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) %>% 
  dplyr::group_by(Item) %>% 
  dplyr::arrange(-H, .by_group=TRUE) %>% 
  dplyr::slice(1, .preserve = TRUE) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(operator="=~") %>% 
  dplyr::arrange(Priority,Loading)




it <- unique(items$Item)

for(i in it){
  cload[i,] <- items %>% 
    slice_max(H)
  items <- anti_join(items,cload,by="rhs")
}

it <- unique(items$Item)
result <- numeric(length(it))

for(i in seq_along(it)){
  tmp <- items$H[items$Item==it[i]]
  result[i] <- max(tmp[!tmp %in% result])
}


items$H[items$Item==it[2]]