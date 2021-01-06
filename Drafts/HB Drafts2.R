model <- "F1 =~ .705*x1 + .445*x2 + .515*x3 + .373*x4 + .497*x5
F2 =~ .489*x4 + .595*x6 + .507*x7 + .559*x8 + .532*x9 + .638*x10
F3 =~ .386*x9 + .546*x11 + .542*x12 + .445*x13 + .570*x14 + .628*x15

F1 ~~ .485*F2
F1 ~~ .657*F3
F2 ~~ .196*F3"

num_fact <- number_factor(model)

itemoptions <- multi_factor_num(model)

facts <- unique(itemoptions$lhs)

results_f <- numeric(length(facts))

for(i in seq_along(it)){
  tmp <- items$H[items$Item==it[i]]
  result[i] <- max(tmp[!tmp %in% result])
}

for (i in seq_along(facts)){
  tmp_f <- itemoptions$Loading[itemoptions$lhs==facts[i]]
  results_f[i] <- min(tmp_f[!tmp_f %in% results_f])
  
}


cbind(unique(itemoptions$lhs)[1:3],results_f)


######

itemoptions

crosses <- itemoptions %>% 
  distinct_at(vars(lhs,Loading), .keep_all = T) %>% 
  group_by(lhs) %>% 
  slice_min(Loading) %>%
  ungroup() %>% 
  arrange(Loading) %>% 
  slice(1:(num_fact-1))

#######


items <- factcor1 %>%
  dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>% 
  dplyr::full_join(crosses,by="lhs") %>% 
  dplyr::full_join(Coef_H,by="rhs") %>% 
  dplyr::filter(Item != "NA") %>% 
  dplyr::arrange(abs(Loading)) %>% 
  distinct_at(vars(rhs,H), .keep_all=T) %>%
  group_by(Item) %>% 
  slice_min(H) %>% 
  ungroup() %>% 
  mutate(operator="=~") %>% 
  arrange(Priority,Loading,-H)


