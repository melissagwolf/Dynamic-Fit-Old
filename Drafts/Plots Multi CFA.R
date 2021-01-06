##Plots

Misspec_dat <- map(results,~select(.,SRMR_M:Type_M) %>% 
      `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

True_dat <- map(results,~select(.,SRMR_T:Type_T) %>% 
      `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

plot <- lapply(seq(length(Misspec_dat)),function(x) bind_rows(Misspec_dat[x],True_dat[x]))

######

SRMR_plot <- map2(plot,misspec_sum,~ggplot(data=.x,aes(x=SRMR,fill=Model))+
                    geom_histogram(position="identity",
                                   alpha=.5)+
                    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                    geom_vline(aes(xintercept=.y$SRMR_M[1],
                                   linetype="misspec_sum$SRMR_M[1]",color="misspec_sum$SRMR_M[1]"),
                               size=.6)+
                    geom_vline(aes(xintercept=.08,
                                   linetype=".08",color=".08"),
                               size=.75)+
                    scale_color_manual(name="Cutoff Values",
                                       labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                       values=c("misspec_sum$SRMR_M[1]"="black",
                                                ".08"="black"))+
                    scale_linetype_manual(name="Cutoff Values",
                                          labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                          values=c("misspec_sum$SRMR_M[1]"="longdash",
                                                   ".08"="dotted"))+
                    theme(axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(color="black"),
                          legend.position = "none",
                          legend.title = element_blank(),
                          legend.box = "vertical"))


RMSEA_plot <- map2(plot,misspec_sum,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
                    geom_histogram(position="identity",
                                   alpha=.5)+
                    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                    geom_vline(aes(xintercept=.y$RMSEA_M[1],
                                   linetype="misspec_sum$RMSEA_M[1]",color="misspec_sum$RMSEA_M[1]"),
                               size=.6)+
                    geom_vline(aes(xintercept=.06,
                                   linetype=".06",color=".06"),
                               size=.75)+
                    scale_color_manual(name="Cutoff Values",
                                       labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                       values=c("misspec_sum$RMSEA_M[1]"="black",
                                                ".06"="black"))+
                    scale_linetype_manual(name="Cutoff Values",
                                          labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                          values=c("misspec_sum$RMSEA_M[1]"="longdash",
                                                   ".06"="dotted"))+
                    theme(axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          panel.background = element_blank(),
                          axis.line = element_line(color="black"),
                          legend.position = "none",
                          legend.title = element_blank(),
                          legend.box = "vertical"))

CFI_plot <- map2(plot,misspec_sum,~ggplot(data=.x,aes(x=CFI,fill=Model))+
                     geom_histogram(position="identity",
                                    alpha=.5)+
                     scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                     geom_vline(aes(xintercept=.y$CFI_M[1],
                                    linetype="misspec_sum$CFI_M[1]",color="misspec_sum$CFI_M[1]"),
                                size=.6)+
                     geom_vline(aes(xintercept=.95,
                                    linetype=".95",color=".95"),
                                size=.75)+
                     scale_color_manual(name="Cutoff Values",
                                        labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                        values=c("misspec_sum$CFI_M[1]"="black",
                                                 ".95"="black"))+
                     scale_linetype_manual(name="Cutoff Values",
                                           labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                           values=c("misspec_sum$CFI_M[1]"="longdash",
                                                    ".95"="dotted"))+
                     theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           panel.background = element_blank(),
                           axis.line = element_line(color="black"),
                           legend.position = "none",
                           legend.title = element_blank(),
                           legend.box = "vertical"))


plots_combo <- lapply(seq(length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))

plots_combo[[1]]



###########

for (i in seq_along(plots_combo)){
  a[[i]] <- plots_combo[[i]]
}


##########


plots_combo[1]


Plot_Final <- lapply(seq(length(plots_combo)),function(x) grid_arrange_shared_legend(plots_combo[x],ncol=3,nrow=1))


grid_arrange_shared_legend(plots_combo[1],ncol=3,nrow=1)


plots_combo[1]



####
testing <- c(SRMR_plot[1],RMSEA_plot[1],CFI_plot[1])






#######

map(plot,~ggplot(data=.,aes(x=SRMR,fill=Model))+
      geom_histogram(position="identity",
                     alpha=.5)+
      scale_fill_manual(values=c("#E9798C","#66C2F5")))

plot

x <- list(1, 1, 1)
y <- list(10, 20, 30)
z <- list(100, 200, 300)

map2(x, y, ~ .x + .y)


######

p <- plot[[1]]

ggplot(p,(aes(x=SRMR)))+
         geom_histogram()


name <- names(plot[[1]])[1:3]
name

plotting <- function(index){
  ggplot(p, aes_string(x=index,fill="Model"))+
           geom_histogram()
}

plotting("SRMR")

for(i in seq_along(plot)){
  um <- map(name,~plotting(.x))
}

um

########


plotting <- function(index){
  ggplot(plot[[i]], aes_string(x=index,fill="Model"))+
    geom_histogram()
}

seq_along(plot)

for(i in seq_along(plot)){
  um2 <- ggplot(plot[[i]], aes(x=SRMR, fill=Model))+
    geom_histogram()
}

um2

UM3 <- map(plot, function(x) map(x, function(y) ggplot(y, aes(x=SRMR, fill=Model))+
                                  geom_histogram()))



class(plot[[1]])





length(plot)O

seq(length(plot))


ii <- map(name,~plotting(.x))

####

map(name,~plotting(.x))



####



map(misspec_cfa, function(x) map_dfr(x, function(y) fitMeasures(y, c("srmr","rmsea","cfi"))))

lapply(combined, function(x) lapply(x, function(y) length(y)))

lapply(plot, function(x) lapply(x, function(y) name,~plotting(y)))






       
       
       
       
class(p$SRMR)

p <- plot[[1]]

p %>% ggplot(aes(x=SRMR))+
  geom_histogram()

plot[[1]]$SRMR







#####

lapply(1, funcion(x) misspec_sum[[x]]$SRMR_M[1])

misspec_sum[[1]]$SRMR_M[1]

#######


SRMR_vline <- lapply(seq(length(misspec_sum)), function(x) misspec_sum[[x, drop=FALSE]][1,][1])


test <- sapply(seq(length(misspec_sum)),function(x) SRMR_vline[[1]] %>% 
                 pull("SRMR_M"))

test

t2 <- unlist(test)
t3 <- as.data.frame(t2)


first(SRMR_vline)

misspec_sum[[1]]$SRMR_M[1]

lapply(combined, function(x) lapply(x, function(y) length(y)))
lapply(misspec_sum,function(x) lapply(x,function(y) y=x[[1]]), lapply(y, function(z) y[[1]]))

misspec_sum[1]

first(misspec_sum[[1]][,1])


sapply(misspec_sum, "[[", 1, "[",1)

SRMR_vline <- lapply(seq(length(misspec_sum)), function(x) misspec_sum[[x]][1,][1])

SRMR_vline[2]

lapply(misspec_sum,function(x) lapply(x,function(y) x[[1]][1]))

lapply(misspec_sum,function(x) lapply(x,function(y) y=x[[1]]), lapply(y, function(z) y[1])) 


x = list(list(1,2), list(3,4), list(5,6))
x1 = lapply(x, function(l) l[[1]])

###########

SRMR_plot <- plot %>%
  ggplot(aes(x=SRMR,fill=Model))+
  geom_histogram(position="identity",
                 alpha=.5)+
  scale_fill_manual(values=c("#E9798C","#66C2F5"))+
  geom_vline(aes(xintercept=misspec_sum$SRMR_M[1],
                 linetype="misspec_sum$SRMR_M[1]",color="misspec_sum$SRMR_M[1]"),
             size=.6)+
  geom_vline(aes(xintercept=.08,
                 linetype=".08",color=".08"),
             size=.75)+
  scale_color_manual(name="Cutoff Values",
                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                     values=c("misspec_sum$SRMR_M[1]"="black",
                              ".08"="black"))+
  scale_linetype_manual(name="Cutoff Values",
                        labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                        values=c("misspec_sum$SRMR_M[1]"="longdash",
                                 ".08"="dotted"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.box = "vertical")

RMSEA_plot <- plot %>%
  ggplot(aes(x=RMSEA,fill=Model))+
  geom_histogram(position="identity",
                 alpha=.5)+
  scale_fill_manual(values=c("#E9798C","#66C2F5"))+
  geom_vline(aes(xintercept=misspec_sum$RMSEA_M[1],
                 linetype="misspec_sum$RMSEA_M[1]",color="misspec_sum$RMSEA_M[1]"),
             size=.6)+
  geom_vline(aes(xintercept=.06,
                 linetype=".06",color=".06"),
             size=.75)+
  scale_color_manual(name="Cutoff Values",
                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                     values=c("misspec_sum$RMSEA_M[1]"="black",
                              ".06"="black"))+
  scale_linetype_manual(name="Cutoff Values",
                        labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                        values=c("misspec_sum$RMSEA_M[1]"="longdash",
                                 ".06"="dotted"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.box = "vertical")

CFI_plot <- plot %>%
  ggplot(aes(x=CFI,fill=Model))+
  geom_histogram(position="identity",
                 alpha=.5)+
  scale_fill_manual(values=c("#E9798C","#66C2F5"))+
  geom_vline(aes(xintercept=misspec_sum$CFI_M[1],
                 linetype="misspec_sum$CFI_M[1]",color="misspec_sum$CFI_M[1]"),
             size=.6)+
  geom_vline(aes(xintercept=.95,
                 linetype=".95",color=".95"),
             size=.75)+
  scale_color_manual(name="Cutoff Values",
                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                     values=c("misspec_sum$CFI_M[1]"="black",
                              ".95"="black"))+
  scale_linetype_manual(name="Cutoff Values",
                        labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                        values=c("misspec_sum$CFI_M[1]"="longdash",
                                 ".95"="dotted"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.box = "vertical")