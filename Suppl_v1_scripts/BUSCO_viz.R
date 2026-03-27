
library(stringr)
library(openxlsx)
library(data.table)
library(magrittr)
#library(GGally)
#library(propr)
library(data.table)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(dplyr)
library(plyr)
library(microshades)
library(scales)
library(ggbreak) 
library(patchwork)
library(ggpubr)
library(gplots)
library(PerformanceAnalytics)

library(patchwork)
library(microshades)
library(ggtext)
library(gridExtra)
library(ggpmisc)



list <- list.files(path = "/home/zavadska/Downloads/", pattern = "BUSCO_out.csv")

big_df <- data.frame()

for (i in list) {
  
  name <- str_replace(i, "_BUSCO_out.csv","")
  
  BUSCO_eugl <- fread(paste0("/home/zavadska/Downloads/",i))
  
  BUSCO_eugl_m <- melt(BUSCO_eugl)
  
  n <- unique(BUSCO_eugl_m$value[which(BUSCO_eugl_m$variable=="n")])
  

  
  BUSCO_eugl_m$setting <- name
  big_df <- rbind(BUSCO_eugl_m, big_df)
  
}



## Plots for the preprint



big_df_plot <- big_df[!big_df$variable%in% c("C","n"),]

big_df_plot$redundancy <- "NR90"
big_df_plot$busco_db <- str_split(big_df_plot$setting, "_", simplify = T)[,1]


db_list_eugl <- c("KUCKIN", "D44", "10D", "13G", "G65133", "RhMon", "OB23", "Gemkin", "BABKIN", "16Ckin", "INUKIN")
sp_to_label_eugl <- paste0( big_df_plot$species[big_df_plot$species %in% db_list_eugl ] , ".", big_df_plot$group[big_df_plot$species %in% db_list_eugl ] )



p1 <- ggplot(big_df_plot[big_df_plot$busco_db=="EUGL",], aes(value/130, interaction(species,group), fill =  as.factor(variable)) )+
  geom_bar(stat = "identity", width = 0.8)+
  scale_y_discrete(labels = ~ if_else( .x %in% sp_to_label_eugl, paste0("<span style='color: chartreuse4'><b>", .x, "</b></span>"),  .x  )) +
  #geom_text(aes(label = value), hjust=-0.5, size = 2) +
  theme_bw()+ theme(axis.text.y = element_markdown())+xlab("BUSCO completeness")+ylab("species") + labs(fill = "BUSCO category")+
  theme(legend.position = "bottom")+
  scale_fill_brewer(palette = "Set3")+
  ggtitle("Euglenozoa")


p2 <- ggplot(big_df_plot[big_df_plot$busco_db=="EUK",], aes(value/255, interaction(species,group), fill =  as.factor(variable)) )+
  geom_bar(stat = "identity", width = 0.8)+
  scale_y_discrete(labels = ~ if_else( .x %in% sp_to_label_eugl, paste0("<span style='color: chartreuse4'><b>", .x, "</b></span>"),  .x  )) +
  #geom_text(aes(label = value), hjust=-0.5, size = 2) +
  theme_bw()+ theme(axis.text.y = element_markdown() )+xlab("BUSCO completeness")+ylab("species") + labs(fill = "BUSCO category")+
  theme(legend.position = "bottom")+
  scale_fill_brewer(palette = "Set3")+
  ggtitle("Eukaryota")


gg <- (p1 + theme(legend.position="none")) | p2  

grid.arrange(p1,p2, nrow = 1)

ggsave("BAD_BUSCO.png", gg, width = 15, height = 10)
