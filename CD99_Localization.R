# Code for visualizing CD99 localization heatmap
library(readxl)
library(tidyverse)
library(ggplot2)

d = read_excel("~/Downloads/intensity.xlsx")
d1 = d[1:11,]
d1$Group = "WT"
d2 = d[15:24,]
d2$Group = "Cyto"
d3 = d[28:38,]
d3$Group = "TM"

dt = rbind.data.frame(d1, d2, d3)
dt = dt[,c(5,2:4)]
dt$name = paste(dt$Group, rownames(dt), sep = "_")
colnames(dt) = c("Group", "cSMAC", "pSMAC", "dSMAC", "name")

## Ordering data.frame
dt.g = dt %>% gather("cSMAC":"dSMAC", key = "Group2", value = "value")
dt.g$value = as.numeric(dt.g$value)
dt.g$Group2 = factor(dt.g$Group2, levels=c("cSMAC","pSMAC","dSMAC"))
dt.g$Group = factor(dt.g$Group, levels = c("WT", "Cyto", "TM"))
dt.g = dt.g[order(dt.g$Group),]
dt.g$rank = as.numeric(rownames(dt.g))
dt.g = dt.g[!is.na(dt.g$value),]

## ggplot
p <- ggplot(dt.g, aes(Group2 , reorder(name, desc(rank)))) +
  geom_tile(aes(fill = value), colour = "white", size = 0) + 
  theme_minimal()+
  theme(legend.position = "right") +
#  scale_fill_gradient(low = "white", high = "#EA3323", na.value = "grey80")+
 # scale_fill_gradient2(low = "white", mid = "yellow", high = "#EA3323", midpoint = 1000,na.value = "grey80")+
                      # labels = scales::rescale(c(min(dt.g[!is.na(dt.g$value),]$value), max(dt.g[!is.na(dt.g$value),]$value)))) +
  scale_fill_gradientn(colors = c("white", "#FFFFC6", "#EA3323"), values=rescale(c(min(dt.g[!is.na(dt.g$value),]$value), 600, max(dt.g[!is.na(dt.g$value),]$value))), 
                       limits = c(0, 1600), oob=squish, na.value= "grey85") +
  labs(title = 'CD99 Localization', x = '', y = '') +
  theme(axis.ticks=element_blank(), 
    #    panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=10, colour="black"),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,vjust = 1, hjust = 0.5))  
ggsave("~/Dropbox/NGR_SNU_2019/CD99_SYK_202201/Figure_localization_noline.pdf",p,  height = 4.5, width = 5.5)


# ggplot(dt.g, (aes(x = Group2, y =  reorder(name, desc(rank)), fill = value))) + 
#  geom_col(position = "fill", width = 1) + 
#  scale_fill_gradientn(colors = c("white", "#FFFFC6", "#EA3323"), values=rescale(c(min(dt.g[!is.na(dt.g$value),]$value), 600, max(dt.g[!is.na(dt.g$value),]$value))), 
#                       limits = c(0, 1600), oob=squish) +
#  facet_grid(Group ~ ., as.table = FALSE, switch = "y") 


