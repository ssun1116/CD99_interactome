######### New CD99 Heatmap_started on.0109 #########
library(readxl)
library(writexl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## Main table
DF = read_excel("Data/IP_raw_data_all.xlsx", sheet = 8)
DF$...26 <- tolower(DF$...26)
DF$...26 <- firstup(DF$...26)

##### How to filter?? #####
# ## list for cytoskeleton group - actomyosin, tubulin
# protein_list <- read.delim('Tables/table.cyto_dup.stringDB_700.CD99_20200628.txt')
# acto <- protein_list[grepl("Actomyosin", protein_list$AttrA), ][,1]
# tub <- protein_list[grepl("Microtubule", protein_list$AttrA), ][,1]

## Load pathways - GO information
go = read.delim('Data/gene_association.mgi', skip=24, header = F)
go = go[,1:7]
colnames(go) <- c('db', 'db_id', 'gene_name', 'qualifier', 'go_id', 'db_ref', 'evidence_code')

## Classify GO IDs into groups 
p_actin = unique(go[go$go_id=='GO:0015629',]$gene_name)
actin_supp = c('Myl6', 'Actbl2', 'Actc1', 'Iqgap1') 
p_actin = c(p_actin, actin_supp) # actin cytoskeleton
p_tubulin = unique(go[go$go_id %in% c('GO:0005874', 'GO:0015630'),]$gene_name) # microtubule
p_myosin = unique(go[go$go_id=='GO:0016459',]$gene_name)
p_myosin = c(p_myosin, c('Myl12a')) # myosin complex
p_actomyo = c(p_actin, p_myosin)
##########################

## 1. Cho_20180514
d1 = DF[,1:7]
colnames(d1) = d1[1,]
d1= d1[2:nrow(d1),c(2,4:ncol(d1))]
d1$`WT(Ctrl IP)_1st` <- round(as.numeric(d1$`WT(Ctrl IP)_1st`), 2)
d1$`WT(Ctrl IP)_2st` <- round(as.numeric(d1$`WT(Ctrl IP)_2st`), 2)
d1$`WT(GFP IP)_1st` <- round(as.numeric(d1$`WT(GFP IP)_1st`), 2)
d1$`WT(GFP IP)_2st` <- round(as.numeric(d1$`WT(GFP IP)_2st`), 2)

### Step 1. Divide d1 by actomyosin & microtubule group
d1_acto = d1[d1$`Gene name` %in% p_actomyo,]
d1_tub = d1[d1$`Gene name` %in% p_tubulin,]

## Step 2. Divide d1_acto and d1_tub by Num & OnOff
d1_acto.Num = d1_acto[d1_acto$`WT(Ctrl IP)_1st` != 0 & d1_acto$`WT(Ctrl IP)_2st` != 0 & d1_acto$`WT(GFP IP)_1st` != 0 & d1_acto$`WT(GFP IP)_2st` != 0, ]
d1_acto.OnOff = d1_acto[d1_acto$`WT(Ctrl IP)_1st` == 0 | d1_acto$`WT(Ctrl IP)_2st` == 0 | d1_acto$`WT(GFP IP)_1st` == 0 | d1_acto$`WT(GFP IP)_2st` == 0, ]
d1_tub.Num = d1_tub[d1_tub$`WT(Ctrl IP)_1st` != 0 & d1_tub$`WT(Ctrl IP)_2st` != 0 & d1_tub$`WT(GFP IP)_1st` != 0 & d1_tub$`WT(GFP IP)_2st` != 0, ]
d1_tub.OnOff = d1_tub[d1_tub$`WT(Ctrl IP)_1st` == 0 | d1_tub$`WT(Ctrl IP)_2st` == 0 | d1_tub$`WT(GFP IP)_1st` == 0 | d1_tub$`WT(GFP IP)_2st` == 0, ]

## d1_acto.Num
d1_acto.Num$WT_1st = as.numeric(round((d1_acto.Num$`WT(GFP IP)_1st`)/(d1_acto.Num$`WT(Ctrl IP)_1st`),2))
d1_acto.Num$WT_2st = as.numeric(round((d1_acto.Num$`WT(GFP IP)_2st`)/(d1_acto.Num$`WT(Ctrl IP)_2st`),2))
d1_acto.Num$CTL_1st = 1
d1_acto.Num$CTL_2st = 1
d1_acto.Num = d1_acto.Num[,c(1,8, 6, 9, 7)]
d1_acto.Num$mean.WT = as.numeric(round(((d1_acto.Num$WT_1st + d1_acto.Num$WT_2st)/2),2))

## d1_acto.OnOff
d1_acto.OnOff$CTL_1st = ifelse(d1_acto.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d1_acto.OnOff$`WT(Ctrl IP)_1st`))
d1_acto.OnOff$CTL_2st = ifelse(d1_acto.OnOff$`WT(Ctrl IP)_2st`== 0, "Off", as.numeric(d1_acto.OnOff$`WT(Ctrl IP)_2st`))
d1_acto.OnOff$WT_1st = ifelse(d1_acto.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                              ifelse(d1_acto.OnOff$CTL_1st == "Off", "On", as.numeric(round((d1_acto.OnOff$`WT(GFP IP)_1st`)/(d1_acto.OnOff$`WT(Ctrl IP)_1st`),2))))
d1_acto.OnOff$WT_2st = ifelse(d1_acto.OnOff$`WT(GFP IP)_2st` == 0, "Off",
                              ifelse(d1_acto.OnOff$CTL_2st == "Off", "On", as.numeric(round((d1_acto.OnOff$`WT(GFP IP)_2st`)/(d1_acto.OnOff$`WT(Ctrl IP)_2st`),2))))
d1_acto.OnOff$CTL_1st = ifelse(d1_acto.OnOff$CTL_1st != "Off", 1, d1_acto.OnOff$CTL_1st)
d1_acto.OnOff$CTL_2st = ifelse(d1_acto.OnOff$CTL_2st != "Off", 1, d1_acto.OnOff$CTL_2st)
d1_acto.OnOff = d1_acto.OnOff[,c(1,6, 8, 7, 9)]
d1_acto.OnOff$mean.WT = ifelse(d1_acto.OnOff$WT_1st %in% c("On", "Off"), d1_acto.OnOff$WT_2st, d1_acto.OnOff$WT_1st)
d1_acto.OnOff$mean.WT = ifelse(d1_acto.OnOff$mean.WT %in% c("On", "Off"),
                               ifelse(d1_acto.OnOff$WT_1st == "On" & d1_acto.OnOff$WT_2st == "On", "On",
                                      ifelse(d1_acto.OnOff$WT_1st == "Off" & d1_acto.OnOff$WT_2st == "Off", "Off", "On/Off")), d1_acto.OnOff$mean.WT)


## Bind Num & OnOff -> Order
d1_acto = rbind.data.frame(d1_acto.Num, d1_acto.OnOff)
d1_acto.1 = d1_acto[d1_acto$mean.WT %in% c("On", "Off", "On/Off"),]
d1_acto.1 = d1_acto.1[order(d1_acto.1$mean.WT),]
d1_acto.2 = d1_acto[!(d1_acto$mean.WT %in% c("On", "Off", "On/Off")),]
d1_acto.2$mean.WT = as.numeric(d1_acto.2$mean.WT)
d1_acto.2= d1_acto.2[order(d1_acto.2$mean.WT, decreasing = TRUE), ]
d1_acto = rbind.data.frame(d1_acto.2, d1_acto.1)
d1_acto = d1_acto[c(1:9,11:14,10),]

d1_acto.tb = d1_acto[,c(2:(ncol(d1_acto)-1),1)]
write_xlsx(d1_acto.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180514.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

d1_acto.tmp = d1_acto
d1_acto.tmp$rank = rownames(d1_acto)
dt.d1_acto <- d1_acto.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:WT_2st)
dt.d1_acto$rank = as.numeric(dt.d1_acto$rank)
dt.d1_acto = dt.d1_acto[order(dt.d1_acto$rank),]
dt.d1_acto$FC.color = ifelse(dt.d1_acto$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d1_acto$FC.color = ifelse(dt.d1_acto$Group %in% c("WT_1st", "WT_2st") & dt.d1_acto$FC == "Off", 0, dt.d1_acto$FC.color)
dt.d1_acto$FC.color = ifelse(dt.d1_acto$Group %in% c("WT_1st", "WT_2st") & dt.d1_acto$FC == "On", 5, dt.d1_acto$FC.color)
dt.d1_acto$FC.color = ifelse(dt.d1_acto$Group %in% c("WT_1st", "WT_2st") & dt.d1_acto$FC < 1, 0, dt.d1_acto$FC.color)
dt.d1_acto$FC.color = ifelse(is.na(dt.d1_acto$FC.color), dt.d1_acto$FC, dt.d1_acto$FC.color)
dt.d1_acto$FC.color = as.numeric(dt.d1_acto$FC.color)

p.d1_acto <- ggplot(dt.d1_acto, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st", "CTL_2st", "WT_2st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d1_acto, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.8) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  
#"#7b9fd4", mid = 'white', high="#Cd1836"
p.d1_acto

## d1_tub.Num
d1_tub.Num$WT_1st = as.numeric(round((d1_tub.Num$`WT(GFP IP)_1st`)/(d1_tub.Num$`WT(Ctrl IP)_1st`),2))
d1_tub.Num$WT_2st = as.numeric(round((d1_tub.Num$`WT(GFP IP)_2st`)/(d1_tub.Num$`WT(Ctrl IP)_2st`),2))
d1_tub.Num$CTL_1st = 1
d1_tub.Num$CTL_2st = 1
d1_tub.Num = d1_tub.Num[,c(1,8, 6, 9, 7)]
d1_tub.Num$mean.WT = as.numeric(round(((d1_tub.Num$WT_1st + d1_tub.Num$WT_2st)/2),2))

## d1_tub.OnOff
d1_tub.OnOff$CTL_1st = ifelse(d1_tub.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d1_tub.OnOff$`WT(Ctrl IP)_1st`))
d1_tub.OnOff$CTL_2st = ifelse(d1_tub.OnOff$`WT(Ctrl IP)_2st`== 0, "Off", as.numeric(d1_tub.OnOff$`WT(Ctrl IP)_2st`))
d1_tub.OnOff$WT_1st = ifelse(d1_tub.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                              ifelse(d1_tub.OnOff$CTL_1st == "Off", "On", as.numeric(round((d1_tub.OnOff$`WT(GFP IP)_1st`)/(d1_tub.OnOff$`WT(Ctrl IP)_1st`),2))))
d1_tub.OnOff$WT_2st = ifelse(d1_tub.OnOff$`WT(GFP IP)_2st` == 0, "Off",
                              ifelse(d1_tub.OnOff$CTL_2st == "Off", "On", as.numeric(round((d1_tub.OnOff$`WT(GFP IP)_2st`)/(d1_tub.OnOff$`WT(Ctrl IP)_2st`),2))))
d1_tub.OnOff$CTL_1st = ifelse(d1_tub.OnOff$CTL_1st != "Off", 1, d1_tub.OnOff$CTL_1st)
d1_tub.OnOff$CTL_2st = ifelse(d1_tub.OnOff$CTL_2st != "Off", 1, d1_tub.OnOff$CTL_2st)
d1_tub.OnOff = d1_tub.OnOff[,c(1,6, 8, 7, 9)]
d1_tub.OnOff$mean.WT = ifelse(d1_tub.OnOff$WT_1st %in% c("On", "Off"), d1_tub.OnOff$WT_2st, d1_tub.OnOff$WT_1st)
d1_tub.OnOff$mean.WT = ifelse(d1_tub.OnOff$mean.WT %in% c("On", "Off"),
                               ifelse(d1_tub.OnOff$WT_1st == "On" & d1_tub.OnOff$WT_2st == "On", "On",
                                      ifelse(d1_tub.OnOff$WT_1st == "Off" & d1_tub.OnOff$WT_2st == "Off", "Off", "On/Off")), d1_tub.OnOff$mean.WT)


## Bind Num & OnOff -> Order
d1_tub = rbind.data.frame(d1_tub.Num, d1_tub.OnOff)
d1_tub.1 = d1_tub[d1_tub$mean.WT %in% c("On", "Off", "On/Off"),]
d1_tub.1 = d1_tub.1[order(d1_tub.1$mean.WT),]
d1_tub.2 = d1_tub[!(d1_tub$mean.WT %in% c("On", "Off", "On/Off")),]
d1_tub.2$mean.WT = as.numeric(d1_tub.2$mean.WT)
d1_tub.2= d1_tub.2[order(d1_tub.2$mean.WT, decreasing = TRUE), ]
d1_tub = rbind.data.frame(d1_tub.2, d1_tub.1)
d1_tub = d1_tub[c(1:13,17,14),]

d1_tub.tb = d1_tub[,c(2:(ncol(d1_tub)-1),1)]
write_xlsx(d1_tub.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180514.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

d1_tub.tmp = d1_tub
d1_tub.tmp$rank = rownames(d1_tub)
dt.d1_tub <- d1_tub.tmp %>%
  gather(Group, FC, CTL_1st:WT_2st)
dt.d1_tub$rank = as.numeric(dt.d1_tub$rank)
dt.d1_tub = dt.d1_tub[order(dt.d1_tub$rank),]
dt.d1_tub$FC.color = ifelse(dt.d1_tub$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d1_tub$FC.color = ifelse(dt.d1_tub$Group %in% c("WT_1st", "WT_2st") & dt.d1_tub$FC == "Off", 0, dt.d1_tub$FC.color)
dt.d1_tub$FC.color = ifelse(dt.d1_tub$Group %in% c("WT_1st", "WT_2st") & dt.d1_tub$FC == "On", 5, dt.d1_tub$FC.color)
dt.d1_tub$FC.color = ifelse(dt.d1_tub$Group %in% c("WT_1st", "WT_2st") & dt.d1_tub$FC < 1, 0, dt.d1_tub$FC.color)
dt.d1_tub$FC.color = ifelse(is.na(dt.d1_tub$FC.color), dt.d1_tub$FC, dt.d1_tub$FC.color)
dt.d1_tub$FC.color = as.numeric(dt.d1_tub$FC.color)

p.d1_tub <- ggplot(dt.d1_tub, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st", "CTL_2st", "WT_2st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d1_tub, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        # legend.position="none",
        plot.margin = unit(c(0.2, 0.5, 0.5, 0.2), "cm")) +
        labs(fill='  FC') 
legend <- get_legend(p.d1_tub)
p.d1_tub <- p.d1_tub + theme(legend.position = "none")
#plot_grid(p.d1_acto, p.d1_tub, legend, ncol = 3)  
p = plot_grid(p.d1_acto, p.d1_tub, ncol = 2)  
ggsave('Figures/Cho_20180514_OnOff.FC.Heatmap_0111.pdf', p, width =5, height = 5, useDingbats=F)


## 2. Cho_20180912
d2 = DF[,9:15]
colnames(d2) = d2[1,]
d2= d2[2:nrow(d2),c(2,4:ncol(d2))]
d2$WT_Con_1st <- round(as.numeric(d2$WT_Con_1st), 2)
d2$WT_Con_2nd <- round(as.numeric(d2$WT_Con_2nd), 2)
d2$WT_GFP_1st <- round(as.numeric(d2$WT_GFP_1st), 2)
d2$WT_GFP_2nd <- round(as.numeric(d2$WT_GFP_2nd), 2)
# d2$gene_name = toupper(d2$`Gene name`)
# d2$`Gene name` <- NULL
colnames(d2) = c("Gene name", "WT(Ctrl IP)_1st", "WT(Ctrl IP)_2st", "WT(GFP IP)_1st", "WT(GFP IP)_2st")

### Step 1. Divide by actomyosin & microtubule group
d2_acto = d2[d2$`Gene name` %in% p_actomyo,]
d2_tub = d2[d2$`Gene name` %in% p_tubulin,]

## Step 2. Divide d2_acto and d2_tub by Num & OnOff
d2_acto.Num = d2_acto[d2_acto$`WT(Ctrl IP)_1st` != 0 & d2_acto$`WT(Ctrl IP)_2st` != 0 & d2_acto$`WT(GFP IP)_1st` != 0 & d2_acto$`WT(GFP IP)_2st` != 0, ]
d2_acto.OnOff = d2_acto[d2_acto$`WT(Ctrl IP)_1st` == 0 | d2_acto$`WT(Ctrl IP)_2st` == 0 | d2_acto$`WT(GFP IP)_1st` == 0 | d2_acto$`WT(GFP IP)_2st` == 0, ]
d2_tub.Num = d2_tub[d2_tub$`WT(Ctrl IP)_1st` != 0 & d2_tub$`WT(Ctrl IP)_2st` != 0 & d2_tub$`WT(GFP IP)_1st` != 0 & d2_tub$`WT(GFP IP)_2st` != 0, ]
d2_tub.OnOff = d2_tub[d2_tub$`WT(Ctrl IP)_1st` == 0 | d2_tub$`WT(Ctrl IP)_2st` == 0 | d2_tub$`WT(GFP IP)_1st` == 0 | d2_tub$`WT(GFP IP)_2st` == 0, ]


## d2_acto.Num
d2_acto.Num$WT_1st = as.numeric(round((d2_acto.Num$`WT(GFP IP)_1st`)/(d2_acto.Num$`WT(Ctrl IP)_1st`),2))
d2_acto.Num$WT_2st = as.numeric(round((d2_acto.Num$`WT(GFP IP)_2st`)/(d2_acto.Num$`WT(Ctrl IP)_2st`),2))
d2_acto.Num$CTL_1st = 1
d2_acto.Num$CTL_2st = 1
d2_acto.Num = d2_acto.Num[,c(1,8, 6, 9, 7)]
d2_acto.Num$mean.WT = as.numeric(round(((d2_acto.Num$WT_1st + d2_acto.Num$WT_2st)/2),2))

## d2_acto.OnOff
d2_acto.OnOff$CTL_1st = ifelse(d2_acto.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d2_acto.OnOff$`WT(Ctrl IP)_1st`))
d2_acto.OnOff$CTL_2st = ifelse(d2_acto.OnOff$`WT(Ctrl IP)_2st`== 0, "Off", as.numeric(d2_acto.OnOff$`WT(Ctrl IP)_2st`))
d2_acto.OnOff$WT_1st = ifelse(d2_acto.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                              ifelse(d2_acto.OnOff$CTL_1st == "Off", "On", as.numeric(round((d2_acto.OnOff$`WT(GFP IP)_1st`)/(d2_acto.OnOff$`WT(Ctrl IP)_1st`),2))))
d2_acto.OnOff$WT_2st = ifelse(d2_acto.OnOff$`WT(GFP IP)_2st` == 0, "Off",
                              ifelse(d2_acto.OnOff$CTL_2st == "Off", "On", as.numeric(round((d2_acto.OnOff$`WT(GFP IP)_2st`)/(d2_acto.OnOff$`WT(Ctrl IP)_2st`),2))))
d2_acto.OnOff$CTL_1st = ifelse(d2_acto.OnOff$CTL_1st != "Off", 1, d2_acto.OnOff$CTL_1st)
d2_acto.OnOff$CTL_2st = ifelse(d2_acto.OnOff$CTL_2st != "Off", 1, d2_acto.OnOff$CTL_2st)
d2_acto.OnOff = d2_acto.OnOff[,c(1,6, 8, 7, 9)]
d2_acto.OnOff$mean.WT = ifelse(d2_acto.OnOff$WT_1st %in% c("On", "Off"), d2_acto.OnOff$WT_2st, d2_acto.OnOff$WT_1st)
d2_acto.OnOff$mean.WT = ifelse(d2_acto.OnOff$mean.WT %in% c("On", "Off"),
                               ifelse(d2_acto.OnOff$WT_1st == "On" & d2_acto.OnOff$WT_2st == "On", "On",
                                      ifelse(d2_acto.OnOff$WT_1st == "Off" & d2_acto.OnOff$WT_2st == "Off", "Off", "On/Off")), d2_acto.OnOff$mean.WT)


## Bind Num & OnOff -> Order
d2_acto = rbind.data.frame(d2_acto.Num, d2_acto.OnOff)
d2_acto.1 = d2_acto[d2_acto$mean.WT %in% c("On", "Off", "On/Off"),]
d2_acto.1 = d2_acto.1[order(d2_acto.1$mean.WT),]
d2_acto.2 = d2_acto[!(d2_acto$mean.WT %in% c("On", "Off", "On/Off")),]
d2_acto.2$mean.WT = as.numeric(d2_acto.2$mean.WT)
d2_acto.2= d2_acto.2[order(d2_acto.2$mean.WT, decreasing = TRUE), ]
d2_acto = rbind.data.frame(d2_acto.2, d2_acto.1)
d2_acto = d2_acto[c(1:10,15:24),]

d2_acto.tb = d2_acto[,c(2:(ncol(d2_acto)-1),1)]
write_xlsx(d2_acto.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180912.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

d2_acto.tmp = d2_acto
d2_acto.tmp$rank = rownames(d2_acto)
dt.d2_acto <- d2_acto.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:WT_2st)
dt.d2_acto$rank = as.numeric(dt.d2_acto$rank)
dt.d2_acto = dt.d2_acto[order(dt.d2_acto$rank),]
dt.d2_acto$FC.color = ifelse(dt.d2_acto$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d2_acto$FC.color = ifelse(dt.d2_acto$Group %in% c("WT_1st", "WT_2st") & dt.d2_acto$FC == "Off", 0, dt.d2_acto$FC.color)
dt.d2_acto$FC.color = ifelse(dt.d2_acto$Group %in% c("WT_1st", "WT_2st") & dt.d2_acto$FC == "On", 5, dt.d2_acto$FC.color)
dt.d2_acto$FC.color = ifelse(dt.d2_acto$Group %in% c("WT_1st", "WT_2st") & dt.d2_acto$FC < 1, 0, dt.d2_acto$FC.color)
dt.d2_acto$FC.color = ifelse(is.na(dt.d2_acto$FC.color), dt.d2_acto$FC, dt.d2_acto$FC.color)
dt.d2_acto$FC.color = as.numeric(dt.d2_acto$FC.color)

p.d2_acto <- ggplot(dt.d2_acto, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st", "CTL_2st", "WT_2st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d2_acto, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.8) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  
#"#7b9fd4", mid = 'white', high="#CD2836"
p.d2_acto

## d2_tub.Num
d2_tub.Num$WT_1st = as.numeric(round((d2_tub.Num$`WT(GFP IP)_1st`)/(d2_tub.Num$`WT(Ctrl IP)_1st`),2))
d2_tub.Num$WT_2st = as.numeric(round((d2_tub.Num$`WT(GFP IP)_2st`)/(d2_tub.Num$`WT(Ctrl IP)_2st`),2))
d2_tub.Num$CTL_1st = 1
d2_tub.Num$CTL_2st = 1
d2_tub.Num = d2_tub.Num[,c(1,8, 6, 9, 7)]
d2_tub.Num$mean.WT = as.numeric(round(((d2_tub.Num$WT_1st + d2_tub.Num$WT_2st)/2),2))

## d2_tub.OnOff
d2_tub.OnOff$CTL_1st = ifelse(d2_tub.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d2_tub.OnOff$`WT(Ctrl IP)_1st`))
d2_tub.OnOff$CTL_2st = ifelse(d2_tub.OnOff$`WT(Ctrl IP)_2st`== 0, "Off", as.numeric(d2_tub.OnOff$`WT(Ctrl IP)_2st`))
d2_tub.OnOff$WT_1st = ifelse(d2_tub.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                             ifelse(d2_tub.OnOff$CTL_1st == "Off", "On", as.numeric(round((d2_tub.OnOff$`WT(GFP IP)_1st`)/(d2_tub.OnOff$`WT(Ctrl IP)_1st`),2))))
d2_tub.OnOff$WT_2st = ifelse(d2_tub.OnOff$`WT(GFP IP)_2st` == 0, "Off",
                             ifelse(d2_tub.OnOff$CTL_2st == "Off", "On", as.numeric(round((d2_tub.OnOff$`WT(GFP IP)_2st`)/(d2_tub.OnOff$`WT(Ctrl IP)_2st`),2))))
d2_tub.OnOff$CTL_1st = ifelse(d2_tub.OnOff$CTL_1st != "Off", 1, d2_tub.OnOff$CTL_1st)
d2_tub.OnOff$CTL_2st = ifelse(d2_tub.OnOff$CTL_2st != "Off", 1, d2_tub.OnOff$CTL_2st)
d2_tub.OnOff = d2_tub.OnOff[,c(1,6, 8, 7, 9)]
d2_tub.OnOff$mean.WT = ifelse(d2_tub.OnOff$WT_1st %in% c("On", "Off"), d2_tub.OnOff$WT_2st, d2_tub.OnOff$WT_1st)
d2_tub.OnOff$mean.WT = ifelse(d2_tub.OnOff$mean.WT %in% c("On", "Off"),
                              ifelse(d2_tub.OnOff$WT_1st == "On" & d2_tub.OnOff$WT_2st == "On", "On",
                                     ifelse(d2_tub.OnOff$WT_1st == "Off" & d2_tub.OnOff$WT_2st == "Off", "Off", "On/Off")), d2_tub.OnOff$mean.WT)


## Bind Num & OnOff -> Order
d2_tub = rbind.data.frame(d2_tub.Num, d2_tub.OnOff)
d2_tub.1 = d2_tub[d2_tub$mean.WT %in% c("On", "Off", "On/Off"),]
d2_tub.1 = d2_tub.1[order(d2_tub.1$mean.WT),]
d2_tub.2 = d2_tub[!(d2_tub$mean.WT %in% c("On", "Off", "On/Off")),]
d2_tub.2$mean.WT = as.numeric(d2_tub.2$mean.WT)
d2_tub.2= d2_tub.2[order(d2_tub.2$mean.WT, decreasing = TRUE), ]
d2_tub = rbind.data.frame(d2_tub.2, d2_tub.1)
d2_tub = d2_tub[c(1:21,27:37),]
d2_tub= d2_tub[c(1:18,22:32,19:21),]

d2_tub.tb = d2_tub[,c(2:(ncol(d2_tub)-1),1)]
write_xlsx(d2_tub.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180912.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

d2_tub.tmp = d2_tub
d2_tub.tmp$rank = rownames(d2_tub)
dt.d2_tub <- d2_tub.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:WT_2st)
dt.d2_tub$rank = as.numeric(dt.d2_tub$rank)
dt.d2_tub = dt.d2_tub[order(dt.d2_tub$rank),]
dt.d2_tub$FC.color = ifelse(dt.d2_tub$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d2_tub$FC.color = ifelse(dt.d2_tub$Group %in% c("WT_1st", "WT_2st") & dt.d2_tub$FC == "Off", 0, dt.d2_tub$FC.color)
dt.d2_tub$FC.color = ifelse(dt.d2_tub$Group %in% c("WT_1st", "WT_2st") & dt.d2_tub$FC == "On", 5, dt.d2_tub$FC.color)
dt.d2_tub$FC.color = ifelse(dt.d2_tub$Group %in% c("WT_1st", "WT_2st") & dt.d2_tub$FC < 1, 0, dt.d2_tub$FC.color)
dt.d2_tub$FC.color = ifelse(is.na(dt.d2_tub$FC.color), dt.d2_tub$FC, dt.d2_tub$FC.color)
dt.d2_tub$FC.color = as.numeric(dt.d2_tub$FC.color)

p.d2_tub <- ggplot(dt.d2_tub, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st", "CTL_2st", "WT_2st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d2_tub, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.6) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        # legend.position="none",
        plot.margin = unit(c(0.2, 0.5, 0.5, 0.2), "cm")) +
  labs(fill='  FC') 
legend <- get_legend(p.d2_tub)
p.d2_tub <- p.d2_tub + theme(legend.position = "none")
#plot_grid(p.d2_acto, p.d2_tub, legend, ncol = 3)  
p = plot_grid(p.d2_acto, p.d2_tub, ncol = 2)  
p
ggsave('Figures/Cho_20180912_OnOff.FC.Heatmap_0111.pdf', p, width =5, height = 7.5, useDingbats=F)


## 3. Cho_20180625
d3 = DF[,17:19]
colnames(d3) = d3[1,] ## No Symbol. Merge from Cho 0912
d2 = DF[,9:15]
colnames(d2) = d2[1,]
d2.symbol = d2[2:nrow(d2),1:2]
d3 = as.data.frame(merge(d3, d2.symbol, by = "Protein name"))
d3= d3[,c(4,2:3)]
colnames(d3) = c("Gene name", "WT(Ctrl IP)_1st", "WT(GFP IP)_1st")
d3$`WT(Ctrl IP)_1st` <- round(as.numeric(d3$`WT(Ctrl IP)_1st`), 2)
d3$`WT(GFP IP)_1st` <- round(as.numeric(d3$`WT(GFP IP)_1st`), 2)

### Step 1. Divide by actomyosin & microtubule group
d3_acto = d3[d3$`Gene name` %in% p_actomyo,]
d3_tub = d3[d3$`Gene name` %in% p_tubulin,]

## Step 2. Divide d3_acto and d3_tub by Num & OnOff
d3_acto.Num = d3_acto[d3_acto$`WT(Ctrl IP)_1st` != 0 & d3_acto$`WT(GFP IP)_1st` != 0, ]
d3_acto.OnOff = d3_acto[d3_acto$`WT(Ctrl IP)_1st` == 0 | d3_acto$`WT(GFP IP)_1st` == 0, ]
d3_tub.Num = d3_tub[d3_tub$`WT(Ctrl IP)_1st` != 0 & d3_tub$`WT(GFP IP)_1st` != 0, ]
d3_tub.OnOff = d3_tub[d3_tub$`WT(Ctrl IP)_1st` == 0 | d3_tub$`WT(GFP IP)_1st` == 0, ]

## d3_acto.Num
d3_acto.Num$WT_1st = as.numeric(round((d3_acto.Num$`WT(GFP IP)_1st`)/(d3_acto.Num$`WT(Ctrl IP)_1st`),2))
d3_acto.Num$CTL_1st = 1
d3_acto.Num = d3_acto.Num[,c(1,5,4)]
d3_acto.Num$mean.WT = as.numeric(round(d3_acto.Num$WT_1st,2))

## d3_acto.OnOff
d3_acto.OnOff$CTL_1st = ifelse(d3_acto.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d3_acto.OnOff$`WT(Ctrl IP)_1st`))
d3_acto.OnOff$WT_1st = ifelse(d3_acto.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                              ifelse(d3_acto.OnOff$CTL_1st == "Off", "On", as.numeric(round((d3_acto.OnOff$`WT(GFP IP)_1st`)/(d3_acto.OnOff$`WT(Ctrl IP)_1st`),2))))
d3_acto.OnOff = d3_acto.OnOff[,c(1,4,5)]
d3_acto.OnOff$mean.WT = "Off"
d3_acto.OnOff$mean.WT = ifelse(d3_acto.OnOff$WT_1st == "On", "On", "Off")

## Bind Num & OnOff -> Order
d3_acto = rbind.data.frame(d3_acto.Num, d3_acto.OnOff)
d3_acto.1 = d3_acto[d3_acto$mean.WT %in% c("On", "Off", "On/Off"),]
d3_acto.1 = d3_acto.1[order(d3_acto.1$mean.WT),]
d3_acto.2 = d3_acto[!(d3_acto$mean.WT %in% c("On", "Off", "On/Off")),]
d3_acto.2$mean.WT = as.numeric(d3_acto.2$mean.WT)
d3_acto.2= d3_acto.2[order(d3_acto.2$mean.WT, decreasing = TRUE), ]
d3_acto = rbind.data.frame(d3_acto.2, d3_acto.1)
rownames(d3_acto) = 1:nrow(d3_acto)
d3_acto = d3_acto[c(1:2,4:8),]

d3_acto.tb = d3_acto[,c(2:(ncol(d3_acto)-1),1)]
write_xlsx(d3_acto.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180625.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

d3_acto.tmp = d3_acto
rownames(d3_acto.tmp) = 1:nrow(d3_acto.tmp)
d3_acto.tmp$rank = rownames(d3_acto.tmp)
#d3_acto.tmp = d3_acto.tmp[d3_acto.tmp$`Gene name` != "Flna",]
dt.d3_acto <- d3_acto.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:WT_1st)
dt.d3_acto$rank = as.numeric(dt.d3_acto$rank)
dt.d3_acto = dt.d3_acto[order(dt.d3_acto$rank),]
dt.d3_acto$FC.color = ifelse(dt.d3_acto$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d3_acto$FC.color = ifelse(dt.d3_acto$Group %in% c("WT_1st", "WT_2st") & dt.d3_acto$FC == "Off", 0, dt.d3_acto$FC.color)
dt.d3_acto$FC.color = ifelse(dt.d3_acto$Group %in% c("WT_1st", "WT_2st") & dt.d3_acto$FC == "On", 5, dt.d3_acto$FC.color)
dt.d3_acto$FC.color = ifelse(dt.d3_acto$Group %in% c("WT_1st", "WT_2st") & dt.d3_acto$FC < 1, 0, dt.d3_acto$FC.color)
dt.d3_acto$FC.color = ifelse(is.na(dt.d3_acto$FC.color), dt.d3_acto$FC, dt.d3_acto$FC.color)
dt.d3_acto$FC.color = as.numeric(dt.d3_acto$FC.color)

p.d3_acto <- ggplot(dt.d3_acto, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d3_acto, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.9) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  
#"#7b9fd4", mid = 'white', high="#CD2836"
p.d3_acto

## d3_tub.Num
d3_tub.Num$WT_1st = as.numeric(round((d3_tub.Num$`WT(GFP IP)_1st`)/(d3_tub.Num$`WT(Ctrl IP)_1st`),2))
d3_tub.Num$CTL_1st = 1
d3_tub.Num = d3_tub.Num[,c(1,5,4)]
d3_tub.Num$mean.WT = as.numeric(round(d3_tub.Num$WT_1st,2))

## d3_tub.OnOff
d3_tub.OnOff$CTL_1st = ifelse(d3_tub.OnOff$`WT(Ctrl IP)_1st`== 0, "Off", as.numeric(d3_tub.OnOff$`WT(Ctrl IP)_1st`))
d3_tub.OnOff$WT_1st = ifelse(d3_tub.OnOff$`WT(GFP IP)_1st` == 0, "Off",
                              ifelse(d3_tub.OnOff$CTL_1st == "Off", "On", as.numeric(round((d3_tub.OnOff$`WT(GFP IP)_1st`)/(d3_tub.OnOff$`WT(Ctrl IP)_1st`),2))))
d3_tub.OnOff = d3_tub.OnOff[,c(1,4,5)]
d3_tub.OnOff$mean.WT = "Off"
d3_tub.OnOff$mean.WT = ifelse(d3_tub.OnOff$WT_1st == "On", "On", "Off")

## Bind Num & OnOff -> Order
d3_tub = rbind.data.frame(d3_tub.Num, d3_tub.OnOff)
d3_tub.1 = d3_tub[d3_tub$mean.WT %in% c("On", "Off", "On/Off"),]
d3_tub.1 = d3_tub.1[order(d3_tub.1$mean.WT),]
d3_tub.2 = d3_tub[!(d3_tub$mean.WT %in% c("On", "Off", "On/Off")),]
d3_tub.2$mean.WT = as.numeric(d3_tub.2$mean.WT)
d3_tub.2= d3_tub.2[order(d3_tub.2$mean.WT, decreasing = TRUE), ]
d3_tub = rbind.data.frame(d3_tub.2, d3_tub.1)

d3_tub.tb = d3_tub[,c(2:(ncol(d3_tub)-1),1)]
write_xlsx(d3_tub.tb, path = "../CD99_SYK_202201/OnOff.Cho_20180625.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

d3_tub.tmp = d3_tub
rownames(d3_tub.tmp) = 1:nrow(d3_tub.tmp)
d3_tub.tmp$rank = rownames(d3_tub.tmp)
#d3_tub.tmp = d3_tub.tmp[d3_tub.tmp$`Gene name` != "Flna",]
dt.d3_tub <- d3_tub.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:WT_1st)
dt.d3_tub$rank = as.numeric(dt.d3_tub$rank)
dt.d3_tub = dt.d3_tub[order(dt.d3_tub$rank),]
dt.d3_tub$FC.color = ifelse(dt.d3_tub$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d3_tub$FC.color = ifelse(dt.d3_tub$Group %in% c("WT_1st", "WT_2st") & dt.d3_tub$FC == "Off", 0, dt.d3_tub$FC.color)
dt.d3_tub$FC.color = ifelse(dt.d3_tub$Group %in% c("WT_1st", "WT_2st") & dt.d3_tub$FC == "On", 5, dt.d3_tub$FC.color)
dt.d3_tub$FC.color = ifelse(dt.d3_tub$Group %in% c("WT_1st", "WT_2st") & dt.d3_tub$FC < 1, 0, dt.d3_tub$FC.color)
dt.d3_tub$FC.color = ifelse(is.na(dt.d3_tub$FC.color), dt.d3_tub$FC, dt.d3_tub$FC.color)
dt.d3_tub$FC.color = as.numeric(dt.d3_tub$FC.color)

p.d3_tub <- ggplot(dt.d3_tub, aes(x=factor(Group, levels = c("CTL_1st", "WT_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d3_tub, aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        #legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  
#"#7b9fd4", mid = 'white', high="#CD2836"
p.d3_tub
legend <- get_legend(p.d3_tub)
p.d3_tub <- p.d3_tub + theme(legend.position = "none")
#plot_grid(p.d3_acto, p.d3_tub, legend, ncol = 3)  
p = plot_grid(p.d3_acto, p.d3_tub, ncol = 2)  
p
ggsave('Figures/Cho_20180625_OnOff.FC.Heatmap_0111.pdf', p, width =3.5, height = 4, useDingbats=F)


## 4. Kang_20140814
d4 = DF[,21:26]
colnames(d4) = d4[1,] 
d4= d4[2:nrow(d4),c(6,4:5)]
colnames(d4) = c("Gene name", "Myc(Ctrl IP)_1st", "Myc(GFP IP)_1st")
d4$`Myc(Ctrl IP)_1st` <- round(as.numeric(d4$`Myc(Ctrl IP)_1st`), 2)
d4$`Myc(GFP IP)_1st` <- round(as.numeric(d4$`Myc(GFP IP)_1st`), 2)

### Step 1. Divide by actomyosin & microtubule group
d4_acto = d4[d4$`Gene name` %in% p_actomyo,]
d4_tub = d4[d4$`Gene name` %in% p_tubulin,]

## Step 2. Divide d4_acto and d4_tub by Num & OnOff
d4_acto.Num = d4_acto[d4_acto$`Myc(Ctrl IP)_1st` != 0 & d4_acto$`Myc(GFP IP)_1st` != 0, ]
d4_acto.OnOff = d4_acto[d4_acto$`Myc(Ctrl IP)_1st` == 0 | d4_acto$`Myc(GFP IP)_1st` == 0, ]
d4_tub.Num = d4_tub[d4_tub$`Myc(Ctrl IP)_1st` != 0 & d4_tub$`Myc(GFP IP)_1st` != 0, ]
d4_tub.OnOff = d4_tub[d4_tub$`Myc(Ctrl IP)_1st` == 0 | d4_tub$`Myc(GFP IP)_1st` == 0, ]

## d4_acto.Num
d4_acto.Num$Myc_1st = as.numeric(round((d4_acto.Num$`Myc(GFP IP)_1st`)/(d4_acto.Num$`Myc(Ctrl IP)_1st`),2))
d4_acto.Num$CTL_1st = 1
d4_acto.Num = d4_acto.Num[,c(1,5,4)]
d4_acto.Num$mean.Myc = as.numeric(round(d4_acto.Num$Myc_1st,2))

## d4_acto.OnOff
d4_acto.OnOff$CTL_1st = ifelse(d4_acto.OnOff$`Myc(Ctrl IP)_1st`== 0, "Off", as.numeric(d4_acto.OnOff$`Myc(Ctrl IP)_1st`))
d4_acto.OnOff$Myc_1st = ifelse(d4_acto.OnOff$`Myc(GFP IP)_1st` == 0, "Off",
                              ifelse(d4_acto.OnOff$CTL_1st == "Off", "On", as.numeric(round((d4_acto.OnOff$`Myc(GFP IP)_1st`)/(d4_acto.OnOff$`Myc(Ctrl IP)_1st`),2))))
d4_acto.OnOff = d4_acto.OnOff[,c(1,4,5)]
d4_acto.OnOff$mean.Myc = "Off"
d4_acto.OnOff$mean.Myc = ifelse(d4_acto.OnOff$Myc_1st == "On", "On", "Off")

## Bind Num & OnOff -> Order
d4_acto = rbind.data.frame(d4_acto.Num, d4_acto.OnOff)
d4_acto.1 = d4_acto[d4_acto$mean.Myc %in% c("On", "Off", "On/Off"),]
d4_acto.1 = d4_acto.1[order(d4_acto.1$mean.Myc),]
d4_acto.2 = d4_acto[!(d4_acto$mean.Myc %in% c("On", "Off", "On/Off")),]
d4_acto.2$mean.Myc = as.numeric(d4_acto.2$mean.Myc)
d4_acto.2= d4_acto.2[order(d4_acto.2$mean.Myc, decreasing = TRUE), ]
d4_acto = rbind.data.frame(d4_acto.2, d4_acto.1)
d4_acto = d4_acto[!(d4_acto$CTL_1st == "Off" & d4_acto$Myc_1st == "Off"),]

d4_acto.tb = d4_acto[,c(2:(ncol(d4_acto)-1),1)]
write_xlsx(d4_acto.tb, path = "../CD99_SYK_202201/OnOff.Kang_20140814.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

d4_acto.tmp = d4_acto
d4_acto.tmp$rank = rownames(d4_acto.tmp)
#d4_acto.tmp = d4_acto.tmp[d4_acto.tmp$`Gene name` != "Flna",]
dt.d4_acto <- d4_acto.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:Myc_1st)
dt.d4_acto$rank = as.numeric(dt.d4_acto$rank)
dt.d4_acto = dt.d4_acto[order(dt.d4_acto$rank),]
dt.d4_acto$FC.color = ifelse(dt.d4_acto$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d4_acto$FC.color = ifelse(dt.d4_acto$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_acto$FC == "Off", 0, dt.d4_acto$FC.color)
dt.d4_acto$FC.color = ifelse(dt.d4_acto$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_acto$FC == "On", 5, dt.d4_acto$FC.color)
dt.d4_acto$FC.color = ifelse(dt.d4_acto$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_acto$FC < 1, 0, dt.d4_acto$FC.color)
dt.d4_acto$FC.color = ifelse(is.na(dt.d4_acto$FC.color), dt.d4_acto$FC, dt.d4_acto$FC.color)
dt.d4_acto$FC.color = as.numeric(dt.d4_acto$FC.color)

p.d4_acto1 <- ggplot(dt.d4_acto[1:56,], aes(x=factor(Group, levels = c("CTL_1st", "Myc_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_acto[1:56,], aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none")
#        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  
#"#7b9fd4", mid = 'white', high="#CD2836"
p.d4_acto2 <- ggplot(dt.d4_acto[57:118,], aes(x=factor(Group, levels = c("CTL_1st", "Myc_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_acto[57:118,], aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none")
#        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

library(patchwork)
p.d4_acto <- p.d4_acto1 + p.d4_acto2
p.d4_acto <- p.d4_acto +
  labs(title = 'Actomyosin',x = '', y = '') + 
  theme(plot.title = element_text(size = 12, face = "bold", vjust = 1, hjust = 4.5))
p.d4_acto

## d4_tub.Num
d4_tub.Num$Myc_1st = as.numeric(round((d4_tub.Num$`Myc(GFP IP)_1st`)/(d4_tub.Num$`Myc(Ctrl IP)_1st`),2))
d4_tub.Num$CTL_1st = 1
d4_tub.Num = d4_tub.Num[,c(1,5,4)]
d4_tub.Num$mean.Myc = as.numeric(round(d4_tub.Num$Myc_1st,2))

## d4_tub.OnOff
d4_tub.OnOff$CTL_1st = ifelse(d4_tub.OnOff$`Myc(Ctrl IP)_1st`== 0, "Off", as.numeric(d4_tub.OnOff$`Myc(Ctrl IP)_1st`))
d4_tub.OnOff$Myc_1st = ifelse(d4_tub.OnOff$`Myc(GFP IP)_1st` == 0, "Off",
                             ifelse(d4_tub.OnOff$CTL_1st == "Off", "On", as.numeric(round((d4_tub.OnOff$`Myc(GFP IP)_1st`)/(d4_tub.OnOff$`Myc(Ctrl IP)_1st`),2))))
d4_tub.OnOff = d4_tub.OnOff[,c(1,4,5)]
d4_tub.OnOff$mean.Myc = "Off"
d4_tub.OnOff$mean.Myc = ifelse(d4_tub.OnOff$Myc_1st == "On", "On", "Off")

## Bind Num & OnOff -> Order
d4_tub = rbind.data.frame(d4_tub.Num, d4_tub.OnOff)
d4_tub.1 = d4_tub[d4_tub$mean.Myc %in% c("On", "Off", "On/Off"),]
d4_tub.1 = d4_tub.1[order(d4_tub.1$mean.Myc),]
d4_tub.2 = d4_tub[!(d4_tub$mean.Myc %in% c("On", "Off", "On/Off")),]
d4_tub.2$mean.Myc = as.numeric(d4_tub.2$mean.Myc)
d4_tub.2= d4_tub.2[order(d4_tub.2$mean.Myc, decreasing = TRUE), ]
d4_tub = rbind.data.frame(d4_tub.2, d4_tub.1)
d4_tub = d4_tub[c(1:11, 18:65),]

d4_tub.tb = d4_tub[,c(2:(ncol(d4_tub)-1),1)]
write_xlsx(d4_tub.tb, path = "../CD99_SYK_202201/OnOff.Kang_20140814.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

d4_tub.tmp = d4_tub
rownames(d4_tub.tmp) = 1:nrow(d4_tub.tmp)
d4_tub.tmp$rank = rownames(d4_tub.tmp)
#d4_tub.tmp = d4_tub.tmp[d4_tub.tmp$`Gene name` != "Flna",]
dt.d4_tub <- d4_tub.tmp %>%
  #    rownames_to_column() %>%
  gather(Group, FC, CTL_1st:Myc_1st)
dt.d4_tub$rank = as.numeric(dt.d4_tub$rank)
dt.d4_tub = dt.d4_tub[order(dt.d4_tub$rank),]
dt.d4_tub$FC.color = ifelse(dt.d4_tub$Group %in% c("CTL_1st", "CTL_2st"), 0, NA)
dt.d4_tub$FC.color = ifelse(dt.d4_tub$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_tub$FC == "Off", 0, dt.d4_tub$FC.color)
dt.d4_tub$FC.color = ifelse(dt.d4_tub$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_tub$FC == "On", 5, dt.d4_tub$FC.color)
dt.d4_tub$FC.color = ifelse(dt.d4_tub$Group %in% c("Myc_1st", "Myc_2st") & dt.d4_tub$FC < 1, 0, dt.d4_tub$FC.color)
dt.d4_tub$FC.color = ifelse(is.na(dt.d4_tub$FC.color), dt.d4_tub$FC, dt.d4_tub$FC.color)
dt.d4_tub$FC.color = as.numeric(dt.d4_tub$FC.color)

p.d4_tub.1 <- ggplot(dt.d4_tub[1:56,], aes(x=factor(Group, levels = c("CTL_1st", "Myc_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_tub[1:56,], aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

p.d4_tub.2 <- ggplot(dt.d4_tub[57:118,], aes(x=factor(Group, levels = c("CTL_1st", "Myc_1st")), y=reorder(`Gene name`, desc(rank)), fill=FC.color)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  #  scale_fill_gradientn(colors = c("#5779c9", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  scale_fill_gradientn(colors = c("#5387d4", "white", "#cc333f"), values=rescale(c(0,1,5)), limits = c(0, 5), oob=squish) +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_tub[57:118,], aes(label =FC), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.7) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.position="none",
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

p.d4_tub <- p.d4_tub.1 + p.d4_tub.2
p.d4_tub <- p.d4_tub +
  labs(title = 'Microtubule',x = '', y = '') + 
  theme(plot.title = element_text(size = 12, face = "bold", vjust = 1, hjust = 4))
p.d4_tub

p = plot_grid(p.d4_acto, p.d4_tub, ncol = 2)  
p
ggsave('Figures/Kang_20140814_OnOff.FC.Heatmap_0111.pdf', p, width =7, height = 8, useDingbats=F)

