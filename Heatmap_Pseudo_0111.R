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

d1_acto = d1[d1$`Gene name` %in% p_actomyo,]
d1_acto$WT_1st = (d1_acto$`WT(GFP IP)_1st`+1)/(d1_acto$`WT(Ctrl IP)_1st`+1)  
d1_acto$WT_2nd = (d1_acto$`WT(GFP IP)_2st`+1)/(d1_acto$`WT(Ctrl IP)_2st`+1)  
d1_acto2 = d1_acto[,c(1,6,7)]
d1_acto2$meanWT = (d1_acto2$WT_1st + d1_acto2$WT_2nd)/2
d1_acto2 = d1_acto2[order(d1_acto2$meanWT, decreasing = T),]
d1_acto2 = d1_acto2[,c(2,3,1)]

write_xlsx(d1_acto2, path = "../CD99_SYK_202201/Pseudo.Cho_20180514.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)


dt.d1_acto <- d1_acto2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st:WT_2nd)

p.d1_acto <- ggplot(dt.d1_acto, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
#  scale_fill_gradient(low="white", high="#Cd1836") +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, max(dt.d1_acto$FC))), 
                       limits = c(1, max(dt.d1_acto$FC)), oob=squish) +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d1_acto, aes(label = round(FC,2)), color = "black", size = 2.3) +
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

d1_tub = d1[d1$`Gene name` %in% p_tubulin,]
d1_tub$WT_1st = (d1_tub$`WT(GFP IP)_1st`+1)/(d1_tub$`WT(Ctrl IP)_1st`+1)  
d1_tub$WT_2nd = (d1_tub$`WT(GFP IP)_2st`+1)/(d1_tub$`WT(Ctrl IP)_2st`+1)  
d1_tub2 = d1_tub[,c(1,6,7)]
d1_tub2$meanWT = (d1_tub2$WT_1st + d1_tub2$WT_2nd)/2
d1_tub2 = d1_tub2[order(d1_tub2$meanWT, decreasing = T),]
d1_tub2 = d1_tub2[,c(1:3)]
d1_tub2 = d1_tub2[!(d1_tub2$WT_1st == 1 & d1_tub2$WT_2nd == 1),]
d1_tub2 = d1_tub2[,c(2:3,1)]
write_xlsx(d1_tub2, path = "../CD99_SYK_202201/Pseudo.Cho_20180514.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d1_tub <- d1_tub2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st:WT_2nd)

p.d1_tub <- ggplot(dt.d1_tub, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
#  scale_fill_gradient(low="white", high="#Cd1836") +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, max(dt.d1_tub$FC))), 
                       limits = c(1, max(dt.d1_tub$FC)), oob=squish) +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d1_tub, aes(label = round(FC,2)), color = "black", size = 2.3) +
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
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

legend <- get_legend(p.d1_tub)
p.d1_tub <- p.d1_tub + theme(legend.position = "none")
#plot_grid(p.d1_acto, p.d1_tub, legend, ncol = 3)  
p = plot_grid(p.d1_acto, p.d1_tub, ncol = 2)  

ggsave('Figures/Cho_20180514_pseudocount.FC.Heatmap.pdf', p, width =3.5, height = 5, useDingbats=F)


## 2. Cho_20180912
d2 = DF[,9:15]
colnames(d2) = d2[1,]
d2= d2[2:nrow(d2),c(2,4:ncol(d2))]
d2$WT_Con_1st <- round(as.numeric(d2$WT_Con_1st), 2)
d2$WT_Con_2nd <- round(as.numeric(d2$WT_Con_2nd), 2)
d2$WT_GFP_1st <- round(as.numeric(d2$WT_GFP_1st), 2)
d2$WT_GFP_2nd <- round(as.numeric(d2$WT_GFP_2nd), 2)

d2_acto = d2[d2$`Gene name` %in% p_actomyo,]
d2_acto$WT_1st = (d2_acto$`WT_GFP_1st`+1)/(d2_acto$`WT_Con_1st`+1)  
d2_acto$WT_2nd = (d2_acto$`WT_GFP_2nd`+1)/(d2_acto$`WT_Con_2nd`+1)  
d2_acto2 = d2_acto[,c(1,6,7)]
d2_acto2 = d2_acto2[!(d2_acto2$WT_1st == "1" & d2_acto2$WT_2nd == "1"),]
d2_acto2$meanWT = (d2_acto2$WT_1st + d2_acto2$WT_2nd)/2
d2_acto2 = d2_acto2[order(d2_acto2$meanWT, decreasing = T),]
d2_acto2 = d2_acto2[,c(1:3)]
d2_acto2 = d2_acto2[,c(2:3,1)]

write_xlsx(d2_acto2, path = "../CD99_SYK_202201/Pseudo.Cho_20180912.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d2_acto <- d2_acto2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st:WT_2nd)

p.d2_acto <- ggplot(dt.d2_acto, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, max(dt.d2_acto$FC))), 
                       limits = c(1, max(dt.d2_acto$FC)), oob=squish) +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d2_acto, aes(label = round(FC,2)), color = "black", size = 2.3) +
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

d2_tub = d2[d2$`Gene name` %in% p_tubulin,]
d2_tub$WT_1st = (d2_tub$`WT_GFP_1st`+1)/(d2_tub$`WT_Con_1st`+1)  
d2_tub$WT_2nd = (d2_tub$`WT_GFP_2nd`+1)/(d2_tub$`WT_Con_2nd`+1)  
d2_tub2 = d2_tub[,c(1,6,7)]
d2_tub2 = d2_tub2[!(d2_tub2$WT_1st == "1" & d2_tub2$WT_2nd == "1"),]
d2_tub2$meanWT = (d2_tub2$WT_1st + d2_tub2$WT_2nd)/2
d2_tub2 = d2_tub2[order(d2_tub2$meanWT, decreasing = T),]
d2_tub2 = d2_tub2[,c(1:3)]
d2_tub2 = d2_tub2[,c(2:3,1)]

write_xlsx(d2_tub2, path = "../CD99_SYK_202201/Pseudo.Cho_20180912.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)


dt.d2_tub <- d2_tub2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st:WT_2nd)

p.d2_tub <- ggplot(dt.d2_tub, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, max(dt.d2_tub$FC))), 
                       limits = c(1, max(dt.d2_tub$FC)), oob=squish) +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d2_tub, aes(label = round(FC,2)), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.5) +
  theme_minimal(base_size=7) + 
  theme(axis.ticks=element_blank(), 
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(size=8, angle=320, hjust = 0, colour="black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 0.5),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

legend <- get_legend(p.d2_tub)
p.d2_tub <- p.d2_tub + theme(legend.position = "none")
#plot_grid(p.d2_acto, p.d2_tub, legend, ncol = 3)  
p = plot_grid(p.d2_acto, p.d2_tub, ncol = 2)  

#ggsave('Figures/Cho_20180912_pseudocount.FC.Heatmap.pdf', p, width =4, height = 7.5, useDingbats=F)



## 3. Cho_20180625
d3 = DF[,17:19]
colnames(d3) = d3[1,] ## No Symbol. Merge from Cho 0912
d2 = DF[,9:15]
colnames(d2) = d2[1,]
d2.symbol = d2[2:nrow(d2),1:2]

d3 = merge(d3, d2.symbol, by = "Protein name")
d3= d3[,c(4,2:3)]
d3$WT_Ctrl_1st <- round(as.numeric(d3$WT_Ctrl_1st), 2)
d3$WT_GFP_1st <- round(as.numeric(d3$WT_GFP_1st), 2)
# d3$gene_name = toupper(d3$`Gene name`)
# d3$`Gene name` <- NULL

d3_acto = d3[d3$`Gene name` %in% p_actomyo,]
d3_acto$WT_1st = (d3_acto$`WT_GFP_1st`+1)/(d3_acto$`WT_Ctrl_1st`+1)  
d3_acto2 = d3_acto[d3_acto$`Gene name` != "Flna",c(1,4)]
d3_acto2 = d3_acto2[order(d3_acto2$WT_1st, decreasing = T),]
d3_acto2 = d3_acto2[,c(2,1)]

write_xlsx(d3_acto2, path = "../CD99_SYK_202201/Pseudo.Cho_20180625.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d3_acto <- d3_acto2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st)

p.d3_acto <- ggplot(dt.d3_acto, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient(low="white", high="#Cd3836") +
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d3_acto, aes(label = round(FC,2)), color = "black", size = 2.3) +
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

d3_tub = d3[d3$`Gene name` %in% p_tubulin,]
d3_tub$WT_1st = (d3_tub$`WT_GFP_1st`+1)/(d3_tub$`WT_Ctrl_1st`+1)  
d3_tub2 = d3_tub[,c(1,4)]
d3_tub2 = d3_tub2[order(d3_tub2$WT_1st, decreasing = T),]
d3_tub2 = d3_tub2[,c(2,1)]
write_xlsx(d3_tub2, path = "../CD99_SYK_202201/Pseudo.Cho_20180625.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d3_tub <- d3_tub2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, WT_1st)

p.d3_tub <- ggplot(dt.d3_tub, aes(x=Group, y=reorder(`Gene name`,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient(low="white", high="#Cd3836") +
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d3_tub, aes(label = round(FC,2)), color = "black", size = 2.3) +
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
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0.2, 0, 0.5, 0.2), "cm"))  

legend <- get_legend(p.d3_tub)
p.d3_tub <- p.d3_tub + theme(legend.position = "none")
#plot_grid(p.d3_acto, p.d3_tub, legend, ncol = 3)  
p = plot_grid(p.d3_acto, p.d3_tub, ncol = 2)  
p
#ggsave('Figures/Cho_20180625_pseudocount.FC.Heatmap.pdf', p, width =3, height = 4, useDingbats=F)



## 4. Kang_20140814
d4 = DF[,21:26]
colnames(d4) = d4[1,] ## No Symbol. Merge from Cho 0912
d4= d4[2:nrow(d4),c(6,4:5)]
d4 = d4[!is.na(d4$Symbol),]
d4$`Myc CD99 control` <- round(as.numeric(d4$`Myc CD99 control`), 2)
d4$`Myc CD99 IP` <- round(as.numeric(d4$`Myc CD99 IP`), 2)
# d4$gene_name = toupper(d4$`Gene name`)
# d4$`Gene name` <- NULL

d4_acto = d4[d4$Symbol %in% p_actomyo,]
d4_acto$Myc_1st = (d4_acto$`Myc CD99 IP`+1)/(d4_acto$`Myc CD99 control`+1)  
d4_acto2 = d4_acto[!(d4_acto$`Myc CD99 control` == "0" & d4_acto$`Myc CD99 IP` == "0"),]
d4_acto2 = d4_acto2[,c(1,4)]
d4_acto2 = d4_acto2[order(d4_acto2$Myc_1st, decreasing = T),]
d4_acto2 = d4_acto2[,c(2,1)]
write_xlsx(d4_acto2, path = "../CD99_SYK_202201/Pseudo.Kang_20140814.Actomyosin_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d4_acto <- d4_acto2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, Myc_1st)
dt.d4_acto = dt.d4_acto[order(dt.d4_acto$FC,decreasing = TRUE),]

p.d4_acto1 <- ggplot(dt.d4_acto[1:30,], aes(x=Group, y=reorder(Symbol,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
                       limits = c(1, 15), oob=squish) +  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_acto[1:30,], aes(label = round(FC,2)), color = "black", size = 2.3) +
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

p.d4_acto2 <- ggplot(dt.d4_acto[31:59,], aes(x=Group, y=reorder(Symbol,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
                       limits = c(1, 15), oob=squish) +  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_acto[31:59,], aes(label = round(FC,2)), color = "black", size = 2.3) +
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

library(patchwork)
p.d4_acto <- p.d4_acto1 + p.d4_acto2
p.d4_acto <- p.d4_acto +
  labs(title = 'Actomyosin',x = '', y = '') + 
  theme(plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 1.5))
p.d4_acto


d4_tub = d4[d4$Symbol %in% p_tubulin,]
d4_tub$Myc_1st = (d4_tub$`Myc CD99 IP`+1)/(d4_tub$`Myc CD99 control`+1)  
d4_tub2 = d4_tub[!(d4_tub$`Myc CD99 control` == "0" & d4_tub$`Myc CD99 IP` == "0"),]
d4_tub2 = d4_tub2[,c(1,4)]
d4_tub2 = d4_tub2[order(d4_tub2$Myc_1st, decreasing = TRUE),]

d4_tub2 = d4_tub2[order(d4_tub2$Myc_1st, decreasing = T),]
d4_tub2 = d4_tub2[,c(2,1)]
write_xlsx(d4_tub2, path = "../CD99_SYK_202201/Pseudo.Kang_20140814.Microtubule_FoldChange_Replicates.xlsx", col_names = TRUE)

dt.d4_tub <- d4_tub2 %>%
  #    rownames_to_column() %>%
  gather(Group, FC, Myc_1st)

p.d4_tub1 <- ggplot(dt.d4_tub[1:30,], aes(x=Group, y=reorder(Symbol,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
                       limits = c(1, 15), oob=squish) +  labs(title = '', x = '', y = '') +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_tub[1:30,], aes(label = round(FC,2)), color = "black", size = 2.3) +
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

p.d4_tub2 <- ggplot(dt.d4_tub[31:59,], aes(x=Group, y=reorder(Symbol,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradientn(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
                       limits = c(1, 15), oob=squish) +  labs(title = '', x = '', y = '') +
  labs(title = '', x = '', y = '') +
  geom_text(data= dt.d4_tub[31:59,], aes(label = round(FC,2)), color = "black", size = 2.3) +
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

p.d4_tub <- p.d4_tub1 + p.d4_tub2
p.d4_tub <- p.d4_tub +
  labs(title = 'Microtubule',x = '', y = '') + 
  theme(plot.title = element_text(size = 11, face = "bold", vjust = 1, hjust = 1.5))
p.d4_tub

legend <- get_legend(p.d4_tub)
p.d4_tub <- p.d4_tub + theme(legend.position = "none")
#plot_grid(p.d4_acto, p.d4_tub, legend, ncol = 3)  
p = plot_grid(p.d4_acto, p.d4_tub, ncol = 2)  

#ggsave('Figures/Kang_20140814_pseudocount.FC.Heatmap.pdf', p, width =5, height = 8, useDingbats=F)


## 20180912 Cyto, TM Mutant

d = read.delim("Tables/table.edgeR_FC_Mutant.CD99_20180912.txt")

d_acto = d[d$gene_name %in% p_actomyo,c(16:18)]
dt.d_acto <- d_acto %>%
  #    rownames_to_column() %>%
  gather(key = "Group", value= "FC", FC_cyto, FC_tm)
dt.d_acto = dt.d_acto[order(dt.d_acto$FC,decreasing = TRUE),]

d_tub = d[d$gene_name %in% p_tubulin,c(16:18)]
dt.d_tub <- d_tub %>%
  #    rownames_to_column() %>%
  gather(key = "Group", value= "FC", FC_cyto, FC_tm)
dt.d_tub = dt.d_tub[order(dt.d_tub$FC,decreasing = TRUE),]

p.d3_acto <- ggplot(dt.d_acto, aes(x=Group, y=reorder(gene_name,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient2(low = "#5387d4", mid = "white" ,high =  "#cc333f", midpoint = 0) +
  # scale_fill_gradient2(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
  #                      limits = c(min(dt.d_acto$FC), max(dt.d_acto$FC)), oob=squish) +  
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d_acto, aes(label = round(FC,2)), color = "black", size = 2.3) +
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

p.d3_tub <- ggplot(dt.d_tub, aes(x=Group, y=reorder(gene_name,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient2(low = "#5387d4", mid = "white" ,high =  "#cc333f", midpoint = 0) +
  # scale_fill_gradient2(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
  #                      limits = c(min(dt.d_tub$FC), max(dt.d_tub$FC)), oob=squish) +  
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d_tub, aes(label = round(FC,2)), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.5) +
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

p = plot_grid(p.d3_acto, p.d3_tub, ncol = 2)  
p
ggsave('Figures/Cho_20180912_FC.Mutant.Heatmap.pdf', p, width =4, height = 8, useDingbats=F)



## 20180514 Cyto, TM Mutant

d = read.delim("Tables/table.edgeR_FC_Mutant.CD99_20180514.txt")

d_acto = d[d$gene_name %in% p_actomyo,c(15:17)]
dt.d_acto <- d_acto %>%
  #    rownames_to_column() %>%
  gather(key = "Group", value= "FC", FC_cyto, FC_tm)
dt.d_acto = dt.d_acto[order(dt.d_acto$FC,decreasing = TRUE),]

d_tub = d[d$gene_name %in% p_tubulin,c(16:18)]
dt.d_tub <- d_tub %>%
  #    rownames_to_column() %>%
  gather(key = "Group", value= "FC", FC_cyto, FC_tm)
dt.d_tub = dt.d_tub[order(dt.d_tub$FC,decreasing = TRUE),]

p.d3_acto <- ggplot(dt.d_acto, aes(x=Group, y=reorder(gene_name,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient2(low = "#5387d4", mid = "white" ,high =  "#cc333f", midpoint = 0) +
  # scale_fill_gradient2(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
  #                      limits = c(min(dt.d_acto$FC), max(dt.d_acto$FC)), oob=squish) +  
  labs(title = 'Actomyosin', x = '', y = '') +
  geom_text(data= dt.d_acto, aes(label = round(FC,2)), color = "black", size = 2.3) +
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

p.d3_tub <- ggplot(dt.d_tub, aes(x=Group, y=reorder(gene_name,FC), fill=FC)) +
  geom_tile(color = 'white', size = 0.5, aes(width = 1, height = 1)) +
  scale_fill_gradient2(low = "#5387d4", mid = "white" ,high =  "#cc333f", midpoint = 0) +
  # scale_fill_gradient2(colors = c("white", "#Cd1836"), values=rescale(c(1, 15)), 
  #                      limits = c(min(dt.d_tub$FC), max(dt.d_tub$FC)), oob=squish) +  
  labs(title = 'Microtubule', x = '', y = '') +
  geom_text(data= dt.d_tub, aes(label = round(FC,2)), color = "black", size = 2.3) +
  coord_fixed(ratio = 0.5) +
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

p = plot_grid(p.d3_acto, p.d3_tub, ncol = 2)  
p
ggsave('Figures/Cho_20180514_FC.Mutant.Heatmap.pdf', p, width =5, height = 8, useDingbats=F)




