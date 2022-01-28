library(edgeR)
library(readxl)
library(ggrepel)
library(dplyr)
library(tidyr)

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

########################################################

# Cho 20180912
## Load data (20180912) - It has TM and Cyto mutants
d2018 = read_excel("Data/IP_raw_data_all.xlsx", sheet = 6)
d2018$gene_name <- d2018$`Gene name`
d2018$`Gene name` <- NULL
d2018$uniprot_id <- d2018$`Uniprot ID`
d2018$`Uniprot ID` <- NULL
d2018$Mock_GFP_1st <- NULL
d2018$Mock_GFP_2nd <- NULL

## To compare a WT GFP network, we apply differentially expressed proteins using edgeR
colnames(d2018)
counts = d2018[,grepl('1st|2nd', colnames(d2018))]
rownames(counts) <- d2018$uniprot_id
# d0 <- calcNormFactors(DGEList(counts))
d <- DGEList(counts)
# cutoff <- 0
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
dim(d) 
snames <- colnames(counts) # Sample names
group <- gsub('_1st|_2nd', '', snames)
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~ 0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupWT_GFP - groupWT_Con, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$uniprot_id <- rownames(top.table)

m <- merge(top.table, d2018, by='uniprot_id', all.x = T)

## Create background interacting proteins for CD99
#pm = m[m$logFC > 0 & m$adj.P.Val <= 0.05,]
pm = m

pm$nm_WT_GFP_1 <- (pm$WT_GFP_1st + 1) / (pm$WT_Con_1st + 1)
pm$nm_Cyto_GFP_1 <- (pm$Cyto_GFP_1st + 1)/ (pm$WT_Con_1st + 1)
pm$nm_TM_GFP_1 <- (pm$TM_GFP_1st + 1) / (pm$WT_Con_1st + 1)

pm$nm_WT_GFP_2 <- (pm$WT_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)
pm$nm_Cyto_GFP_2 <- (pm$Cyto_GFP_2nd + 1)/ (pm$WT_Con_2nd + 1)
pm$nm_TM_GFP_2 <- (pm$TM_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)

# pm$isCytoskeleton <- ifelse(toupper(pm$gene_name) %in% toupper(c(p_myosin, p_actin, p_tubulin, p_cyto)), 
#                             'Y', 'N')
# pm$isCytoplasm <- ifelse( toupper(pm$gene_name) %in% toupper(p_ctp), 'Y', 'N')
# pm$isCortex <- ifelse( toupper(pm$gene_name) %in% toupper(p_ccrt), 'Y', 'N')
# pm$isTCR <- ifelse( toupper(pm$gene_name) %in% toupper(p_tcr), 'Y', 'N')

## Fold Change compared to WT
pm$fc_Cyto_GFP_1 <- (pm$nm_Cyto_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_Cyto_GFP_2 <- (pm$nm_Cyto_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
pm$fc_TM_GFP_1 <- (pm$nm_TM_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_TM_GFP_2 <- (pm$nm_TM_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
write.table(pm, 'Tables/table.edgeR_WT_GFP.CD99_20190815.txt', sep='\t', quote = F, row.names = F, col.names=T)

## ------------------------------ vv
## To create CT and TM proteins, we apply differentially expressed proteins using edgeR
contr <- makeContrasts(groupCyto_GFP - groupWT_GFP, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_cyto <- topTable(tmp, sort.by = "P", n = Inf)
top.table_cyto$uniprot_id <- rownames(top.table_cyto)  
m_cyto <- merge(top.table_cyto, d2018, by='uniprot_id', all.x = T)

contr <- makeContrasts(groupTM_GFP - groupWT_GFP, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_tm <- topTable(tmp, sort.by = "P", n = Inf)
top.table_tm$uniprot_id <- rownames(top.table_tm)  
m_tm <- merge(top.table_tm, d2018, by='uniprot_id', all.x = T)

m <- merge(setNames(top.table_cyto[,c('logFC', 'adj.P.Val', 'uniprot_id')], 
                    c('logFC_cyto', 'padj_cyto', 'uniprot_id')), 
           setNames(top.table_tm[,c('logFC', 'adj.P.Val', 'uniprot_id')], 
                    c('logFC_tm', 'padj_tm', 'uniprot_id')), by='uniprot_id') 
m1 = merge(m[m$uniprot_id %in% pm$uniprot_id,], d2018, by='uniprot_id', all.x = T)

m1$FC_cyto <- log2(( ((m1$Cyto_GFP_1st + 1)/ (m1$WT_GFP_1st + 1)) + ((m1$Cyto_GFP_2nd + 1)/ (m1$WT_GFP_2nd + 1)))/2)
m1$FC_tm <- log2(( ((m1$TM_GFP_1st + 1)/ (m1$WT_GFP_1st + 1)) + ((m1$TM_GFP_2nd + 1)/ (m1$WT_GFP_2nd + 1)))/2)
write.table(m1, 'Tables/table.edgeR_FC_Mutant.CD99_20180912.txt', sep='\t', quote = F, row.names = F, col.names=T)


m1$Sig <- ifelse(m1$padj_cyto <= 0.05  & m1$padj_tm <= 0.05, 'Both', 
                 ifelse(m1$padj_cyto <= 0.05, 'Cyto', ifelse(m1$padj_tm <= 0.05, 'TM', 'None')))

m1$Sig2 <- 'None'
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm > 0, 'Both Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm < 0, 'Both Down', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm < 0, 'Cyto Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm > 0, 'Cyto Down', m1$Sig2)

m1$Sig2 <- ifelse(m1$Sig %in% c('Cyto') & m1$Sig2=='None', ifelse(m1$logFC_cyto > 0, 'Cyto Up', 'Cyto Down'), m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('TM') & m1$Sig2=='None', ifelse(m1$logFC_tm > 0, 'TM Up', 'TM Down'), m1$Sig2)

m1$isSig_cyto <- ifelse(m1$padj_cyto <= 0.05, ifelse( m1$logFC_cyto >0 , 'Up', 'Down'), 'None')
m1$isSig_tm <- ifelse(m1$padj_tm <= 0.05, ifelse( m1$logFC_tm >0 , 'Up', 'Down'), 'None')

m1$Group <- ifelse(m1$gene_name %in% p_actomyo, 'Acto', 
                    ifelse(m1$gene_name %in% p_tubulin, 'Tubule','other'))


## Plot: Compare the fold change between CT and TM mutants
mm <- gather(m1[m1$Group %in% c('Acto', 'Tubule'), c('FC_cyto', 'FC_tm', 'gene_name')], 
             group, val, FC_cyto:FC_tm, -gene_name, factor_key=TRUE)
mm = mm[order(mm$val, decreasing = T),]
mm$gene_name = factor(mm$gene_name, levels=unique(mm$gene_name))
mm$group2 <- ifelse(grepl('tm', mm$group), 'TM', 'Cyto')

####

mm1 <- mm[mm$group2 == "Cyto",]
mm1 = mm1[order(mm1$val),]
mm1$rank = 1:nrow(mm1)
mm2 <- mm[mm$group2 == "TM",]
mm2 = mm2[order(mm2$val),]
mm2$rank = 1:nrow(mm2)
mm = rbind.data.frame(mm1, mm2)

gs1 = c('Iqgap1') 
gs2 = c('Myh9')
p4 <- ggplot(mm, aes(x=rank, y=val, group = group2)) + 
#  geom_violin(aes(fill = group2), alpha=0.3, scale = "width") +
#  geom_point(color = 'white', fill = 'black', alpha = 3) + 
  geom_point(color = 'black', 
             fill = 'white', size = 0.75,stroke = 1) +
  facet_wrap(~group2, strip.position = "bottom") +
  # scale_x_continuous(breaks = 1:2, labels = unique(counts_df$line)) +
#  geom_line(data = mm[mm$gene_name %in% gs1,], aes(group=gene_name), color='red', alpha=0.3) + 
#  geom_line(data = mm[mm$gene_name %in% gs2,], aes(group=gene_name), color='blue', alpha=0.3) + 
  geom_point(data = mm[mm$gene_name %in% gs1,], color='red', size = 3) + 
  geom_point(data = mm[mm$gene_name %in% gs2,], color='blue', size = 3) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(data = mm[mm$gene_name %in% gs1,], aes(label=gene_name), 
                  color='red',
                  size = 4, nudge_x = -15) +
  geom_text_repel(data = mm[mm$gene_name %in% gs2,], aes(label=gene_name), 
                  color='blue',
                  size = 4, nudge_x=20) +
  theme_linedraw(base_size = 10) + 
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"),
        strip.background = element_rect(colour="black",fill="white"),
        panel.grid = element_line(colour = "gray", size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = '', y='Log fold change compared to WT')
p4


###########

# Cho 20180514
## Load data (20180912) - It has TM and Cyto mutants
d2018 = read_excel("Data/IP_raw_data_all.xlsx", sheet = 4)
d2018 = d2018[,c(1:11, 13:14)]
d2018$gene_name <- d2018$`Gene name`
d2018$`Gene name` <- NULL
d2018$uniprot_id <- d2018$`Uniprot ID`
d2018$`Uniprot ID` <- NULL
colnames(d2018)[2:11] = c("WT_Con_1st", "WT_Con_2nd", "WT_GFP_1st", "WT_GFP_2nd", "Cyto_GFP_1st", 
                          "Cyto_GFP_2nd", "TM_GFP_1st", "TM_GFP_2nd", "Mock_GFP_1st", "Mock_GFP_2nd")
d2018$Mock_GFP_1st <- NULL
d2018$Mock_GFP_2nd <- NULL

## To compare a WT GFP network, we apply differentially expressed proteins using edgeR
colnames(d2018)
counts = d2018[,grepl('1st|2nd', colnames(d2018))]
rownames(counts) <- d2018$uniprot_id
# d0 <- calcNormFactors(DGEList(counts))
d <- DGEList(counts)
# cutoff <- 0
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
dim(d) 
snames <- colnames(counts) # Sample names
group <- gsub('_1st|_2nd', '', snames)
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~ 0 + group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupWT_GFP - groupWT_Con, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))
top.table$uniprot_id <- rownames(top.table)

m <- merge(top.table, d2018, by='uniprot_id', all.x = T)

## Create background interacting proteins for CD99
#pm = m[m$logFC > 0 & m$adj.P.Val <= 0.05,]
#pm = m[m$logFC > 0,]
pm = m

pm$nm_WT_GFP_1 <- (pm$WT_GFP_1st + 1) / (pm$WT_Con_1st + 1)
pm$nm_Cyto_GFP_1 <- (pm$Cyto_GFP_1st + 1)/ (pm$WT_Con_1st + 1)
pm$nm_TM_GFP_1 <- (pm$TM_GFP_1st + 1) / (pm$WT_Con_1st + 1)

pm$nm_WT_GFP_2 <- (pm$WT_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)
pm$nm_Cyto_GFP_2 <- (pm$Cyto_GFP_2nd + 1)/ (pm$WT_Con_2nd + 1)
pm$nm_TM_GFP_2 <- (pm$TM_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)

# pm$isCytoskeleton <- ifelse(toupper(pm$gene_name) %in% toupper(c(p_myosin, p_actin, p_tubulin, p_cyto)), 
#                             'Y', 'N')
# pm$isCytoplasm <- ifelse( toupper(pm$gene_name) %in% toupper(p_ctp), 'Y', 'N')
# pm$isCortex <- ifelse( toupper(pm$gene_name) %in% toupper(p_ccrt), 'Y', 'N')
# pm$isTCR <- ifelse( toupper(pm$gene_name) %in% toupper(p_tcr), 'Y', 'N')

## Fold Change compared to WT
pm$fc_Cyto_GFP_1 <- (pm$nm_Cyto_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_Cyto_GFP_2 <- (pm$nm_Cyto_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
pm$fc_TM_GFP_1 <- (pm$nm_TM_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_TM_GFP_2 <- (pm$nm_TM_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
#write.table(pm, 'Tables/table.edgeR_WT_GFP.CD99_20180514.txt', sep='\t', quote = F, row.names = F, col.names=T)

## ------------------------------ vv
## To create CT and TM proteins, we apply differentially expressed proteins using edgeR
contr <- makeContrasts(groupCyto_GFP - groupWT_GFP, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_cyto <- topTable(tmp, sort.by = "P", n = Inf)
top.table_cyto$uniprot_id <- rownames(top.table_cyto)  
m_cyto <- merge(top.table_cyto, d2018, by='uniprot_id', all.x = T)

contr <- makeContrasts(groupTM_GFP - groupWT_GFP, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_tm <- topTable(tmp, sort.by = "P", n = Inf)
top.table_tm$uniprot_id <- rownames(top.table_tm)  
m_tm <- merge(top.table_tm, d2018, by='uniprot_id', all.x = T)

m <- merge(setNames(top.table_cyto[,c('logFC', 'adj.P.Val', 'uniprot_id')], 
                    c('logFC_cyto', 'padj_cyto', 'uniprot_id')), 
           setNames(top.table_tm[,c('logFC', 'adj.P.Val', 'uniprot_id')], 
                    c('logFC_tm', 'padj_tm', 'uniprot_id')), by='uniprot_id') 
m1 = merge(m[m$uniprot_id %in% pm$uniprot_id,], d2018, by='uniprot_id', all.x = T)

m1$FC_cyto <- log2(( ((m1$Cyto_GFP_1st + 1)/ (m1$WT_GFP_1st + 1)) + ((m1$Cyto_GFP_2nd + 1)/ (m1$WT_GFP_2nd + 1)))/2)
m1$FC_tm <- log2(( ((m1$TM_GFP_1st + 1)/ (m1$WT_GFP_1st + 1)) + ((m1$TM_GFP_2nd + 1)/ (m1$WT_GFP_2nd + 1)))/2)
write.table(m1, 'Tables/table.edgeR_FC_Mutant.CD99_20180514.txt', sep='\t', quote = F, row.names = F, col.names=T)


m1$Sig <- ifelse(m1$padj_cyto <= 0.05  & m1$padj_tm <= 0.05, 'Both', 
                 ifelse(m1$padj_cyto <= 0.05, 'Cyto', ifelse(m1$padj_tm <= 0.05, 'TM', 'None')))

m1$Sig2 <- 'None'
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm > 0, 'Both Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm < 0, 'Both Down', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm < 0, 'Cyto Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm > 0, 'Cyto Down', m1$Sig2)

m1$Sig2 <- ifelse(m1$Sig %in% c('Cyto') & m1$Sig2=='None', ifelse(m1$logFC_cyto > 0, 'Cyto Up', 'Cyto Down'), m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('TM') & m1$Sig2=='None', ifelse(m1$logFC_tm > 0, 'TM Up', 'TM Down'), m1$Sig2)

m1$isSig_cyto <- ifelse(m1$padj_cyto <= 0.05, ifelse( m1$logFC_cyto >0 , 'Up', 'Down'), 'None')
m1$isSig_tm <- ifelse(m1$padj_tm <= 0.05, ifelse( m1$logFC_tm >0 , 'Up', 'Down'), 'None')

m1$Group <- ifelse(m1$gene_name %in% p_actomyo, 'Acto', 
                   ifelse(m1$gene_name %in% p_tubulin, 'Tubule','other'))


## Plot: Compare the fold change between CT and TM mutants
mm <- gather(m1[m1$Group %in% c('Acto', 'Tubule'), c('FC_cyto', 'FC_tm', 'gene_name')], 
             group, val, FC_cyto:FC_tm, -gene_name, factor_key=TRUE)
mm = mm[order(mm$val, decreasing = T),]
mm$gene_name = factor(mm$gene_name, levels=unique(mm$gene_name))
mm$group2 <- ifelse(grepl('tm', mm$group), 'TM', 'Cyto')

####

mm1 <- mm[mm$group2 == "Cyto",]
mm1 = mm1[order(mm1$val),]
mm1$rank = 1:nrow(mm1)
mm2 <- mm[mm$group2 == "TM",]
mm2 = mm2[order(mm2$val),]
mm2$rank = 1:nrow(mm2)
mm = rbind.data.frame(mm1, mm2)

gs1 = c('Iqgap1') 
gs2 = c('Myh9')
p4 <- ggplot(mm, aes(x=rank, y=val, group = group2)) + 
  #  geom_violin(aes(fill = group2), alpha=0.3, scale = "width") +
  #  geom_point(color = 'white', fill = 'black', alpha = 3) + 
  geom_point(color = 'black', 
             fill = 'white', size = 0.75,stroke = 1) +
  facet_wrap(~group2, strip.position = "bottom") +
  # scale_x_continuous(breaks = 1:2, labels = unique(counts_df$line)) +
  #  geom_line(data = mm[mm$gene_name %in% gs1,], aes(group=gene_name), color='red', alpha=0.3) + 
  #  geom_line(data = mm[mm$gene_name %in% gs2,], aes(group=gene_name), color='blue', alpha=0.3) + 
  geom_point(data = mm[mm$gene_name %in% gs1,], color='red', size = 3) + 
  geom_point(data = mm[mm$gene_name %in% gs2,], color='blue', size = 3) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(data = mm[mm$gene_name %in% gs1,], aes(label=gene_name), 
                  color='red',
                  size = 4, nudge_x = -15) +
  geom_text_repel(data = mm[mm$gene_name %in% gs2,], aes(label=gene_name), 
                  color='blue',
                  size = 4, nudge_x=20) +
  theme_linedraw(base_size = 10) + 
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 10, color = "black"),
        strip.background = element_rect(colour="black",fill="white"),
        panel.grid = element_line(colour = "gray", size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = '', y='Log fold change compared to WT')
p4




















