options(stringsAsFactors = F)
library(readxl)
library(cowplot)
library(dplyr)
library(tidyr) 
library(readr)
library(edgeR) # https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
library(ggrepel)
setwd('~/Dropbox/NGR_SNU_2019/wgcna_t_cell_201912/')

## Load data (2018) - It has TM and Cyto mutants
d2018 = read_excel('Data/20180912_CD99_Total ID protein.xlsx')
d2018 = as.data.frame(d2018)
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
pm = m[m$logFC > 0 & m$adj.P.Val <= 0.05,]

pm$nm_WT_GFP_1 <- (pm$WT_GFP_1st + 1) / (pm$WT_Con_1st + 1)
pm$nm_Cyto_GFP_1 <- (pm$Cyto_GFP_1st + 1)/ (pm$WT_Con_1st + 1)
pm$nm_TM_GFP_1 <- (pm$TM_GFP_1st + 1) / (pm$WT_Con_1st + 1)

pm$nm_WT_GFP_2 <- (pm$WT_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)
pm$nm_Cyto_GFP_2 <- (pm$Cyto_GFP_2nd + 1)/ (pm$WT_Con_2nd + 1)
pm$nm_TM_GFP_2 <- (pm$TM_GFP_2nd + 1) / (pm$WT_Con_2nd + 1)

pm$isCytoskeleton <- ifelse(toupper(pm$gene_name) %in% toupper(c(p_myosin, p_actin, p_tubulin, p_cyto)), 
                            'Y', 'N')
pm$isCytoplasm <- ifelse( toupper(pm$gene_name) %in% toupper(p_ctp), 'Y', 'N')
pm$isCortex <- ifelse( toupper(pm$gene_name) %in% toupper(p_ccrt), 'Y', 'N')
pm$isTCR <- ifelse( toupper(pm$gene_name) %in% toupper(p_tcr), 'Y', 'N')

## Fold Change compared to WT
pm$fc_Cyto_GFP_1 <- (pm$nm_Cyto_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_Cyto_GFP_2 <- (pm$nm_Cyto_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
pm$fc_TM_GFP_1 <- (pm$nm_TM_GFP_1 + 1)/ (pm$nm_WT_GFP_1 + 1)
pm$fc_TM_GFP_2 <- (pm$nm_TM_GFP_2 + 1)/ (pm$nm_WT_GFP_2 + 1)
#write.table(pm, 'Tables/table.edgeR_WT_GFP.CD99_20190815.txt', sep='\t', quote = F, row.names = F, col.names=T)

## Create a gene network
pp = pp0[pp0$combined_score >= 700, 1:2]
pp1 = pp[toupper(pp$NodeA) %in% toupper(pm$gene_name) & toupper(pp$NodeB) %in% toupper(pm$gene_name), ]
pp1 = pp1[pp1$NodeA!=pp1$NodeB,]

pp1$AttrA <- ifelse(toupper(pp1$NodeA) %in% p_cyto2, 'Cytoskeleton', 
                    ifelse(toupper(pp1$NodeA) %in% p_cellcycle, 'Cell Cycle', 
                           ifelse(toupper(pp1$NodeA) %in% toupper(p_ccrt), 'Cortex', 
                                  ifelse(toupper(pp1$NodeA) %in% toupper(p_tcr), 'TCR', 
                                         ifelse(toupper(pp1$NodeA) %in% toupper(p_atp), 'ATP', 
                                                ifelse(toupper(pp1$NodeA) %in% toupper(p_ribo), 'Ribo', 'Other'))))))

gs = pm$gene_name
gs = gs[!is.na(gs)]
writeLines(sort(gs), 'Tables/list.edgeR_WT_GFP.CD99_20190815.md')
table(pp1$AttrA)
write.table(pp1, 'Tables/table.stringDB_700.CD99_20190815.txt', sep='\t', quote = F, row.names = F, col.names = T)




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

m1$Sig <- ifelse(m1$padj_cyto <= 0.05  & m1$padj_tm <= 0.05, 'Both', 
                 ifelse(m1$padj_cyto <= 0.05, 'Cyto', ifelse(m1$padj_tm <= 0.05, 'TM', 'None')))

m1$Sig2 <- 'None'
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm > 0, 'Both Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm < 0, 'Both Down', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto > 0 & m1$FC_tm < 0, 'Cyto Up', m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('Both') & m1$FC_cyto < 0 & m1$FC_tm > 0, 'Cyto Down', m1$Sig2)

m1$Sig2 <- ifelse(m1$Sig %in% c('Cyto') & m1$Sig2=='None', ifelse(m1$logFC_cyto > 0, 'Cyto Up', 'Cyto Down'), m1$Sig2)
m1$Sig2 <- ifelse(m1$Sig %in% c('TM') & m1$Sig2=='None', ifelse(m1$logFC_tm > 0, 'TM Up', 'TM Down'), m1$Sig2)
m1$Group <- ifelse( toupper(m1$gene_name) %in% p_cyto2, 'cyto', 
                    ifelse( toupper(m1$gene_name) %in% p_cyto2, 'cyto','other'))

m1$Group <- ifelse(toupper(m1$gene_name) %in% p_myosin, 'Myosin',
                   ifelse(toupper(m1$gene_name) %in% p_actin, 'Actin Cytoskeleton', 
                          ifelse(toupper(m1$gene_name) %in% toupper(p_tubulin), 'Tubulin',
                                 ifelse(toupper(m1$gene_name) %in% toupper(p_traf), 'Molecular trafficking',
                                        ifelse(toupper(m1$gene_name) %in% toupper(p_tcr), 'TCR', 
                                               ifelse(toupper(m1$gene_name) %in% toupper(p_atp), 'ATP', 
                                                      ifelse(toupper(m1$gene_name) %in% toupper(p_ribo), 
                                                             'Ribo', 'Other')))))))
m1$Group <- factor(m1$Group, levels=c('Myosin','Actin Cytoskeleton', 'Tubulin', 'Molecular trafficking', 'TCR', 'ATP', 'Ribo', 'Other'))

m1$isSig_cyto <- ifelse(m1$padj_cyto <= 0.05, ifelse( m1$logFC_cyto >0 , 'Up', 'Down'), 'None')
m1$isSig_tm <- ifelse(m1$padj_tm <= 0.05, ifelse( m1$logFC_tm >0 , 'Up', 'Down'), 'None')
col_dep = c('Up'='firebrick', 'Down'='dodgerblue', 'None'='grey')

## Plot: Compare the fold change between CT and TM mutants
mm <- gather(m1[m1$Group %in% c('Myosin', 'Cytoskeleton'), c('FC_cyto', 'FC_tm', 'gene_name')], 
             group, val, FC_cyto:FC_tm, -gene_name, factor_key=TRUE)
mm = mm[order(mm$val, decreasing = T),]
mm$gene_name = factor(mm$gene_name, levels=unique(mm$gene_name))
mm$group2 <- ifelse(grepl('tm', mm$group), 'TM', 'Cyto')

gs1 = c('Myh9', 'Myo1g', 'Myl6') 
gs2 = c('Tubb5', 'Tuba1b', 'Tubb4b', 'Tubb4a')

p4 <- ggplot(mm, aes(group2, val )) +
  geom_violin(aes(fill = group2), alpha=0.3, scale = "width") +
  geom_point() + 
  geom_point(data = mm[mm$gene_name %in% gs1,], color='red') + 
  geom_point(data = mm[mm$gene_name %in% gs2,], color='blue') + 
  scale_fill_manual(values =  c('#D55E00', '#56B4E9')) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_text_repel(data = mm[mm$gene_name %in% gs1,], aes(label=gene_name), 
                  color='red',
                  size = 2, nudge_x=0.2) +
  geom_text_repel(data = mm[mm$gene_name %in% gs2,], aes(label=gene_name), 
                  color='blue',
                  size = 2, nudge_x=0.2) +
  theme_minimal(base_size = 10) + 
  theme(legend.position="none") +
  labs(x = '', y='Log fold change compared to WT')
p4