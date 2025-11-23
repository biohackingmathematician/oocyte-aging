

counts <- read.table("all_samples.featureCounts", sep="\t", header=T, row.names = 1)
counts <- counts[,-c(1:5)]
colnames(counts) <- sub("\\.","",sub("^X","",sub("_.*","",colnames(counts))))

counts_new <- read.table("13_18_48_55_85.featureCounts", sep="\t", header=T, row.names = 1)
counts_new <- counts_new[,-c(1:5)]
colnames(counts_new) <- sub("\\.","",sub("^X","",sub("_.*","",colnames(counts_new))))

counts <- merge(counts,counts_new, by="row.names")
rownames(counts) <- counts$Row.names
counts <- counts[,-1]

metadata <- data.frame(sample=colnames(counts),group="MI")
rownames(metadata)<-metadata$sample
metadata$group[grep("MII",colnames(counts))]<-"MII"
metadata$group[grep("VG",colnames(counts))]<-"VG"

metadata$group <- factor(metadata$group, levels=c("VG","MI","MII"))
metadata


# Remove 9MII since it was a failed sample
counts <- counts[,! colnames(counts) %in% c("9MII")]
metadata <- metadata[! metadata$sample %in% c("9MII"),]


library(edgeR)
y <- DGEList(counts = counts, group=metadata$group)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)
y$samples

cpms <- cpm(y, normalized.lib.sizes = TRUE, log=T)

library(pheatmap)
library(RColorBrewer)
col = colorRampPalette(c("red", "yellow"))(100)
pheatmap(as.matrix(dist(t(cpms))), color = col, 
         filename = "fig1c_heatmap.pdf",
         annotation_col = metadata[,-1, drop=FALSE],
         annotation_colors = list(group=c(VG="#F8766D",MI="#39B600",MII="#00B0F6")),
         width = 8,
         height = 7)


library(PCAtools)
p <- PCAtools::pca(cpms, removeVar = 0.1, metadata = metadata[colnames(cpms),-1, drop=FALSE])
pdf("fig1b_pca.pdf", width=8, height=6)
biplot(p, 
       showLoadings = F,
       colby = 'group',
       hline = 0, vline = 0,
       legendPosition = 'right',
       xlim=c(-200,400))
dev.off()

tr_gene_names <- read.table("e106_biomart_HomoSapiens_Transcript_Gene_Names.txt", sep="\t", header=T, quote = "")
cpms_named <- merge(cpms, unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
cpms_named <- cpms_named[!duplicated(cpms_named$Gene.name),]
rownames(cpms_named) <- cpms_named$Gene.name 
cpms_named <- cpms_named[,colnames(cpms)]

mm <- model.matrix( ~ 0 + metadata$group)
y <- estimateDisp(y, design = mm, robust = T)

fit <- glmQLFit(y, mm, robust = T)

# MII vs MI
qlf <- glmQLFTest(fit, contrast=c(0,-1,1))
mii_vs_mi_featureCounts <- topTags(qlf, sort.by = "P", n = Inf)

mii_vs_mi_featureCounts_res <- merge(mii_vs_mi_featureCounts$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
mii_vs_mi_featureCounts_res <- merge(mii_vs_mi_featureCounts_res,cpms_named[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MII"])], by.y="row.names", by.x="Gene.name")
mii_vs_mi_featureCounts_res <- mii_vs_mi_featureCounts_res[order(mii_vs_mi_featureCounts_res$PValue, decreasing = F),]

write.table(mii_vs_mi_featureCounts_res,"Mii_vs_Mi_FeatureCounts_diff_analysis.tab", sep="\t", row.names = F)

# MI vs VG
qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
vg_vs_mi_featureCounts <- topTags(qlf, sort.by = "P", n = Inf)

vg_vs_mi_featureCounts_res <- merge(vg_vs_mi_featureCounts$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
vg_vs_mi_featureCounts_res <- merge(vg_vs_mi_featureCounts_res,cpms_named[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MI"])], by.y="row.names", by.x="Gene.name")
vg_vs_mi_featureCounts_res <- vg_vs_mi_featureCounts_res[order(vg_vs_mi_featureCounts_res$PValue, decreasing = F),]

write.table(vg_vs_mi_featureCounts_res,"Vg_vs_Mi_featureCounts_diff_analysis.tab", sep="\t", row.names = F)


# MII vs VG
qlf <- glmQLFTest(fit, contrast=c(-1,0,1))
vg_vs_mii_featureCounts <- topTags(qlf, sort.by = "P", n = Inf)

vg_vs_mii_featureCounts_res <- merge(vg_vs_mii_featureCounts$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
vg_vs_mii_featureCounts_res <- merge(vg_vs_mii_featureCounts_res,cpms_named[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MII"])], by.y="row.names", by.x="Gene.name")
vg_vs_mii_featureCounts_res <- vg_vs_mii_featureCounts_res[order(vg_vs_mii_featureCounts_res$PValue, decreasing = F),]

write.table(vg_vs_mii_featureCounts_res,"Vg_vs_Mii_featureCounts_diff_analysis.tab", sep="\t", row.names = F)



library(EnhancedVolcano)


pdf("Fig3_a_b_c_volcano_plots.pdf", width=10,height=6)


EnhancedVolcano(mii_vs_mi_featureCounts_res,
                subtitle = "Mii vs Mi",
                lab = mii_vs_mi_featureCounts_res$Gene.name,
                pCutoff = 0.05,
                ylim =  c(0,2),
                x = 'logFC',
                y = 'FDR'
                )

EnhancedVolcano(vg_vs_mi_featureCounts_res,
                subtitle = "MI vs VG",
                lab = vg_vs_mi_featureCounts_res$Gene.name,
                pCutoff = 0.05,
                ylim =  c(0,2),
                x = 'logFC',
                y = 'FDR')


EnhancedVolcano(vg_vs_mii_featureCounts_res,
                subtitle = "MII vs VG",
                lab = vg_vs_mii_featureCounts_res$Gene.name,
                pCutoff = 0.05,
                ylim =  c(0,3.5),
                x = 'logFC',
                y = 'FDR')


dev.off()


# violin_plot

genes_to_plot <- c("DNMT1", "DPPA3", "MED30", "ZFAND2A", "TMEM216", "PANX1", 
                   "CDC20", "PATL2", "TRIP13", "WEE2", "ZNF738", "ZAR1", 
                   "CDK1", "CCNB1")


metadata_extended <- metadata
metadata_extended <- cbind(metadata_extended,
                           t(cpms_named[genes_to_plot,metadata$sample]))

metadata_extended$group <- factor(metadata_extended$group, levels=c("VG","MI", "MII"))

pdf("fig4_violin_plots.pdf", width=8, height=6)
for(gene in genes_to_plot){
p <- ggplot(metadata_extended, aes(x=group, y=get(gene), fill=group)) +
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center',
                               position=position_dodge(1)) +
  labs(title=paste(gene, "expression by group"),x="Group", y = "Normalized CPM")+
  theme_classic()+ theme(text = element_text(size = 20))
print(p)
}
dev.off()


library(readxl)
library(pheatmap)
library(RColorBrewer)

col = colorRampPalette(c("steelblue", "coral"))(100)

sheet_names <- c("HDAC","SIRT","KDM","Histone acetiltransferase","SETD",
                 "KMT","Imprinted genes","DNMT, TET, UHRF1 and DPPA3")
for (sheet in sheet_names){
  sheet_content <- read_xlsx("log2_CPM_EPIGENETICS.xlsx", sheet = sheet, skip = 1)
  pheatmap(sheet_content[,-1], color = col, 
           gaps_col = c(6,12,19),
           show_rownames = T, labels_row = sheet_content$GENE,
           filename = paste("fig5_heatmap_",gsub(" ","_",sheet) ,".pdf", sep=""),
           #annotation_col = metadata[sub("-","",colnames(sheet_content)),-1, drop=FALSE],
           cluster_rows = F, cluster_cols = F, width = 8, height = 7)
  
}


# voom
library(edgeR)
group <- metadata$group
y <- DGEList(counts = counts, group=group)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

mm <- model.matrix( ~ 0 + group)
y <- voom(y, mm, plot = T)

fit <- lmFit(y, mm)

#plotMDS(y, col = group)

contr <- makeContrasts(groupMII - groupMI, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_featurecounts_voom_mii_vs_mi <- topTable(tmp, sort.by = "P", n = Inf)

contr <- makeContrasts(groupMI - groupVG, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_featurecounts_voom_vg_vs_mi <- topTable(tmp, sort.by = "P", n = Inf)

contr <- makeContrasts(groupMII - groupVG, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_featurecounts_voom_vg_vs_mii <- topTable(tmp, sort.by = "P", n = Inf)

write.table(top.table_featurecounts_voom_mii_vs_mi,"Mi_vs_Mii_featureCounts_diff_analysis_voom.tab", sep="\t", row.names = F)
write.table(top.table_featurecounts_voom_vg_vs_mi,"Vg_vs_Mi_featureCounts_diff_analysis_voom.tab", sep="\t", row.names = F)
write.table(top.table_featurecounts_voom_vg_vs_mii,"Vg_vs_Mii_featureCounts_diff_analysis_voom.tab", sep="\t", row.names = F)



# DESEq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design= ~ 0 + group)
keep <- rowSums(counts(dds) >= 10) > 3
dds <- dds[keep,]
dds <- DESeq(dds)
res_featureCounts_vg_mi <- results(dds, contrast=c("group", "VG", "MI"))
res_featureCounts_vg_mii <- results(dds, contrast=c("group", "VG", "MII"))
res_featureCounts_mi_mii <- results(dds, contrast=c("group", "MI", "MII"))

write.table(res_featureCounts_mi_mii,"Mi_vs_Mii_featureCounts_diff_analysis_deseq2.tab", sep="\t", row.names = F)
write.table(res_featureCounts_vg_mi,"Vg_vs_Mi_featureCounts_diff_analysis_deseq2.tab", sep="\t", row.names = F)
write.table(res_featureCounts_vg_mii,"Vg_vs_Mii_featureCounts_diff_analysis_deseq2.tab", sep="\t", row.names = F)


#res <- res[order(res$padj, decreasing = F),]
#
#full_res <- merge(results_table,res, by="Row.names")
#full_res <- full_res[order(full_res$FDR, decreasing = F),]

library(VennDiagram)
venn.diagram(
  x = list(gfp_vs_yfp, gfp_vs_dn, yfp_vs_dn),
  category.names = c("GFP VS YFP" , "GFP VS DN " , "YFP VS DN"),
  fill = c("blue", "green", "red"),
  alpha = 0.50,
  filename = 'venn_diagramm_original.tif',
  output=TRUE
)


#### Kallisto analysis

library(tximport)

dirs <- grep("_kallisto$",list.dirs("kallisto/",recursive = F), value=T)
files <- paste(dirs, "/abundance.h5", sep="")
names(files) <- sub("_.*","",sub(".*\\/\\/","",sub("_S.*","",files)))

txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene = tr_gene_names[,c(1,2)])
#head(txi.kallisto$counts)
#head(txi.kallisto$abundance)

kallisto_counts <- txi.kallisto$counts

metadata <- data.frame(sample=colnames(txi.kallisto$abundance),group="MI")
rownames(metadata)<-metadata$sample
metadata[c("11VG","37VG","38VG","50VG","54VG", "55-VG"),"group"] <- "VG"
metadata[c("12MII","13-MII","16MII","17-MII","53MII","7-MII","85-MII","9-MII"),"group"] <- "MII"

# Remove 9-MII since it was a failed sample
kallisto_counts <- kallisto_counts[,! colnames(kallisto_counts) %in% c("9-MII")]
metadata <- metadata[! metadata$sample %in% c("9-MII"),]


library(edgeR)
y <- DGEList(counts = kallisto_counts, group=metadata$group)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)
y$samples

cpms_kallisto <- cpm(y, normalized.lib.sizes = TRUE, log=T)

library(pheatmap)
library(RColorBrewer)
col = colorRampPalette(c("red", "yellow"))(100)
pheatmap(as.matrix(dist(t(cpms_kallisto))), color = col, 
         filename = "fig1c_heatmap_kallisto.pdf",
         annotation_col = metadata[,-1, drop=FALSE],
         annotation_colors = list(group=c(VG="#F8766D",MI="#39B600",MII="#00B0F6")),
         width = 8,
         height = 7)


library(PCAtools)
p <- PCAtools::pca(cpms_kallisto, removeVar = 0.1, metadata = metadata[colnames(cpms),-1, drop=FALSE])
pdf("fig1b_pca_kallisto.pdf", width=8, height=6)
biplot(p, 
       showLoadings = F,
       colby = 'group',
       hline = 0, vline = 0,
       legendPosition = 'right',
       xlim=c(-200,400))
dev.off()

tr_gene_names <- read.table("e106_biomart_HomoSapiens_Transcript_Gene_Names.txt", sep="\t", header=T, quote = "")
cpms_named_kallisto <- merge(cpms_kallisto, unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
cpms_named_kallisto <- cpms_named_kallisto[!duplicated(cpms_named_kallisto$Gene.name),]
rownames(cpms_named_kallisto) <- cpms_named_kallisto$Gene.name 
cpms_named_kallisto <- cpms_named_kallisto[,colnames(cpms_kallisto)]

mm <- model.matrix( ~ 0 + metadata$group)
y <- estimateDisp(y, design = mm, robust = T)

fit <- glmQLFit(y, mm, robust = T)

# MII vs MI
qlf <- glmQLFTest(fit, contrast=c(-1,1,0))
mii_vs_mi_kallisto <- topTags(qlf, sort.by = "P", n = Inf)

mii_vs_mi_kallisto_res <- merge(mii_vs_mi_kallisto$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
mii_vs_mi_kallisto_res <- merge(mii_vs_mi_kallisto_res,cpms_named_kallisto[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MII"])], by.y="row.names", by.x="Gene.name")
mii_vs_mi_kallisto_res <- mii_vs_mi_kallisto_res[order(mii_vs_mi_kallisto_res$PValue, decreasing = F),]

write.table(mii_vs_mi_kallisto_res,"Mii_vs_Mi_kallisto_diff_analysis.tab", sep="\t", row.names = F)

# MI vs VG
qlf <- glmQLFTest(fit, contrast=c(-1,0,1))
vg_vs_mi_kallisto <- topTags(qlf, sort.by = "P", n = Inf)

vg_vs_mi_kallisto_res <- merge(vg_vs_mi_kallisto$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
vg_vs_mi_kallisto_res <- merge(vg_vs_mi_kallisto_res,cpms_named_kallisto[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MI"])], by.y="row.names", by.x="Gene.name")
vg_vs_mi_kallisto_res <- vg_vs_mi_kallisto_res[order(vg_vs_mi_kallisto_res$PValue, decreasing = F),]

write.table(vg_vs_mi_kallisto_res,"Vg_vs_Mi_kallisto_diff_analysis.tab", sep="\t", row.names = F)


# MII vs VG
qlf <- glmQLFTest(fit, contrast=c(0,-1,1))
vg_vs_mii_kallisto <- topTags(qlf, sort.by = "P", n = Inf)

vg_vs_mii_kallisto_res <- merge(vg_vs_mii_kallisto$table,unique(tr_gene_names[,c("Gene.stable.ID","Gene.name")]), by.x="row.names", by.y="Gene.stable.ID")
vg_vs_mii_kallisto_res <- merge(vg_vs_mii_kallisto_res,cpms_named_kallisto[c(metadata$sample[metadata$group=="VG"],metadata$sample[metadata$group=="MII"])], by.y="row.names", by.x="Gene.name")
vg_vs_mii_kallisto_res <- vg_vs_mii_kallisto_res[order(vg_vs_mii_kallisto_res$PValue, decreasing = F),]

write.table(vg_vs_mii_kallisto_res,"Vg_vs_Mii_kallisto_diff_analysis.tab", sep="\t", row.names = F)

sum(mii_vs_mi_kallisto$table$FDR<0.05)


# voom
library(edgeR)
group <- metadata$group
y <- DGEList(counts = kallisto_counts, group=group)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

mm <- model.matrix( ~ 0 + group)
y <- voom(y, mm, plot = T)

fit <- lmFit(y, mm)

#plotMDS(y, col = group)

contr <- makeContrasts(groupMII - groupMI, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_kallisto_voom_mii_vs_mi <- topTable(tmp, sort.by = "P", n = Inf)


contr <- makeContrasts(groupVG - groupMI, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_kallisto_voom_vg_vs_mi <- topTable(tmp, sort.by = "P", n = Inf)

contr <- makeContrasts(groupVG - groupMII, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_kallisto_voom_vg_vs_mii <- topTable(tmp, sort.by = "P", n = Inf)

write.table(top.table_kallisto_voom_mii_vs_mi,"Mi_vs_Mii_kallisto_diff_analysis_voom.tab", sep="\t", row.names = F)
write.table(top.table_kallisto_voom_vg_vs_mi,"Vg_vs_Mi_kallisto_diff_analysis_voom.tab", sep="\t", row.names = F)
write.table(top.table_kallisto_voom_vg_vs_mii,"Vg_vs_Mii_kallisto_diff_analysis_voom.tab", sep="\t", row.names = F)




# DESEq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = round(kallisto_counts),
                              colData = metadata,
                              design= ~ 0 + group)
keep <- rowSums(counts(dds) >= 10) > 3
dds <- dds[keep,]
dds <- DESeq(dds)
res_kallisto_vg_mi <- results(dds, contrast=c("group", "VG", "MI"))
res_kallisto_vg_mii <- results(dds, contrast=c("group", "VG", "MII"))
res_kallisto_mi_mii <- results(dds, contrast=c("group", "MI", "MII"))

write.table(res_kallisto_mi_mii,"Mi_vs_Mii_kallisto_diff_analysis_deseq2.tab", sep="\t", row.names = F)
write.table(res_kallisto_vg_mi,"Vg_vs_Mi_kallisto_diff_analysis_deseq2.tab", sep="\t", row.names = F)
write.table(res_kallisto_vg_mii,"Vg_vs_Mii_kallisto_diff_analysis_deseq2.tab", sep="\t", row.names = F)



library(VennDiagram)
venn.diagram(
  x = list(gfp_vs_yfp, gfp_vs_dn, yfp_vs_dn),
  category.names = c("GFP VS YFP" , "GFP VS DN " , "YFP VS DN"),
  fill = c("blue", "green", "red"),
  alpha = 0.50,
  filename = 'venn_diagramm_original.tif',
  output=TRUE
)

edgeR <- vg_vs_mi_featureCounts_res$Row.names[vg_vs_mi_featureCounts_res$FDR<0.05]
voom <- rownames(top.table_featurecounts_voom_vg_vs_mi)[top.table_featurecounts_voom_vg_vs_mi$adj.P.Val<0.05]
deseq <- rownames(res_featureCounts_vg_mi)[res_featureCounts_vg_mi$padj<0.05]
deseq <- deseq[!is.na(deseq)]

category_names <- c("FCounts EdgeR" , "FCounts Voom" , "FCounts DEseq2")
filename <- 'venn_diagramm_FC_vg_vs_mi.tif'

library(VennDiagram)
venn.diagram(
  x = list(edgeR, 
           voom, 
           deseq),
  category.names = category_names,
  fill = c("blue", "green", "red"),
  alpha = 0.50,
  filename = filename,
  output=TRUE
)

edgeR <- vg_vs_mii_featureCounts_res$Row.names[vg_vs_mii_featureCounts_res$FDR<0.05]
voom <- rownames(top.table_featurecounts_voom_vg_vs_mii)[top.table_featurecounts_voom_vg_vs_mii$adj.P.Val<0.05]
deseq <- rownames(res_featureCounts_vg_mii)[res_featureCounts_vg_mii$padj<0.05]
deseq <- deseq[!is.na(deseq)]

category_names <- c("FCounts EdgeR" , "FCounts Voom" , "FCounts DEseq2")
filename <- 'venn_diagramm_FC_vg_vs_mii.tif'

library(VennDiagram)
venn.diagram(
  x = list(edgeR, 
           voom, 
           deseq),
  category.names = category_names,
  fill = c("blue", "green", "red"),
  alpha = 0.50,
  filename = 'venn_diagramm_FC_vg_vs_mi.tif',
  output=TRUE
)

