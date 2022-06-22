### This analysis is performed based on gene counts table

require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)
require(preprocessCore)
require(edgeR)
require(writexl)


samples <- data.frame(c('sample_BD1', 'sample_BD2', 'sample_BD3', 'sample_BD4', 'sample_BD5', 'sample_BW1', 'sample_BW2', 'sample_BW3', 'sample_BW4',
             'sample_BW5', 'sample_MD1', 'sample_MD2', 'sample_MD3', 'sample_MD4', 'sample_MD5', 'sample_MW1', 'sample_MW2', 'sample_MW3',
             'sample_MW4', 'sample_MW5'), c(rep('wt',10), rep('mdx', 10)),
             c(rep('10days', 5), rep('10weeks', 5), rep('10days', 5), rep('10weeks', 5)),
             c(rep('WT10D', 5), rep('WT10W', 5), rep('MDX10D', 5), rep('MDX10W', 5)))
colnames(samples) <- c('id', 'genotype', 'age', 'group')
rownames(samples) <- samples$id

##### construct a results table #####
counts <- read_tsv('/projects/portsmouth-darek-cerebellum/results/seq-results/cerebellum_gene-level-counts-annotated.tsv')

write_xlsx(counts,'/projects/portsmouth-darek-cerebellum/results/seq-results/gene-counts-table.xlsx')
counts <- data.frame(counts)

results <- data.frame(counts$`Transcript stable ID`, counts$`Gene name`, counts$`Gene Synonym`, counts$`Gene stable ID`)
colnames(results) <- c('transcript_stable_id', 'gene_name', 'gene_synonyms', 'gene_stable_id')

counts[,11:30] %>% data.matrix %>% normalize.quantiles() -> normalised
colnames(normalised) <- colnames(counts[,11:30])
rownames(normalised) <- counts$`Gene stable ID`
log_tpms <- log2(normalised + 1)

rownames(log_tpms) <- counts$`Gene stable ID`
write.csv(log_tpms, '/projects/portsmouth-darek-cerebellum/results/seq-results/log_counts.tsv')

log_tpms <- read_csv('/projects/portsmouth-darek-cerebellum/results/seq-results/log_counts.tsv') %>%
  column_to_rownames(var="...1")

##### EdgeR analysis #####

design <- model.matrix(~0+group, data=samples)
colnames(design) <- levels(as.factor(samples$group))

counts_edger <- DGEList(counts=counts[,samples$id], group=samples$group)
keep <- filterByExpr(counts_edger)
counts_edger <- counts_edger[keep, , keep.lib.sizes=FALSE]
counts_edger <- calcNormFactors(counts_edger)

counts_edger <- estimateDisp(counts_edger, design)
fit <- glmQLFit(counts_edger, design)

my.contrasts <- makeContrasts(
  #mdx10d vs mdx10w but controlled for wt effects
  mdx_age_corrected_wt=(MDX10D-MDX10W-WT10D+WT10W),
  #mdx10d vs mdx10w not controlled
  mdx_age=MDX10W-MDX10D,
  #genotype effect 10 days
  mdx_vs_wt_10d=MDX10D-WT10D,
  #genotype effect 10 weeks
  mdx_vs_wt_10w=MDX10W-WT10W,
  #age effect overall
  age=WT10W+MDX10W-WT10D-MDX10D,
  #wt age
  wt_age=WT10W-WT10D,
  #overall genotype effect
  genotype=(WT10D+WT10W)-(MDX10D+MDX10W),
  levels=design)

mdx_age_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_age"])
mdx_age_corrected_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_age_corrected_wt"])
mdx_vs_wt_10d_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_vs_wt_10d"])
mdx_vs_wt_10w_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_vs_wt_10w"])
wt_age_fit <- glmQLFTest(fit, contrast=my.contrasts[,"wt_age"])
age_fit <- glmQLFTest(fit, contrast=my.contrasts[,"age"])
genotype_fit <- glmQLFTest(fit, contrast=my.contrasts[,"genotype"])

mdx_age_test <- topTags(mdx_age_fit, n="inf")
mdx_age_corrected_test <- topTags(mdx_age_corrected_fit, n="inf")
mdx_vs_wt_10d_test <- topTags(mdx_vs_wt_10d_fit, n="inf")
mdx_vs_wt_10w_test <- topTags(mdx_vs_wt_10w_fit, n="inf")
wt_age_test <- topTags(wt_age_fit, n='inf')
age_test <- topTags(age_fit, n='inf')
genotype_test <- topTags(genotype_fit, n='inf')

mdx_age_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(mdx_age_test$table), rownames(counts))]
mdx_age_corrected_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(mdx_age_corrected_test$table), rownames(counts))]
mdx_vs_wt_10d_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(mdx_vs_wt_10d_test$table), rownames(counts))]
mdx_vs_wt_10w_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(mdx_vs_wt_10w_test$table), rownames(counts))]
wt_age_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(wt_age_test$table), rownames(counts))]
age_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(age_test$table), rownames(counts))]
genotype_test$table$gene_id <- counts$`Gene stable ID`[match(rownames(genotype_test$table), rownames(counts))]

mdx_age_test$table$FDR[match(results$gene_stable_id, mdx_age_test$table$gene_id)] -> results$FDR_mdx_age
mdx_age_test$table$logCPM[match(results$gene_stable_id, mdx_age_test$table$gene_id)] -> results$logCPM_mdx_age
mdx_age_test$table$logFC[match(results$gene_stable_id, mdx_age_test$table$gene_id)] -> results$logFC_mdx_age
mdx_age_test$table$PValue[match(results$gene_stable_id, mdx_age_test$table$gene_id)] -> results$pvalue_mdx_age
mdx_age_test$table$`F`[match(results$gene_stable_id, mdx_age_test$table$gene_id)] -> results$Fvalue_mdx_age

mdx_age_corrected_test$table$FDR[match(results$gene_stable_id, mdx_age_corrected_test$table$gene_id)] -> results$FDR_age_interaction
mdx_age_corrected_test$table$logCPM[match(results$gene_stable_id, mdx_age_corrected_test$table$gene_id)] -> results$logCPM_age_interaction
mdx_age_corrected_test$table$logFC[match(results$gene_stable_id, mdx_age_corrected_test$table$gene_id)] -> results$logFC_age_interaction
mdx_age_corrected_test$table$PValue[match(results$gene_stable_id, mdx_age_corrected_test$table$gene_id)] -> results$pvalue_age_interaction
mdx_age_corrected_test$table$`F`[match(results$gene_stable_id, mdx_age_corrected_test$table$gene_id)] -> results$Fvalue_age_interaction

mdx_vs_wt_10d_test$table$FDR[match(results$gene_stable_id, mdx_vs_wt_10d_test$table$gene_id)] -> results$FDR_mdx_vs_wt_10d
mdx_vs_wt_10d_test$table$logCPM[match(results$gene_stable_id, mdx_vs_wt_10d_test$table$gene_id)] -> results$logCPM_mdx_vs_wt_10d
mdx_vs_wt_10d_test$table$logFC[match(results$gene_stable_id, mdx_vs_wt_10d_test$table$gene_id)] -> results$logFC_mdx_vs_wt_10d
mdx_vs_wt_10d_test$table$PValue[match(results$gene_stable_id, mdx_vs_wt_10d_test$table$gene_id)] -> results$pvalue_mdx_vs_wt_10d
mdx_vs_wt_10d_test$table$`F`[match(results$gene_stable_id, mdx_vs_wt_10d_test$table$gene_id)] -> results$Fvalue_mdx_vs_wt_10d

mdx_vs_wt_10w_test$table$FDR[match(results$gene_stable_id, mdx_vs_wt_10w_test$table$gene_id)] -> results$FDR_mdx_vs_wt_10w
mdx_vs_wt_10w_test$table$logCPM[match(results$gene_stable_id, mdx_vs_wt_10w_test$table$gene_id)] -> results$logCPM_mdx_vs_wt_10w
mdx_vs_wt_10w_test$table$logFC[match(results$gene_stable_id, mdx_vs_wt_10w_test$table$gene_id)] -> results$logFC_mdx_vs_wt_10w
mdx_vs_wt_10w_test$table$PValue[match(results$gene_stable_id, mdx_vs_wt_10w_test$table$gene_id)] -> results$pvalue_mdx_vs_wt_10w
mdx_vs_wt_10w_test$table$`F`[match(results$gene_stable_id, mdx_vs_wt_10w_test$table$gene_id)] -> results$Fvalue_mdx_vs_wt_10w

wt_age_test$table$FDR[match(results$gene_stable_id, wt_age_test$table$gene_id)] -> results$FDR_wt_age
wt_age_test$table$logCPM[match(results$gene_stable_id, wt_age_test$table$gene_id)] -> results$logCPM_wt_age
wt_age_test$table$logFC[match(results$gene_stable_id, wt_age_test$table$gene_id)] -> results$logFC_wt_age
wt_age_test$table$PValue[match(results$gene_stable_id, wt_age_test$table$gene_id)] -> results$pvalue_wt_age
wt_age_test$table$`F`[match(results$gene_stable_id, wt_age_test$table$gene_id)] -> results$Fvalue_wt_age

age_test$table$FDR[match(results$gene_stable_id, age_test$table$gene_id)] -> results$FDR_age
age_test$table$logCPM[match(results$gene_stable_id, age_test$table$gene_id)] -> results$logCPM_age
age_test$table$logFC[match(results$gene_stable_id, age_test$table$gene_id)] -> results$logFC_age
age_test$table$PValue[match(results$gene_stable_id, age_test$table$gene_id)] -> results$pvalue_age
age_test$table$`F`[match(results$gene_stable_id, age_test$table$gene_id)] -> results$Fvalue_age

genotype_test$table$FDR[match(results$gene_stable_id, genotype_test$table$gene_id)] -> results$FDR_genotype
genotype_test$table$logCPM[match(results$gene_stable_id, genotype_test$table$gene_id)] -> results$logCPM_genotype
genotype_test$table$logFC[match(results$gene_stable_id, genotype_test$table$gene_id)] -> results$logFC_genotype
genotype_test$table$PValue[match(results$gene_stable_id, genotype_test$table$gene_id)] -> results$pvalue_genotype
genotype_test$table$`F`[match(results$gene_stable_id, genotype_test$table$gene_id)] -> results$Fvalue_genotype

#na.omit(results) -> results
results %>% mutate(gene_name = coalesce(gene_name,transcript_stable_id)) -> results
results <- cbind(results, log_tpms)
write_xlsx(results,'/projects/portsmouth-darek-cerebellum/results/edger-results-table.xlsx')


### prepare for plotting ###
wt_age_test$table %>% filter(FDR_age < 0.00000000001) -> to_plot

log_tpms[genes_to_plot$gene_id,samples$id] -> to_plot
rownames(to_plot) <- results$gene_name[match(rownames(to_plot), results$gene_stable_id)]




##### plotting #####
mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

#group.names <- unique(as.character(samples.info$group[order(samples.info$diets)]))
#col.labels <- c(rep("", 5), group.names[1], rep(" ", 9), 
#                group.names[2], rep(" ", 9),
#                group.names[3], rep(" ", 9),
#                group.names[4], rep(" ", 9),
#                group.names[5], rep(" ", 9),
#                group.names[6], rep(" ", 9),
#                group.names[7], rep(" ", 9),
#                group.names[8], rep(" ", 4)
#)

cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


to_plot %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to_plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    main="",
    key=TRUE,
    keysize = 0.3,
    Colv = FALSE,
    scale="row",
    colsep = c(5,10,15),
    sepwidth = c(0.3,0.3),
    labRow=rownames(to_plot),
    #labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.7,
    offsetCol = -0.6,
    margins = c(5, 5)
  )




