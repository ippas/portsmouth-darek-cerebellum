### This analysis is performed based on Callisto transcripts table - ANOVA on TPM and EgdeR on counts
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)
require(preprocessCore)
require(edgeR)


samples <- data.frame(c('sample_BD1', 'sample_BD2', 'sample_BD3', 'sample_BD4', 'sample_BD5', 'sample_BW1', 'sample_BW2', 'sample_BW3', 'sample_BW4',
             'sample_BW5', 'sample_MD1', 'sample_MD2', 'sample_MD3', 'sample_MD4', 'sample_MD5', 'sample_MW1', 'sample_MW2', 'sample_MW3',
             'sample_MW4', 'sample_MW5'), c(rep('wt',10), rep('mdx', 10)),
             c(rep('10days', 5), rep('10weeks', 5), rep('10days', 5), rep('10weeks', 5)),
             c(rep('WT10D', 5), rep('WT10W', 5), rep('MDX10D', 5), rep('MDX10W', 5)))
colnames(samples) <- c('id', 'genotype', 'age', 'group')
rownames(samples) <- samples$id


##### PART 1 ANOVA analysis #####

tpms <- read_tsv('/projects/portsmouth-darek-cerebellum/results/seq-results/cerebellum_transcript-level-tpm-annotated.tsv')

##### transform tmps - this requires an older version of R (done on R 3.6)####

#tpms[,12:31] %>% data.matrix %>% normalize.quantiles() -> normalised
#colnames(normalised) <- colnames(tpms[,12:31])
#rownames(normalised) <- tpms$`Gene name`
#normalised <- data.frame(normalised)

#write_tsv(normalised, '/projects/portsmouth-darek-cerebellum/data/quantile-normalised-tpm.tsv')

read_tsv('/projects/portsmouth-darek-cerebellum/data/quantile-normalised-tpm.tsv') %>% data.matrix() -> normalised
rownames(normalised) <- tpms$`Transcript stable ID`
log_tpms <- log2(normalised + 1)


##### construct a results table #####
results <- data.frame(tpms$`Transcript stable ID`, tpms$`Gene name`, tpms$`Gene Synonym`)
colnames(results) <- c('transcript_stable_id', 'gene_name', 'gene_synonyms')

##### one-way anova to check the overall group effect #####

stat_one_way <- function(x) {
  y <- unlist(x, use.names = FALSE)
  summary(aov(y ~ samples$group)) %>%
    unlist %>%
    extract(c("Pr(>F)1"))
}


apply(log_tpms[,samples$id], 1, stat_one_way) %>% p.adjust(., method="fdr") -> results$one_way_fdr

# prepare to plot
results %>% filter(one_way_fdr < 0.0001) %>% select(transcript_stable_id) -> transcripts_to_plot
log_tpms[transcripts_to_plot$transcript_stable_id,samples$id] -> to_plot
rownames(to_plot) <- results$gene_name[match(rownames(to_plot), results$transcript_stable_id)]



##### two-way anova ####

stat <- function(x) {
  y <- unlist(x, use.names = FALSE)
  summary(aov(y ~ samples$genotype * samples$age)) %>%
    unlist %>%
    extract(c("Pr(>F)1", "Pr(>F)2", "Pr(>F)3" )) # genotype, age, interaction
}


apply(log_tpms[,samples$id], 1, stat)[1,] %>% p.adjust(., method="fdr") -> results$genotype_fdr
apply(log_tpms[,samples$id], 1, stat)[2,] %>% p.adjust(., method="fdr") -> results$age_fdr
apply(log_tpms[,samples$id], 1, stat)[3,] %>% p.adjust(., method="fdr") -> results$interaction_fdr


# prepare to plot
results %>% filter(age_fdr < 0.0000000001) %>% select(transcript_stable_id) -> transcripts_to_plot
log_tpms[transcripts_to_plot$transcript_stable_id,samples$id] -> to_plot
rownames(to_plot) <- results$gene_name[match(rownames(to_plot), results$transcript_stable_id)]


##### pairwise t-tests for requested comparisons #####

stat_paired_t <- function(x) {
                  if (is.na(x[1])) { c(NA, NA, NA, NA)
                  } else { 
                    pairwise.t.test(
                      x, samples$group, p.adjust.method = 'bonferroni') %>% 
                      unlist() %>% 
                      extract(c("p.value1", "p.value2", "p.value9")) #mdx vs wt 10 days, mdx vs mdx age, 10weeks mdx vs wt
                  }}
                

apply(log_tpms[,samples$id], 1, stat_paired_t)[1,] -> results$t_test_mdx_10d_vs_wt_10d
apply(log_tpms[,samples$id], 1, stat_paired_t)[2,] -> results$t_test_mdx_10d_vs_mdx_10w
apply(log_tpms[,samples$id], 1, stat_paired_t)[3,] -> results$t_test_mdx_10w_vs_wt_10w

# prepare to plot
results %>% filter(t_test_mdx_10d_vs_wt_10d < 0.001) %>% select(transcript_stable_id) -> transcripts_to_plot
log_tpms[transcripts_to_plot$transcript_stable_id,samples$id] -> to_plot
rownames(to_plot) <- results$gene_name[match(rownames(to_plot), results$transcript_stable_id)]



##### PART 2 EdgeR analysis #####

counts <- read_tsv('/projects/portsmouth-darek-cerebellum/results/seq-results/cerebellum_transcript-level-counts-annotated.tsv')
  
design <- model.matrix(~0+group, data=samples)
colnames(design) <- levels(as.factor(samples$group))


counts_edger <- DGEList(counts=counts[,12:31], group=samples$group)
keep <- filterByExpr(counts_edger)
counts_edger <- counts_edger[keep, , keep.lib.sizes=FALSE]
counts_edger <- calcNormFactors(counts_edger)


counts_edger <- estimateDisp(counts_edger, design)
fit <- glmQLFit(counts_edger, design)

my.contrasts <- makeContrasts(mdx_age_corrected_wt=(MDX10D-MDX10W-WT10D-WT10W), #mdx10d vs mdx10w but controlled for wt effects
  #mdx10d vs mdx10w not controlled
  mdx_age=MDX10W-MDX10D,
  #genotype effect 10 days
  mdx_vs_wt_10d=MDX10D-WT10D,
  #genotype effect 10 weeks
  mdx_vs_wt_10w=MDX10W-WT10W,
  wt_ctrl=WT10W+MDX10W-WT10D-MDX10D,
  levels=design)

mdx_age_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_age"])
mdx_age_corrected_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_age_corrected_wt"])
mdx_vs_wt_10d_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_vs_wt_10d"])
mdx_vs_wt_10w_fit <- glmQLFTest(fit, contrast=my.contrasts[,"mdx_vs_wt_10w"])
wt_fit <- glmQLFTest(fit, contrast=my.contrasts[,"wt_ctrl"])

mdx_age_test <- topTags(mdx_age_fit, n="inf")
mdx_age_corrected_test <- topTags(mdx_age_corrected_fit, n="inf")
mdx_vs_wt_10d_test <- topTags(mdx_vs_wt_10d_fit, n="inf")
mdx_vs_wt_10w_test <- topTags(mdx_vs_wt_10w_fit, n="inf")
wt_test <- topTags(wt_fit, n='inf')

mdx_age_test$table$transcript_id <- tpms$`Transcript stable ID`[match(rownames(mdx_age_test$table), rownames(tpms))]
mdx_age_corrected_test$table$transcript_id <- tpms$`Transcript stable ID`[match(rownames(mdx_age_corrected_test$table), rownames(tpms))]
mdx_vs_wt_10d_test$table$transcript_id <- tpms$`Transcript stable ID`[match(rownames(mdx_vs_wt_10d_test$table), rownames(tpms))]
mdx_vs_wt_10w_test$table$transcript_id <- tpms$`Transcript stable ID`[match(rownames(mdx_vs_wt_10w_test$table), rownames(tpms))]
wt_test$table$transcript_id <- tpms$`Transcript stable ID`[match(rownames(wt_test$table), rownames(tpms))]


### prepare for plotting ###
wt_test$table %>% filter(PValue < 0.01) %>% select(transcript_id) -> transcripts_to_plot
log_tpms[transcripts_to_plot$transcript_id,samples$id] -> to_plot
rownames(to_plot) <- results$gene_name[match(rownames(to_plot), results$transcript_stable_id)]




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
    Colv = FALSE,
    scale="row",
    colsep = c(5,10,15),
    sepwidth = c(0.3,0.3),
    labRow=rownames(to_plot),
    #labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.6,
    offsetCol = 0.1
  )


#### export results ####

results <- cbind(results, log_tpms)
write.table(results, ('/projects/portsmouth-darek-cerebellum/results/aov-all.csv')
)
