### This analysis is performed based on transcript counts table

require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)
require(writexl)
require(reshape2)

samples <- data.frame(c('sample_BD1', 'sample_BD2', 'sample_BD3', 'sample_BD4', 'sample_BD5', 'sample_BW1', 'sample_BW2', 'sample_BW3', 'sample_BW4',
                        'sample_BW5', 'sample_MD1', 'sample_MD2', 'sample_MD3', 'sample_MD4', 'sample_MD5', 'sample_MW1', 'sample_MW2', 'sample_MW3',
                        'sample_MW4', 'sample_MW5'), c(rep('wt',10), rep('mdx', 10)),
                      c(rep('10days', 5), rep('10weeks', 5), rep('10days', 5), rep('10weeks', 5)),
                      c(rep('WT10D', 5), rep('WT10W', 5), rep('MDX10D', 5), rep('MDX10W', 5)))
colnames(samples) <- c('id', 'genotype', 'age', 'group')
rownames(samples) <- samples$id


### load and filter different transcripts
transcripts_kallisto <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/results/seq-results/cerebellum_transcript-level-tpm-annotated.tsv'
  )

dmd_kallisto <- transcripts_kallisto[transcripts_kallisto$`Gene stable ID` == "ENSMUSG00000045103",]
all_dmd_transcripts <- dmd_kallisto$`Transcript stable ID`
dmd_kallisto <- dmd_kallisto[,c(2,12:31)]

transcripts_cuff_1 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-1/isoforms.fpkm_table.tsv'
  )
dmd_cuff_1 <- transcripts_cuff_1[transcripts_cuff_1$tracking_id %in% all_dmd_transcripts,]

transcripts_cuff_2 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-2/isoforms.fpkm_table.tsv'
)
dmd_cuff_2 <- transcripts_cuff_2[transcripts_cuff_2$tracking_id %in% all_dmd_transcripts,]

transcripts_cuff_3 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-3/isoforms.fpkm_table.tsv'
)
dmd_cuff_3 <- transcripts_cuff_3[transcripts_cuff_3$tracking_id %in% all_dmd_transcripts,]

transcripts_cuff_4 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-4/isoforms.fpkm_table.tsv'
)
dmd_cuff_4 <- transcripts_cuff_4[transcripts_cuff_4$tracking_id %in% all_dmd_transcripts,]


transcripts_cuff_5 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-5/isoforms.fpkm_table.tsv'
)
dmd_cuff_5 <- transcripts_cuff_5[transcripts_cuff_5$tracking_id %in% all_dmd_transcripts,]

transcripts_cuff_6 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-6/isoforms.fpkm_table.tsv'
)
dmd_cuff_6 <- transcripts_cuff_6[transcripts_cuff_6$tracking_id %in% all_dmd_transcripts,]


transcripts_cuff_7 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-7/isoforms.fpkm_table.tsv'
)
dmd_cuff_7 <- transcripts_cuff_7[transcripts_cuff_7$tracking_id %in% all_dmd_transcripts,]


transcripts_cuff_8 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-8/isoforms.fpkm_table.tsv'
)
dmd_cuff_8 <- transcripts_cuff_8[transcripts_cuff_8$tracking_id %in% all_dmd_transcripts,]

transcripts_cuff_9 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-9/isoforms.fpkm_table.tsv'
)
dmd_cuff_9 <- transcripts_cuff_9[transcripts_cuff_9$tracking_id %in% all_dmd_transcripts,]


### plot different transcript levels:

dmd_kallisto %>% melt() -> dmd_kallisto
colnames(dmd_kallisto) <- c('id', 'sample', 'fpkm')
dmd_cuff_1 %>% melt() -> dmd_cuff_1
colnames(dmd_cuff_1) <- c('id', 'sample', 'fpkm')
dmd_cuff_2 %>% melt() -> dmd_cuff_2
colnames(dmd_cuff_2) <- c('id', 'sample', 'fpkm')
dmd_cuff_3 %>% melt() -> dmd_cuff_3
colnames(dmd_cuff_3) <- c('id', 'sample', 'fpkm')
dmd_cuff_4 %>% melt() -> dmd_cuff_4
colnames(dmd_cuff_4) <- c('id', 'sample', 'fpkm')
dmd_cuff_5 %>% melt() -> dmd_cuff_5
colnames(dmd_cuff_5) <- c('id', 'sample', 'fpkm')
dmd_cuff_6 %>% melt() -> dmd_cuff_6
colnames(dmd_cuff_6) <- c('id', 'sample', 'fpkm')
dmd_cuff_7 %>% melt() -> dmd_cuff_7
colnames(dmd_cuff_7) <- c('id', 'sample', 'fpkm')
dmd_cuff_8 %>% melt() -> dmd_cuff_8
colnames(dmd_cuff_8) <- c('id', 'sample', 'fpkm')
dmd_cuff_9 %>% melt() -> dmd_cuff_9
colnames(dmd_cuff_9) <- c('id', 'sample', 'fpkm')


dmd_kallisto %>% mutate(
  option = rep('kallisto', nrow(dmd_kallisto))) %>%
  mutate(group = samples$group[match(dmd_kallisto$sample, samples$id)]) -> dmd_kallisto

dmd_cuff_1 %>% mutate(
  option = rep('dmd_cuff_1', nrow(dmd_cuff_1)), group = samples$group[match(dmd_cuff_1$sample, samples$id)]
  ) -> dmd_cuff_1

dmd_cuff_2 %>% mutate(
  option = rep('dmd_cuff_2', nrow(dmd_cuff_2)), 
  group = samples$group[match(dmd_cuff_2$sample, samples$id)]
  ) -> dmd_cuff_2

dmd_cuff_3 %>% mutate(
  option = rep('dmd_cuff_3', nrow(dmd_cuff_3)), 
  group = samples$group[match(dmd_cuff_3$sample, samples$id)]
  ) -> dmd_cuff_3

dmd_cuff_4 %>% mutate(
  option = rep('dmd_cuff_4', nrow(dmd_cuff_4)), 
  group = samples$group[match(dmd_cuff_4$sample, samples$id)]
  ) -> dmd_cuff_4

dmd_cuff_5 %>% mutate(
  option = rep('dmd_cuff_5', nrow(dmd_cuff_5)), 
               group = samples$group[match(dmd_cuff_5$sample, samples$id)]
  ) -> dmd_cuff_5

dmd_cuff_6 %>% mutate(
  option = rep('dmd_cuff_6', nrow(dmd_cuff_6)), 
               group = samples$group[match(dmd_cuff_6$sample, samples$id)]
  ) -> dmd_cuff_6

dmd_cuff_7 %>% mutate(
  option = rep('dmd_cuff_7', nrow(dmd_cuff_7)), 
  group = samples$group[match(dmd_cuff_7$sample, samples$id)]
) -> dmd_cuff_7

dmd_cuff_8 %>% mutate(
  option = rep('dmd_cuff_8', nrow(dmd_cuff_8)), 
  group = samples$group[match(dmd_cuff_8$sample, samples$id)]
) -> dmd_cuff_8

dmd_cuff_9 %>% mutate(
  option = rep('dmd_cuff_9', nrow(dmd_cuff_9)), 
  group = samples$group[match(dmd_cuff_9$sample, samples$id)]
) -> dmd_cuff_9

dmds <- bind_rows(
  list(
    dmd_cuff_1, dmd_cuff_3, dmd_cuff_4, dmd_cuff_5, dmd_cuff_6,
  dmd_cuff_7, dmd_cuff_8, dmd_cuff_9
  ))
dmds$id_short <- str_sub(dmds$id, start= -4)


ggplot(data=dmds, aes(x=id_short, y=fpkm, color=option)) +
  geom_point(position=position_dodge(width=0.3)) +facet_grid(~group)


### dmd option 9 - with all transcripts looks all right, thus we will take it, transform it and calculate statistics
transcripts_cuff_9 <- read_tsv(
  '/projects/portsmouth-darek-cerebellum/data/cuffcuant-options/cuffnorm-9/isoforms.fpkm_table.tsv'
)
dmd <- transcripts_cuff_9[transcripts_cuff_9$tracking_id %in% all_dmd_transcripts,]

dmd[,2:21] %>% data.matrix() -> dmd.to.transform
log2(dmd.to.transform + 1) -> dmd[,2:21]



#### aov

stat <- function(x) {
  #y <- unlist(x, use.names = FALSE)
  summary(aov(as.numeric(x) ~ samples$genotype * samples$age)) %>%
    unlist -> res
    return(res[c("Pr(>F)1", "Pr(>F)2", "Pr(>F)3" )]) # genotype, age, interaction
}


apply(data.matrix(dmd[,samples$id]), 1, stat)[1,] -> dmd$genotype_fdr
apply(data.matrix(dmd[,samples$id]), 1, stat)[2,] -> dmd$age_fdr
apply(data.matrix(dmd[,samples$id]), 1, stat)[3,] -> dmd$interaction_fdr


write.csv(dmd, '/projects/portsmouth-darek-cerebellum/results/seq-results/dmd-quant.csv')

