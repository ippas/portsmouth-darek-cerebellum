### This analysis is performed based on transcript counts table

require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)
require(stringr)
require(writexl)


transcripts <- read_tsv('/projects/portsmouth-darek-cerebellum/results/seq-results/cerebellum_transcript-level-tpm-annotated.tsv')

transcripts[transcripts$`Gene stable ID` == "ENSMUSG00000045103",] -> dmd

dmd$`Transcript stable ID`

dmd[,colnames(dmd)[12:31]] %>% rowMeans()

anova_results <- read_tsv('/projects/portsmouth-darek-cerebellum/results/)
