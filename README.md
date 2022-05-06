# Portsmouth-darek-cerebellum

#### RNA seq of mouse cerebellum
5 samples/group - mdx mice and wt mice 10 days vs 10 weeks old


## Methods
This sections should be a description of preprocessin and analysis ready to be included in the publication


## Preprocessing
Intelliseq flow RNA-seq pipeline was used to generate bam files

- transcripts were quantified with Kallisto
- genes were quantified with feature count


## Analysis

#### Differential Gene Expression analysis
Run in R with EdgeR library


#### DMD transcript analysis

##### Various DMD transcript levels were quantified with multiple methods:

1. All transcripts with Kallisto
2. All transcripts with cufflinks
3. Three main transcripts with cufflinks (compare a few options)
4. Two main transcripts with cufflinks (compare a few options)

##### Transcript levels were then compared with two-way ANOVA

- are there differences between the algorithms?

