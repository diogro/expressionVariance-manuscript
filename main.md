---
title: The landscape of gene expression variation in human populations
bibliography: [./references.bib]
csl: [./evolution.csl]
mainfont: Skolar PE TEST Regular
sansfont: Skolar Sans PE TEST
mainfontoptions:
- Numbers=Lowercase
- Numbers=Proportional
geometry:
- top=20mm
- left=25mm
- right=25mm
- bottom=20mm
---

# Intro

Gene expression variation drives phenotypic difference.

Here we use public gene expression data sets to evaluate how the differences in gene expression variation is structures across several independent studies. We collected and compared the gene expression variation across several studies, and used the similarities across these studies to create a gene expression variation ranking, which orders genes from least variable to most variable. We then explore the expected drivers of this gene expression ranking, showing that both cis and trans regulation are involved with the determination of gene expression variance. 


# Methods

## Data sources

We selected 60 studies with large sample sizes from public gene expression repositories recount3 [@Wilks2021-uj] and Expression Atlas [@Papatheodorou2020-dn].

## Data processing

Filtering by min cpm and mean cpm. Variance stabilizing transformation from DESeq2 [@Love2014-mp]. Fixed effects correction. Outlier removal using [@Chen2020-fy]. Gene expression variance is measured in the residuals after fixed effect correction and outlier removal.

## Gene connectivity



## Variance correlation



# Results


![Standard deviation correlation PCoA](figures/sd_PCoA_plot.png){#fig:sd_pcoa}

\newpage

# Discussion

Gene expression variance is reasonably conserved across studies.
Gene expression variance is predictive of biological function.
Gene expression variance can be partially explained by genetic variation and genetic associations between gene expression.


# References





