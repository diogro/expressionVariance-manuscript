---
title: "Characterizing the landscape of gene expression variance in humans"
# author:
  # - Scott Wolf, Diogo Melo, Kristina M. Garske, Luisa Pallares, Julien Ayroles
author:
  - Scott Wolf:
      institute: lsi
      email: swwolf@princeton.edu
      orcid: 0000-0003-4397-1395
      equal_contributor: "yes"
  - Diogo Melo:
      institute:
        - lsi
        - eeb
      email: damelo@princeton.edu
      orcid: 0000-0002-7603-0092
      equal_contributor: "yes"
  - Kristina M. Garske:
      institute: lsi
      orcid: 0000-0002-7228-8125
  - Luisa F. Pallares:
      institute: fml
      orcid: 0000-0001-6547-1901
  - Julien F. Ayroles:
      institute:
        - lsi
        - eeb
      email: jayroles@princeton.edu
      correspondence: "yes"
      orcid: 0000-0001-8729-0511
institute:
  - lsi:
      name: Lewis-Sigler Institute for Integrative Genomics, Princeton University
  - eeb:
      name: Department of Ecology and Evolutionary Biology, Princeton University
  - fml:
      name: Friedrich Miescher Laboratory, Max Planck Society
# classoption: twocolumn
output: pdf_document
geometry:
- top=20mm
- left=25mm
- right=25mm
- bottom=20mm
header-includes:
- \usepackage[left]{lineno}
- \linenumbers
- \usepackage[backref=true,style=authoryear]{biblatex}
- \DefineBibliographyStrings{english}{backrefpage = {page}, backrefpages = {pages}}
- \usepackage{multicol}
- \usepackage{setspace}
- \usepackage{float}
- \usepackage{afterpage}
- \usepackage{stfloats}
- \usepackage{graphicx}
- \newcommand{\hideFromPandoc}[1]{#1}
- \hideFromPandoc{ \let\Begin\begin \let\End\end}
- \addbibresource{references.bib}
link-citations: yes
mainfont: Skolar PE TEST Regular
mainfontoptions:
- Numbers=Lowercase
- Numbers=Proportional
csl: ./pandoc/cse.csl
sansfont: Skolar Sans PE TEST
bibliography: ./references.bib
abstract: Gene expression variance has been linked to organismal function and fitness but remains a commonly neglected aspect of gene expression research. We lack a comprehensive view of the patterns of variance across genes, and how this variance is linked to context-specific gene regulation and gene function. Here, we use large publicly available human RNA-seq data sets to investigate the landscape of gene expression variance. In particular, we ask if there are consistently more or less variable genes across tissues and across data sets and what mechanisms drive these patterns. We show that gene expression variance is broadly similar across tissues and studies. We use this similarity to create both global and within-tissue rankings of variation, which we use to show that function, sequence variation, and gene regulatory signatures contribute to gene expression variance. Gene expression variance is strongly predictive of gene function, with low-variance genes being associated with fundamental cell processes, and high-variance genes being linked to responding to the environment. Our results show differences in the regulatory mechanisms of high and low gene expression variance, in addition to a clear link between function and gene expression variance, suggesting that these differences are adaptive. We expect these results will help to place the pattern of variation at the center of our understanding of molecular phenotypes.
---
<!--
# Author Summary

Required for PLOS Genetics... -->

# Introduction

Molecular phenotypes such as gene expression are powerful tools for understanding physiology, disease, and evolutionary adaptations.
In this context, average trait values are usually the focus of investigation, while variation is treated as a nuisance [@De_Jong2019-po].
However, variability is ubiquitous in nature and is, alongside robustness, a fundamental feature of most complex systems. 
Gene expression variance can be directly involved in determining fitness [@Fraser2004-sv; @Wang2011-ts], and changes in the associations between gene expression can be indicative of disease, even in the absence of changes in mean expression [@Lea2019-pq].
Yet, few studies have investigated variance itself as a regulated trait, and the consequences of transcriptional variance on complex traits and diseases [@Li2010-qs; @Mar2011-dr; @Hagai2018-fu].
From an evolutionary perspective, the availability of transcriptional variance can drive phenotypic variation - the substrate for evolutionary change [@Hansen2021-zo].
Previous work also demonstrate that  the genetic architecture of variance it self can evolve [@Bruijning2020-bf].
Focusing on the landscape of gene expression variance, and how variable it is across genes and across human populations is a neglected avenue for understanding biological evolution, disease traits, and our relation to the environment. 
In particular, we still have a very poor understanding of whether some genes are inherently more variable than others, or if the pattern of transcriptional variance is consistent across populations. 
Nor do we know how the degree of variance is linked to context-specific gene regulation and gene function.

From a mechanistic perspective, several competing forces act to shape gene expression variance [@Houle1998-mj; @Bruijning2020-bf], and the outcome of the interaction between these processes is still poorly understood [@Hansen2011-es].
At the genomic sequence level, we expect the influx of new mutations to increase the variance, while the selective removal of these polymorphisms, via purifying selection or selective sweeps, would decrease variation. 
The extent to which mutations and polymorphisms can contribute to variance depends on aspects of genetic architecture, like mutational target size and the presence of epistatic interactions.
From a quantitative trait perspective, stabilizing selection should decrease variation around an optimal value, and directional selection can lead to a transient increase in variance while selected alleles sweep to fixation, followed by a reduction in variance as these alleles become fixed.
Pleiotropic effects are also important, as they allow selection on one trait to influence the variance of other traits [@Wagner1997-hw; @Pavlicev2011-xm].
Both indirect effects of directional selection on variance opens the possibility that the main driver of gene expression variance is not direct selection on variance but indirect effects due to selection on trait means [@Hansen2011-es].
To what extent these processes shape gene expression variance is an open question.
If homogeneous selection across groups is the main driver of gene expression variance, we would expect to have consistently more or less variable genes.
If idiosyncratic selection patterns and context-specific environmental interactions are more important, we could observe large differences in gene expression variance across groups.

Even within individuals, gene expression is also variable across tissues [@GTEx2017-xb].
To what extent differential expression (i.e., differences in mean expression level) translate into differences in expression variance is not clear.
Higher expression could lead to higher variation, but other processes could also affect transcriptional variance.
For example, if a gene is expressed in more than one tissues, stabilizing selection on gene expression could be more intense depending on the role of that gene in a particular tissue, leading to a local reduction in variation that causes differences in variance across tissues that is not necessarily correlated with mean expression.
Alternatively, expression variation across tissues could be tightly coupled, and in this example, selection in one tissue would lead to a reduction in variance across tissues, resulting in a consistent pattern of variation.
Alemu et al. [@Alemu2014-jo] used microarray data from several human tissues to show that epigenetic markers were linked to gene expression variation and that these markers were variable across tissues and between high- and low-variance genes.

Here, we use publicly available human gene expression data sets to evaluate how the differences in gene expression variance are structured across independent samples.
By comparing the gene expression variance measured across many studies, we show that the patterns of gene expression variance are broadly similar across studies and tissues.
We used the observed similarities across these studies to create an across-study gene expression variance ranking, which orders genes from least variable to most variable.
We then integrate various functional annotations as well as sequence variation to probe the drivers of this across-study ranking.
Finally, we explore the link between gene expression variance and biological function by leveraging gene ontology and disease annotations.


# Results


![Overview of the distribution of transcriptional variance across studies. A. Heatmap showing the correlation in transcriptional variance across studies (as the Spearman correlation of standard deviations). Pairs of studies with more similar patterns of gene expression variance have higher correlations. Studies are shown in the same order as in @fig:corr_model, panel A; B. Histogram of the correlations shown in the previous panel; C. Standard deviation correlation PCoA. There is no clear structuring of the studies with respect to their source, which is indicated by the colors; D. Density plot of standard deviations after z-normalization. The inset plot shows the distribution of mean-centered standard deviations grouped by study without normalization. The corresponding rug plots show the location of the highest-ranking gene in standard deviation rank (*HBB*) (right, blue) and lowest (*WDR33*) (left, red).](figures/fig1.png){#fig:sd_corr}

## Datasets

We use 57 publicly available human gene expression RNA-seq data sets which are derived from the publications listed in table \ref{tab1} of the [Methods](#Methods) section. We only use dataset from population samples, so we exclude studies using single-cell data. We also did not use studies those with no corresponding publication and studies without sample-level metadata. Several data sets were derived from two large consortiums: GTEx [@GTEx2017-xb] and TCGA [@tcga2013-gx], and we note the origin of the data sets in the figures.
We refer to data sets and studies interchangeably, and so each tissue in GTEx is referred to as a different study.

## Gene expression standard deviations

For each study, gene expression standard deviations (SDs) were calculated using a unified pipeline that normalized the mean-variance relation in read-count data, controlled for batch effects, and removed outliers (see [Methods](#Methods) for details).
Spearman correlations ($\rho_s$) between gene expression SDs reveal a broadly similar rank of gene expression variance, such that genes that are most variable in one study tend to be most variable in all studies (@fig:sd_corr A and B).
A principal coordinate analysis [@Gower1966-dk] using $|1 - \rho_s|$ as an between study distance measure does not show clearly delineated groups, but GTEx  and TCGA studies are clustered among themselves and close together (@fig:sd_corr C).
This clustering indicates some effect of study source on the similarity between gene expression SD across studies, which we explore in detail below.
The observed range of gene expression SD across genes is variable across studies but can be normalized so that the distributions are comparable (@fig:sd_corr D).
Given that the correlations across studies are mostly positive and high (75% of correlations are between 0.45 and 0.9), indicating similar ordering of the genes, we seek to summarize the differences in variance across genes by using a single cross-study rank, averaging the ordering across all studies.
To create this rank, we used the score of each gene in the first principal component of the Spearman correlation matrix.
Ordering genes using these scores generate a ranked list of genes, with the most variable genes having the highest rank.
We create a similar across-study rank for mean expression.
The red and blue ticks at the bottom of @fig:sd_corr D show the positions on the SD distributions of the least and most variable genes in our variance rank. The position of these highlighted genes in the SD distributions illustrates how the extremes of the rank are indeed some of the least and most variable genes across all studies.
We also create a set of tissue-specific SD ranks, which use the same procedure outlined above but using only studies that were performed on the same tissue.
This creates a series of gene ranks, one for each sampled tissue, which describes the gene expression SD rank in that particular tissue.
Both tissue-specific and across-study ranks are available in the Supporting Information.

![Modeling the correlations between transcriptional variance across studies. The panels show coefficient estimates from a linear model using the among studies Spearman correlations between gene expression SDs as the response variable. These correlations are shown in @fig:sd_corr A and B. In the linear model (see Methods for model equation), correlations are Fisher z-transformed. Study source and tissue are added as fixed effects. Coefficient estimates are shown with 50% and 95% credibility intervals. Panel A: The per-study random effect which accounts for the non-independence between the pairwise correlation values and estimates the characteristic contribution of each study to these correlations. For example, the lowest estimate among these parameters, which corresponds to the study \textsc{bone marrow} (from GTEx), indicates that correlations involving this study tend to be lower than the others. Panels B and C: Fixed effect estimates for the effects of tissue congruence and study-source effect. In (B) we see that correlations among studies that use the same tissue are slightly higher; and (C) correlations involving studies in the "Misc." category (non-GTEx and non-TCGA) tend to be lower, while comparisons involving GTEx and TCGA are higher.](figures/correlationModeling.png){#fig:corr_model}



## What drives differences in gene expression variance?

To characterize what factors may explain differences in across-study similarity, we directly modeled the correlations across-study using a mixed-effect linear model designed to account for the non-independence in pairwise correlation data [@Dias2021-wk; @Dias2021-hb].
In this model (see [Methods](#Methods)), we use a random effect for individual study ID, a fixed effect for pairwise tissue congruence (whether a comparison is within the same tissue or different tissue), and a fixed effect for pairwise study source (which pair of sources among GTEx, TCGA, and miscellaneous is involved in a comparison) as predictors of the correlations (see [Methods](#Methods)).
This modeling (@fig:corr_model) shows that comparisons of studies within GTEx and TCGA have on average higher values of $\rho_s$, but also that comparing studies across GTEx and TCGA also shows a mild increase in the average correlation (@fig:corr_model C).
Correlations that do not involve studies from TCGA and GTEx (marked as "Misc.") are on average lower (@fig:corr_model C).
Since TCGA and GTEx are independent, this mild effect on the similarities could be due to the quality of the data coming from these two large projects.
Tissue also affects the similarity between gene expression SD, with studies using the same tissue being, on average, more similar (@fig:corr_model B).
However, all of these pairwise effects are mild, and the largest effects on the correlations are those associated with individual studies, in particular some specific tissues, i.e., comparisons involving \textsc{bone marrow} (from GTEx) and study \textsc{srp057500} (which used platelets) are on average lower (@fig:corr_model A).
The only negative correlation we observe is between these two studies, which also appear further away in the PCoA plot in @fig:sd_corr C.

## Does biological function explain variance in expression?

As a first step toward explaining the factors that may drive variation in variability between transcripts, we focused on the top 5% most variable and the bottom 5% least variable genes in our ranking (560 genes in each group) and performed a Gene Ontology (GO) enrichment analysis within each group.
This analysis allowed us to establish the representative functions of these consistently high and low-variance genes.
In total, we found 59 enriched terms in the low variance genes, and 738 enriched terms in the high-variance genes (using a hypergeometric test and Benjamini-Hochberg (BH) adjusted p-value threshold of 10\textsuperscript{-3}; see S2 Table for a complete listing).
Among the 5% most variable genes we observe enrichment for biological processes such as immune function, response to stimulus, maintenance of homeostasis, and tissue morphogenesis (@fig:go_tails A).
Notably, we see a 7.7-fold enrichment for genes that encode secreted proteins in the top 5% most variable genes, relative to all other genes (hypergeometric test, p < 10\textsuperscript{-3}).
Among the 5% least variable genes we see enrichment for housekeeping functions such as mRNA processing, cell cycle regulation, methylation, histone modification, translation, transcription, and DNA repair (@fig:go_tails B); and accordingly, we find that previously characterized human housekeeping genes [@Hounkpe2020-yq] are enriched within the 5% least variable genes 2.0-fold relative to all other genes (hypergeometric test, p < 10\textsuperscript{-3}).
The genes exhibiting the lowest variance (lowest 5%) are also enriched for those that have been previously shown to have a high probability of being loss-of-function intolerant (pLI) [@lek2016analysis] (1.2-fold enrichment, hypergeometric test, p < 10\textsuperscript{-3}).
Genes with a high pLI have been shown to be important in housekeeping functions, and have higher mean expression values across a broad set of tissues and cell types [@lek2016analysis]. The observation that genes with low variance are enriched for both housekeeping genes and genes with high pLI is consistent with this previous report; and we further see that the mean expression of genes positively correlates with pLI (Partial Spearman correlation $\rho_s$ = 0.32, p < 10\textsuperscript{-3}), showing the opposite relationship between variance and mean expression when considering pLI.

In the previous analysis, we explored the relation between transcriptional variance and function by starting from the extremes in the variational distribution and searching for GO enrichment among these high- and low-variance genes.
We also approach the problem from the opposite direction, starting from the genes associated with each GO term and searching for enrichment for high- or low-variance genes.
To do this, we gathered all biological process GO terms in level 3 (i.e. terms that are at a distance of 3 from the top of the GO hierarchy).
Using only the set of genes that are associated with at least one of these level-3 terms, we separated the genes into expression variance deciles, with the first decile having the lowest variance.
We then counted how many genes in each decile have been associated with each term.
If variance rank is not associated with the GO annotations, terms should have an equal proportion of genes in each decile.
We measured how far from this uniform allocation each term is by measuring the Shannon entropy of the proportion of genes in each decile.
Higher entropy is associated with a more uniform distribution of genes across deciles.
GO terms with low entropy indicate some deciles are over-represented in the genes associated with that term.
We also measured skewness for each term, which should be zero if no decile is over-represented, negative if high-variance terms are over-represented, and positive if low-variance deciles are over-represented.
The relation between skewness and entropy for each GO term can be seen in @fig:skew_entropy .
Positive-skew low-entropy terms, those enriched with low-variance genes, are associated with housekeeping functions, like RNA localization, translation initiation, methylation and chromosome segregation (@fig:go_skewness  A).
Likewise, terms with negative skew and low entropy, enriched for high-variance genes, are related to immune response, tissue morphogenesis, chemotaxis---all dynamic biological functions related to interacting with the environment (@fig:go_skewness  B).

Both GO analyses suggest a strong association between biological function and the degree of transcriptional variance.
Genes associated with baseline fundamental functions, expected to be under strong stabilizing selection, are also low-variance; high-variance genes are associated with responding to external stimuli (i.e., tissue reorganization and immune response).


![Gene set enrichment analyses testing for over-representation of gene ontology categories in the upper and lower 5\% quantiles of the gene variance rank. (A) High-variance genes are enriched for terms related to immune function, response to wounding, blood vessel morphogenesis, and inflammatory response. In contrast, (B) low-variance genes are associated with translation, control of methylation, RNA processing, chromosome separation, and other cell housekeeping functions. All displayed terms are significant with a 5% FDR corrected p-value below 10\textsuperscript{-3}.](figures/local_go_lowerUpper.png){#fig:go_tails}

![Relationship between skew and entropy of rank decile distributions for each GO term. High entropy terms, to the right of the plot, are associated with a more egalitarian proportion of genes in each of the SD rank deciles. Terms on the left of the plot are associated with more genes in some particular decile. The skewness in the y-axis measures if the high- or low-variance deciles are more represented for a particular term. Terms on the positive side of the y-axis are associated with low-variance genes, and terms on the negative side of the y-axis are associated with high-variance genes. The GO terms are filtered for gene counts greater than 100, as in @fig:go_skewness. Some of the top high- and low-skewness terms are labeled for illustration.](figures/GOterm_entropy_by_skewness.png){#fig:skew_entropy}


![Distributions of decile ranks of level-3 GO terms. Each plot shows the count of genes in each decile of the rank. We only use GO terms that are associated with at least 100 genes, and sort these terms by the skewness of the distribution. The top panel shows the 5 most positively skewed terms and the bottom panel shows the 5 most negatively skewed terms.](figures/GOterm_decile_barplot.png){#fig:go_skewness}



## Sequence variation and gene expression connectivity

Next, we explore the role evolutionary forces may play in driving variation in transcriptional variance across genes.
To that end, we use gene-level summary statistics focused on three measures: nucleotide diversity ($\pi$), gene expression connectivity, and the estimated proportion of substitutions that are adaptive ($\alpha$).
For all the correlations in this section we use partial Spearman correlations that include the mean gene expression rank as a covariate, which accounts for any residual mean-variance correlation.
Nucleotide diversity is used as a proxy for the effect of cis-regulatory variation on transcriptional variance. We expect variance to increase with nucleotide diversity and indeed, we find a positive correlation of 0.184 between $\pi$ and the gene expression rank (p < 10\textsuperscript{-3}).
Gene-gene connectivity, a proxy for gene regulatory interactions and selective constraints [@Mahler2017-bb], in turn, should be negatively correlated with transcriptional variance, as highly connected genes are expected to be more constrained in their variation.
Consistent with this expectation, we find a negative partial Spearman correlation of -0.024 between connectivity and the gene expression rank (p $\approx$ 6 \times 10\textsuperscript{-3}).
Finally, we also find a negative partial Spearman correlation of -0.046 (p $\approx$ 10\textsuperscript{-3}) for $\alpha$.
Although all associations are significant and in the expected direction, their effect sizes are very small, suggesting a weak link between these broad measures and transcriptional variance.

## How do molecular signatures of gene regulation relate to gene expression variance?

We assess how local epigenetic features relate to gene expression variance. We use each gene, including the surrounding 10 kb on both ends, to calculate the proportion of gene regions that correspond to epigenetic signatures of gene regulation defined through ChromHMM [@ernst2012chromhmm] chromatin states. Chromatin states associated with distal (i.e., non-promoter) gene regulation are positively correlated with the across-study variance rank, regardless of whether the regulatory effect on gene expression is positive or negative (@fig:lineplot; see also "across-study" correlations in supplementary fig. 1A). For example, both the proportion of gene regions made up of enhancers and repressed genomic states are positively correlated with gene expression variance (BH adjusted spearman correlation p < 0.05). Histone modifications associated with active promoters, as well as transcribed states, are inversely correlated with gene expression variance (supplementary fig. 1A), whereas they are positively correlated with the mean rank (supplementary fig. 1B). Taken together, these results are compatible with gene expression variance being more associated with distal (i.e., non-promoter) gene regulation, rather than the overall active transcriptional state of a gene region, as is the case with mean gene expression.

![Proportion of gene regions made up of ChromHMM chromatin states for genes in the top and bottom 5% of the across-study variance rank metric. Line plot contrasts the proportion of gene regions made up of the indicated chromatin states for genes in the top and bottom 5% of the across-study variance rank metric. Ends denote the median proportion of gene regions made up of the chromatin state, and error bars represent the standard error of the mean. States colored black are not significant, all others exhibit significant differences in gene region made up of the chromatin state for genes in the top and bottom 5% of the variance rank metric (BH adjusted Wilcoxon signed-rank test, p < 0.05). Het indicates heterochromatin; TSS, transcription start sites; znf, zinc finger genes.](figures/top_bottom_5pVarrank_fxnlGen_lineplot_10kb_simplified.png){#fig:lineplot}

We also explore the relationship between tissue-specific ChromHMM chromatin states and SD rank and contrast these tissue-level analyses to the across-study analysis outlined above. Many of the across-study correlations are recapitulated at the tissue-specific level, including a strong and highly consistent positive correlation between the proportion of gene regions made up of enhancer states and that gene’s expression variance, and an inverse relationship between gene expression variance and histone marks associated with gene transcription  (supplementary fig. 1A). Two blood associations stand out as being different from the consistent effects across the other tissue-level and across-study associations. First, the weak (i.e., histone marks associated with both activating and repressive functions) promoter state is positively correlated with gene expression variance in all comparisons except blood. Second, the consistent inverse correlation of gene expression variance with weak transcription is reversed in blood, such that there is a positive correlation between histone marks associated with weak transcription and blood gene expression variance (supplementary fig. 1A). Taken together, these results suggest that, rather than genes with a bivalent promoter state (i.e., poised genes) exhibiting more expression variance, blood high-variance genes are more likely already expressed at basal levels (i.e., weakly transcribed), as discussed previously [@Rogatsky2014-bi].

Immediate early genes (IEGs) respond quickly to external signals without requiring _de novo_ protein synthesis, and a bivalent state has been reported to be associated with IEG promoters [reviewed in @Bahrami2016-dx]. Given our results that genes with high expression variance are enriched for cellular signaling and response mechanisms (@fig:go_tails A), and bivalent promoter states are correlated with the gene expression variance rank (supplementary fig. 1A), we hypothesized that IEGs would be enriched within genes in the top expression variance ranks. This was the case for all tissue-level gene expression variance ranks (enrichment ratios range from 3.3-8.8, Bonferroni-adjusted hypergeometric test, p < 0.05), except for in blood (enrichment ratio = 1.2, hypergeometric test, p = 0.3). Thus, once again blood stands out when attempting to understand genomic regulatory drivers of expression variance. In all, while high-variance genes are generally shared across tissues and enriched for immune and environmental signaling pathways, it seems that the gene regulatory mechanisms governing their expression are distinct between immune cell types and other tissues studied here.


## Linking expression variance and disease

To explore the link between transcription variance and genes known to be associated with human diseases, we used a dataset designed to provide causal relationships between gene expressions and complex traits/diseases (based on a probabilistic transcriptome-wide association study (PTWAS) [@Zhang2020-cl]).
Using the list of significant gene-disease pairs at 5% FDR provided by Zhang et al. [-@Zhang2020-cl], we performed a hypergeometric enrichment test for the top 5% high- and low-variance genes in our across-study rank and in all tissue-specific gene variance ranks.
We use both across-study and tissue-specific ranks because some genes only appear in the tissue-specific rank due to their limited tissue-specific gene expression.
In the high-variance group, we find no enrichment in the across-study rank, but we do find enrichment of genes annotated for allergy, immune disease, and endocrine system disease among the high-variance genes in several tissue-specific variance ranks.
For example, among high-variance genes in the colon rank, we see enrichment for endocrine system disease (1.77-fold, hypergeometric test, p < 10\textsuperscript{-4}).
Among high-variance genes in the immune cell tissue rank, we see enrichment for endocrine system disease (1.67-fold, hypergeometric test, p < 10\textsuperscript{-3}),
allergy (1.7-fold, hypergeometric test, p < 10\textsuperscript{-3}), and immune disease (1.32-fold, hypergeometric test, p < 10\textsuperscript{-2}).
Among high-variance genes in the thyroid rank, we see enrichment for endocrine system disease (1.9-fold, hypergeometric test, p < 10\textsuperscript{-5}),
allergy (1.85-fold, hypergeometric test, p < 10\textsuperscript{-4}), and immune disease (1.45-fold, hypergeometric test, p < 10\textsuperscript{-4}).
These are all rather similar and suggest a stable pattern of high-variance gene expression across these tissues, with enrichment for these three classes of diseases.
The link with immune diseases is expected given the high enrichment for immune-related genes in the high-variance group [@Hagai2018-fu].
As for the low-variance group, we found strong enrichment for genes associated with psychiatric and neurological disorders in the across-study rank and in some tissue-specific ranks (breast, liver, and stomach; ~1.2-fold enrichment, hypergeometric test, p < 0.05, for all cases).
The psychiatric disease link is consistent with previous work [@Mar2011-dr] and is discussed below; however, the enrichment among the low-variance genes is weaker.

# Discussion

Using large publicly available data sets allowed us to probe the landscape of transcriptional variance in several human tissues.
We find a broadly similar pattern of gene expression variance across studies, with high correlations between gene expression variance across most studies, consistent with measurements of expression variance in single cells and in populations of cells for various tissues [@Li2010-qs; @Dong2011-sa; @Alemu2014-jo].
Leveraging this similarity between gene expression variance, we used a multivariate strategy to create a single rank of expression variance, which allowed us to order almost 13k genes according to their expression variance.
Using this rank, we were able to study the general properties associated with high and low variance gene as well as factors driving variation in variance between genes.

Some differences in gene expression variance were driven by technical aspects of gene expression measurement (with data derived from large consortia showing more similar patterns of variance across genes); and by tissue (with studies using the same tissues also showing higher similarities).
This would suggest that careful consideration of sample sizes and experimental design are fundamental to the study of gene expression variance, and the usual small samples of RNA-seq studies might be underpowered for the study of this particular aspect of gene expression.
However, both the effects of study origin and tissue were small, and the largest drivers of differences across studies were idiosyncratic differences related to single data sets, with tissues known to have divergent gene expression patterns (i.e. bone marrow, blood, testis, and platelets) also showing the largest differences in gene expression variance.
Understanding the consequences of these differences in variance for specific tissues is still an open field.
It is clear, however, that differences in variance are informative beyond the differences in mean expression.
Even after we account for differences in mean expression, differences in gene expression variance carry information about tissue origin and function.

Functional analysis using GO enrichment indicated a clear link between function and gene expression variance.
First, genes with high gene expression variance were enriched for biological functions related to response to environmental stimulus, such as immune function and tissue reconstruction.
Likewise, low-variance genes were enriched for basic cell functions, (e.g.  RNA processing, translation, DNA methylation, and cell duplication).
These results are consistent with previous analysis of gene expression variance on a tissue-by-tissue basis [@Alemu2014-jo].
This pattern of enrichment is also observed when we look at enrichment for high- or low-variance genes within the genes associated with each term in the GO hierarchy.
Basic cell function terms are enriched for low-variance genes, and terms involved in response to external stimulus are enriched for high-variance genes.

While indirect, all these patterns point to a selective structuring of gene expression variance.
Stabilizing and purifying selection are consistent: genes expected to be under strong stabilizing selection, those linked with fundamental baseline biological processes, are indeed overrepresented in the least variable genes.
These same genes are also expected to be under strong purifying selection and to show low levels of polymorphisms, which we observe.
Likewise, genes whose function is constrained by myriad interactions with several other genes, those with high connectivity, are less variable.
Furthermore, genes involved with direct interaction with the environment, which must change their pattern of expression depending on external conditions, are expected to be more variable, and again we see a strong enrichment of genes related to interacting with the environment among the most variable.
Given this strong functional link between function and variance, it is not surprising that the gene variance ranking is similar across studies.

One interesting aspect of the GO term analysis shown in @fig:skew_entropy and @fig:go_skewness is that there is no biological process term associated with enrichment for intermediate variance genes: the low-entropy terms have either positive or negative skew, never zero skew.
In other words, there is no annotated biological process for which the associated genes are kept at some intermediary level of variation.
For the GO terms we used, either there is no relation between the transcriptional variance and the biological process, or there is a strong bias toward high or low-variance genes.
This suggests that selective shaping of gene expression has two modes, corresponding with (1) biological processes under strong stabilizing selection (i.e., variance-reducing selection) or (2) biological processes under disruptive selection (i.e., variance-increasing selection).
In short, we find strong support for the idea that there are genes with consistently more (or less) variable expression levels, and that these differences in variance are the result of different patterns of selection.

Following Alemu et al. [@Alemu2014-jo], we observe that epigenetic signatures of gene regulation, such as enhancer histone marks, make up a higher proportion of the surrounding genomic regions of genes that exhibit higher variance in expression.
In contrast, an accumulation of strong promoter elements and overall transcriptional activation is associated with genes with lower expression variance.
These results suggest the presence of distinct modes of regulation for genes with high vs. low variance.
Combined, the differences in the types of genomic regulatory features surrounding the high- and low-variance genes and their distinct functional annotations suggest different mechanisms of regulation of their gene expression variance [@Alemu2014-jo].
This heterogeneity could lead to detectable differences in selection signatures between distal regulatory elements and promoters depending on the gene expression variance.
This heterogeneity in regulation for high and low-variance genes is also notable due to the usual focus on gene expression robustness, in the sense of reducing variation [@Siegal2014-dv; @Payne2015-wn; @Macneil2011-ax; @Denby2012-as]. For example, Siegal and Leu [@Siegal2014-dv] provide several examples of known regulatory mechanisms for reducing gene expression variance, but no examples for the maintenance of high gene expression variance.
We posit that it should be possible to go beyond the usual characterization of strategies of gene expression robustness, in the sense of reducing variation, and to explore mechanisms for the _robustness of plasticity_, that is, the maintenance of high levels of gene expression variation given environmental cues.

Given the broad consistency of gene expression variance in healthy tissues, a natural question is how do these well-regulated levels of variation behave in perturbed or disease conditions.
We find some suggestive links between tissue-specific variance ranks and disease, but these links need to be better explored using more specific methods.
Comparing two HapMap populations, Li et al. [-@Li2010-qs] showed that gene expression variance was similar in both populations and that high-variance genes were enriched for genes related to HIV susceptibility, consistent with our observation of enrichment for immune-related genes among those with more variable expression.
In a case-control experiment, Mar et al. [-@Mar2011-dr] showed that expression variance was related to disease status in Schizophrenia and Parkinson's disease patients, with altered genes being non-randomly distributed across signaling networks.
These authors also find a link between gene network connectivity and expression variance, in agreement with the effect we find using the gene expression variance rank.
The pattern of variance alteration differed across diseases, with Parkinson's patients showing increased expression variance, and Schizophrenia patients showing more constrained patterns of expression.
The authors hypothesize that the reduced variance in Schizophrenia patients reduces the robustness of their gene expression networks, what we refer to as a loss of plasticity.
This suggests several types of shifts in gene expression variation are possible, with different outcomes.
We highlight three distinct possibilities:
First, low-variance genes, under strong stabilizing selection, could become more variable under stress, indicating a reduced capacity for maintaining homeostasis.
Second, high-variance genes, expected to be reactive to changes in the environment, could become less variable, indicating a reduced capacity to respond to external stimuli.
Third, the covariance between different genes could be altered, leading to decoherence between interdependent genes [@Lea2019-pq].
Any one of these changes in expression variance patterns could have physiological consequences, and exploring these differences should be a major part of linking gene expression to cell phenotypes and function (see Hagai et al. [-@Hagai2018-fu] for example).
Genes are also expected to differ in their capacity to maintain an optimal level of gene expression variance [@Macneil2011-ax].
Variation in robustness is linked to gene regulatory networks and epigenetic gene expression regulation [@Payne2015-wn; @Chalancon2012-ul], and therefore, should differ across high- and low-variance genes.
Our results suggest that low- and high-variance genes could use different strategies in order to maintain their optimal levels of variation and that this variability in strategies is the result of different patterns of selection.

\footnotesize

# Methods

## Data sources

We selected 57 human RNA-seq studies with large sample sizes downloaded from the public gene expression repositories recount3 [@Wilks2021-uj] and Expression Atlas [@Papatheodorou2020-dn]. Metadata and details on the included data sets can be found in the supporting information. Because we are interested in population-level variation of gene expression, we exclude single-cell studies. We only used studies for which raw read count data was available, and for which we could parse the metadata for batch effects.
We use the word "studies" to refer to independent data sets, which could have been generated by the same consortium.
For example, the GTEx data are separated by tissue, and we refer to each tissue as a separate study.
We divide our data sets into three categories depending on their origin: GTEx, TCGA, and Miscellaneous.

Table: Data Source Table \label{tab1}

|Study ID                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Citation                                                             |
|:----------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------|
|ADIPOSE_TISSUE (Fat), ADRENAL_GLAND (Adrenal), BLOOD (Blood), BLOOD_VESSEL (Blood_vessel), BONE_MARROW (Marrow), BRAIN (Neuron), HEART (Heart), BREAST (Breast), SALIVARY_GLAND (Salivary), COLON (Colon), LIVER (Liver), NERVE (Neuron), LUNG (Lung), PANCREAS (Pancreas), MUSCLE (Muscle), THYROID (Thyroid), OVARY (Ovary), STOMACH (Stomach), ESOPHAGUS (Esophagus), SPLEEN (Spleen), PROSTATE (Prostate), SKIN (Skin), PITUITARY (Pituitary), TESTIS (Testis) |The GTEx Consortium, 2020 - @GTEx_Consortium2020-xl                  |
||
|LUSC (Lung), STAD (Stomach), COAD (Colon), LUAD (Lung), BRCA (Breast), KIRC (Kidney), KIRP (Kidney), LIHC (Liver), THCA (Thyroid), PRAD (Prostate), UCEC (Uterus)                                                                                                                                                                                                                                                                                                  |The Cancer Genome Atlas Research Network et al., 2013 - @tcga2013-gx |
||
|SRP150552 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Altman et al., 2019 - @Altman2019-wf                                 |
|SRP101294 (Fat)                                                                                                                                                                                                                                                                                                                                                                                                                                                    |Armenise et al., 2017 - @Armenise2017-qn                             |
|SRP057500 (Platelets)                                                                                                                                                                                                                                                                                                                                                                                                                                              |Best et al., 2015 - @Best2015-ha                                     |
|SRP051848 (Immune)                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Breen et al., 2015 - @Breen2015-fm                                   |
|SRP187978 (Liver)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Çalışkan et al., 2019 - @Caliskan2019-od                  |
|E-ENAD-34 (Immune)                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Chen et al., 2016 - @Chen2016-dl                                     |
|SRP059039 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |DeBerg et al., 2018 - @DeBerg2018-ea                                 |
|SRP174638 (Immune)                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Dufort et al., 2019 - @Dufort2019-nw                                 |
|E-GEOD-57945 (Colon)                                                                                                                                                                                                                                                                                                                                                                                                                                               |Haberman et al., 2014 - @Haberman2014-fw                             |
|SRP162654 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Harrison et al., 2019 - @Harrison2019-ld                             |
|SRP095272 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Jadhav et al., 2019 - @Jadhav2019-fo                                 |
|SRP102999 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Kuan et al., 2017 - @Kuan2017-hp                                     |
|SRP145493 (Immune)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Kuan et al., 2019 - @Kuan2019-ju                                     |
|E-GEUV-1 (Immune)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Lappalainen et al., 2013 - @Lappalainen2013-io                       |
|SRP035988 (Skin)                                                                                                                                                                                                                                                                                                                                                                                                                                                   |Li et al., 2014 - @Li2014-lw                                         |
|SRP192714 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Michlmayr et al., 2020 - @Michlmayr2020-dp                           |
|ERP115010 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Roe et al., 2020 - @Roe2020-di                                       |
|E-ENAD-33 (Neuron)                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Schwartzentruber et al., 2018 - @Schwartzentruber2018-fa             |
|SRP181886 (Neuron)                                                                                                                                                                                                                                                                                                                                                                                                                                                 |Srinivasan et al., 2020 - @Srinivasan2020-es                         |
|SRP098758 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Suliman et al., 2018 - @Suliman2018-ga                               |
|SRP032775 (Blood)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Tran et al., 2016 - @Tran2016-ro                                     |
|SRP069212 (Liver)                                                                                                                                                                                                                                                                                                                                                                                                                                                  |Yang et al., 2017 - @Yang2017-lu                                     |


__Processing pipeline__: We use a standardized pipeline to measure gene expression variance while removing extraneous sources of variation.
Data from case-control studies were filtered to keep only control samples.
Technical replicates were summed.
For each study, we filtered genes that did not achieve a minimum of 1 count per million (cpm) reads in all samples and a mean of 5 cpm reads across samples.
To account for the mean-variance relation in RNA-seq count data, we applied a variance stabilizing transformation implemented in DESeq2 [@Love2014-mp] to the genes passing the read-count filters.
This mean-variance correction was verified by plotting mean-variance relations before and after correction, and these plots can be seen in the supporting information.
Various technical covariates (like experimental batch, sex, etc.) were manually curated from the metadata for each study and accounted for using a linear fixed-effects model. 
A list of covariates used for each study is available in the supporting information.
Outlier individuals in the residual distribution were removed using a robust Principal Component Analysis (PCA) approach of automatic outlier detection described in @Chen2020-fy. 
This procedure first estimates a robust Principal Components for each study and then measures the Mahalanobis distance between each sample and the robust mean. 
Samples that are above the 0.99 percentile in Mahalanobis distance to the mean are marked as outliers and removed.
We verify that the batch effect correction and outlier removal are reasonable by using PCA scatter plots after each step of the pipeline to check the result for residual problems like groupings or other artifacts.
These PCA plots before and after batch correction and outlier removal are also included in the supporting information.
Gene expression standard deviation is measured as the residual standard deviation after fixed effect correction and outlier removal.
We choose standard deviation as a measure of variation to have a statistic on a linear scale, and we do not use coefficient of variation because we have already corrected for mean differences and for the mean-variance relation inherent to RNA-seq count data [@De_Jong2019]. 
The full annotated pipeline is available at [the github repository ayroles-lab/ExpressionVariance](https://github.com/ayroles-lab/ExpressionVariance).

## Gene expression variance across-study correlation

We assessed the similarity in gene expression variance across studies by using a between-study Spearman correlation matrix of the measured SDs.
Only genes present in all studies were used to calculate the Spearman correlation matrix, ~4200 genes in total.
Using Spearman correlations avoids problems related to overall scaling or coverage differences, and allows us to assess if the same genes are usually more- or less-variable across studies.
To investigate the factors involved in determining correlations between studies, we used a Bayesian varying effects model to investigate the effect of study origin and tissue on the correlations across studies.
This model is designed to take the non-independent nature of a set of correlations into account when modeling the correlation between gene expression SDs.
This is accomplished by adding a per-study random effect, see [@Dias2021-hb] for details.
The Fisher z-transformed Spearman correlations across studies ($z(\rho_{ij})$) are modeled as:

$$\begin{aligned}
z(\rho_{ij}) &\sim N(\mu_{ij}, \sigma) \\
\mu_{ij} &= \mu_0 + \alpha_i + \alpha_j + \beta X \\
\alpha_i &\sim N(0, \sigma_{\alpha})
\end{aligned}$$

The $\alpha$ terms account for the non-independence between the pairs of correlations and estimate the idiosyncratic contribution of each study to all the correlations it is involved in. The fixed effects encoded in the design matrix $X$ measure the effects of tissue congruence and study-origin congruence. All fixed effect parameters ($\beta$) and per-study parameters ($\alpha$) receive weakly informative normal priors with a standard deviation of one quarter. For the overall variance ($\sigma$) we use a unit exponential prior, and for the intercept ($\mu_0$) a unit normal prior. This model was fit in Stan [@carpenter2017stan] via the rethinking R package [@mcelreath2020statistical], using eight chains, with 4000 warm-up iterations and 2000 sampling iterations. Convergence was assessed using R-hat diagnostics [@Gelman2013-ae], and we observed no warnings or divergent transitions.

__Gene expression SD rank:__ Given that most of the variation in the Spearman correlation across studies is explained by a single principal component (PC1 accounts for 62% of the variation in the across-study Spearman correlation matrix, while PC2 accounts for only 5%; see SI fig. 3), we use the ranked projections of gene expression SDs in this principal component (PC1) to create an across-study rank of gene variation.
The higher the rank, the higher the expression SD of a given gene.
Genes that were expressed in at least 50% of the studies were included in the rank.
In order to project a particular gene onto the PC1 of the between-study correlation matrix, we impute missing values using a PCA-based imputation [@Husson2019-sl].
The imputation procedure has minimal effect on the ranking, and imputing missing SD ranks at the beginning or at the end of the ranks produces similar results.
We also create a tissue-specific variance ranking, using the same ranking procedure but joining studies done in the same tissue type.
For this tissue-level ranking, we only use genes that are expressed in all studies of a given tissue.
For tissues that are represented by a single study, we use the SD ranking for that study as the tissue rank.
We further investigate the tissue-level expression variance ranks as they relate to genomic regulation.

__Gene expression mean rank:__ We also use the same strategy to create a mean gene expression rank, repeating the process but using mean expression instead of standard deviation. All ranks are available in the supporting information.

## Gene level statistics

__Genetic variation__: Genetic variation measures were obtained from the PopHuman project, which provides a comprehensive set of genomic information for human populations derived from the 1000 Genomes Project.
Gene-level metrics were used when available.
If only window-based metrics are available, we assembled gene-level information from 10 kb window tracks where each window that overlaps with a given gene was assigned to the gene and the mean metric value is reported.
In parallel, we use the PopHumanScan data set, which expands PopHuman by compiling and annotating regions under selection.
Similarly, we used gene-level information when possible, and for tracks with only window-based metrics, gene-level information was assembled from the 10 kb windows using the same assignment method described above.
Nucleotide diversity ($\pi$), the average pairwise number of differences per site among the chromosomes in a population [@Nei1979-hg], provides insight into the genetic diversity within a population, in this case, the CEU population within 1000 genomes.
The nucleotide diversity can also be used as an estimator of the central population genetic parameter, normally given as $\theta$.

__Gene connectivity__: We calculated the average weighted connectivity for all genes by creating a fully connected gene-by-gene graph in which each edge is weighted by the Spearman correlation between gene expression levels.
We then trimmed this graph by keeping only edges for which the Spearman correlation is significant at a BH false discovery rate of 1%.
In this trimmed network, we then took the average of the Spearman correlation of all remaining edges for each gene.
So, for each study, we have a measure of the average correlation of each gene with every other gene.
The average connectivity for each gene is the average across all studies in which that gene is expressed.

__Cross-tissue vs. tissue-level chromatin states__: We use the universal [@vu2022universal] and tissue-specific [@Ernst2015-zk] ChromHMM [@ernst2012chromhmm] chromatin states to compare the non-overlapping genome segmentation to cross-tissue and tissue-level gene expression variance metrics. We use the proportion of the gene regions (gene +/- 10 kb) made up of each of the chromHMM chromatin states.

__Correlations:__ We use the ppcor R package v1.1 [@kim2015ppcor] to run the pairwise partial Spearman correlations between gene-level statistics and the gene expression variance rank while controlling for the mean expression rank. P-values are corrected using the Benjamini-Hochberg procedure and comparisons with an adjusted p<0.05 are considered significant.

## Gene function assessment

__GO term enrichment__: All gene ontology (GO) analyses were done using the clusterProfiler R package v4.2.2 [@Wu2021-db] and the Org.Hs.eg.db database package v3.14.0 [@godb]. GO and all further enrichment analysis used the hypergeometric test to assess the significance of the enrichment.

__Secreted genes__: We use The Protein Atlas [@uhlen2015tissue] to extract information on which proteins are secreted [@uhlen2019human] and test for enrichment of genes with secreted products in the genes within the highest and lowest 5% of gene expression variance rank.

__Housekeeping genes__: Human housekeeping genes were identified as genes that are expressed with low variance in all 52 human cell and tissue types, assessed in over 10,000 samples [@Hounkpe2020-yq]. We test for enrichment of housekeeping genes in the genes within the highest and lowest 5% of gene expression variance rank.

__Immediate early genes (IEGs):__ Human IEGs were curated from the literature in @Arner2015-be as genes that respond to experimental stimulation through up-regulation within the first 60 minutes of the experiment. We use the hypergeometric test to assess the significance of the enrichment. Immediate early genes (IEGs): Human IEGs were curated from the literature in  @Arner2015-be as genes that respond to experimental stimulation through up-regulation within the first 60 minutes of the experiment.

__Probability of being loss-of-function intolerant (pLI)__: Genes that are likely haploinsufficient (i.e., intolerant of heterozygous loss-of-function variants) were detected as those with fewer than expected protein-truncating variants (PTVs) in ExAC [@Lek2016-xw]. We use genes with a pLI > 0.9 to test for the enrichment of loss-of-function intolerant genes in the genes exhibiting the highest and lowest 5% gene expression variance estimates.

__Disease annotations__: We use the gene annotations for involvement with diseases provided by the supporting information Table S2 from Zhang et al. [-@Zhang2020-cl].

# Code availability

All code for reproducing all analyses and figures, along with a walk-through, is available at [github.com/ayroles-lab/ExpressionVariance](https://github.com/ayroles-lab/ExpressionVariance).

# Supporting information

Supporting information is available at [github.com/diogro/expVarManuscript](https://github.com/diogro/expVarManuscript).

1. SI figure 1 - Across-study and tissue-specific gene expression variance and mean correlations
with non-overlapping chromatin states through ChromHMM.

1. SI figure 2 - Proportion of gene regions made up of ChromHMM chromatin states for genes in
the top and bottom 5% of the across-study mean rank metric.

1. SI table 1 - Variance and mean rank metrics and the corresponding ChromHMM annotations
used.

1. SI data - Study metadata -  Metadata file describing the data used in the study as well as some intermediate processing information.

1. SI data - Gene ranks - Gene expression mean and variance ranks, across-study and tissue-specific.

1. SI data - GO enrichment - Combined table describing gene ontology enrichment in the top 5% and bottom 5% of genes as ranked by variance.


# References

<!-- %% \printbibliography -->
