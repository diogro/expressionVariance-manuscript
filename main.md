---
title: "Characterizing the landscape of gene expression variance in humans"
short_title: "The landscape of transcriptional variance"
date: Nov 15th, 2022
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
      correspondence: "yes"
      orcid: 0000-0002-7603-0092
      equal_contributor: "yes"
  - Kristina M. Garske:
      institute: lsi
      orcid: 0000-0002-7228-8125
  - Luisa F. Pallares:
      institute: lsi, fml
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
      name: "Current address: Friedrich Miescher Laboratory, Max Planck Society"
# classoption: twocolumn
output: pdf_document
geometry:
- top=20mm
- left=25mm
- right=25mm
- bottom=20mm
header-includes:
- \usepackage[left]{lineno}
- \modulolinenumbers[5]
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
---
<!--
# Author Summary

Required for PLOS Genetics... -->
# Abstract
Gene expression variance has been linked to organismal function and fitness but remains a commonly neglected aspect of gene expression research. As a result, we lack a comprehensive understanding of the patterns of variance across genes, and how this variance is linked to context-specific gene regulation and gene function. Here, we use 57 large publicly available RNA-seq data sets to investigate the landscape of gene expression variance. These studies cover a wide range of tissues and allowed us to assess if there are consistently more or less variable genes across tissues and data sets and what mechanisms drive these patterns. We show that gene expression variance is broadly similar across tissues and studies, indicating that the pattern of transcriptional variance is consistent. We use this similarity to create both global and within-tissue rankings of variation, which we use to show that function, sequence variation, and gene regulatory signatures contribute to gene expression variance. Low-variance genes are associated with fundamental cell processes and have lower levels of genetic polymorphisms, have higher gene-gene connectivity, and tend to be associated with chromatin states associated with transcription. In contrast, high-variance genes are enriched for genes involved in immune response, environmentally responsive genes, immediate early genes, and are associated with higher levels of polymorphisms. These results show that the pattern of transcriptional variance is not noise. Instead, it is a consistent gene trait that seems to be functionally constrained in human populations. Furthermore, this commonly neglected aspect of molecular phenotypic variation harbors important information to understand complex traits and disease.

# Introduction
\linenumbers
Molecular phenotypes such as gene expression are powerful tools for understanding physiology, disease, and evolutionary adaptations.
In this context, average trait values are usually the focus of investigation, while variation around the average is often considered a nuisance and treated as noise [@De_Jong2019-po].
However, gene expression variance can be directly involved in determining fitness [@Fraser2004-sv; @Wang2011-ts], can drive phenotypic variation [@Hansen2021-zo], and the genetic architecture of variance itself can evolve [@Bruijning2020-bf].
This suggests that studying gene expression variance as a bona fide trait, its genetic architecture, and the evolutionary mechanisms shaping and maintaining gene-specific patterns of variance has the potential to further our understanding of complex traits and disease [@Li2010-qs; @Mar2011-dr; @Hagai2018-fu].

Variability is ubiquitous in nature and is, alongside robustness, a fundamental feature of most complex systems.
But, at the same time, the degree of variability seems to differ between genes [@De_Jong2019-po] suggesting that it might be associated with biological function and therefore be constrained by selection.
From a mechanistic perspective, several competing forces act to shape transcriptional variance [@Houle1998-mj; @Bruijning2020-bf], and the outcome of the interaction between these processes is still poorly understood [@Hansen2011-es].
For example, we expect the influx of new mutations to increase the variance, while the selective removal of these polymorphisms, via purifying selection or selective sweeps, to decrease it [@Schmutzer2020-cd; @Pettersson2012-om].
From a quantitative trait perspective, stabilizing selection should decrease variation around an optimal value, and directional selection can lead to a transient increase in variance while selected alleles sweep to fixation, followed by a reduction in variance as these alleles become fixed.
Pleiotropic effects are also important, as they allow selection on one trait to influence the variance of other traits [@Wagner1997-hw; @Pavlicev2011-xm].
Both indirect effects of directional selection on variance open the possibility that the main driver of gene expression variance is not direct selection on variance but indirect effects due to selection on trait means [@Hansen2011-es].
How the interaction of these processes shape gene expression variance is an open question.
However, some general predictions can be made.
If a homogeneous pattern of stabilizing selection is the main driver of gene expression variance, we would expect transcriptional variance to be consistent regardless of the population, tissue, or environmental context.
If idiosyncratic selection patterns and context-specific environmental interactions are more important, we could observe large differences in gene
expression variance.

A key difficulty in addressing these questions is that the constraints on gene expression variance might also be dependent on the gene tissue specificity.
Mean expression is known to differ across tissues [@GTEx2017-xb], however, to what extent differential expression (i.e., differences in mean expression level) translate into differences in expression variance is not clear.
Higher mean expression could lead to higher variance, but other processes can also affect transcriptional variance.
For example, if a gene is expressed in more than one tissue and variance regulation is independent across tissues, stabilizing selection on gene expression could be more intense depending on the role of that gene in a particular tissue, causing a local reduction in variation that leads to differences in variance across tissues (@fig:variance_control A).
These across-tissue differences would not necessarily follow mean expression.
Alternatively, expression variation across tissues could be tightly coupled and, in this example, selection in one tissue would lead to a reduction in variance across tissues, resulting in a consistent pattern of variation (@fig:variance_control B).
While we lack a clear picture of how tissue-specific gene expression variation is regulated, Alemu et al. [@Alemu2014-jo] used microarray data from several human tissues to show that epigenetic markers were linked to gene expression variation and that these markers were variable across tissues and between high- and low-variance genes.

![Example of how differences in the regulation of transcriptional variance can drive changes in the correlations between gene expression variance measures. In (A), independent regulation causes the reduction in variation to be restricted to context 1 (context here can refer to different tissues, environments, populations, studies, etc.). On the right side of panel A, independent regulation results in low correlation across contexts. In (B), a shared regulatory architecture maintains consistent variance across both conditions, leading to high similarity in transcriptional variance across contexts. In (C), we see how the similarity seen in panel B can be leveraged to create an across-context rank of gene expression variance. When transcriptional variance ranks are highly correlated, the rank of the projection onto the first principal component (PC1) allows us to summarize the across-context pattern of transcriptional variance.](figures/variance_control_horizontal.png){#fig:variance_control}

To explore the landscape of gene expression variance and the association between transcriptional variance and biological function, we use 57 publicly available human gene expression data sets spanning a wide range of experimental contexts and tissues.
By comparing the gene expression variance measured across such heterogeneous data sets, we show that the degree of expression variance is indeed consistent across studies and tissues.
We use the observed similarities to create an across-study gene expression variance ranking, which orders genes from least variable to most variable.
We then integrate various genomic-level functional annotations as well as sequence variation to probe the drivers of this variance ranking.
Finally, we explore the link between gene expression variance and biological function by leveraging gene ontology and other gene annotations.

# Results


![Overview of the distribution of transcriptional variance across studies. (A) Heatmap showing the correlation in transcriptional variance across studies (as the Spearman correlation of gene expression standard deviations). Pairs of studies with more similar patterns of gene expression variance have higher correlations. Studies are shown in the same order as in SI fig. 1, panel A. (B) Distribution of the pairwise Spearman correlations between studies shown in the previous panel. (C) PCoA using the distance between studies derived from the pairwise correlations. There is no clear structuring of the studies with respect to their source, which is indicated by the colors. (D) Density plot of standard deviations after z-normalization. The inset plot shows the distribution of mean-centered standard deviations grouped by study without normalization. The corresponding rug plots show the location of the highest-ranking gene in standard deviation rank ([_hbb_]{.smallcaps}) (right, blue) and lowest ([_wdr33_]{.smallcaps}) (left, red).](figures/fig1.png){#fig:sd_corr}

## Data sets

We use 57 publicly available human gene expression RNA-seq data sets which were derived from the publications listed in table \ref{tab1} of the [Methods](#Methods) section, and a complete metadata table for each study is available in the supporting information (SI data 1).
We only use data sets that fulfilled the following conditions: samples came from bulk RNA-seq (and no single cell approaches), data sets were associated with a publication, sample-level metadata was available, and the post-filtering sample size was greater than 10.
These data sets span 13 different tissue types and the post-filtering mean sample size we used for each data set was 390, with a median of 251, and ranged from 12 to 2931 samples.
Several data sets were derived from two large consortia: GTEx [@GTEx2017-xb] and TCGA [@tcga2013-gx], and we note the origin of the data sets in the figures where appropriate.
We refer to data sets and studies interchangeably, and so each tissue in GTEx is referred to as a different study.
The final list of genes used from each study can be found in SI data 2.

## Gene expression variance

For each study, transcriptional variance per gene was measured as the standard deviation (SD) of the distribution of gene expression values for all individuals in a particular study.
Mean and variance are known to be correlated in RNA-seq data, both due to the nature of count data and the expectation that more highly expressed genes should have more variation.
As our focus here is on variance, we control for both of these expected drivers of transcriptional variation.
To achieve this, SD was calculated using a unified pipeline that normalized the mean-variance relation in read-count data, controlled for batch effects, and removed outliers (see [Methods](#Methods) for details, and the calculated values for means and standard deviations are available in SI data 3).
The observed range of gene expression SDs across genes is variable but can be normalized so that the distributions are comparable (@fig:sd_corr D).
This comparison reveals differences in the range of gene expression SDs that can be due to any number of methodological or biological differences between the data sets.
We avoid having to deal with these global differences in the range of variation by using only the ranking of the genes according to their gene expression SD in each study.
Therefore, patterns of transcriptional variance were compared across studies using Spearman correlations (ρ~s~) between gene expression SDs.
This comparison reveals a broadly similar rank of gene expression variance as the correlations across studies are mostly positive and high (75% of correlations are between 0.45 and 0.9, @fig:sd_corr A and B), indicating that genes that are most variable in one study tend to be most variable in all studies.
A principal coordinate analysis [@Gower1966-dk] using |1 -- ρ~s~| as a between-study distance measure does not show clearly delineated groups, but GTEx and TCGA studies are clustered among themselves and close together (@fig:sd_corr C).
This clustering indicates some effect of study source on the similarity between gene expression SD across studies, which we explore in detail below.

To characterize what factors may explain differences in across-study similarity, we directly modeled the across-study correlations using a mixed-effect linear model designed to account for the non-independence in pairwise correlation data [@Dias2021-wk; @Dias2022-au].
In this model (see [Methods](#Methods)), we use a random effect for individual study ID, a fixed effect for pairwise tissue congruence (whether a comparison is within the same tissue or different tissue), and a fixed effect for pairwise study source (which pair of sources among GTEx, TCGA, and miscellaneous is involved in a comparison) as predictors of the correlations (see [Methods](#Methods)).
This model (SI fig. 1) shows that comparisons of studies within GTEx and TCGA have on average higher values of $\rho_s$, but also that comparing studies across GTEx and TCGA also shows a mild increase in the average correlation (SI fig. 1 C).
Correlations that do not involve studies from TCGA and GTEx (marked as "Misc.") are on average lower (SI fig. 1 C).
While we do not have a clear explanation for this pattern, since TCGA and GTEx are independent, this mild effect on the similarities could be due to the level of standardization of the data coming from these two large consortia.
Tissue type also affects the degree of similarity in transcriptional variance, with studies using the same tissue being, on average, more similar (SI fig. 1 B).
However, all these pairwise effects are mild, and the largest effects on the correlations are those associated with individual studies, in particular some specific tissues, i.e., comparisons involving [bone marrow]{.smallcaps} (from GTEx) and study [srp057500]{.smallcaps} (which used platelets) are on average lower (SI fig. 1 A).
The only negative correlation we observe is between these two studies, which also appear further away in the PCoA plot in @fig:sd_corr C.

## Transcriptional variance rank

The strong correlations between transcriptional variance across studies suggest that variance rank is indeed a property of genes that can be robustly estimated.
To estimate this gene-level rank, we devised an across-study approach that allowed us to rank individual genes according to their degree of transcriptional variance by averaging the ordering across all studies.
We do this by calculating the score of each gene on the first principal component of the across-study Spearman correlation matrix shown in @fig:sd_corr A. This procedure is illustrated in @fig:variance_control C.
Ordering genes using these scores generate a ranked list of genes, with the most variable genes having the highest rank.
The position in the SD distributions shown in @fig:sd_corr D of the most and least variable genes in this rank illustrates how the extremes of the rank are indeed some of the least and most variable genes across all studies.
In addition, to be able to account for any residual effect of mean expression on the variance we also created a similar across-study rank for mean expression.
To explore tissue-specific divers or transcriptional variation, we also create a set of tissue-specific SD ranks.
To that end, we used the same procedure outlined above but using only studies that were performed on  the same tissue.
Both tissue-specific and across-study ranks are available in the Supporting Information (SI data 4).

<!-- ![Modeling the correlations between transcriptional variance across studies. The panels show coefficient estimates from a linear model using the among studies Spearman correlations between gene expression SDs as the response variable. These correlations are shown in @fig:sd_corr, panels A and B. In the linear model (see Methods for model equation), correlations are Fisher z-transformed. Study source and tissue are added as fixed effects. Coefficient estimates are shown with 50% and 95% credibility intervals. (A) The per-study random effect which accounts for the non-independence between the pairwise correlation values and estimates the characteristic contribution of each study to these correlations. For example, the lowest estimate among these parameters, which corresponds to the study [bone marrow]{.smallcaps} (from GTEx), indicates that correlations involving this study tend to be lower than the others. (B) and (C) Fixed effect estimates for the effects of tissue congruence and study-source effect. In (B) we see that correlations among studies that use the same tissue are slightly higher; and (C) correlations involving studies in the "Misc." category (non-GTEx and non-TCGA) tend to be lower, while comparisons involving GTEx and TCGA are higher.](figures/correlationModeling.png){#fig:corr_model width=550px} -->


## Biological function explains gene-level transcriptional variance

As a first step toward explaining the factors that drive variation in variability between transcripts, we focused on the top 5% most variable and the bottom 5% least variable genes in the across-study ranking (560 genes in each group).
A Gene Ontology (GO) enrichment analysis shows 59 enriched terms in the low-variance genes, and 738 enriched terms in the high-variance genes (using a hypergeometric test and a conservative Benjamini-Hochberg (BH) adjusted p-value threshold of 10^-3^; see supporting information SI data 5 for a complete listing).

Among the most variable genes, we observe enrichment for biological processes such as immune function, response to stimulus, maintenance of homeostasis, and tissue morphogenesis (SI fig. 2 A).
Furthermore, we see a 7.7-fold enrichment for genes that encode secreted proteins in the most variable genes, relative to all other genes (hypergeometric test, p < 10^-3^).
Given that the GO enrichment suggests high-variance genes are involved in responding to stimulus, we compare them to a recently generated catalog of environmentally responsive genes. This catalog was generated using 11 environmental exposures in 544 immortalized Lymphoblastoid Cell Lines (LCL) from the 1000 Genomes Project [@Lea2021-hf]. We find a strong enrichment of high-variance genes among environmentally responsive genes across 7 out of the 10 environmental exposures we used (hypergeometric test, p < 10^-3^). We do not find any enrichment among the least variable genes (SI table 2).

Among the least variable genes, we see enrichment for housekeeping functions such as mRNA processing, cell cycle regulation, methylation, histone modification, translation, transcription, and DNA repair (SI fig. 2 B); accordingly, we also find a 2.0-fold enrichment in previously characterized human housekeeping genes [@Hounkpe2020-yq] (hypergeometric test, p < 10^-3^).
The genes exhibiting the lowest variance are also enriched for genes that have been previously shown to have a high probability of being loss-of-function intolerant (pLI) [@lek2016analysis] (1.2-fold enrichment, hypergeometric test, p < 10^-3^).
Genes with a high pLI have been shown to be important in housekeeping functions and have higher mean expression values across a broad set of tissues and cell types [@lek2016analysis]. The observation that genes with low variance are enriched for both housekeeping genes and genes with high pLI is consistent with this previous report; and we further see that the mean expression of genes positively correlates with pLI (partial Spearman correlation $\rho_s$ = 0.32, p < 10^-3^), showing the opposite relationship between variance and mean expression when considering pLI.

![Relationship between skew and entropy of rank decile distributions for each GO term. High entropy terms, to the right of the plot, are associated with a more egalitarian proportion of genes in each of the SD rank deciles. The terms on the left of the plot are associated with more genes in some particular decile. The skewness in the y-axis measures if the high- or low-variance deciles are more represented for a particular term. Terms on the positive side of the y-axis are associated with low-variance genes, and terms on the negative side of the y-axis are associated with high-variance genes. The GO terms are filtered for gene counts greater than 100, as in @fig:go_skewness. Some of the top high- and low-skewness terms are labeled for illustration.](figures/GOterm_entropy_by_skewness.png){#fig:skew_entropy width=500px}

In the previous analysis, we explored the relationship between transcriptional variance and function by starting from the extremes of the variance distribution and searching for GO enrichment among these high- and low-variance genes.
We also approach the problem from the opposite direction, starting from the genes associated with each GO term and searching for enrichment for high- or low-variance genes among them.
To this end, we gathered all biological process GO terms in level 3 (i.e., terms that are at a distance of 3 from the top of the GO hierarchy).
Using level-3 terms gives us a good balance between number of terms and genes per term.
We separated the genes associated with at least one of these level-3 terms into expression variance deciles, with the first decile having the lowest variance.
We then counted how many genes in each decile have been associated with each term.
If variance rank is not associated with the GO annotations, terms should have an equal proportion of genes in each decile.
We measured how far from this uniform allocation each term is by measuring the Shannon entropy of the proportion of genes in each decile.
Higher entropy is associated with a more uniform distribution of genes across deciles.
GO terms with low entropy indicate some deciles are over-represented in the genes associated with that term.
We also measured skewness for each term, which should be zero if no decile is over-represented, negative if high-variance terms are over-represented, and positive if low-variance deciles are over-represented.
The relation between skewness and entropy for each GO term can be seen in @fig:skew_entropy.
Positive-skew low-entropy terms, those enriched with low-variance genes, are associated with housekeeping functions, like RNA localization, translation initiation, methylation, and chromosome segregation (@fig:go_skewness A).
Likewise, terms with negative skew and low entropy, enriched for high-variance genes, are related to immune response, tissue morphogenesis, chemotaxis---all dynamic biological functions related to interacting with the environment (@fig:go_skewness B).

Both GO analyses suggest a strong association between biological function and the degree of transcriptional variance.
Genes associated with baseline fundamental functions, expected to be under strong stabilizing selection, are also low-variance; high-variance genes are associated with responding to external stimuli (i.e., tissue reorganization and immune response).


<!-- ![Gene set enrichment analyses testing for over-representation of gene ontology categories in the upper and lower 5\% quantiles of the gene variance rank. (A) High-variance genes are enriched for terms related to immune function, response to wounding, blood vessel morphogenesis, and inflammatory response. In contrast, (B) low-variance genes are associated with translation, control of methylation, RNA processing, chromosome separation, and other cell housekeeping functions. All displayed terms are significant with a 5% FDR corrected p-value below 10^-3^.](figures/local_go_lowerUpper.png){#fig:go_tails width=500px} -->




![Distributions of decile ranks of level-3 GO terms. Each plot shows the count of genes in each decile of the rank. Only GO terms that are associated with at least 100 genes are used. We sort these terms by the skewness of the distribution. The top panel (A) shows the 5 most positively skewed terms, and the bottom panel (B) shows the 5 most negatively skewed terms.](figures/GOterm_decile_barplot.png){#fig:go_skewness width=500px}



## Evolutionary forces at play in shaping transcriptional variance

We use three gene-level summary statistics, nucleotide diversity (π), gene expression connectivity, and the rate of adaptive substitutions (α), as a proxy to assess whether selection might be involved in shaping gene expression variance.
For all the correlations in this section, we use partial Spearman correlations that include the mean gene expression rank as a covariate, which accounts for any residual mean-variance correlation.
Nucleotide diversity in the gene region is used as a proxy for the impact of cis-regulatory genetic variation on transcriptional variance.
As expected, low-variance genes tend to have lower levels of polymorphisms (partial Spearman correlation, ρ~s~ = 0.184, p < 10^-10^).
Gene-gene connectivity, a proxy for gene regulatory interactions and selective constraints [@Mahler2017-bb], is, in turn, negatively correlated with the expression variance (partial Spearman correlation, ρ~s~ = -0.024, p < 10^-2^), supporting the expectation that highly connected genes are more constrained in their variation.
Finally, we also find that low-variance genes tend to have fewer substitutions by comparing the across-study rank with $\alpha$ (partial Spearman correlation, ρ~s~ = -0.044, p < 10^-2^), in line with the expectation that genes under stronger selection should be less variable.
Despite all associations being significant and in the expected direction, their effect sizes are very small, suggesting a weak link between these broad measures and transcriptional variance.

## Specific gene regulatory signatures are associated with transcriptional variance

To assess how local epigenetic features relate to gene expression variance we calculate the proportion of the gene (±10 kb) that corresponds to epigenetic signatures of gene regulation defined through ChromHMM [@ernst2012chromhmm] chromatin states. Chromatin states associated with distal (i.e., non-promoter) gene regulation are positively correlated with the across-study variance rank, regardless of whether the regulatory effect on gene expression is positive or negative (@fig:lineplot; see across-study correlations in SI fig. 3A). For example, both the proportion of gene regions made up of enhancers and repressed genomic states are positively correlated with gene expression variance (BH adjusted Spearman correlation, p < 0.05). In contrast, histone modifications associated with active promoters, as well as transcribed states, are inversely correlated with gene expression variance (SI fig. 3A), whereas they are positively correlated with the mean rank (SI fig. 3B). Taken together, these results are compatible with gene expression variance being regulated through distal (i.e., non-promoter) gene regulatory mechanisms, rather than the overall active transcriptional state of a gene region, as is the case with mean gene expression.

![Proportion of gene regions made up of ChromHMM chromatin states for low- and high-variance genes. The line plot contrasts the proportion of gene regions made up of the indicated chromatin states for genes in the top and bottom 5% of the across-study variance rank metric. Ends denote the median proportion of gene regions made up of the chromatin state, and error bars represent the standard error of the mean. States colored black are not significant, all others exhibit significant differences between low- and high-variance genes (BH adjusted Wilcoxon signed-rank test, p < 0.05). Het indicates heterochromatin; TSS, transcription start sites; znf, zinc finger genes. The mean rank version of this analysis is shown in SI fig. 4.](figures/top_bottom_5pVarrank_fxnlGen_lineplot_10kb.png){#fig:lineplot width=350px}

Given that ChromHMM chromatin states are available for specific tissues, we asked whether the regulatory signatures associated with the across-study variance rank are recapitulated at the tissue level.
Many of the across-study correlations are recapitulated at the tissue-specific level (with two exceptions noted below), including a strong and highly consistent positive correlation between the proportion of gene regions made up of enhancer states and that gene’s expression variance, and an inverse relationship between gene expression variance and histone marks associated with gene transcription (SI fig. 3A).
Two blood associations stand out as being different from the consistent effects across the other tissue-level and across-study associations.
First, the weak (i.e., histone marks associated with both activating and repressive functions) promoter state is positively correlated with transcriptional variance in all comparisons except blood.
Second, the consistent inverse correlation of gene expression variance with weak transcription is reversed in blood, such that there is a positive correlation between histone marks associated with weak transcription and blood gene expression variance (SI fig. 3A).
Taken together, these results suggest that, rather than genes with a bivalent promoter state (i.e., poised genes) exhibiting more expression variance, blood high-variance genes are more likely already expressed at basal levels (i.e., weakly transcribed), as discussed previously [@Rogatsky2014-bi].

Immediate early genes (IEGs) respond quickly to external signals without requiring _de novo_ protein synthesis, and a bivalent state has been reported to be associated with IEG promoters [reviewed in @Bahrami2016-dx]. Given our results that genes with high expression variance are enriched for cellular signaling and response mechanisms (SI fig. 2 A), and bivalent promoter states are correlated with the gene expression variance rank (SI fig. 3A), we hypothesized that IEGs would be enriched within genes in the top expression variance ranks. This was the case for all tissue-level gene expression variance ranks (enrichment ratios range from 3.3-8.8, Bonferroni-adjusted hypergeometric test, p < 0.05), except for blood (enrichment ratio = 1.2, hypergeometric test, p = 0.3). Thus, once again blood stands out when attempting to understand genomic regulatory drivers of expression variance. In all, while high-variance genes are generally shared across tissues and enriched for immune and environmental signaling pathways, it seems that the gene regulatory mechanisms governing their expression are distinct between immune cell types and other tissues studied here.


## Linking expression variance and disease

To explore the link between transcription variance and genes known to be associated with human diseases, we used a data set designed to provide causal relationships between gene expressions and complex traits/diseases (based on a probabilistic transcriptome-wide association study (PTWAS) [@Zhang2020-cl]).
Using the list of significant gene-disease pairs at 5% FDR provided by Zhang et al. [-@Zhang2020-cl], we performed a hypergeometric enrichment test for the top 5% high- and low-variance genes in our across-study rank and in all tissue-specific gene variance ranks.
We use both across-study and tissue-specific ranks because some genes only appear in the tissue-specific rank due to their limited tissue-specific gene expression.
In the high-variance group, we find no enrichment in the across-study rank, but we do find enrichment of genes annotated for allergy, immune disease, and endocrine system disease among the high-variance genes in several tissue-specific variance ranks.
For example, among high-variance genes in the colon, we see enrichment for endocrine system disease (1.77-fold, hypergeometric test, p < 10^-4^).
Among high-variance genes in the immune cells, we see enrichment for endocrine system disease (1.67-fold, hypergeometric test, p < 10^-3^),
allergy (1.7-fold, hypergeometric test, p < 10^-3^), and immune disease (1.32-fold, hypergeometric test, p < 10^-2^).
Among high-variance genes in the thyroid, we see enrichment for endocrine system disease (1.9-fold, hypergeometric test, p < 10^-5^),
allergy (1.85-fold, hypergeometric test, p < 10^-4^), and immune disease (1.45-fold, hypergeometric test, p < 10^-4^).
These are all rather similar and suggest a stable pattern of high-variance gene expression across these tissues, with enrichment for these three classes of diseases.
The link with immune diseases is expected given the high enrichment for immune-related genes in the high-variance group [@Hagai2018-fu].
As for the low-variance group, we found strong enrichment for genes associated with psychiatric and neurological disorders in the across-study rank and in some tissue-specific ranks (breast, liver, and stomach; ~1.2-fold enrichment, hypergeometric test, p < 0.05, for all cases).
The psychiatric disease link is consistent with previous work [@Mar2011-dr] and is discussed below; however, the enrichment among the low-variance genes is weaker.

# Discussion

Using large publicly available data sets allowed us to probe the landscape of transcriptional variance in humans.
We find a broadly similar pattern of transcriptional variance, evidenced by the high correlations between gene expression variance across most studies, consistent with measurements of expression variance in single cells and in populations of cells for various tissues [@Li2010-qs; @Dong2011-sa; @Alemu2014-jo].
Leveraging this similarity between gene expression variance across tissues and contexts, we developed a multivariate strategy to create a single rank of expression variance, which allowed us to order almost 13k genes (~65% of the genes expressed in humans) according to their transcriptional variance.
Using this rank, we were able to study the general properties associated with high- and low-variance genes as well as factors driving variation in variance across genes.

Some differences in gene expression variance were driven by technical aspects of gene expression measurement (with data derived from large consortia showing more similar patterns of variance across genes), and by tissue (with studies using the same tissues also showing higher similarities).
This suggests that careful consideration of sample sizes and experimental design are fundamental to the study of gene expression variance, and the usual small samples of RNA-seq studies might be underpowered for the study of this particular aspect of gene expression.
However, both the effects of study origin and tissue were small, and the largest drivers of differences across studies were idiosyncratic differences related to single data sets, with tissues known to have divergent gene expression patterns (i.e., bone marrow, blood, testis, and platelets) also showing the largest differences in gene expression variance.
Understanding the consequences of these differences in variance for specific tissues is still an open field.
It is clear, however, that differences in variance are informative beyond the differences in mean expression.
Even after we account for differences in mean expression, differences in gene expression variance carry information about tissue origin and function.

Functional analyses using GO enrichment indicated a clear link between function and gene expression variance.
On the one hand, genes with high transcriptional variance were enriched for biological functions related to response to environmental stimuli, such as immune function and tissue reconstruction.
On the other hand, low-variance genes were enriched for basic cell functions, (e.g., RNA processing, translation, DNA methylation, and cell duplication).
These results are consistent with previous analyses of gene expression variance on a tissue-by-tissue basis [@Alemu2014-jo].
This pattern of enrichment is also observed when we look at enrichment for high- or low-variance genes within the genes associated with each term in the GO hierarchy.
Basic cell function terms are enriched for low-variance genes, and terms involved in response to external stimulus are enriched for high-variance genes.

While indirect, all these patterns point to a selective structuring of gene expression variance.
Stabilizing and purifying selection are consistent: genes expected to be under strong stabilizing selection, those linked with fundamental baseline biological processes, are indeed overrepresented in the least variable genes.
These same genes are also expected to be under strong purifying selection and to show low levels of polymorphisms, which we observe.
Likewise, genes whose function is constrained by myriad interactions with several other genes, those with high connectivity, are less variable.
Furthermore, genes involved with direct interaction with the environment, which must change their pattern of expression depending on external conditions, are expected to be more variable, and again we see a strong enrichment of environmentally responsive genes among the most variable.
Given this strong link between function and variance, it is not surprising that the gene variance ranking is similar across data sets.

One interesting aspect of the GO term analysis shown in @fig:skew_entropy and @fig:go_skewness is that there is no GO biological process term associated with enrichment for intermediate variance genes: the low-entropy terms have either positive or negative skew, never zero skew.
In other words, there is no annotated biological process for which the associated genes are kept at some intermediary level of variation.
For the GO terms we used, either there is no relation between the transcriptional variance and the biological process, or there is a strong bias toward high or low-variance genes.
This suggests that selective shaping of gene expression has two modes, corresponding with (1) biological processes under strong stabilizing selection (i.e., variance-reducing selection) or (2) biological processes under disruptive selection (i.e., variance-increasing selection).
In short, we find strong support for the idea that there are genes with consistently more (or less) variable expression levels, and that these differences in variance are the result of different patterns of selection.

Following Alemu et al. [@Alemu2014-jo], we observe that epigenetic signatures of gene regulation, such as enhancer histone marks, make up a higher proportion of the surrounding genomic regions of genes that exhibit higher variance in expression.
In contrast, an accumulation of strong promoter elements and overall transcriptional activation is associated with genes with lower expression variance.
These results suggest the presence of distinct modes of regulation for genes with high vs. low variance.
Combined, the differences in the types of genomic regulatory features surrounding the high- and low-variance genes and their distinct functional annotations suggest different mechanisms of regulation of their gene expression variance [@Alemu2014-jo].
This heterogeneity could lead to detectable differences in selection signatures between distal regulatory elements and promoters depending on the transcriptional variance.
This heterogeneity in regulation for high and low-variance genes suggests that important biological information has been overlooked given the focus that the field has placed on understanding gene expression robustness, in the sense of reducing variation [@Siegal2014-dv; @Payne2015-wn; @Macneil2011-ax; @Denby2012-as]. For example, Siegal and Leu [@Siegal2014-dv] provide several examples of known regulatory mechanisms for reducing gene expression variance, but no examples for the maintenance of high gene expression variance.
We posit that it should be possible to go beyond the usual characterization of mechanisms of gene expression robustness, in the sense of reducing variation, and to explore mechanisms for the _robustness of plasticity_, that is, the maintenance of high levels of gene expression variation given environmental cues.

Given the broad consistency of gene expression variance in healthy tissues, a natural question is how do these well-regulated levels of variation behave in perturbed or disease conditions.
We find some suggestive links between tissue-specific variance ranks and disease, but these links need to be further explored using more specific methods.
Comparing two HapMap populations, Li et al. [-@Li2010-qs] showed that gene expression variance was similar in both populations and that high-variance genes were enriched for genes related to HIV susceptibility, consistent with our observation of enrichment for immune-related genes among those with more variable expression.
In a case-control experiment, Mar et al. [-@Mar2011-dr] showed that expression variance was related to disease status in Schizophrenia and Parkinson's disease patients, with altered genes being non-randomly distributed across signaling networks.
These authors also find a link between gene network connectivity and expression variance, in agreement with the effect we find using the gene expression variance rank.
The pattern of variance alteration differed across diseases, with Parkinson's patients showing increased expression variance, and Schizophrenia patients showing more constrained patterns of expression.
The authors hypothesize that the reduced variance in Schizophrenia patients reduces the robustness of their gene expression networks, what we refer to as a loss of plasticity.
This suggests that several types of shifts in gene expression variation are possible, each with different outcomes.
We highlight three distinct possibilities:
First, low-variance genes, under strong stabilizing selection, could become more variable under stress, indicating a reduced capacity for maintaining homeostasis.
Second, high-variance genes, expected to be reactive to changes in the environment, could become less variable, indicating a reduced capacity to respond to external stimuli.
Third, the covariance between different genes could be altered, leading to decoherence between interdependent genes [@Lea2019-pq].
Any one of these changes in expression variance patterns could have physiological consequences and exploring these differences should be a major part of linking gene expression to cell phenotypes and function (see Hagai et al. [-@Hagai2018-fu] for example).
Genes are also expected to differ in their capacity to maintain an optimal level of gene expression variance [@Macneil2011-ax].
Variation in robustness is linked to gene regulatory networks and epigenetic gene expression regulation [@Payne2015-wn; @Chalancon2012-ul] and, therefore, should differ across high- and low-variance genes.
Our results suggest that the mechanisms responsible for maintaining optimal levels of variation in high- and low-variance could differ and that this variability is the result of different patterns of selection.

\footnotesize

# Methods

## Data sources

We selected 57 human RNA-seq data sets from the public gene expression repositories recount3 [@Wilks2021-uj] and Expression Atlas [@Papatheodorou2020-dn].
We only used data sets with an associated publication, for which raw read count and sample-level metadata were available.
Because we are interested in individual-level variation of gene expression, we exclude single-cell studies.
Metadata and details on the included data sets can be found in the supporting information.
We use the word "studies" to refer to independent data sets, which could have been generated by the same consortium.
For example, the GTEx data are separated by tissue, and we refer to each tissue as a separate study.
We divide our data sets into three categories depending on their origin: GTEx, TCGA, and Miscellaneous.

Table: Data set source references. Columns show the study ID, with the corresponding tissue in parenthesis, and the source publication. \label{tab1}

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


## Processing pipeline

We use a standardized pipeline to measure gene expression variance while removing extraneous sources of variation.
Because we are interested in variation under non-perturbed conditions, data from case-control studies were filtered to keep only control samples.
Technical replicates were summed.
For each study, we filtered genes that did not achieve a minimum of 1 count per million (cpm) reads in all samples and a mean of 5 cpm reads across samples.
To account for library size and the mean-variance relation in RNA-seq count data, we applied a variance stabilizing transformation implemented in the function vst from the DESeq2 R package [@Love2014-mp] to the genes passing the read-count filters.
This mean-variance correction was verified by plotting mean-variance relations before and after correction, and these plots can be seen in the supporting information (SI appendix 1).
Various technical covariates (like experimental batch, sex, etc.) were manually curated from the metadata associated with each study and accounted for using an independent linear fixed-effects model for each study.
A list of covariates used for each study is available in the supporting information (SI data 1).
Outlier individuals in the residual distribution were removed using a robust Principal Component Analysis (PCA) approach of automatic outlier detection described in @Chen2020-fy.
This procedure first estimates robust Principal Components for each study and then measures the Mahalanobis distance between each sample and the robust mean.
Samples that are above the 0.99 percentile in Mahalanobis distance to the mean are marked as outliers and removed.
We verify that the batch effect correction and outlier removal are reasonable by using PCA scatter plots after each step of the pipeline to check the result for residual problems like groupings or other artifacts.
These PCA plots before and after batch correction and outlier removal are also included in SI appendix 1.
After all sample filtering, the mean sample size we used for each data set was 390, with a median of 251, and ranged from 12 to 2931 samples.
Gene expression standard deviations (SDs) are measured as the residual standard deviations after fixed effect correction and outlier removal.
We choose standard deviation as a measure of variation to have a statistic on a linear scale, and we do not use the coefficient of variation because we have already corrected for mean differences and for the mean-variance relation inherent to RNA-seq count data [@De_Jong2019-po].
The full annotated pipeline is available on GitHub at [github.com/ayroles-lab/expressionVariance-code](https://github.com/ayroles-lab/expressionVariance-code).

## Correlations in transcriptional variance

We assessed the similarity in gene expression variance across studies by using a across-study Spearman correlation matrix of the measured SDs.
Only genes present in all studies were used to calculate the Spearman correlation matrix, ~4200 genes in total.
Using Spearman correlations avoids problems related to overall scaling or coverage differences, and allows us to assess if the same genes are usually more or less variable across studies.
To investigate the factors involved in determining correlations between studies, we used a Bayesian varying effects model to investigate the effect of study origin and tissue on the correlations across studies.
This model is designed to take the non-independent nature of a set of correlations into account when modeling the correlation between gene expression SDs.
This is accomplished by adding a per-study random effect, see [@Dias2022-au] for details.
The Fisher z-transformed Spearman correlations across studies ($z(\rho_{ij})$) are modeled as:

$$\begin{aligned}
z(\rho_{ij}) &\sim N(\mu_{ij}, \sigma) \\
\mu_{ij} &= \mu_0 + \alpha_i + \alpha_j + \beta X \\
\alpha_i &\sim N(0, \sigma_{\alpha})
\end{aligned}$$

The $\alpha_i$ terms account for the non-independence between the pairs of correlations and estimate the idiosyncratic contribution of each study to all the correlations it is involved in. The fixed effects encoded in the design matrix $X$ measure the effects of tissue congruence and study-origin congruence. We also explored a version of this model that included the effect of sample size on the pairwise correlations, but sample size did not have a relevant effect and so was dropped in the final model. All fixed effect parameters ($\beta$) and per-study parameters ($\alpha_i$) receive weakly informative normal priors with a standard deviation of one quarter. For the overall variance ($\sigma$) we use a unit exponential prior, and for the intercept ($\mu_0$) a unit normal prior. This model was fit in Stan [@Carpenter2017-pn] via the _rethinking_ R package [@mcelreath2020statistical], using eight chains, with 4000 warm-up iterations and 2000 sampling iterations per chain. Convergence was assessed using R-hat diagnostics [@Gelman2013-ae], and we observed no warnings or divergent transitions.

__Gene expression SD rank:__ Given that most of the variation in the Spearman correlation across studies is explained by a single principal component (PC1 accounts for 62% of the variation in the across-study Spearman correlation matrix, while PC2 accounts for only 5%; see SI fig. 5), we use the ranked projections of gene expression SDs in this principal component (PC1) to create an across-study rank of gene variation.
The higher the rank, the higher the expression SD of a given gene.
Genes that were expressed in at least 50% of the studies were included in the rank.
To project a particular gene onto the PC1 of the across-study correlation matrix, we impute missing values using a PCA-based imputation [@Husson2019-sl].
The imputation procedure has minimal effect on the ranking and imputing missing SD ranks at the beginning or at the end of the ranks produces similar results.
We also create a tissue-specific variance ranking, using the same ranking procedure but joining studies done in the same tissue type.
For this tissue-level ranking, we only use genes that are expressed in all studies of a given tissue, and in this case, no imputation is required.
For tissues that are represented by a single study, we use the SD ranking for that study as the tissue rank.

__Gene expression mean rank:__ We also use the same strategy to create a mean gene expression rank, repeating the process but using mean expression instead of standard deviation. All ranks are available in the supporting information.

## Gene level statistics

__Genetic variation__: Genetic variation measures were obtained from the PopHuman project, which provides a comprehensive set of genomic information for human populations derived from the 1000 Genomes Project.
Gene-level metrics were used when available.
If only window-based metrics are available, we assembled gene-level information from 10 kb window tracks where each window that overlaps with a given gene was assigned to the gene and the mean metric value is reported.
In parallel, we use the PopHumanScan data set, which expands PopHuman by compiling and annotating regions under selection.
Similarly, we used gene-level information when possible, and for tracks with only window-based metrics, gene-level information was assembled from the 10 kb windows using the same assignment method described above.
Nucleotide diversity ($\pi$), the average pairwise number of differences per site among the chromosomes in a population [@Nei1979-hg], provides insight into the genetic diversity within a population, in this case, the CEU population within 1000 genomes.

__Gene connectivity__: For each data set, we calculated the average weighted connectivity for all genes by creating a fully connected gene-by-gene graph in which each edge is weighted by the Spearman correlation between gene expression levels across samples.
We then trimmed this graph by keeping only edges for which the Spearman correlation is significant at a BH false discovery rate of 1%.
In this trimmed network, we then took the average of the Spearman correlation of all remaining edges for each gene.
So, for each study, we have a measure of the average correlation of each gene with every other gene.
The average connectivity for each gene is the average across all studies in which that gene is expressed.

__Cross-tissue vs. tissue-level chromatin states__: We use the universal [@vu2022universal] and tissue-specific [@Ernst2015-zk] ChromHMM [@ernst2012chromhmm] chromatin states to compare the non-overlapping genome segmentation to cross-tissue and tissue-level gene expression variance metrics. We use the proportion of the gene regions (gene ± 10 kb) made up of each of the ChromHMM chromatin states.

__Correlations:__ We use the ppcor R package v1.1 [@kim2015ppcor] to run the pairwise partial Spearman correlations between gene-level statistics and the gene expression variance rank while controlling for the mean expression rank. P-values are corrected using the Benjamini-Hochberg procedure and comparisons with an adjusted p<0.05 are considered significant.

## Gene function assessment

__GO term enrichment__: All gene ontology (GO) analyses were done using the clusterProfiler R package v4.2.2 [@Wu2021-db] and the Org.Hs.eg.db database package v3.14.0 [@godb]. GO and all further enrichment analyses used the hypergeometric test to assess the significance of the enrichment.

__Environmentally responsive genes__: We use the list of environmentally responsive genes available in the supporting information from Lea et al. [@Lea2021-hf] and test for enrichment of environmentally responsive genes in the genes within the highest and lowest 5% of gene expression variance rank.

__Housekeeping genes__: Human housekeeping genes were identified as genes that are expressed with low variance in all 52 human cell and tissue types, assessed in over 10,000 samples [@Hounkpe2020-yq]. We test for enrichment of housekeeping genes in the genes within the highest and lowest 5% of gene expression variance rank.

__Probability of being loss-of-function intolerant (pLI)__: Genes that are likely haploinsufficient (i.e., intolerant of heterozygous loss-of-function variants) were detected as those with fewer than expected protein-truncating variants (PTVs) in ExAC [@Lek2016-xw]. We use genes with a pLI > 0.9 to test for the enrichment of loss-of-function intolerant genes in the genes exhibiting the highest and lowest 5% gene expression variance estimates.

__Secreted genes__: We use The Protein Atlas [@uhlen2015tissue] to extract information on which proteins are secreted [@uhlen2019human] and test for enrichment of genes with secreted products in the genes within the highest and lowest 5% of gene expression variance rank.

__Immediate early genes (IEGs):__ Human IEGs were curated from the literature in @Arner2015-be as genes that respond to experimental stimulation through up-regulation within the first 60 minutes of the experiment. We use the hypergeometric test to assess the significance of the enrichment. Immediate early genes (IEGs): Human IEGs were curated from the literature in @Arner2015-be as genes that respond to experimental stimulation through up-regulation within the first 60 minutes of the experiment.

__Disease annotations__: We use the gene annotations for involvement with diseases provided by the supporting information Table S2 from Zhang et al. [-@Zhang2020-cl] and test for enrichment for disease annotations in the genes within the highest and lowest 5% of gene expression variance rank.

## Code availability

Code for reproducing all analyses and figures, along with a walk-through, is available at [github.com/ayroles-lab/ExpressionVariance](https://github.com/ayroles-lab/ExpressionVariance).

## Supporting information

Supporting information is available at [github.com/ayroles-lab/expressionVariance-manuscript](https://github.com/ayroles-lab/expressionVariance-manuscript).

1. SI figure 1 - Modeling the correlations between transcriptional variance across studies.

1. SI figure 2 - GO enrichment analysis of the most and least variable genes.

2. SI figure 3 - Across-study and tissue-specific gene expression variance and mean correlations
with non-overlapping chromatin states through ChromHMM.

1. SI figure 4 - Proportion of gene regions made up of ChromHMM chromatin states for genes in
the top and bottom 5% of the across-study mean rank metric.

1. SI figure 5 - Scree plot showing variance explained by each PC of the across-study Spearman correlation matrix of gene expression standard deviations.

1. SI table 1 - Variance and mean rank metrics and the corresponding ChromHMM annotations
used.

1. SI table 2 - Enrichment analysis of environmentally responsive genes in LCLs.

1. SI appendix 1 - Diagnostics plots for processing pipeline.

1. SI data 1 - Study metadata - Metadata file describing the data used in the study as well as some intermediate processing information.

1. SI data 2 - Study gene lists - List of genes included in each study after filtering.

1. SI data 3 - Gene expression means and standard deviations - Tables with final calculated means and standard deviations.

2. SI data 4 - Gene ranks - Gene expression mean and variance ranks, across-study and tissue-specific.

3. SI data 5 - GO enrichment - Combined table describing gene ontology enrichment in the top 5% and bottom 5% of genes as ranked by variance.

# Author Contributions

**Conceptualization**: S.W., D.M., L.P., and J.A.
**Analysis**: S.W., D.M., and K.G.
**Draft**: D.M.
**Review and Editing**: S.W., D.M., K.G., L.P., and J.A.
**Funding Acquisition**: S.W., D.M., K.G., L.P., and J.A.

# Acknowledgments

We thank the Ayroles lab members for their support, and Noah Rose and Cara Weisman for their thoughtful comments. We thank Pedro Madrigal for help with the Expression Atlas interface. S.W. is supported by the NSF Graduate Research Fellowship Program (DGE-2039656). D.M. was funded by a fellowship from the Princeton Presidential Postdoctoral Research Fellows Program. K.G. is funded by NIH grant F32ES034668. L.P. was funded by a Long-Term Postdoctoral Fellowship from the Human Frontiers Science Program. J.A. is funded by NIH/NIEHS grant R01-ES029929. This study was supported in part by the Lewis-Sigler Institute for Integrative Genomics at Princeton University.

The author(s) are pleased to acknowledge that the work reported on in this paper was substantially performed using the Princeton Research Computing resources at Princeton University which is consortium of groups led by the Princeton Institute for Computational Science and Engineering (PICSciE) and Office of Information Technology's Research Computing.

# References

<!-- %% \printbibliography -->
