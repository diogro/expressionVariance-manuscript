---
title: "Characterizing the landscape of gene expression variance in humans"
author:
  - Scott Wolf:
      institute:
        - lsi
      email: swwolf@princeton.edu
      orcid: 0000-0000-0000-0000
      equal_contributor: "yes"
  - Diogo Melo:
      institute: eeb
      equal_contributor: "yes"
  - Kristina M. Garske:
      institute: lsi
  - Luisa Pallares:
      institute: fml
  - Julien Ayroles:
      institute:
        - lsi
        - eeb
      email: jayroles@princeton.edu
      correspondence: "yes"
institute:
  - lsi:
      name: Lewis-Sigler Institute for Integrative Genomics, Princeton University
  - eeb:
      name: Department of Ecology and Evolutionary Biology, Princeton University
  - fml:
      name: Friedrich Miescher Laboratory, Max Planck Society
classoption: twocolumn
output: pdf_document
geometry:
- top=20mm
- left=25mm
- right=25mm
- bottom=20mm
header-includes:
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
# csl: ./cse.csl
sansfont: Skolar Sans PE TEST
bibliography: ./references.bib
---

<!-- https://tex.stackexchange.com/questions/536353/biblatex-colors-and-links-only-the-year-not-the-rest-of-the-citation -->
\makeatletter
\renewbibmacro*{cite:plabelyear+extradate}{%
  \iffieldundef{labelyear}{}
    {\clearfield{labelmonth}% don't want months in citations
     \clearfield{labelday}% don't want days in citations
     \clearfield{labelendmonth}% don't want months in citations
     \clearfield{labelendday}% don't want days in citations
     \iffieldsequal{labelyear}{labelendyear}% Don't want no-op year ranges
       {\clearfield{labelendyear}}
       {}%
     \iffieldundef{origyear}
       {}
       {\printorigdate%
        \setunit*{\addslash}}%
     \iffieldundef{related}
       {}
       {\iffieldequalstr{relatedtype}{reprintfrom}
         {\entrydata*{\thefield{related}}{\printlabeldateextra}%
          \setunit*{\addslash}}
         {}}%
     \printlabeldateextra}}

\renewbibmacro*{cite}{%
  \iffieldequals{fullhash}{\cbx@lasthash}
   {\setunit{\compcitedelim}%
    \printtext[bibhyperref]{%
      \usebibmacro{cite:plabelyear+extradate}}}%
   {\printtext[bibhyperref]{%
      \ifnameundef{labelname}
       {\usebibmacro{cite:noname}%
         \setunit{\printdelim{nameyeardelim}}%
         \usebibmacro{cite:plabelyear+extradate}%
         \savefield{fullhash}{\cbx@lasthash}}
       {\ifnameundef{shortauthor}
         {\printnames{labelname}}%
         {\cbx@apa@ifnamesaved
           {\printnames{shortauthor}}
           {\ifnameundef{groupauthor}
             {\printnames[labelname]{author}}
             {\printnames[labelname]{groupauthor}}%
            \addspace\printnames[sabrackets]{shortauthor}}}%
         \setunit{\printdelim{nameyeardelim}}%
        \usebibmacro{cite:plabelyear+extradate}%
        \savefield{fullhash}{\cbx@lasthash}}}}%
   \setunit{\multicitedelim}}

\renewbibmacro*{textcite}{%
  \iffieldequals{fullhash}{\cbx@lasthash}
    {\setunit{\compcitedelim}%
     \printtext[bibhyperref]{%
       \usebibmacro{cite:plabelyear+extradate}}}
    {%
    \ifbool{cbx:parens}
      {\bibcloseparen\global\boolfalse{cbx:parens}}
      {}%
      \setunit{\compcitedelim}%
      \ifnameundef{labelname}
       {\iffieldundef{shorthand}%
         {\printtext[bibhyperref]{%
            \usebibmacro{cite:noname}}%
          \setunit{\ifbool{cbx:np}%
                   {\printdelim{nameyeardelim}}%
                   {\global\booltrue{cbx:parens}\addspace\bibopenparen}}%
          \printtext[bibhyperref]{%
            \usebibmacro{cite:plabelyear+extradate}}}
         {\printtext[bibhyperref]{%
            \usebibmacro{cite:shorthand}}}}
       {\printtext[bibhyperref]{%
          \ifnameundef{shortauthor}%
           {\printnames{labelname}}
           {\cbx@apa@ifnamesaved
             {\printnames{shortauthor}}
             {\ifnameundef{groupauthor}
               {\printnames[labelname]{author}}
               {\printnames[labelname]{groupauthor}}}}}%
        \setunit{\ifbool{cbx:np}
                  {\printdelim{nameyeardelim}}
                  {\global\booltrue{cbx:parens}\addspace\bibopenparen}}%
        \printtext[bibhyperref]{%
          \ifnameundef{shortauthor}
           {}
           {\cbx@apa@ifnamesaved
             {}
             {\printnames{shortauthor}\setunit{\printdelim{nameyeardelim}}}}%
          \usebibmacro{cite:plabelyear+extradate}}%
        \savefield{fullhash}{\cbx@lasthash}}}}
\makeatother
<!-- # Abstract -->

# Intro

Molecular phenotypes such as gene expression are a powerful tool for understanding physiology, disease, and evolutionary adaptations.
In this context, average trait values are usually the focus of investigation, while variation is treated as a nuisance [@De_Jong2019-po].
However, gene expression variance can be directly involved in determining fitness [@Fraser2004-sv; @Wang2011-ts], and changes in the associations between gene expression can be indicative of disease, even in the absence of changes in mean expression [@Lea2019-pq].
From an evolutionary perspective, the availability of gene expression variance is what allows evolutionary change, and the genetic architecture of gene expression variance can also evolve [@Bruijning2020-bf].
Understanding the landscape of gene expression variance, and how variable it is across genes and across human populations is then a neglected avenue to understand biological evolution and our relation to the environment.
In particular, we lack a clear picture of which genes show more gene expression variance, or even if the pattern of gene expression variance is consistent across populations.

Several competing forces act to shape gene expression variance [@Houle1998-mj; @Bruijning2020-bf], and the outcome of the interaction between these processes is still poorly understood [@Hansen2011-es].
From a genomic perspective, we expect the influx of new mutations to increase observed variation, while the selective removal of polymorphisms, via purifying selection or selective sweeps, would decrease variation.
From a trait-centric perspective, stabilizing selection should decrease variation around an optimal value, and directional selection can lead to transient increase in variance while selected alleles sweep to fixation, followed by a reduction in variance as these alleles become fixed.
This simple picture is complicated by epistatic interactions between loci and other aspects of genetic architecture.
For example, pleiotropic effects allow selection on one trait to influence the variance of other traits, potentially limiting the direct response to selection [@Wagner1997-hw; @Pavlicev2011-xm].
The indirect effect of directional selection on variance opens the possibility that the main driver of gene expression variance is not direct selection on variance but indirect effects due to selection on trait means [@Hansen2011-es].
Furthermore, gene by environment (GxE) interactions can also lead to changes in the observed phenotypic variance of gene expression, further complicating the landscape of variation.
To what extent these different processes shape gene expression variance is an open question.
If consistent selection across populations is the main driver of gene expression variance, we would expect to have consistently more or less variable genes.
If idiosyncratic selection patterns and context specific environmental interactions are more important, we could observe large differences in gene expression variance across populations.

Even within individuals, gene expression is also variable across tissues [@GTEx2017-xb].
To what extent differences in mean expression level translate to differences in expression variance is not clear.
Of course, genes that are exclusively expressed in a single cell type or tissue are necessarily more variable in that particular tissue, but differentially expressed genes could also be more variable in a particular context.
For example, stabilizing selection on gene expression could be more intense depending on the role of that gene in a particular tissue, leading to a local reduced variation and differences in variation across tissue.
Alternatively, expression variation across tissues could be tightly coupled, and in this example, selection in one tissue would lead to a reduction in variance across tissues, resulting in a consistent pattern of variation.

Here, we use public gene expression data sets to evaluate how the differences in gene expression variance is structured across independent samples.
We collected and compared the gene expression variance across many studies and used the similarities across these studies to create a gene expression variance ranking, which orders genes from least variable to most variable.
We then explore the expected drivers of this gene expression ranking, showing that both cis and trans regulation are involved with the determination of gene expression variance.
Finally, we explored the link between gene expression variance and biological function by leveraging gene ontology annotations.

# Results


\begin{figure*}[t!]
    \centering
     \includegraphics[width=\dimexpr 1\textwidth]{figures/fig1.png}
    \caption{A. Correlation heatmap showing the across study Spearman correlation of standard deviations. Pairs of studies with more similar patterns of gene expression variance have higher correlations; B. Histogram of the correlations shown in the previous panel; C. Standard deviation correlation PCoA, with colors ; D. Density plot of standard deviations after z-normalization. Inset plot shows distribution of mean centered standard deviations grouped by study without normalization. The corresponding rug plots show the location of the highest ranking gene in standard deviation rank (right, blue) and lowest (left, red).}
    \label{fig:sd_corr}
\end{figure*}

Gene expression standard deviations (SDs) were calculated for each data set using a unified pipeline that normalized the mean-variance relation in count data, controlled for batch effects, and removed outliers (see methods for details).
Spearman correlations ($\rho_s$) between gene expression SDs reveals a broadly similar rank of gene expression variance, so genes that are most variable in one study tend to be most variable in all studies ([@fig:sd_corr]A and B).
A principal coordinate analysis [@Gower1966-dk] using $|1 - \rho_s|$ as a distance measure does not show clearly delineated groups, but GTEx and TCGA studies are clustered among themselves and close together ([@fig:sd_corr]C).
This indicates some effect of study source on the similarity between gene expression SD across studies, which we explore in detail below.
Observed range of gene expression SD across genes is variable across studies, but can be normalized so that the distributions are comparable ([@fig:sd_corr]D).
Given that the correlations across studies are broadly high, indicating similar ordering of the genes, we seek to summarize the differences in variance across genes by using a single cross-study rank, averaging the ordering across all studies.
To create this rank, we use the score of each gene in the first principal component of the Spearman correlation matrix.
This generates a ranked list of genes, with most variable genes having highest rank.
The red and blue ticks at the bottom of [@fig:sd_corr]D show the positions on the SD distributions of the least and most variable gene in our variance rank.

\begin{figure*}[t!]
    \centering
    \includegraphics[width=\linewidth]{figures/correlationModeling.png}
    \caption{Coefficients estimates from a linear model using the among studies Spearman correlations as the response variable. These correlations are shown in fig.~\ref{fig:sd_corr}A and B. In the linear model, correlations are Fisher z-transformed. Study source and tissue are added as fixed effects. Coefficient estimates are shown with 50\% and 95\% credibility intervals. Panel A: The per-study random effect captures the non-independence of the correlation values and estimates the characteristic contribution of each study to the correlation. For example: comparisons involving bone marrow (from GTEx) tend to be lower than the others. Panels B and C: Fixed effect estimates: correlations among studies that use the same tissue are higher, and correlations involving studies in the ``Other'' category (non GTEx and TCGA) tend to be lower, while comparison involving GTEx and TCGA are higher.} <!-- `` -->
    \label{fig:corr_model}
\end{figure*} 


## What drives differences in gene expression variance?

To characterize the drivers of across study similarity, we directly model the correlations across studies using a mixed effect linear model [@Dias2021-wk; @Dias2021-hb].
In this mode, we use study, sampled tissue, and study origin as predictors of the pairwise correlations (see Methods).
This modeling ([@fig:corr_model]) shows that comparisons of studies within GTEx and TCGA have on average higher values of $\rho_s$, but also that comparing studies across GTEx and TCGA also shows a mild increase in the average correlation ([@fig:corr_model]C).
Correlation involving studies that are not from TCGA and GTEx (marked as "Other") are on average lower ([@fig:corr_model]C).
Since these two sources are independent, this effect on the similarities could be due to the quality of the data coming from these two large projects.
Tissue also affects the similarity between gene expression SD, with studies using the same tissue being, on average, more similar ([@fig:corr_model]B).
The largest effects on the correlations are those associated with individual studies, in particular some specific tissues, i.e., comparisons involving bone marrow (from GTEx) and study SRP057500 (which used platelets) are on average lower ([@fig:corr_model]A).
These studies also show up further away in the PCoA plot in [@fig:sd_corr]C.

## Does biological function explain variance in expression?

To explore the relationship between variance and function, we took the top 5% most variable and the bottom 5% least variable genes in our ranking (about ~560 genes in each group) and performed a Gene Ontology (GO) enrichment analysis within each group.
This allows us to establish the representative functions of these consistently high and low-variance genes.
In total, using a Benjamini-Hochberg adjusted p-value threshold of $10^{-3}$, we found 59 enriched terms in the low variance genes, and 738 enriched terms in the high variance genes (see supporting table 1 for a complete listing).
Among the 5% most variable genes we observe enrichment for biological processes like immune function, response to stimulus, maintenance of homeostasis, and tissue morphogenesis (@fig:go_tails, left).
In line with this GO term enrichment, the top 5% most variable genes are enriched 7.7-fold for genes that encode secreted proteins, relative to all other genes ($p < 10^{-3}$).
Among the 5% least variable genes we see enrichment for housekeeping functions like mRNA processing, cell cycle regulation, methylation, histone modification, translation, transcription, and DNA repair (@fig:go_tails, right); and accordingly we find that previously characterized human housekeeping genes [@Hounkpe2020-yq] are enriched within the 5% least variable genes 2.0-fold relative to all other genes ($p < 10^{-3}$).
The genes exhibiting the lowest variance (lowest 5%) are also enriched for those that have been previously shown to have a high probability of being loss-of-function intolerant (pLI) [@lek2016analysis] (1.2-fold enrichment, $p < 10^{-3}$).
Genes with a high pLI have been shown to be important in housekeeping functions, and have higher mean expression values across a broad set of tissues and cell types [@lek2016analysis]. Our result that genes with low variance are enriched for both housekeeping genes and genes with high pLI is consistent with this previous report; and we further see that the mean expression of genes positively correlates with pLI (Partial Spearman correlation $\rho_s$ = 0.32, $p < 10^{-3}$), showing the opposite relationship between variance and mean expression when considering pLI.

We also explore the distribution of expression variance among the genes associated with GO terms.
For this, we gather all biological process GO terms in level 3 (i.e. terms that are at a distance of 3 for the top of the GO hierarchy).
Using only the set of genes that are associated with at least one of these level-3 terms, we separate the genes into expression variance deciles, with the first decile having the lowest variance.
We then count how many genes in each decile has been associated with each term. If variance rank is not linked to the GO annotations, terms should have an equal proportion of genes in each decile.
We measure how far from this uniform allocation each term is by measuring the Shannon entropy of the proportion of genes in each decile.
Higher entropy is associated with more uniform distribution of genes across deciles.
GO terms with low entropy indicate some decile is over-represented in the genes associated with that term. We also measure skewness for each term, which should be zero if no decile is over-represented, negative if high-variance terms are over-represented, and positive if low-variance deciles are over-represented.
Skewness by entropy for each GO term can be seen in @fig:skew_entropy. Positive-skew low-entropy terms, those enriched with low-variance genes, are associated with house keeping functions, like RNA localization, translation initiation, methylation and chromosome segregation (@fig:go_skewness A).
Likewise, terms with negative skew and low entropy, enriched for high-variance genes, are related to immune response, tissue morphogenesis, chemotaxis---all dynamic biological functions related to interacting with the environment (@fig:go_skewness B).

Both GO analyses suggests a strong influence of biological function in determining gene expression variance. Genes associated with baseline fundamental functions, expected to be under strong stabilizing selection, are also low-variance; high-variance genes are associated with responding to external stimuli (i.e., tissue reorganization and immune response).

\begin{figure}
    \centering
    \includegraphics[width=\linewidth]{figures/local_go_lowerUpper.png}
    \caption{Gene set enrichment analyses testing for over representation of gene ontology categories in the upper and lower 5\% quantiles of the gene variance rank. High-variance gene are enriched for terms related to immune function, response to wounding, blood vessel morphogenesis and inflammatory response. In contrast, low-variance genes are associated with translation, control of methylation, RNA processing, chromosome separation, and other cell housekeeping functions. All displayed terms are significant with a 5\% FDR corrected p-value below $10^{-3}$.}
    \label{fig:go_tails}
\end{figure}

\begin{figure*}[t!]
    \centering
    \includegraphics[width=\linewidth]{figures/skew_entropy.png}
    \caption{Relationship between skew and entropy of rank decile distributions for each GO term. The GO terms are filtered for gene counts greater than 100 as in fig. \ref{fig:go_skewness}.}
    \label{fig:skew_entropy}
\end{figure*}

\begin{figure*}[t!]
    \centering
    \includegraphics[width=\linewidth]{figures/GOterm_decile_barplot.png}
    \caption{Distributions of decile ranks of second level GO terms. Each plot shows the count of genes in each decile of the rank. These GO terms are filtered for gene counts greater than 100 and sorted by the skewness of the distribution. The top panel shows the top 5 and the bottom panel shows the bottom 5.}
    \label{fig:go_skewness}
\end{figure*}



## Sequence variation and gene expression connectivity

We use gene-level statistics capturing evolutionary and population variation to link processes that potentially influence variation in gene expression to the observed variance rank.
We focus on three gene-level measures: nucleotide diversity($\pi$), gene expression connectivity, and the proportion of substitutions that are adaptive ($\alpha$).
Nucleotide diversity is used as a proxy for cis-regulation sites, and we expect variation to increase with diversity. Here, we find a partial Spearman's correlation of 0.184 ($p < 10^{-3}$). Connectivity, a proxy for regulatory interactions with other genes, in turn, should be negatively correlated with variation, as highly connected genes are expected to be more constrained in their variability. The resulting partial Spearman's correlation is -0.024 ($p \approx 6 \times 10^{-3}$). Finally, we find a partial Spearman's correlation of -0.046 ($p \approx 1 \times 10^{-3}$) for the proportion of substitutions that are adaptive.

<!-- Finally, $d_{XY}$....

We also use linear models to measure the association between rank and these statistics while accounting for the effect of mean expression.
The strongest effect was of... (stats), followed by ...(stats), and  ... (stats). -->

<!-- \begin{table}[]
\resizebox{\linewidth}{!}{%
\begin{tabular}{|l|l|l|}
\hline
Covariate         & P-value      & Partial Spearman Correlation \\ \hline
pi                & $9.57 \times 10^{-85}$ & 0.184              \\ \hline
mean connectivity & $5.87 \times 10^{-3}$ & -0.024              \\ \hline
alpha             & $1.18 \times 10^{-3}$ & -0.046              \\ \hline
\end{tabular}%
}
\end{table} -->

## How do molecular signatures of gene regulation relate to gene expression variance?

We assess how local chromatin state relates to gene expression variance.
We use each gene, including the surrounding 10 kb on both ends, to calculate the proportion of gene regions that corresponds to functional chromatin states and annotations previously used to stratify the genome into interpretable functional categories [@finucane2015partitioning], including promoter and enhancer regions, open chromatin (assayed through DNase hypersensitivity (DHS)), and transcription factor binding sites (TFBS).
Biochemical features associated with cis gene regulation are positively correlated with the gene expression variance rank metric, regardless of whether the regulatory effect on gene expression is positive or negative [KG supp fig 1].
For example, both the proportion of gene regions made up of enhancers and repressed genomic states are positively correlated with gene expression variance (Benjamini-Hochberg adjP<0.05) [KG supp fig 1].
As expected, the proportion of gene regions made up of repressed genomic states is inversely correlated with mean expression of the gene [KG supp fig 1].
The magnitude of the correlation with general RefSeq gene features, such as promoter and coding sequence, is lower for both the variance and mean, and we see that this coincides with an overal positive (in the case of the mean) and negative (in the case of the variance) associations with gene density in the expanded gene regions (gene +/- 250 kb) [KG supp table 1; KG supp fig 1].
Furthermore, the biochemical properties associated with promoter flanking regions, as well as transcribed states, are inversely correlated with gene expression variance [KG supp fig 1].
Taken together, these results are in line with gene expression variance being more associated with distal (i.e., non-promoter) gene regulation, rather than overall active transcriptional state of a gene region, as is the case with mean gene expression.

<!-- tissue-level results incoming... -->

# Discussion

By using large publicly available data sets, we were able to probe the landscape of gene expression variance in several human tissues.
Differences in gene expression variance were driven by technical aspects of gene expression measurement, with data derived from large consortia showing more similar patters of variance across genes; and by tissue, with studies using the same tissues also showing higher similarities.
This would suggest that careful consideration of sample sizes and experiment design are fundamental to the study of gene expression variance, and the usual small samples of RNAseq studies might be underpowered for the study of this particular aspect of gene expression.
Furthermore, the largest driver of differences across studies was idiosyncratic differences related to single data sets, with tissues know to have divergent gene expression patterns (i.e. bone marrow, blood, testis, and platelets) also showing the largest differences in gene expression variance.
Understanding the consequences of these differences in variance for specific tissues is still an open field.
It is clear, however, that differences in variance are informative beyond the differences in mean expression.
Even after we account for differences in mean expression, differences in gene expression variance carry information about tissue origin and function.

While these observed differences are notable, we also find a broadly similar pattern of gene expression variance across studies, with high correlations between gene expression variance across most studies (75% of correlations are between 0.45 and 0.9), consistent with measurements of expression variance in single cells and in populations of cells for various tissues [@Li2010-qs; @Dong2011-sa; @Alemu2014-jo].
Leveraging this similarity between gene expression variance, we used a multivariate strategy to create a single rank of expression variance, which allowed us to order almost 13k genes according to their expression variance.
This rank is associated with within-gene sequence variation, with more polymorphic genes being more variable.
Furthermore, genes with high connectivity, those with higher levels of gene expression correlations with other genes, are less variable.

Functional analysis using GO enrichment indicated a clear link between function and gene expression variance.
First, genes with high gene expression variance were enriched for biological functions related to reacting to environmental pressures, like immune function and tissue reconstruction.
Likewise, low variance genes were enriched for basic cell function, like RNA processing, translation, DNA methylation, and cell duplication.
This pattern of enrichment is also observed when we look at enrichment for high or low variance genes within the genes associated with each terms in the GO hierarchy.
Basic cell function terms are enriched for low variance genes, and terms involved in response to external stimulus are enriched for high variance genes.

While indirect, all these patterns point to a selective structuring of gene expression variance.
Stabilizing and purifying selection are consistent, genes expected to be under strong variance reducing stabilizing selection, those linked with fundamental baseline biological processes, are indeed over represented in the least variable genes.
These same genes are also expected to be under strong purifying selection and show low levels of substitution and polymorphisms, which we observe.
Likewise, genes whose function is constrained by myriad interactions with several other genes, those with high connectivity, are less variable.
Furthermore, genes involved with direct interaction to the environment, which must change their pattern of expression depending on external conditions, are expected to be more variable, and again we see a strong enrichment of genes related to interacting with the environment among the most variable.
Given this strong functional linkage between function and variance, it is not surprising that the gene variance ranking be somewhat consistent across studies, allowing us to create our ranking in the first place.
We find strong support for the idea that there are indeed genes with consistently more (or less) variable expression levels, and that these differences in variance are the result of different patterns of selection.

Given this consistency, the natural question is then how do these well regulated levels of gene expression variance behave in perturbed or disease conditions.
Comparing two HapMap populations, Li et al [-@Li2010-qs] showed that gene expression variance was similar in both populations, and that high variance genes were enriched for genes related to HIV susceptibility, consistent with our observation of enrichment for immune related genes among those with more variable expression.
In a case-control experiment, Mar et al. [-@Mar2011-dr] showed that expression variance was related to disease status in Schizophrenia and Parkinson's disease patients, with altered genes being non randomly distributed across signaling networks.
These authors also find a link between gene network connectivity and expression variance, consistent with the effect we find using the gene expression variance rank.
Also, the pattern of variance alteration differed across diseases, with Parkinson's patients showing increased expression variance, and Schizophrenia patients showing more constrained patters of expression.
The authors hypothesizes that the reduced variance in Schizophrenia patients reduces the robustness of their gene expression networks.
This suggests several types of shifts in gene expression variation are possible, with different outcomes. 
We highlight three different possibilities:
First, low variance genes, under strong stabilizing selection, could become more variable under stress, indicating a reduced capacity for maintaining homeostasis. 
Second, high variance genes, expected to be reactive to changes in the environment, could become less variable, indicating reduced capacity to responding to external stimuli.
Third, the covariance between different genes could be altered, leading to decoherence between interdependent genes [@Lea2019-pq].
Any one of these changes in expression variance pattern could have major physiological consequences.

Presumably genes will differ in their capacity to maintain their baseline variation levels.



__Drafts:__

<!-- - Differences in gene expression variance can be driven by experimental features, so care must be taken when designing experiments focused on finding gene expression differences.
- Tissue differences in gene expression variance are an unexplored field. -->
- Variation in perturbed conditions
- Differences in robustness
<!-- - Did we look into tissue specificity and rank?  -->


\footnotesize

# Methods

## Data sources

We selected 60 RNA-seq studies with large sample sizes from public gene expression repositories recount3 [@Wilks2021-uj] and Expression Atlas [@Papatheodorou2020-dn]. Because we are interested in population variation of gene expression, we exclude single-cell studies and focus only on studies derived from populational samples. We only used studies for which raw read count data was available, and for which we could parse the metadata for batch effects.
We use studies to refer to independent data sets, which could have been generated by the same consortium.
For example, the GTEx data are separated by tissue, and we refer to each tissue as a separate study.


## Data processing pipeline

We use a standardized pipeline to measure gene expression variance while removing extraneous sources of variation.
Data from case-control studies was filtered to keep only control samples.


For each study, we filtered genes that did not achieve a minimum of 1 count per million (cpm) reads in all samples and a mean 5 cpm reads.
To account for the mean variance relation in count data, remaining genes were subjected to the variance stabilizing transformation implemented in DESeq2 [@Love2014-mp].
Fixed effects were manually curated from the metadata for all studies and removed using a linear fixed effect model.
Outlier individuals in the residual distribution were removed using a robust Principal Component Analysis (PCA) approach of automatic outlier detection [@Chen2020-fy].
Gene expression standard deviation is measured as the residual standard deviation after fixed effect correction and outlier removal.

## Variance correlation

We assessed the similarity in gene expression variance across studies by using a between-study Spearman correlation matrix of the measured SDs.
Only genes present in all studies were used to calculate the Spearman correlation matrix, 4300 genes in total.
Using Spearman correlations avoids problems related to overall scaling or coverage differences, and allows us to assess if the same genes are usually more or less variable across studies.
To investigate the factors involved in determining correlations between studies, we used a varying effects model to investigate the effect of study origin and tissue on the correlations across studies.
This model is designed to take the non-independent nature of a set of correlations into account when modeling the correlation between gene expression variance.
This is accomplished by adding a per-study random effect, see [@Dias2021-hb] for details.
Given that most of the variation in the Spearman correlation across studies is explained by a single principal component, we use the ranked projections of gene expression variance in this principal component (PC1) to create an across study rank of gene variation.
The higher the rank, the higher the gene SD of a given gene.
Genes that were expressed in at least 50% of the studies were included in the rank.
In order to project a particular gene onto the PC1 of the between study correlation matrix, we impute missing values using a PCA based imputation [@Husson2019-sl].
The imputation procedure has minimal effect on the ranking, and imputing missing SD ranks at the beginning or at the end of the ranks produces similar results.

## Gene level statistics

__Genetic variation__: Genetic variation measures were obtained from the PopHuman project, which provides a comprehensive set of genomic information for human populations derived from the 1000 Genomes Project.
Gene level metrics were used when available.
If only window based metrics are available, we assembled gene level information from 5kb window tracks where each window that overlaps with a given gene was assigned to the gene and the mean metric value is reported.
In parallel, we use the PopHumanScan data set, which expands PopHuman by compiling and annotating regions under selection.
Similarly, we used gene level information when possible, and for tracks with only window based metrics, gene level information was assembled from the 10kb windows using the same assignment method described above.
Nucleotide diversity ($\pi$), the average pairwise number of differences per site among the chromosomes in a population [@Nei1979-hg], provides insight in the genetic diversity within a population, in this case CEU population within 1000 genomes.
The nucleotide diversity can also be used as an estimator of the central population genetic parameter, normally given as $\theta$.

__Gene connectivity__: As a proxy for the degree of trans regulation that each gene is subjected to, we calculate the average weighted connectivity for all genes.
To do this, for each study, we create a fully connected gene-by-gene graph in which each edge is weighted by the Spearman correlation between gene expression.
We then trim this graph by keeping only edges for which the Spearman correlation is significant at a false discovery rate of 1%.
In this trimmed network, we then take the average of the Spearman correlation of all remaining edges for each gene.
So, for each study we have a measure of the average correlation of each gene with every other gene.
The average connectivity for each gene is the average across all studies in which that gene is expressed.
As a proxy for the degree of trans regulation that each gene is subjected to, we calculate the average weighted connectivity for all genes.
To do this, for each study, we create a fully connected gene-by-gene graph in which each edge is weighted by the Spearman correlation between gene expression.

__Chromatin state correlates of gene expression variance__: We first obtain various annotations previously used to stratify the genome into interpretable functional categories [@finucane2015partitioning].
A subset of these annotations are used to quantify functional and molecular correlates of the gene expression variance metric: 1) promoter, coding, and 3' and 5' UTR are annotations from the RefSeq gene model; 2) CTCF, promoter flanking, transcribed, transcription start site, and enhancer categories were defined as the union [@finucane2015partitioning] of these annotations derived from ChromHMM/Segway across 6 cell types [@hoffman2013integrative]; 3) the repressed category was defined as the intersection [@finucane2015partitioning] of these annotations derived from ChromHMM/Segway across 6 cell types [@hoffman2013integrative]; 4) conserved elements were identified across 29 mammalian species [@lindblad2011high; @ward2012evidence]; 5) TFBS were identified from digital genomic footprinting of DNase hypersensitive sites in 57 cell lines [@gusev2014partitioning; @encode2012integrated]; super-enhancers were defined as the union [@finucane2015partitioning] of all super-enhancers identified in 86 human cell and tissue types [@hnisz2013super]; 6) DHS sites were defined as the union [@finucane2015partitioning] of DHSs identified across 13 cell lines [@encode2012integrated; @trynka2013chromatin].


## Code availability

All code for reproducing all analysis and figures, along with a walthrough, is available at [github.com/Wolfffff/exp_var](https://github.com/Wolfffff/exp_var).

\printbibliography
