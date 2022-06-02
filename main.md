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
      institute: eeb
      email: damelo@princeton.edu
      orcid: 0000-0002-7603-0092
      equal_contributor: "yes"
  - Kristina M. Garske:
      institute: lsi
      orcid: 0000-0002-7228-8125
  - Luisa Pallares:
      institute: fml
      orcid: 0000-0001-6547-1901
  - Julien Ayroles:
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
csl: ./cse.csl
sansfont: Skolar Sans PE TEST
bibliography: ./references.bib
abstract: Gene expression provides a basis for understanding physiology, disease, and evolution. Transcriptional profiling has revealed how differences in mean expression across groups can drive phenotypic variation in human populations. Recent work has expanded towards understanding the role variation in expression plays in shaping this phenotypic variation.  However, the precise landscape in which this variance exists remains unknown and the factors affecting variation across the landscape remain understudied. Here we show the landscape of expression variation over 20,000 samples across 60 studies and 13 tissues. Using both within study rankings of variation and a cross-study variance score,  we show that gene function, sequence variation, and molecular signatures are key regulators of gene expression variance. Our results serve both a baseline for understanding the landscape of gene expression variance and for study the drivers of the variance. We anticipate that this work will serve as starting point for understanding the factors driving variation in gene expression further push the field torwards a comprehensive understanding of gene expression variation.
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
If homogeneous selection across populations is the main driver of gene expression variance, we would expect to have consistently more or less variable genes.
If idiosyncratic selection patterns and context specific environmental interactions are more important, we could observe large differences in gene expression variance across populations.

Even within individuals, gene expression is also variable across tissues [@GTEx2017-xb].
To what extent differences in mean expression level translate to differences in expression variance is not clear.
Of course, genes that are exclusively expressed in a single cell type or tissue are necessarily more variable in that particular tissue, but differentially expressed genes could also be more variable in a particular context.
For example, stabilizing selection on gene expression could be more intense depending on the role of that gene in a particular tissue, leading to a local reduced variation and differences in variation across tissue.
Alternatively, expression variation across tissues could be tightly coupled, and in this example, selection in one tissue would lead to a reduction in variance across tissues, resulting in a consistent pattern of variation. 
Regarding tissue-specific gene expression variation, Alemu et al. [@Alemu2014-jo] used microarray data for several human tissues to investigate the epigenetic drivers of gene expression variation, showing that several epigenetic markers were indeed linked to gene expression variation, and that these were variable across tissues and between high- and low-variance genes. 
This apparent heterogeneity in regulation for high and low-variance genes is interesting expecially due to the usual focus on gene expression robustness, in the sense of reducing variation [@Siegal2014-dv; @Payne2015-wn; @Macneil2011-ax; @Denby2012-as]. For example, Siegal and Leu [@Siegal2014-dv] provide several examples of known regulatory mechanisms for reducing gene expression variance, but no examples for the maintenance of high gene expression variance.

Here, we use public gene expression data sets to evaluate how the differences in gene expression variance is structured across independent samples.
We collected and compared the gene expression variance across many studies and used the similarities across these studies to create a gene expression variance ranking, which orders genes from least variable to most variable.
We then explored the expected drivers of this gene expression ranking, showing that both cis and trans regulation are involved with the determination of gene expression variance.
Finally, we explored the link between gene expression variance and biological function by leveraging gene ontology annotations.

# Results


\begin{figure*}[t!]
    \centering
     \includegraphics[width=\dimexpr 1\textwidth]{figures/fig1.png}
    \caption{A. Correlation heatmap showing the across study Spearman correlation of standard deviations. Pairs of studies with more similar patterns of gene expression variance have higher correlations. Studies are shown in the same order as in fig. \ref{fig:corr_model}, panel A; B. Histogram of the correlations shown in the previous panel; C. Standard deviation correlation PCoA. There is no clear structuring of the studies with respect to their source, which is indicated by the colors; D. Density plot of standard deviations after z-normalization. Inset plot shows distribution of mean centered standard deviations grouped by study without normalization. The corresponding rug plots show the location of the highest ranking gene in standard deviation rank (right, blue) and lowest (left, red).}
    \label{fig:sd_corr}
\end{figure*}

Gene expression standard deviations (SDs) were calculated for each data set using a unified pipeline that normalized the mean-variance relation in count data, controlled for batch effects, and removed outliers (see methods for details).
Spearman correlations ($\rho_s$) between gene expression SDs reveals a broadly similar rank of gene expression variance, so genes that are most variable in one study tend to be most variable in all studies (fig. \ref{fig:sd_corr}A and B).
Several studies were conducted under the umbrella of two large research projects: GTEx [@GTEx2017-xb] and TCGA [@tcga2013-gx], and we note these study origins in the figures.
A principal coordinate analysis [@Gower1966-dk] using $|1 - \rho_s|$ as a distance measure does not show clearly delineated groups, but GTEx  and TCGA studies are clustered among themselves and close together (fig. \ref{fig:sd_corr}C).
This indicates some effect of study source on the similarity between gene expression SD across studies, which we explore in detail below.
Observed range of gene expression SD across genes is variable across studies, but can be normalized so that the distributions are comparable (fig. \ref{fig:sd_corr}D).
Given that the correlations across studies are broadly high, indicating similar ordering of the genes, we seek to summarize the differences in variance across genes by using a single cross-study rank, averaging the ordering across all studies.
To create this rank, we used the score of each gene in the first principal component of the Spearman correlation matrix.
This generates a ranked list of genes, with most variable genes having highest rank.
The red and blue ticks at the bottom of fig. \ref{fig:sd_corr}D show the positions on the SD distributions of the least and most variable gene in our variance rank.

\begin{figure*}[t!]
    \centering
    \includegraphics[width=\linewidth]{figures/correlationModeling.png}
    \caption{Coefficients estimates from a linear model using the among studies Spearman correlations as the response variable. These correlations are shown in fig. \ref{fig:sd_corr}A and B. In the linear model, correlations are Fisher z-transformed. Study source and tissue are added as fixed effects. Coefficient estimates are shown with 50\% and 95\% credibility intervals. Panel A: The per-study random effect captures the non-independence of the correlation values and estimates the characteristic contribution of each study to the correlation. For example: comparisons involving bone marrow (from GTEx) tend to be lower than the others. Panels B and C: Fixed effect estimates: correlations among studies that use the same tissue are higher, and correlations involving studies in the ``Misc.'' category (non GTEx and TCGA) tend to be lower, while comparison involving GTEx and TCGA are higher.}
    \label{fig:corr_model}
\end{figure*}


## What drives differences in gene expression variance?

To characterize the drivers of across study similarity, we directly modeled the correlations across studies using a mixed effect linear model [@Dias2021-wk; @Dias2021-hb].
In this model, we use individual study, sampled tissue (whether a comparison is between same tissue or different tissue), and study source (GTEx, TCGA and miscellaneous) as predictors of the pairwise correlations (see Methods).
This modeling (fig. \ref{fig:corr_model}) shows that comparisons of studies within GTEx and TCGA have on average higher values of $\rho_s$, but also that comparing studies across GTEx and TCGA also shows a mild increase in the average correlation (fig. \ref{fig:corr_model}C).
Correlation involving studies that are not from TCGA and GTEx (marked as "Misc.") are on average lower (fig. \ref{fig:corr_model}C).
Since these two sources are independent, this mild effect on the similarities could be due to the quality of the data coming from these two large projects.
Tissue also affects the similarity between gene expression SD, with studies using the same tissue being, on average, more similar (fig. \ref{fig:corr_model}B).
The largest effects on the correlations are those associated with individual studies, in particular some specific tissues, i.e., comparisons involving bone marrow (from GTEx) and study SRP057500 (which used platelets) are on average lower (fig. \ref{fig:corr_model}A).
These studies also show up further away in the PCoA plot in fig. \ref{fig:sd_corr}C.

## Does biological function explain variance in expression?

To explore the relationship between variance and function, we took the top 5% most variable and the bottom 5% least variable genes in our ranking (about ~560 genes in each group) and performed a Gene Ontology (GO) enrichment analysis within each group.
This allowed us to establish the representative functions of these consistently high and low-variance genes.
In total, using a Benjamini-Hochberg (BH) adjusted p-value threshold of $10^{-3}$, we found 59 enriched terms in the low variance genes, and 738 enriched terms in the high variance genes (see supporting table 1 for a complete listing).
Among the 5% most variable genes we observe enrichment for biological processes like immune function, response to stimulus, maintenance of homeostasis, and tissue morphogenesis (fig. \ref{fig:go_tails}A).
In line with this GO term enrichment, the top 5% most variable genes are enriched 7.7-fold for genes that encode secreted proteins, relative to all other genes ($p < 10^{-3}$).
Among the 5% least variable genes we see enrichment for housekeeping functions like mRNA processing, cell cycle regulation, methylation, histone modification, translation, transcription, and DNA repair (fig. \ref{fig:go_tails}B); and accordingly we find that previously characterized human housekeeping genes [@Hounkpe2020-yq] are enriched within the 5% least variable genes 2.0-fold relative to all other genes ($p < 10^{-3}$).
The genes exhibiting the lowest variance (lowest 5%) are also enriched for those that have been previously shown to have a high probability of being loss-of-function intolerant (pLI) [@lek2016analysis] (1.2-fold enrichment, $p < 10^{-3}$).
Genes with a high pLI have been shown to be important in housekeeping functions, and have higher mean expression values across a broad set of tissues and cell types [@lek2016analysis]. Our result that genes with low variance are enriched for both housekeeping genes and genes with high pLI is consistent with this previous report; and we further see that the mean expression of genes positively correlates with pLI (Partial Spearman correlation $\rho_s$ = 0.32, $p < 10^{-3}$), showing the opposite relationship between variance and mean expression when considering pLI.

We also explored the distribution of expression variance among the genes associated with GO terms.
For this, we gathered all biological process GO terms in level 3 (i.e. terms that are at a distance of 3 for the top of the GO hierarchy).
Using only the set of genes that are associated with at least one of these level-3 terms, we separated the genes into expression variance deciles, with the first decile having the lowest variance.
We then counted how many genes in each decile has been associated with each term.
If variance rank is not linked to the GO annotations, terms should have an equal proportion of genes in each decile.
We measured how far from this uniform allocation each term is by measuring the Shannon entropy of the proportion of genes in each decile.
Higher entropy is associated with more uniform distribution of genes across deciles.
GO terms with low entropy indicate some decile is over-represented in the genes associated with that term.
We also measured skewness for each term, which should be zero if no decile is over-represented, negative if high-variance terms are over-represented, and positive if low-variance deciles are over-represented.
Skewness by entropy for each GO term can be seen in fig. \ref{fig:skew_entropy}.
Positive-skew low-entropy terms, those enriched with low-variance genes, are associated with house keeping functions, like RNA localization, translation initiation, methylation and chromosome segregation (fig. \ref{fig:go_skewness} A).
Likewise, terms with negative skew and low entropy, enriched for high-variance genes, are related to immune response, tissue morphogenesis, chemotaxis---all dynamic biological functions related to interacting with the environment (fig. \ref{fig:go_skewness} B).

Both GO analyses suggests a strong influence of biological function in determining gene expression variance.
Genes associated with baseline fundamental functions, expected to be under strong stabilizing selection, are also low-variance; high-variance genes are associated with responding to external stimuli (i.e., tissue reorganization and immune response).

\begin{figure}
    \centering
    \includegraphics[width=\linewidth]{figures/local_go_lowerUpper.png}
    \caption{Gene set enrichment analyses testing for over representation of gene ontology categories in the upper and lower 5\% quantiles of the gene variance rank. High-variance gene are enriched for terms related to immune function, response to wounding, blood vessel morphogenesis and inflammatory response. In contrast, low-variance genes are associated with translation, control of methylation, RNA processing, chromosome separation, and other cell housekeeping functions. All displayed terms are significant with a 5\% FDR corrected p-value below $10^{-3}$.}
    \label{fig:go_tails}
\end{figure}

\begin{figure*}[t!]
    \centering
    \includegraphics[width=\linewidth]{figures/skew_entropy.png}
    \caption{Relationship between skew and entropy of rank decile distributions for each GO term. High entropy terms, to the right of the plot, are associated with a more egalitarian proportion of genes in each of the SD rank deciles. Terms on the left of the plot are associated with more genes in some particular decile. The skewness in the y-axis measures if the high- or low-variance deciles are more represented for a particular term. Terms on the positive side of the y-axis are associated with low-variation genes, and terms on the negative side of the y-axis are associated with high variation genes. The GO terms are filtered for gene counts greater than 100, as in fig. \ref{fig:go_skewness}.}
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
Nucleotide diversity is used as a proxy for cis-regulation sites, and we expect variation to increase with diversity.
Here, we find a partial Spearman's correlation of 0.184 ($p < 10^{-3}$).
Connectivity, a proxy for regulatory interactions with other genes and of selective constraints [@Mahler2017-bb], in turn, should be negatively correlated with variation, as highly connected genes are expected to be more constrained in their variability.
The resulting partial Spearman's correlation is -0.024 ($p \approx 6 \times 10^{-3}$).
Finally, we find a partial Spearman's correlation of -0.046 ($p \approx 1 \times 10^{-3}$) for the proportion of substitutions that are adaptive.
In spite of all of these associations being significant and in the expected direction, their effect sizes are very small, suggesting a weak link between these broad measures and gene expression variance.

<!--
\begin{table}[]
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

We assess how local epigenetic features relate to gene expression variance.
We use each gene, including the surrounding 10kb on both ends, to calculate the proportion of gene regions that correspond to epigenetic marks and gene annotations previously used to stratify the genome into interpretable functional categories [@finucane2015partitioning], including promoter and enhancer regions, open chromatin (assayed through DNase hypersensitivity (DHS)), and transcription factor binding sites (TFBS).
Biochemical features associated with cis gene regulation are positively correlated with the gene expression variance rank metric, regardless of whether the regulatory effect on gene expression is positive or negative [KG supp fig 1].
For example, both the proportion of gene regions made up of enhancers and repressed genomic states are positively correlated with gene expression variance (BH adjusted $p < 0.05$) [KG supp fig 1].
As expected, the proportion of gene regions made up of repressed genomic states is inversely correlated with mean expression of the gene, and that made up of enhancers is positively correlated with the mean expression of genes [KG supp fig 1].
This shows that gene expression variance is not simply associated with the same features as mean expression levels.
The magnitude of the correlation with general RefSeq gene features, such as promoter and coding sequence, is lower for both the variance and mean, and we see that this coincides with an overall positive (in the case of the mean) and negative (in the case of the variance) associations with gene density in the expanded gene regions (gene +/- 250kb) [KG supp table 1; KG supp fig 1].
Furthermore, the biochemical properties associated with promoter flanking regions, as well as transcribed states, are inversely correlated with gene expression variance, whereas they are positively correlated with the mean expression [KG supp fig 1].
Taken together, these results are compatible with gene expression variance being more associated with distal (i.e., non-promoter) gene regulation, rather than overall active transcriptional state of a gene region, as is the case with mean gene expression.

These results are largely in line with a previous assessment of human microarray data across 41 tissues to identify gene expression variance correlates with epigenetic marks on a tissue-by-tissue basis [@Alemu2014-jo], with the notable difference that the gene expression variance metric used here is a single expression variance rank for all tissues and studies assessed, and its relationship with global genomic annotations also defined across many tissues or cell types (see Methods and [@finucane2015partitioning]).
The concordance between these two sets of results is consistent with the high correlation of gene expression variance across tissues (fig. \ref{fig:sd_corr}A and B), which implies that a global view into expression variance should, for the most part, recapitulate what is seen at the tissue-specific level.
Two major considerations arise when interpreting these results.
First, there is considerable overlap between the different epigenetic marks either globally or in a given tissue, making it difficult to parse out the differential effects of individual regulatory states on gene expression variance.
For example, open chromatin sites are broadly associated with regions that are available for gene regulation and overlap with enhancers, polycomb-mediated repressive sites, and promoters [ref and analysis], among others.
Second, there has been a massive increase in epigenetic data collected in over 100 human tissues and cell types since these previous microarray and epigenetic mark data were curated and published, meaning we now have more cell-type-specific information and increased understanding of the diversity of regulatory states that can take shape within the nucleus.
To address these points, we investigate both cross-tissue and tissue-specific expression variance relationships with non-overlapping annotations of chromatin states as defined through ChromHMM [@ernst2012chromhmm].
The genome segment annotations were defined using epigenetic data collected through ENCODE [@ENCODE2012-mz] and Roadmap [@Roadmap_Epigenomics_Consortium2015-mq], either at the universal level across 127 cell and tissue types [@vu2022universal] or in each tissue independently [ref].
For eight of the tissue types assessed in the current study, we use the ChromHMM states from the corresponding tissue when available, and we use a representative cell type when the tissue itself is not available [KG supp table 2].

For the cross-tissue gene expression variance comparison with the universal chromatin states, we mostly reproduce the results obtained when using the previously curated gene regulatory feature annotations, such as the positive correlation between gene expression variance and both enhancer and polycomb-mediated repressed chromatin states; and the inverse relationship between gene expression variance and active promoters or transcribed states [KG supp fig 2].
One notable difference is that the strong positive correlation seen between gene expression variance and the union of DHS among cell types [KG supp fig 1] is not seen when using the universal chromatin state for DNase [KG supp fig 2].
This is likely due in part to the aforementioned difference between overlapping [@finucane2015partitioning] and non-overlapping [@vu2022universal] annotations, such that regions that contain both DNase hypersensitive sites and other gene regulatory epigenetic marks are defined as the chromatin state associated with the other epigenetic marks [@vu2022universal].
This suggests that the DNase state represents a distinct form of gene regulation not clearly defined through the histone marks profiled and used to define the universal chromatin states.
Indeed, Vu et al. [-@vu2022universal] find that the DNase chromatin state that is associated with DNase only across all cell types studied is most strongly enriched for CTCF-specific chromatin states.
CTCF is a transcription factor that can function as an activator, repressor, or insulator protein [Dunn2003-cu], and the diverse roles it plays in gene regulation, particularly at the universal level, likely have widespread differential effects on gene expression variance, thus leading to the lack of correlation between the DNase state and gene expression variance. __(need to develop and then refine this more - lit review and any additional analysis)__

## Do tissue-specific chromatin states associate with tissue-level gene expression variance?

We also explore the relation of tissue-specific chromatin states and SD rank, and contrast these local analysis to the global analysis outlined above.
Many of the cross-tissue correlations are recapitulated at the tissue-level, including a strong and highly consistent positive correlation between enhancer states and gene expression variance, and an inverse relationship between gene expression variance and gene transcription or ZNF states [KG supp fig 2].
Two blood associations stand out as being different from the consistent effects across the other tissue-level and cross-tissue associations.
The weak promoter state is positively correlated with gene expression variance in all comparisons except blood, reflective of a likely role of bivalent promoters in context-dependent gene expression __(? lit review of bivalent promoters, development, the role in differentiated cells)__.
Furthermore, the consistent inverse correlation of gene expression variance with weak transcription is reversed in blood [KG supp fig 2].
__(These results may have something to do with the fact that the top 5% most variable genes are enriched for GO terms related to immunity but I need to get everything into context and think of potential other analyses)__

Some notable differences exist between the universal and tissue-specific chromatin state associations with gene expression variance.
First, while the universal heterochromatin state positively correlates with cross-tissue gene expression variance, the tissue-specific heterochromatin states are inversely correlated with the tissue-level gene expression variance [KG supp fig 2].
Heterochromatin states in a given cell type should show reduced variance because there should be drastically reduced gene expression overall, as we see in the inverse correlation of heterochromatin states with mean gene expression [KG supp fig 2].
The reason for a universal heterochromatin state showing a positive association with gene expression variance remains to be determined.
The universal promoter chromatin state is inversely associated with gene expression variance, in line with our result that genes that are ubiquitously expressed and involved in housekeeping functions are enriched in the low variance genes.
However, interestingly in both adipose and liver tissues, the tissue-specific promoter state is positively correlated with the tissue-level gene expression variance [KG supp fig 2].
This could be reflective of the necessity for rapid environmental responses at key expressed genes in these metabolic tissues [ref?].
Taken together, XX...
__(look more closely at the Alemu paper supplement to see if they show similar results - I don't think they do a great job comparing the tissue patterns in their results/dicussion?)__

We ask whether the associations seen between a given tissue's gene expression variance and its corresponding chromatin state is specific to that tissue by comparing the tissue-level expression variance to the universal and other tissue chromatin states [KG supp fig 3].

Maybe do some TF binding site analyses in the different tissue-specific regions to see what might be contributing to gene expression variance at the cross-tissue vs. tissue level.

How to bring back to selection/conservation etc? ConsHMM could also be used?

## Linking expression variance and disease

To explore the link between expression variance and disease, we use the gene annotations derived from a probabilistic transcriptome-wide association study (PTWAS, @Zhang2020-cl).
Using the list of significant gene-trait pairs at 5% FDR provided by Zhang et al. [-@Zhang2020-cl], we performed a hypergeometric enrichment test for the top 5% high- and low-variance genes in our across-tissue rank and in all tissue-specific gene variance ranks.
Despite their overall high similarity, we use both across tissue and tissue specific ranks because some genes only appear in the tissue specific rank due to their limited tissue specific gene expression.
In the high-variance group, we find no enrichment in the across-tissue rank, but we do find enrichment of genes annotated for allergy, immune disease, and endocrine system disease among the high variance genes in several tissue-specific variance ranks.
Among high-variance genes in the colon rank, we see enrichment for endocrine system disease (1.77-fold, BH adjusted $p < 10^{-4}$).
Among high-variance genes in immune cell rank, we see enrichment for endocrine system disease (1.67-fold, BH adjusted $p < 10^{-3}$),
allergy (1.7-fold, BH adjusted $p < 10^{-3}$), and immune disease (1.32-fold, BH adjusted $p < 10^{-2}$).
Among high-variance genes in the thyroid rank, we see enrichment for endocrine system disease (1.9-fold, BH adjusted $p < 10^{-5}$),
allergy (1.85-fold, BH adjusted $p < 10^{-4}$), and immune disease (1.45-fold, BH adjusted $p < 10^{-4}$).
These are all quite similar and suggest a stable pattern of high-variance gene expression across these tissues, with enrichment for these three classes of diseases.
The link with immune diseases is expected given the high enrichment for immune-related genes in the high-variance group [@Hagai2018-fu].
As for the low-variance group, we found strong enrichment for genes associated with psychiatric and neurological disorders in the across-tissue rank and in some tissue specific ranks (breast, liver, and stomach; ~1.2-fold enrichment, $p < 0.05$, for all cases).
The psychiatric disease links is consistent with previous work[@Mar2011-dr] and is discussed below; however, the enrichment among the low-variance genes is weaker.

# Discussion

By using large publicly available data sets, we were able to probe the landscape of gene expression variance in several human tissues.
Differences in gene expression variance were driven by technical aspects of gene expression measurement, with data derived from large consortia showing more similar patters of variance across genes; and by tissue, with studies using the same tissues also showing higher similarities.
This would suggest that careful consideration of sample sizes and experiment design are fundamental to the study of gene expression variance, and the usual small samples of RNA-seq studies might be underpowered for the study of this particular aspect of gene expression.
However, both the effects of study origin and tissue were small, and the largest drivers of differences across studies were idiosyncratic differences related to single data sets, with tissues know to have divergent gene expression patterns (i.e. bone marrow, blood, testis, and platelets) also showing the largest differences in gene expression variance.
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
These results are consistent with previous analysis of gene expression variance on a tissue-by-tissue basis [@Alemu2014-jo].
This pattern of enrichment is also observed when we look at enrichment for high or low variance genes within the genes associated with each terms in the GO hierarchy.
Basic cell function terms are enriched for low variance genes, and terms involved in response to external stimulus are enriched for high variance genes.

While indirect, all these patterns point to a selective structuring of gene expression variance.
Stabilizing and purifying selection are consistent: genes expected to be under strong stabilizing selection, those linked with fundamental baseline biological processes, are indeed over represented in the least variable genes.
These same genes are also expected to be under strong purifying selection and to show low levels of polymorphisms, which we observe.
Likewise, genes whose function is constrained by myriad interactions with several other genes, those with high connectivity, are less variable.
Furthermore, genes involved with direct interaction to the environment, which must change their pattern of expression depending on external conditions, are expected to be more variable, and again we see a strong enrichment of genes related to interacting with the environment among the most variable.
Given this strong functional linkage between function and variance, it is not surprising that the gene variance ranking is similar across studies, allowing us to create our ranking in the first place.

One interesting aspect of the GO term analysis shown in figures \ref{fig:skew_entropy} and \ref{fig:go_skewness} is that there is no biological process term associated with enrichment for intermediate variance genes: the low-entropy terms have either positive or negative skew, never zero skew.
In other words, there is no annotated biological process for which the associated genes are kept at some intermediary level of variation.
Either there is not relation between the gene expression variance and the biological process, or there is a strong bias towards high or low-variance genes.
This suggests that selective shaping of gene expression has two modes, corresponding with (1) biological processes under strong stabilizing selection (i.e., variance reducing selection) or (2) biological processes under disruptive selection (i.e., variance increasing selection).
In short, we find strong support for the idea that there are genes with consistently more (or less) variable expression levels, and that these differences in variance are the result of different patterns of selection.

Following Alemu et al. [@Alemu2014-jo], we observe that epigenetic signatures of gene regulation, such as enhancer histone marks, make up a higher proportion of the surrounding genomic regions of genes that exhibit higher variance in expression.
In contrast, an accumulation of strong promoter elements and overall transcriptional activation is associated with genes with lower expression variance.
These results suggest the presence of distinct modes of regulation for genes with high vs. low variance.
Combined, the differences in the types of genomic regulatory features surrounding the high- and low-variance genes and their distinct functional annotations suggest different mechanisms of regulation of their gene expression variance [@Alemu2014-jo].
This heterogeneity could lead to detectable differences in selection signatures between distal regulatory elements and promoters depending on the gene expression variance.
We posid that it should be possible to go beyond the usual caracterization of strategies of gene expression robustness, and to explore mechanisms for the _robustness of plasticity_, that is, the maintenance of high levels of gene expression variation given environmental cues.

Given the broad consistency of gene expression variance, a natural question is how do these well regulated levels of variation behave in perturbed or disease conditions.
We find some suggestive links between tissue-specific variance ranks and disease, but these links need to be better explored using more specific methods.
Comparing two HapMap populations, Li et al. [-@Li2010-qs] showed that gene expression variance was similar in both populations, and that high variance genes were enriched for genes related to HIV susceptibility, consistent with our observation of enrichment for immune related genes among those with more variable expression.
In a case-control experiment, Mar et al. [-@Mar2011-dr] showed that expression variance was related to disease status in Schizophrenia and Parkinson's disease patients, with altered genes being non randomly distributed across signaling networks.
These authors also find a link between gene network connectivity and expression variance, in agreement with the effect we find using the gene expression variance rank.
Also, the pattern of variance alteration differed across diseases, with Parkinson's patients showing increased expression variance, and Schizophrenia patients showing more constrained patters of expression.
The authors hypothesizes that the reduced variance in Schizophrenia patients reduces the robustness of their gene expression networks.
This suggests several types of shifts in gene expression variation are possible, with different outcomes.
We highlight three distinct possibilities:
First, low variance genes, under strong stabilizing selection, could become more variable under stress, indicating a reduced capacity for maintaining homeostasis.
Second, high variance genes, expected to be reactive to changes in the environment, could become less variable, indicating reduced capacity to responding to external stimuli.
Third, the covariance between different genes could be altered, leading to decoherence between interdependent genes [@Lea2019-pq].
Any one of these changes in expression variance pattern could have physiological consequences, and exploring these differences should be a major part of linking gene expression to cell phenotypes and function (see Hagai et al. -@Hagai2018-fu for example).
Genes are also expected to differ in their capacity to maintain an optimal level of gene expression variance [@Macneil2011-ax].
These differences in robustness are linked to gene regulatory networks and epigenetic gene expression regulation [@Payne2015-wn; @Chalancon2012-ul].
Our results suggest that low- and high-variance genes could use different strategies in order to maintain their optimal levels of variation, and that these are the result of different patterns of selection.

## Draft

- We could use some epigenetics discussion in the introduction. Maybe some comments on the Alemu paper and references from it. 
- Classes of high and low variance genes solve a tension between robustness and plasticity. Literature has plenty of examples of strategies to acheive robustness in the sense of limiting variation [@Siegal2014-dv], but not as many mechanisms to robustness in maintaining high variation, or robustness in maintaining apropriate levels of plasticity.

\footnotesize

# Methods

## Data sources

We selected 60 RNA-seq studies with large sample sizes from public gene expression repositories recount3 [@Wilks2021-uj] and Expression Atlas [@Papatheodorou2020-dn]. Because we are interested in population level variation of gene expression, we exclude single-cell studies and focus only on studies derived from tissue samples. We only used studies for which raw read count data was available, and for which we could parse the metadata for batch effects.
We use studies to refer to independent data sets, which could have been generated by the same consortium.
For example, the GTEx data are separated by tissue, and we refer to each tissue as a separate study.
We divide out studies into three categories depending on their origin: GTEx, TCGA, and Miscellaneous.

__Processing pipeline__: We use a standardized pipeline to measure gene expression variance while removing extraneous sources of variation.
Data from case-control studies was filtered to keep only control samples.
For each study, we filtered genes that did not achieve a minimum of 1 count per million (cpm) reads in all samples and a mean 5 cpm reads.
To account for the mean variance relation in count data, remaining genes were subjected to the variance stabilizing transformation implemented in DESeq2 [@Love2014-mp].
Fixed effects were manually curated from the metadata for all studies and removed using a linear fixed effect model.
Outlier individuals in the residual distribution were removed using a robust Principal Component Analysis (PCA) approach of automatic outlier detection [@Chen2020-fy].
Gene expression standard deviation is measured as the residual standard deviation after fixed effect correction and outlier removal.

## Gene expression variance across-tissue correlation

We assessed the similarity in gene expression variance across studies by using a between-study Spearman correlation matrix of the measured SDs.
Only genes present in all studies were used to calculate the Spearman correlation matrix, 4300 genes in total.
Using Spearman correlations avoids problems related to overall scaling or coverage differences, and allows us to assess if the same genes are usually more or less variable across studies.
To investigate the factors involved in determining correlations between studies, we used a varying effects model to investigate the effect of study origin and tissue on the correlations across studies.
This model is designed to take the non-independent nature of a set of correlations into account when modeling the correlation between gene expression variance.
This is accomplished by adding a per-study random effect, see [@Dias2021-hb] for details.

__Gene expression SD rank:__ Given that most of the variation in the Spearman correlation across studies is explained by a single principal component, we use the ranked projections of gene expression SDs in this principal component (PC1) to create an across study rank of gene variation.
The higher the rank, the higher the expression SD of a given gene.
Genes that were expressed in at least 50% of the studies were included in the rank.
In order to project a particular gene onto the PC1 of the between study correlation matrix, we impute missing values using a PCA based imputation [@Husson2019-sl].
The imputation procedure has minimal effect on the ranking, and imputing missing SD ranks at the beginning or at the end of the ranks produces similar results.
We also create a tissue specific variance ranking, using the same ranking procedure but joining studies done in the same tissue type.
For this tissue level ranking, we only use genes that are expressed in all studies of a given tissue.
For tissues that are represented by a single study, we use the SD ranking for that study as the tissue rank.
We further investigate the tissue-level expression variance ranks as they relate to genomic regulation.

## Gene level statistics

__Genetic variation__: Genetic variation measures were obtained from the PopHuman project, which provides a comprehensive set of genomic information for human populations derived from the 1000 Genomes Project.
Gene level metrics were used when available.
If only window based metrics are available, we assembled gene level information from 10kb window tracks where each window that overlaps with a given gene was assigned to the gene and the mean metric value is reported.
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

## Gene function assessment

__GO term enrichment__: All gene ontology analysis were done using the clusterProfiler R package v4.2.2 [@Wu2021-db] and the Org.Hs.eg.db database package v3.14.0 [@godb].

__Secreted genes__: We use The Protein Atlas [@uhlen2015tissue] to extract information on which proteins are secreted [@uhlen2019human] and test for an enrichment of genes with secreted products in the genes within the highest and lowest 5% of gene expression variance rank. We use the hypergeometric test to assess the significance of the enrichment.

__Housekeeping genes__: Human housekeeping genes were identified as genes that are expressed with low variance in all 52 human cell and tissue types, assessed in over 10,000 samples [@Hounkpe2020-yq]. We test for an enrichment of housekeeping genes in the genes within the highest and lowest 5% of gene expression variance rank. We use the hypergeometric test to assess the significance of the enrichment.

__Probability of being loss-of-function intolerant (pLI)__: Genes that are likely haploinsufficient (i.e., intolerant of heterozygous loss-of-function variants) were detected as those with fewer than expected protein-truncating variants (PTVs) in ExAC [@Lek2016-xw]. We use genes with a pLI > 0.9 to test for the enrichment of loss-of-function intolerant genes in the genes exhibiting the highest and lowest 5% gene expression variance estimates. We use the hypergeometric test to assess the significance of the enrichment.

__Disease annotations__: We use the gene annotations for involvement with diseases provided by the supporting information Table S2 from Zhang et al. [-@Zhang2020-cl]

## Epigenetic marks and gene features

__Data used__: We first obtain various annotations previously used to stratify the genome into interpretable functional categories [@finucane2015partitioning].
A subset of these annotations are used to quantify functional and molecular correlates of the gene expression variance metric: 1) promoter, coding, and 3' and 5' UTR are annotations from the RefSeq gene model; 2) CTCF, promoter flanking, transcribed, transcription start site, and enhancer categories were defined as the union [@finucane2015partitioning] of these annotations derived from ChromHMM/Segway across 6 cell types [@hoffman2013integrative]; 3) the repressed category was defined as the intersection [@finucane2015partitioning] of these annotations derived from ChromHMM/Segway across 6 cell types [@hoffman2013integrative]; 4) conserved elements were identified across 29 mammalian species [@lindblad2011high; @ward2012evidence]; 5) TFBS were identified from digital genomic footprinting of DNase hypersensitive sites in 57 cell lines [@gusev2014partitioning; @ENCODE2012-mz]; super-enhancers were defined as the union [@finucane2015partitioning] of all super-enhancers identified in 86 human cell and tissue types [@hnisz2013super]; 6) DHS sites were defined as the union [@finucane2015partitioning] of DHSs identified across 13 cell lines [@ENCODE2012-mz; @trynka2013chromatin].

__Correlations__: We use the ppcor R package [@kim2015ppcor] v1.1 to run the pairwise partial Spearman correlations for three variables: the gene expression variance and mean ranks and the proportion of the gene regions (gene +/- 10 kb) made up of the various features described above, one at a time.
We extract the partial Spearman correlation rho and p-values for the variance and mean associations with the chromatin and gene features.
P-values are corrected using the Benjamini-Hochberg procedure and comparisons with an adjusted $p < 0.05$ are considered significant.

## Cross-tissue vs. tissue-level chromatin states

__Data used__: We use the universal [@vu2022universal] and tissue-specific [ref] ChromHMM [@ernst2012chromhmm] chromatin states to compare the non-overlapping genome segmentation to cross-tissue and tissue-level gene expression variance metrics.

__Correlations__: Correlations were performed in the same manner as the global assessment (above) and corrected (BH) for all tests (all tissue-level (n = 8 tissues) gene expression variance ranks plus the cross-tissue ranks correlated with all 13 chromatin state categories in each tissue plus the universal annotation (n = 1,053 tests)).

# Code availability

All code for reproducing all analysis and figures, along with a walk-through, is available at [github.com/Wolfffff/exp_var](https://github.com/Wolfffff/exp_var).

# Supporting information

S1 - **Metadata file describing the data used in the study as well as some intermediate processing information.**

S1 - Table

S1 - Column Descriptions

S2 - **Combined table describing gene ontology enrichment in the top 5% and bottom 5% of genes as ranked by variance**

S2 - Table

S2 - Column Descriptions


# References

<!-- %% \printbibliography -->
