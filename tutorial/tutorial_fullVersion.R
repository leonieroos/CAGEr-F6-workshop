---
title: "1. CAGE Data Analysis with CAGEr"
output:
  html_notebook: default
df_print: paged
number_sections: yes
toc: yes
---

<br>
In this part of the workshop we will go through a sample workflow with the [CAGEr](http://bioconductor.org/packages/release/bioc/html/CAGEr.html) R package (Haberle _et al._, 2015). Two samples were chosen at different stages in zebrafish development for this practical: 512 cell stage (~2,7 hours hpf) and Prim 6 (~24 hpf). These stages have been previously published by Nepal _et al_, 2013 and are also available in one of the R-packages that can be found [here](http://promshift.genereg.net/CAGEr/PackageSource/). 
<br>
<br>
We will only go through all the standard functions and analyses of _CAGEr_ in this practical. In the following practicals of this workshop we will explore the data further.

## Summary and goals of this practical
<br>
  _CAGEr Analysis_

* Import publicly available CAGE-seq data
* Creating a CAGEset object 
* Normalization
* Create tag clusters of CTSSs
* Assess promoter width
* Create consensus clusters among samples
* Global expression patterns
* Promoter shifting


## 1.1 Loading in data from the R-package
<br>
  There are various publicly available CAGE-seq data that can be easily accessed and uploaded into R with _CAGEr_. See the vignette (section 5) for all the publicly available data and ways how to import them in _CAGEr_ including:
  
  * ___FANTOM 5:___ 
Human and mouse primary cells, cell lines and tissue

* ___FANTOM 3 and 4:___
Human and mouse tissues as well as several timecourses

* ___ENCODE:___ 
Common human cell lines

* ___Zebrafish Developmental Timecourse:___
Twelve developmental stages

<br>
  As we will be using the two stages (512 cells and Prim6) of zebrafish development, let's start by loading in the R-package for the _zebrafish developmental timecourse_ as well as _CAGEr_:

```{r, message = FALSE, warning = FALSE}
# packages
require(ZebrafishDevelopmentalCAGE)
require(CAGEr)
```
We now can display the samples in this package. First by loading data from the zebrafish samples that will return a dataframe with information about the twelve samples:
```{r}
# loading data
data(ZebrafishSamples)
head(ZebrafishSamples) 

# samples
samples <- as.character(ZebrafishSamples$sample) # to see all the samples
samples
```
We will pick just two developmental stages: "zf_512cells" and "zf_prim6". In the vignette (and manual) the way of loading in these data is by the function `importPublicData()`. 
We need to specify the _source_, _dataset_, _group_, and _samples_. The source is "ZebrafishDevelopment" and the rest are found in the dataframe of the `ZebrafishSamples`: ZebrafishCAGE, development, samplenames.
The avoid spelling everything out we'll use the character vector of samples:
  
  ```{r}
# creating a CAGEset object:
myCAGEset <- importPublicData(source = "ZebrafishDevelopment", 
                              dataset = "ZebrafishCAGE", 
                              group = "development", 
                              sample = samples[c(4,11)]  )
```
<br>
  What does this look like? Let's type `myCAGEset` in R:
```{r}
myCAGEset
```

As you can see the S4 object contains slots to store the data as we go along that is grouped in the summary under the sections:

* Input data information

* CTSS information

* Tag cluster (TC) information

* Consensus cluster information

* Number of consensus clusters

* Expression profiling

* Promoter shifting

> Double check the input data information section for the reference genome and the sample labels! 

At this stage we only have Cage TSS (CTSS) information. This is already the form of the data of this package. Only the (actual) tag count for each CTSS for both samples is stored here. The slot of Normalized tpm is (still) empty as is the rest of the object. <br> 
If for example bam files are used, you would have to determine these CTSS still (see example code below and _CAGEr_ vignette). 


```{r, eval = FALSE}
myCAGEset = new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7", 
inputFiles = file.path(..), inputFilesType = "bam",
sampleLabels = c("zf_512cells","zf_prim6")

getCTSS(myCAGEset)
```

## 1.2 Functions of CAGEr
With the following few standard _CAGEr_ functions we will "fill in" all these slots and check some general features of the two samples.
<br>

### Correlation plots
How similar are the two samples? _CAGEr_ has a function to look at the correlation between the samples and plots a png in the working directory: `plotCorrelation()`. 

```{r, message=FALSE}
plotCorrelation(myCAGEset, samples = "all", method = "pearson") 
```
<img src="images/CTSS_raw_values_pairwise_correlation.png" width="350">
<br>
<br>
To access the correlations displayed in the terminal, just assign the function to `corr.m` to keep a matrix of correlations:

```{r}
corr.m <- plotCorrelation(myCAGEset, samples = "all", method = "pearson") 
```
<br>

### Normalization
At this point we have the actual tag counts per CTSS. Do we have different library sizes for these samples? To check our library sizes use `librarySizes()` on our cageset.

```{r, message = FALSE}
# overview library sizes:
librarySizes(myCAGEset)
```
To make samples comparable, we will need to normalize our data. _CAGEr_ offers both a _tags per million normalization_ and a _power-law based normalization_. Because many CAGE-seq data follow a power-law distrubtion (Balwierz _et al._, 2009), we'll use the power-law based normalization here.
<br>
  First we need to plot the reverse cumulatives with `plotReverseCumulatives()`. On the log-log scale this reverse cumulative distribution will manifest as a monotonically decreasing linear function of which the slope (alpha) is determined per sample.
```{r, message = FALSE}
# normalisation first plot to assess alpha
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE) 
```
From this we can see the alpha for each sample. Select the alpha for the normalization step, the `normalizeTagCount()` function when method = "powerlaw". <br> <br>
  What is the average alpha of the two samples? Is that the "Ref distribution" in the plot?

```{r, message = FALSE}
normalizeTagCount(myCAGEset, method = "powerLaw",fitInRange = c(5, 1000), alpha = 1.20, T = 1*10^6)
```

Let's check our object again, what is added and in which slot?
```{r, eval = FALSE}
myCAGEset
```
<br>

### Tag clusters of CTSSs
CTSS' in close proximity of each other give rise to functionally similar set of transcripts within the same promoter elements. Thus, they are basically larger transcriptional units that correspond to individual promoters. 
<br>
  <br>
  To capture these, tag clusters (TC) of CTSS' can be easily defined within _CAGEr_ with the function: `clusterCTSS()`. Run the code below to see the information in the viewer
```{r, eval=FALSE}
?clusterCTSS
```

We will use a simple distance measurement (method = "distclu"). This means that the true distance in bp between individual CTSS is used. Let's set this at 20 bp (maxDist = 20), so that  individual CTSS can not be more than 20 bp apart to form the TC. 
<br>
  We also set the threshold at 1, which means that each CTSS should have 2 or more tag counts in all the samples prior to clustering so as to exclude low fidelity CTSSs. 
<br>
  Finally, no TC are to be called comprised of a single CTSS if the normalized signal is below 5:
  ```{r, message = FALSE}
clusterCTSS(object = myCAGEset, 
            threshold = 1, 
            thresholdIsTpm = TRUE, 
            nrPassThreshold = 1, 
            method = "distclu", 
            maxDist = 20, 
            removeSingletons = TRUE, 
            keepSingletonsAbove = 5)
```
> Keep in mind that these are determined per sample

Have a look again at the myCAGEset! The Tag cluster information slot has now been filled in.
<br>
  
  ### Promoter width
  Another feature we can determine with _CAGEr_ is the promoter width. That is in this case, the width of each tag cluster per sample. _CAGEr_ has a function that calculates the TC width based on the cumulative distribution of tag signal per TC. See the figure below.
![cumulative promoter](images/PromoterWidth.png)
<br>
  
  Generally, the width of a TC (promoter) is set between the quantiles 0.1 and 0.9 to capture 80% of CTSS within the TC. To determine the width we follow these two steps:
  
  1. Calculate the cumulative distribution per TC (`cumulativeCTSSdistribution()`)
2. Determine the quantile positions for each TC (`quantilePositions()`) 


```{r, message = FALSE}
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, 
                  clusters = "tagClusters", 
                  qLow = 0.1, 
                  qUp = 0.9)

```
<br>
  We can easily plot a histogram of the interquantile widths per sample with _CAGEr_ too:
  ```{r}
plotInterquantileWidth(myCAGEset, 
                       clusters = "tagClusters",
                       tpmThreshold = 3, 
                       qLow = 0.1, 
                       qUp = 0.9)
```
If you look at these two plots, you can see that on avarage, the 512 cell stage has more narrow TCs than prim6 (which is what we would expect). This is great to see general differences in promoter width between samples.
<br>
  
  ### Consensus Clusters
  Until now, we have determined the TC and the width of those per individual sample. However, TCs are often sample-specific. This could mean that they are present in one sample but not in the other. More so, TCs do not coincide perfectly within the same region but overlap each other as well as the possibility of more than one TC in one sample and a less amount of TC in the other but larger in width (we'll visualize an example of this later today). <br> 
                                                                                                                                                                                                                                                                                                                                                                                                                   When we want to compare transcriptional activity across samples, regions that encompass the TCs across the samples are needed. These are called ___consensus clusters___ that are aggregates of TCs from all samples. _CAGEr_ has a function for this too: `aggregateTagClusters()`. 
                                                                                                                                                                                                                                                                                                                                                                                                                   <br> <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   These we will use and extract these later during the workshop. For now we will add these to our myCAGEset. We specify that the maximum distance is a 100 bp between the TCs (with the coordinates of quantile 0.1 and quantile 0.9) across the samples.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   aggregateTagClusters(myCAGEset, 
                                                                                                                                                                                                                                                                                                                                                                                                                   tpmThreshold = 5, 
                                                                                                                                                                                                                                                                                                                                                                                                                   qLow = 0.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                   qUp = 0.9, 
                                                                                                                                                                                                                                                                                                                                                                                                                   maxDist = 100)
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ### Expression - Self organising maps
                                                                                                                                                                                                                                                                                                                                                                                                                   CAGE signals are essentially measuring transcription starting from CTSS and thus the transcription levels can be assessed. This can be done per CTSSs, TCs, or consensus clusters. <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   _CAGEr_ offers two methods to cluster gene expression to identify gene expression patterns:
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   - k-means clustering
                                                                                                                                                                                                                                                                                                                                                                                                                   - self-organising maps (SOM)
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Both require the number of expression clusters to be known _a priori_. <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   Here we will perform the SOM algorithm at consensus cluster level for our two samples to identify expression patterns. We set the threshold to at least 10 tpm (normalized) in at least one sample to make sure we have robustly expressed TCs within the consensus clusters. <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   The function is `getExpressionProfiles()` and in here we define that we expect in this case 6 clusters (`xDim = 3, yDim = 2`).
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   getExpressionProfiles(myCAGEset, 
                                                                                                                                                                                                                                                                                                                                                                                                                   what = "consensusClusters", 
                                                                                                                                                                                                                                                                                                                                                                                                                   tpmThreshold = 10,
                                                                                                                                                                                                                                                                                                                                                                                                                   nrPassThreshold = 1, 
                                                                                                                                                                                                                                                                                                                                                                                                                   method = "som", 
                                                                                                                                                                                                                                                                                                                                                                                                                   xDim = 3, 
                                                                                                                                                                                                                                                                                                                                                                                                                   yDim = 2)
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   What happens if you change the tpmThreshold? Or change the xDim and yDim?
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   <br><br>
                                                                                                                                                                                                                                                                                                                                                                                                                   These results can be plotted:
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   plotExpressionProfiles(myCAGEset, what = "consensusClusters")
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   ![SOM](images/expression_profiles.png)
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Subsequently, a character vector of the expression classes can be extracted from the cageset. This will be the same order as the consensus clusters.
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   expr.clus = expressionClasses(myCAGEset, what = "consensusClusters")
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   <br>
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ### Promoter shifting
                                                                                                                                                                                                                                                                                                                                                                                                                   Promoter shifting is described by Haberle _et al_, 2014. They've shown that the same promoter can be used differently in samples at different times in development. The method from this paper has been implemented in _CAGEr_. Shifting can be detected between two individual samples. If there are more than two samples, they can be merged per group (if possible) and the two groups are then compared. 
                                                                                                                                                                                                                                                                                                                                                                                                                   <br>
                                                                                                                                                                                                                                                                                                                                                                                                                     <br>
                                                                                                                                                                                                                                                                                                                                                                                                                     The shifting uses, similary to determining promoter width, the cumulative distribution per sample of CAGE signal (here: within the consensus clusters). The difference between the cumulative distributions is calculated as a _shifting score_ (see figure below). Additionally, a Kolmogorov-Smirnov test on the cumulative distributions is performed to give a general assessment of differential TSS usage.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   <img src="../images/shiftingScore.png" width="550">
                                                                                                                                                                                                                                                                                                                                                                                                                     <br>
                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                     We will start by determing the cumulative distributions along all consensus clusters:
                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                     ```{r, message = FALSE}
                                                                                                                                                                                                                                                                                                                                                                                                                   cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Determine the shifting score with `scoreShift()`. The variables include our myCAGEset, as well as the sample labels of the samples we want to compare (`groupX, groupY`), and use the normalized tag counts by `useTpmKS = TRUE`.
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r, message = FALSE}
                                                                                                                                                                                                                                                                                                                                                                                                                   scoreShift(myCAGEset, 
                                                                                                                                                                                                                                                                                                                                                                                                                              groupX = samples[4], 
                                                                                                                                                                                                                                                                                                                                                                                                                              groupY = samples[11], 
                                                                                                                                                                                                                                                                                                                                                                                                                              testKS = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                              useTpmKS = TRUE)
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   __Positive values__ are interpreted as the proportion of transcription initiation in the sample with lower expression that is happening "outside" (either upstream or downstream) of the region used for transcription initiation in the other sample.
                                                                                                                                                                                                                                                                                                                                                                                                                   Whereas __negative values__ indicate no physical separation, the region used for transcription initiation in the sample with lower expression is completely contained within the region of other sample.
                                                                                                                                                                                                                                                                                                                                                                                                                   <br> <br>
                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                     The two sample Kolmogorov-Smirnov (K-S) test determines if the two probability distributions are different. The cumulative sums in both samples are scaled to range
                                                                                                                                                                                                                                                                                                                                                                                                                   between 0 and 1. In our example by setting _testKS = TRUE_, the normalized tpm values were used (if FALSE, raw tag counts are used). The p-values are already corrected for multiple testing (Benjamini and Hochenberg) and a false discovery rate (FDR) is determined for each p-value. 
                                                                                                                                                                                                                                                                                                                                                                                                                   <br> <br>
                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                     These values are all again stored in our object. The retrieve only the shifting promoters as a dataframe we can set a few tresholds:
                                                                                                                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                                                                                                     - `tpmThreshold`
                                                                                                                                                                                                                                                                                                                                                                                                                   - `scoreThreshold` (shifting score)
                                                                                                                                                                                                                                                                                                                                                                                                                   - `fdrThreshold`
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Here we will set the threshold of shifting score at 0.6 like the original paper as well as setting the FDR at 0.01:
                                                                                                                                                                                                                                                                                                                                                                                                                     ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   shifting.promoters <- getShiftingPromoters(myCAGEset, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                              tpmThreshold = 5, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                              scoreThreshold = 0.6, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                              fdrThreshold = 0.01)
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   What is stored in this dataframe?
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   head(shifting.promoters)
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   This returns genomic coordinates (of conensusclusters), shifting score, adjusted p-value, and FDR p-value of the consensus clusters as well as the value of CAGE signal and position
                                                                                                                                                                                                                                                                                                                                                                                                                   of the dominant TSS in the two compared (groups of) samples.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ## 1.3 Saving the object 
                                                                                                                                                                                                                                                                                                                                                                                                                   As the CAGEset object nicely stores all the information, we can easily save the one object as an intermediate file (.RData) for revisiting it later on.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   save(myCAGEset, file ="../Data/intermediate/CAGEobject_twoSamples_PowNom_allSlots.RData")
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   ## 1.4 Empty R environment
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r}
                                                                                                                                                                                                                                                                                                                                                                                                                   rm(list = ls())
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   # References
                                                                                                                                                                                                                                                                                                                                                                                                                   Balwierz PJ, Carninci P, Daub CO, et al. Methods for analyzing deep sequencing expression data: constructing the human and mouse promoterome with deepCAGE data. Genome Biology. 2009;10(7):R79. doi:10.1186/gb-2009-10-7-r79.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Haberle V, Li N, Hadzhiev Y, et al. Two independent transcription initiation codes overlap on vertebrate core promoters. Nature. 2014;507(7492):381-385. doi:10.1038/nature12974.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Haberle V, Forrest ARR, Hayashizaki Y, Carninci P, Lenhard B. CAGEr: precise TSS data retrieval and high-resolution promoterome mining for integrative analyses. Nucleic Acids Research. 2015;43(8):e51. doi:10.1093/nar/gkv054.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   Nepal C, Hadzhiev Y, Previti C, et al. Dynamic regulation of the transcription initiation landscape at single nucleotide resolution during vertebrate embryogenesis. Genome Research. 2013;23(11):1938-1950. doi:10.1101/gr.153692.112.
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   # Session Info
                                                                                                                                                                                                                                                                                                                                                                                                                   ```{r, echo=FALSE, error=FALSE, warning=FALSE}
                                                                                                                                                                                                                                                                                                                                                                                                                   sessionInfo()
                                                                                                                                                                                                                                                                                                                                                                                                                   ```
                                                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                                   