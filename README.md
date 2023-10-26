# ImmuMicrobiome
This package is based on the code of Dr. Villette and Dr. Larsen. The package is written and maintained by Dr. Villette. This package is meant to facilitate microbiome exploration and ensuring nice plotting.\
This package covers :

-   The dada2 pipeline with wrapper functions that ease the processing of multiple projects

-   Some plotting functions for beta diversity, heatmap and differential abundance analysis using directly a phyloseq object

-   A pipeline for IgASeq analysis

Most of the functions in this package were used in the following papers : ...
---
output:
  html_document:   
---


# Introduction

<div style="text-align: justify">
This package is based on the code of Dr. Villette and Dr. Larsen. The package is written and maintained by Dr. Villette. This package is meant to facilitate microbiome exploration and ensuring nice plotting.\
This package covers :\

  -   The DADA2 pipeline with wrapper functions that ease the processing of multiple projects

  -   Some plotting functions for beta diversity, heatmap and differential abundance analysis using directly a phyloseq object

  -   A pipeline for **IgASeq** analysis

</div>

## DADA2 pipeline
<div style="text-align: justify">
I strongly encourage you to read the paper from Callahan and also the [tutorial](https://benjjneb.github.io/dada2/tutorial.html). DADA2 (and other pipeline) introduced the concept of Exact Sequence Variant which are very different from the OTU approach.\
Quickly, ESV (ASV for DADA2) considers that every reads with a single base difference is different from another by definition. In this approach the authors came with the idea that until proven otherwise a read is different from another even when only a single base is different.\
To decipher if the read is truly different or not, the pipeline will use the Phred score associated to each base to determine if the base read here is different between a read A and B because of:  

  -   A biological difference, for example to different strains  
  
  -   A technical difference, a misreading from the sequencer


</div>

## Nice plotting.
<div style="text-align: justify">

Microbiomics is a very nice playground for biostatistic and data modelisation. Based on the work done in my lab, I've created a bunch of wrapper functions and new functions to explore microbiomics. A part of the functions are working directly on [phyloseq](https://joey711.github.io/phyloseq/) objects and the others accept directly data.frames.
</div>

## IgASeq analysis.
<div style="text-align: justify">

IgA and more importantly secretory IgA gain a growing interest as it is the main immunological microbial regulator. This immunoglobulin has the particular ability to control both pathogens and commensals, with a negative action on pathogens (pathogen removal) and a positive action on commensals (ecological niche formation).  
In this package we describe a new method to produce the samples and also to make the analysis after the sequencing. This method is actually published in [...].

</div>
## Machine learning
<div style="text-align: justify">
While we are not expert on machine learning and statistical models, more and more authors use them to predict or classify samples according to variables and features, we decided to support machine learning steps in this package. In this package we made wrapper functions to ease the use of models. 

This part is not yet furnished and will available in future updates.

</div>
