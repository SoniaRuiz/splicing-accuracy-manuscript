[![DOI](https://zenodo.org/badge/470977869.svg)](https://zenodo.org/badge/latestdoi/470977869)

## Splicing accuracy varies across human introns, tissues, age and disease

*Sonia Garcia-Ruiz, [David Zhang](https://github.com/dzhang32), [Emil K Gustavsson](https://github.com/egustavsson), Guillermo Rocamora-Perez, [Melissa Grant-Peters](https://github.com/mgrantpeters), Aine Fairbrother-Browne, Regina H Reynolds, Jonathan W Brenton, Ana L Gil-Martinez, Zhongbo Chen, Donald C Rio, Juan A Botia, Sebastian Guelfi, [Leonardo Collado-Torres](https://lcolladotor.github.io/), Mina Ryten*

bioRxiv 2023.03.29.534370;
doi: [https://doi.org/10.1101/2023.03.29.534370](https://doi.org/10.1101/2023.03.29.534370)


## Repository 

This 'splicing-accuracy-manuscript' repository contains the code used to generate all the analyses and figures produced for the manuscript [**Splicing accuracy varies across human introns, tissues and age**](https://doi.org/10.1101/2023.03.29.534370).


### 1. *"Splicing"* Database - Paper figure generation

To produce the main and supplementary figures and tables made using the *'Splicing'* database supplied for the manuscript *"Splicing accuracy varies across introns, tissues and age"*, please use the functions provided within the R script *"splicing_accuracy_manuscript_figures.R"*. 


### 2. *"Age Stratification"* Database - Paper figure generation

To produce the main and supplementary figures and tables made using the *'Age Stratification'* database supplied for the manuscript *"Splicing accuracy varies across introns, tissues and age"*, please use the functions provided within the R script *"splicing_accuracy_manuscript_age_figures.R"*.

To evaluate whether the expression levels of RNA-binding proteins involved in post-transcriptional processes [Van Nostrand et at. 2020](https://www.nature.com/articles/s41586-020-2077-3) are affected by the age of the sample donor across BRAIN tissues, please use the code provided within the R script: *"splicing_accuracy_manuscript_age_figures.R"* (function 'age_stratification_RBPs_affected_age').


### 3. *"ENCODE shRNA"* Database - Paper figure generation

To produce the main and supplementary figures and tables made using the *'ENCODE shRNA'* database supplied for the manuscript *"Splicing accuracy varies across introns, tissues and age"*, please use the functions provided within the R script *"splicing_accuracy_manuscript_ENCODE_figures.R"*.


### 4. Supplementary Tables

All supplementary tables can be accessed through Zenodo using the DOI: [10.5281/zenodo.7732872](https://zenodo.org/record/7732872)


## Environments

The code included within this repository has been successfully tested on:
* Ubuntu version "16.04.7 LTS (Xenial Xerus)"
* Ubuntu version "22.04.2 LTS (Jammy Jellyfish)"

## Acknowledgments

* All [RytenLab](https://rytenlab.com/) members
* [Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/#)
* [UCL Great Ormond Street Institute Of Child Health](https://www.ucl.ac.uk/child-health/great-ormond-street-institute-child-health-0)

![Aligning Science Across Parkinson's (ASAP)](https://parkinsonsroadmap.org/wp-content/uploads/2020/10/cropped-ASAP_Logo_FullColor.png)
