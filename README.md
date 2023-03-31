[![DOI](https://zenodo.org/badge/470977869.svg)](https://zenodo.org/badge/latestdoi/470977869)


# Manuscript

#### **Splicing accuracy varies across human introns, tissues and age**.

*Sonia Garcia-Ruiz, David Zhang, Emil K Gustavsson, Guillermo Rocamora-Perez, Melissa Grant-Peters, Aine Fairbrother-Browne, Regina H Reynolds, Jonathan W Brenton, Ana L Gil-Martinez, Zhongbo Chen, Donald C Rio, Juan A Botia, Sebastian Guelfi, Leonardo Collado-Torres, Mina Ryten*

bioRxiv 2023.03.29.534370;
doi: [https://doi.org/10.1101/2023.03.29.534370](https://doi.org/10.1101/2023.03.29.534370)


# Repository
This 'splicing-accuracy-manuscript' repository contains the code used to generate all the analyses and figures produced for the manuscript:

## Repository Structure

### 1. *"Splicing"* Database Generation

To produce the *"Splicing"* intron database, please follow the file pipeline indicated below:

1. *init.R*. Main file. It starts the exon-exon junction download from recount3. Then, it annotates the split reads in three different categories: annotated intron, novel donor and novel acceptor junction. Once the split reads are annotated, it pairs the split reads from the annotated category with the novel donor and novel acceptor junctions across the samples of each tissue. Finally, after some intermediary calls to multiple auxiliary functions, it starts the SQL commands and creates the "Splicing" intron database.
2. *database_junction_pairing.R*. Auxiliary file. It contains the code to do the junction pairing between the annotated introns and the novel donor and novel acceptor junctions.
3. *database_SQL_helper.R*. Auxiliary file. It contains the SQL helper code to assist in the creation of the "Splicing" intron database.
4. *database_SQL_generation.R*. Auxiliary file. It contains the main SQL code to create the "Splicing" database.

### 2. *"Age Stratification"* Database Generation

To produce the *"Age Stratification"* intron database, please follow the file pipeline indicated below:

1. *database_age_stratification_generation.R*. Main file. It starts the GTEx v8 sample clustering by age within the age groups "20-39", "40-59" and "60-79" years-old. Then, using the previously downloaded and QC'ed exon-exon junction data used for the creation of the Splicing database, this script groups the exon-exon split reads and extracts read counts metrics across the samples of each age category. Then, it pairs the split reads from the annotated category with the split reads from the novel donor and novel acceptor junctions across the samples of each age cluster at the tissue level. Finally, after some intermediary calls to multiple auxiliary functions, it starts the SQL commands and creates the "Age Stratification" intron database.
2. *database_junction_pairing.R*. Auxiliary file. It contains the code to do the junction pairing between the annotated introns and the novel donor and novel acceptor junctions.
3. *database_SQL_helper.R*. Auxiliary file. It contains the SQL helper code to assist in the creation of the "Splicing" intron database.
4. *database_SQL_generation.R*. Auxiliary file. It contains the main SQL code to create the "Splicing" database.

### 3. RBP expression changes with age

To evaluate whether the expression levels of RNA-binding proteins involved in post-transcriptional processes [Van Nostrand et at. 2020](https://www.nature.com/articles/s41586-020-2077-3) are affected by the age of the sample donor across BRAIN tissues, please use the code provided within the R script: *rbp_expression.R*.

### 4. Paper figure generation

To produce the main and supplementary figures and tables supplied with the manuscript *"Splicing accuracy varies across introns, tissues and age"*, please use the functions provided within the R script *"splicing_accuracy_manuscript_figures.R"* and *"splicing_accuracy_manuscript_age_figures.R"*.

### 5. Supplementary Tables

All supplementary tables can be accessed through Zenodo using the DOI: [10.5281/zenodo.7732872](https://zenodo.org/record/7732872)


## Environments

The code included within this repository has been successfully tested on:
* Ubuntu version "16.04.7 LTS (Xenial Xerus)"
* Ubuntu version "22.04.2 LTS (Jammy Jellyfish)"

## Acknowledgments

* All [RytenLab](https://rytenlab.com/) members
* [Aligning Science Across Parkinson’s (ASAP)](https://parkinsonsroadmap.org/#)
* [UCL Great Ormond Street Institute Of Child Health](https://www.ucl.ac.uk/child-health/great-ormond-street-institute-child-health-0)


