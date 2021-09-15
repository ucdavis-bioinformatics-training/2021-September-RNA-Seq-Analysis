
# Create a new RStudio project

Open RStudio and create a new project, for more info see [Using-Projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects)

* File > New Project > New Directory > New Project (name the new directory, Ex. mRNA_Seq_Workshop).

Learn more about [renv](https://rstudio.github.io/renv/articles/renv.html)

## Install the needed R packages

Set some options and make sure the packages edgeR, gplots, RColorBrewer, topGO, KEGGREST, Rgraphviz and org.Mm.eg.db are installed (if not install it), and then load

In the R console run the following commands
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!any(rownames(installed.packages()) == "edgeR")){
  BiocManager::install("edgeR")
}
library(edgeR)

if (!any(rownames(installed.packages()) == "topGO")){
  BiocManager::install("topGO")
}
library(topGO)

if (!any(rownames(installed.packages()) == "KEGGREST")){
  BiocManager::install("KEGGREST")
}
library(KEGGREST)

if (!any(rownames(installed.packages()) == "Rgraphviz")){
  BiocManager::install("Rgraphviz")
}
library(Rgraphviz)

if (!any(rownames(installed.packages()) == "org.Mm.eg.db")){
  BiocManager::install("org.Mm.eg.db")
}
library(org.Mm.eg.db)

if (!any(rownames(installed.packages()) == "gplots")){
    BiocManager::install("gplots")
}
library(gplots)

if (!any(rownames(installed.packages()) == "RColorBrewer")){
    BiocManager::install("RColorBrewer")
}
library(RColorBrewer)
```

## Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-September-RNA-Seq-Analysis/master/data_analysis/DE_Analysis_mm.Rmd", "DE_Analysis_mm.Rmd")
```

## Download the data file for the workshop document and preview/open it

This is the the counts file generated after running [Generating counts tables](https://ucdavis-bioinformatics-training.github.io/2021-September-RNA-Seq-Analysis/data_reduction/counts).

I've also uploaded to the github repo. In the R console run the following command.
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-September-RNA-Seq-Analysis/master/datasets/rnaseq_workshop_counts.txt", "rnaseq_workshop_counts.txt")
```

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-September-RNA-Seq-Analysis/master/datasets/ensembl_mm_104.tsv", "ensembl_mm_104.tsv")
```

#### For the salmon datasets

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-September-RNA-Seq-Analysis/master/datasets/rnaseq_salmon_workshop_counts.txt", "rnaseq_salmon_workshop_counts.txt")
```

### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "RNAseq Data Analysis in R"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>
