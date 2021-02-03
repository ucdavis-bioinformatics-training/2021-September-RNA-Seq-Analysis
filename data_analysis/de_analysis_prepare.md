
# Create a new RStudio project

Open RStudio and create a new project, for more info see (Using-Projects)[https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects]

* File > New Project > New Directory > New Project (name the new directory, Ex. mRNA_Seq_Workshop) and check "use packrat with this project", or "use renv with this project" if your using the devel version.

Learn more about [renv](https://rstudio.github.io/renv/articles/renv.html)

Learn more about [packrat](https://rstudio.github.io/packrat/)

## Install the needed R packages

Set some options and make sure the packages edgeR, gplots, RColorBrewer, topGO, KEGGREST, Rgraphviz and org.Hs.eg.db are installed (if not install it), and then load

In the R console run the following commands
```r
if (!any(rownames(installed.packages()) == "edgeR")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("edgeR")
}
library(edgeR)

if (!any(rownames(installed.packages()) == "topGO")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("topGO")
}
library(topGO)

if (!any(rownames(installed.packages()) == "KEGGREST")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}
library(KEGGREST)

if (!any(rownames(installed.packages()) == "Rgraphviz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Rgraphviz")
}
library(Rgraphviz)

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

library(edgeR)
if (!any(rownames(installed.packages()) == "gplots")){
install.packages("gplots")
}
library(gplots)

if (!any(rownames(installed.packages()) == "RColorBrewer")){
install.packages("RColorBrewer")
}
library(RColorBrewer)
```

## Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/data_analysis/DE_Analysis.Rmd", "DE_Analysis.Rmd")
```

## Download the data file for the workshop document and preview/open it

This is the the counts file generated after running [Generating Summarized Counts](https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/counts).

I've also uploaded to the github repo. In the R console run the following command.
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/datasets/rnaseq_workshop_counts.txt", "rnaseq_workshop_counts.txt")
```

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/datasets/ensembl_hg_100.tsv", "ensembl_hg_100.tsv")
```

### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Data_in_R"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>
