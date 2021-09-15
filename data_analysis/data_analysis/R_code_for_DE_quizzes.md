R code for selected quizzes/questions from the DE Analysis section

## Quiz 1

How many genes are in the counts table?
```
nrow(counts)
```

How many samples are in the counts table?
```
ncol(counts)
```

What is the total count across all genes for sample mouse_110_WT_C?
```
colSums(counts)["sample mouse_110_WT_C"]
```

How many genes have a count of 0 in every sample?
sum(rowSums(counts) == 0)

## Quiz 2
Which sample has the largest normalization factor?
```
rownames(d0$samples)[which.max(d0$samples$norm.factors)]
```

Is the sample with the largest normalization factor the sample with the smallest total counts?
```
rownames(d0$samples)[which.min(d0$samples$lib.size)] == rownames(d0$samples)[which.max(d0$samples$norm.factors)]
```

Make an MDS plot of the unfiltered data.  How does it differ from the MDS plot of the filtered data?
```
plotMDS(d0, col = as.numeric(metadata$group), cex=1)
```

## Quiz 3
Based on the above model, how many genes are significantly differentially expressed between WT C and WT NC?
```
length(which(top.table$adj.P.Val < 0.05))
```

Based on the above model, and without taking significance into account, how many genes have higher expression in WT C than in WT NC?
```
length(which(top.table$logFC > 0))
```

How many genes have an unadjusted p-value less than 0.05 for the comparison of WT C to WT NC in the above model?
```
length(which(top.table$P.Value < 0.05))
```

What is the adjusted p-value for the last gene with unadjusted P < 0.05?
```
top.table$adj.P.Val[max(which(top.table$P.Value < 0.05))]
```

## Quiz 4
For the model ~0 + group + mouse, how many genes are differentially expressed between WT NC and KOMIR150 NC?
```
mm.1 <- model.matrix(~0 + group + mouse)
y.1 <- voom(d, mm.1); fit.1 <- lmFit(y.1, mm.1)
tmp.1 <- contrasts.fit(fit.1, contrasts = makeContrasts(groupWT.NC - groupKOMIR150.NC, levels = colnames(coef(fit.1))))
tmp.1 <- eBayes(tmp.1)
length(which(topTable(tmp.1, n = Inf, sort.by = "P")$adj.P.Val < 0.05))
```
