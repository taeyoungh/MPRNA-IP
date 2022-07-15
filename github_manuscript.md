MPRNA-IP source code and manuscript figures
================
Taeyoung Hwang, Jul 14, 2022

``` r
library(cowplot)
library(ggrepel)
library(ggpubr)
library(GGally)
library(gridExtra)

library(reshape2)
library(GenomicRanges)
```

``` r
theme_set(theme_cowplot())
source("github_functions.R")
```

# 1\. Optimization of MPRNA-IP (using hTR pool)

### Barcode number (Figure 1C)

Load raw count data.

``` r
load("RawData/final_hTR_pipeline.Rdata")
```

Randomly select barcodes for different numbers of barcodes.  
Calculate correlation coefficients of log2FC between two samples.

``` r
# Using the counts of replicate 2 and replicate 3
barcode.sim.corr <- simCorrLog2FC(c("Rep2", "Rep3"), 1:23, 10)
temp <- melt(barcode.sim.corr, varnames=c("barcodeNum", "simulation"), value.name="corr")
temp$barcodeNum <- factor(temp$barcodeNum, levels=1:23)

p <- ggplot(temp, aes(x=barcodeNum, y=corr))
p <- p + geom_boxplot()
p <- p + scale_y_continuous(breaks=seq(0,1,0.1))
p <- p + background_grid()
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Sequencing depth (Suppl. Figure 1C)

Load raw count data.

``` r
load("RawData/final_Assay_subsample.Rdata")
```

Select random oligos.

``` r
temp <- subset(tile.table, type=="Random")$tileID
count.table.random <- subset(count.table, tileID %in% temp)
```

Calculate RPM by using total depth with the sum of depth across random
oligos.

``` r
temp <- dcast(count.table.random, sample~., value.var = "raw", fun.aggregate = sum)
count.table.random$depth <- temp[match(count.table.random$sample, temp$sample), 2]
count.table.random$rpm <- count.table.random$raw/(count.table.random$depth/10^6)
```

Count the number of missing oligos at every sequencing depths.

``` r
# Using the replicate 1 (Input) experiment
temp <- mOligoCount(count.table.random, "rep1in", c(0,1,2))
temp <- melt(temp, id.vars = "sample", variable.name = "threshold", value.name="count")
temp$replicate <- count.table.random$replicate[match(temp$sample, count.table.random$sample)]
temp$sampling <- count.table.random$sampling[match(temp$sample, count.table.random$sample)]
temp$simID <- count.table.random$simID[match(temp$sample, count.table.random$sample)]
temp$depth <- count.table.random$depth[match(temp$sample, count.table.random$sample)]
tileNum <- length(subset(tile.table, type=="Random")$tileID)
temp$avgDepth <- temp$depth/(tileNum*BARCODE_NUM)
temp$countRate <- temp$count/(tileNum*BARCODE_NUM) # Proportion relative to 100*25 random oligos

p <- ggplot(temp, aes(x=avgDepth, y=countRate, col=threshold, shape=threshold))
p <- p + geom_point(size=2) + geom_smooth()
p <- p + xlab("Average depth") + ylab("Proportion of missing oligos")
p <- p + scale_x_continuous(breaks=seq(200,800,by=100))
p <- p + scale_color_manual(values=c("black", "blue", "red"))
p <- p + background_grid() # always place this after the theme
p <- p + theme(axis.text.x=element_text(angle=70, hjust=1, vjust=1))
p <- p + background_grid()
p
```

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](github_manuscript_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# 2\. PUM2 MPRNA-IP

Load raw count data and tile annotations.

``` r
load("RawData/final_PUM2_pipeline.Rdata") # raw data

gene.table <- read.table("pum2_annot_gene.txt", sep="\t", header=T, stringsAsFactors = F) # gene annotation

pre.table <- read.table("pum2_annot_pre.txt", sep="\t", header=T, stringsAsFactors = F) # PRE annotation

eCLIP.k562 <- read.table("encode_in_oligoPool.txt", header=T, sep="\t", stringsAsFactors = F) # eCLIP data
```

### PUM2-binding RNAs according to eCLIP data (Suppl. Figure 2A)

``` r
temp <- gene.table[order(gene.table$encodePeakNum, decreasing=T), ]
temp$gene_name <- factor(temp$gene_name, levels=temp$gene_name)
p <- ggplot(temp, aes(x=gene_name, y=encodePeakNum))
p <- p + geom_bar(stat="identity")
p <- p + geom_vline(xintercept = 15.5, linetype="dashed")
p <- p + xlab("Gene") + ylab("Number of eCLIP peaks")
p <- p + theme(axis.text.x=element_text(angle=70, hjust=1, vjust=1))
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Cumulative distributions of RPM in Input (Suppl. Figure 2B)

``` r
rpm.cdf <- list()
for (s in fastqc$Sample) {
  rpm.cdf[[s]] <- sapply(c(0.1,1,2,3,5,10), function(x) nrow(subset(count.table, sample==s & rpm>=x))) / (nrow(tile.table)*BARCODE_NUM)
}
rpm.cdf <- as.data.frame(rpm.cdf)
rownames(rpm.cdf) <- paste0("RPM>=", c(0.1,1,2,3,5,10))
```

``` r
temp <- melt(as.matrix(rpm.cdf), varnames = c("cutoff", "sample"), value.name="prop")
temp$pulldown <- fastqc$Pulldown[match(temp$sample, fastqc$Sample)]
temp$replicate <- fastqc$Replicate[match(temp$sample, fastqc$Sample)]
temp$cutoff <- factor(temp$cutoff, levels=paste0("RPM>=", c(0.1,1,2,3,5,10)))

p <- ggplot(subset(temp, pulldown=="Input"), aes(x=cutoff, y=prop, col=replicate, group=replicate))
p <- p + geom_point() + geom_line()
p <- p + scale_color_manual(values = c("black", "green", "red"))
p <- p + background_grid()
p <- p + theme(axis.text.x=element_text(angle=70, hjust=1, vjust=1, color = "black"))
p <- p + xlab("RPM threshold") + ylab("Proportion to total oligo #")
p <- p + theme(axis.text.x=element_text(angle=70, hjust=1, vjust=1))
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### RPM distribution (Suppl. Figure 2C)

``` r
p <- ggplot(count.table, aes(log10(rpm+1), fill=pulldown)) + facet_grid(.~replicate)
p <- p + geom_density(alpha=0.3)
p <- p + scale_fill_manual(values=c("black", "red"))
p <- p + background_grid()
p <- p + xlab("Log10 (RPM+1)") + ylab("Density")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### Fold change correlation (Suppl. Figure 3)

``` r
p <- ggpairs(log2fc.table[,c("rep1.mean", "rep2.mean", "rep3.mean")]) 
p <- p + xlab("Log2 (IP/Input)") + ylab("Log2 (IP/Input)")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

### GLM fitting

``` r
glm.table <- runGLMtest(BARCODE_LIST = BARCODE_LIST, REPLICATE_LIST = REPLICATE_LIST, IP = "IP", INPUT_RPM_TH=1, IP_COUNT_TH=0)
```

    ## Loading required package: MASS

### Multiple test correction

``` r
require(fdrtool)
```

    ## Loading required package: fdrtool

``` r
glm.table$poisson.epadj <- NA
idx <- which(!is.na(glm.table$poisson.zval))
temp <- fdrtool(glm.table$poisson.zval[idx], statistic= "normal", plot=F)
```

    ## Step 1... determine cutoff point
    ## Step 2... estimate parameters of null distribution and eta0
    ## Step 3... compute p-values and estimate empirical PDF/CDF
    ## Step 4... compute q-values and local fdr

``` r
glm.table$poisson.epadj[idx] <- temp$qval
```

### Identification of significant tiles

``` r
tile.table$n <- glm.table$n
tile.table$rpm <- glm.table$rpm
tile.table$log2FC <- log2(exp(glm.table$poisson.beta))
tile.table$pval <- glm.table$poisson.pval
tile.table$padj <- glm.table$poisson.epadj
tile.table$sig <- factor(tile.table$padj<0.05, levels=c("TRUE", "FALSE"), labels=c("Sig", "Insig"))
tile.table$enriched <- "Insig"
tile.table$enriched[which(tile.table$sig=="Sig" & tile.table$log2FC>0)] <- "Enriched"
tile.table$enriched[which(tile.table$sig=="Sig" & tile.table$log2FC<0)] <- "Depleted"
table(tile.table$enriched) 
```

    ## 
    ## Depleted Enriched    Insig 
    ##      178      155     1266

### MA plot (Suppl. Figure 4)

``` r
temp <- subset(tile.table, !is.na(log2FC))
temp$log2FC2 <- temp$log2FC
idx <- which(temp$log2FC<(-10))
temp$log2FC2[idx] <- (-10) # for optimal figure
p <- ggplot(temp, aes(x=log10(rpm), y=log2FC2))
p <- p + geom_point(data=subset(temp, enriched=="Insig"), col="gray", size=2, alpha=0.5)
p <- p + geom_point(data=subset(temp, enriched=="Depleted"), col="blue", size=2, alpha=0.8)
p <- p + geom_point(data=subset(temp, enriched=="Enriched"), col="red", size=2, alpha=0.8)
p <- p + geom_hline(yintercept = 0, linetype="dashed")
p <- p + xlab("Log10 (Mean RPM)") + ylab("Log2 (IP/Input)")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### NORAD log2FC (Figure 2A)

``` r
temp <- subset(tile.table, gene_name=="NORAD")

p <- ggplot(temp)
p <- p + geom_rect(data=subset(temp, enriched!="Enriched"), aes(xmin=start, xmax=end, ymin=0, ymax=2^(log2FC)), fill="gray", alpha=0.5)
p <- p + geom_rect(data=subset(temp, enriched=="Enriched"), aes(xmin=start, xmax=end, ymin=0, ymax=2^(log2FC)), fill="red", alpha=0.5)
p <- p + geom_rect(data=subset(eCLIP.k562, gene_name=="NORAD" & width>=8), aes(xmin=mRNA_start, xmax=mRNA_end, ymin=-2, ymax=0), fill="black")
p <- p + geom_rect(data=subset(pre.table, seqnames=="NORAD"), aes(xmin=start, xmax=end, ymin=-3.5, ymax=-2.5), fill="black")
p <- p + scale_x_continuous(breaks = c(1,1000,2000,3000,4000,5000)) # ,5343
p
```

    ## Warning: Removed 2 rows containing missing values (geom_rect).

![](github_manuscript_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Correlation with eCLIP (Figure 2B)

``` r
temp <- dcast(subset(tile.table, enriched=="Enriched"), gene_name~., value.var = "tileID")
```

    ## Aggregation function missing: defaulting to length

``` r
gene.table$eTileNum <- temp[match(gene.table$gene_name, temp$gene_name), 2]
gene.table$eTileNum[which(is.na(gene.table$eTileNum))] <- 0 

temp <- gene.table
temp <- temp[order(temp$encodePeakNum/temp$mRNA_size, decreasing = T),]
temp$label <- ""
temp$label[1:5] <- temp$gene_name[1:5]
temp$y <- temp$eTileNum/temp$tileNum
temp$x <- temp$encodePeakNum/(temp$mRNA_size/1000)

p <- ggplot(temp, aes(x=x, y=y, label=label))
p <- p + geom_point(size=3) + geom_smooth(method="lm", se = T) 
p <- p + geom_text_repel(size=4)
p <- p + background_grid()
p
```

    ## `geom_smooth()` using formula 'y ~ x'

![](github_manuscript_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
cor.test(temp$x, temp$y)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  temp$x and temp$y
    ## t = 4.4546, df = 21, p-value = 0.0002191
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.3996469 0.8616591
    ## sample estimates:
    ##       cor 
    ## 0.6970245

### Correlation with PRE (Figure 2D)

``` r
temp <- subset(tile.table, !is.na(log2FC))
table(temp$pre.num)
```

    ## 
    ##    0    1    2 
    ## 1423   85   10

``` r
idx <- which(temp$log2FC<(-20))
temp$log2FC[idx] <- min(temp$log2FC[-idx])-5
temp$outlier <- FALSE
temp$outlier[idx] <- TRUE

p <- ggviolin(temp, x="pre.num", y="log2FC", color="pre.num", shape="outlier", palette=c("black", "blue", "red"), alpha=1, add="jitter", add.params=list(size=2, alpha=0.7))
p <- p + stat_compare_means(comparisons = list( c("0", "1"), c("1", "2") ))
p <- p + geom_hline(yintercept = 0, linetype="dashed")
p <- p + theme(legend.position="none") # This removes all legends.
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# 3\. MS2 MPRNA-IP

Load raw count data.

``` r
load("RawData/final_4xMS2_pipeline.Rdata") # raw data
```

### FLAG vs. HA (Suppl. Figure 5B)

FLAG

``` r
library(ggcorrplot)
temp <- matrix(rep(NA, 3*3), nrow=3, dimnames = list(c("Rep1", "Rep2", "Rep3"), c("Rep1", "Rep2", "Rep3")))
temp[cbind(1:nrow(temp), 1:nrow(temp))] <- 1
temp[1,2] <- cor.test(log2fc.table.flag$rep1.mean, log2fc.table.flag$rep2.mean)$estimate
temp[1,3] <- cor.test(log2fc.table.flag$rep1.mean, log2fc.table.flag$rep3.mean)$estimate
temp[2,3] <- cor.test(log2fc.table.flag$rep2.mean, log2fc.table.flag$rep3.mean)$estimate

p <- ggcorrplot(temp, type = "upper", lab = TRUE, lab_size = 5, show.diag = T, digits = 3, colors = c("gray", "gray", "gray"), outline.color = "black")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

HA

``` r
temp <- matrix(rep(NA, 3*3), nrow=3, dimnames = list(c("Rep1", "Rep2", "Rep3"), c("Rep1", "Rep2", "Rep3")))
temp[cbind(1:nrow(temp), 1:nrow(temp))] <- 1
temp[1,2] <- cor.test(log2fc.table.ha$rep1.mean, log2fc.table.ha$rep2.mean)$estimate
temp[1,3] <- cor.test(log2fc.table.ha$rep1.mean, log2fc.table.ha$rep3.mean)$estimate
temp[2,3] <- cor.test(log2fc.table.ha$rep2.mean, log2fc.table.ha$rep3.mean)$estimate

p <- ggcorrplot(temp, type = "upper", lab = TRUE, lab_size = 5, show.diag = T, digits = 3, colors = c("gray", "gray", "gray"), outline.color = "black")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

### Calculation of stem score

Load stem design
information.

``` r
stem.design <- read.table("ms2_mt_stem.txt", header=T, sep="\t", stringsAsFactors = F)
stem.design$stemID <- paste(stem.design$idx, stem.design$mutation, sep="_")

# map to tileID
temp <- subset(tile.table, group=="Stem")
temp$stemID <- paste(temp$mutLoc, temp$mutType, sep="_")
stem.design$tileID <- temp$tileID[match(stem.design$stemID, temp$stemID)]

# Fix duplicated sequences in the assignment above.
temp <- split(stem.design, f = stem.design$seq)
temp <- lapply(temp, function(x) {idx <- which(!is.na(x$tileID)); if (length(idx==1)) {x$tileID<-x$tileID[idx]; return(x)} else {cat("Error"); return(0)}})
stem.design <- do.call(rbind, temp)
rownames(stem.design) <- 1:nrow(stem.design)
```

Initialize a main object.

``` r
stem.table <- dcast(stem.design, idx~mutation, value.var = "tileID")
colnames(stem.table) <- c("stemID", "tileID.mut3p", "tileID.mut5p", "tileID.comp")
stem.table <- stem.table[order(as.numeric(gsub("stem:", "", stem.table$stemID))), ]
```

GLM
fitting

``` r
stem.glm.out <- data.frame(tileID=stem.table$stemID, n=NA, rpm=NA, p.seq.beta=NA, p.seq.sd=NA, p.seq.zval=NA, p.seq.pval=NA, p.stem.beta=NA, p.stem.sd=NA, p.stem.zval=NA, p.stem.pval=NA, poisson.overdisp=NA, poisson.overdisp.pval=NA, nb.seq.beta=NA, nb.seq.sd=NA, nb.seq.zval=NA, nb.seq.pval=NA, nb.stem.beta=NA, nb.stem.sd=NA, nb.stem.zval=NA, nb.stem.pval=NA, nb.alpha=NA, resid.diff=NA)
```

``` r
X.wt <- data.frame(mutType = rep("WT", each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(subset(tile.table, group=="WT")$tileID, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(paste0("Barcode:", 1:BARCODE_NUM), each=REPLICATE_NUM*2),
                  replicate = rep(rep(c("Rep1", "Rep2", "Rep3"), each=2), BARCODE_NUM),
                  pulldown = rep(rep(c("Input", "HA"), REPLICATE_NUM), BARCODE_NUM),
                  stringsAsFactors = F)

for (i in 1:nrow(stem.glm.out)) {
  #print(i)
  
  temp <- unlist(stem.table[i, c("tileID.mut5p", "tileID.mut3p", "tileID.comp"), drop=T])
  X.mut <- data.frame(mutType = rep(names(temp), each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(temp, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(rep(paste0("Barcode:", 1:BARCODE_NUM), each=REPLICATE_NUM*2), length(temp)),
                  replicate = rep(rep(rep(c("Rep1", "Rep2", "Rep3"), each=2), BARCODE_NUM), length(temp)),
                  pulldown = rep(rep(rep(c("Input", "HA"), REPLICATE_NUM), BARCODE_NUM), length(temp)),
                  stringsAsFactors = F)
  X <- rbind(X.wt, X.mut)
  X$oligoID <- paste(X$tileID, X$barcodeID, sep="_")
  
  X$seqEffect <- 0
  X$seqEffect[which(X$mutType %in% c("WT"))] <- 1

  X$stemEffect <- 0
  X$stemEffect[which(X$mutType %in% c("WT", "tileID.comp"))] <- 1

  y <- subset(count.table, (tileID %in% unique(X$tileID)) & (pulldown %in% c("Input", "HA")))
  
  # filter barcodes by input rpm
  rpm.th <- 1
  temp <- dcast(subset(y, pulldown=="Input"), barcodeID~., value.var = "rpm", fun.aggregate = mean)
  y <- subset(y, barcodeID %in% temp[temp[, 2]>=rpm.th, 1])
  if (length(unique(y$barcodeID))<3) {next; }

  # input of glm
  idx <- match(paste(y$oligoID, y$replicate, y$pulldown, sep="_"), paste(X$oligoID, X$replicate, X$pulldown, sep="_"))
  cur.data <- data.frame(X[idx, ], count=y$normalized, stringsAsFactors = F)
  cur.data$pulldown <- factor(cur.data$pulldown, levels=c("Input", "HA"))
  cur.data$rep <- factor(cur.data$rep)
  cur.data$seqEffect <- factor(cur.data$seqEffect)
  cur.data$stemEffect <- factor(cur.data$stemEffect)
  #if (min(table(cur.data$rep, cur.data$pulldown)) < 2) {next;}

  stem.glm.out$n[i] = nrow(y)
  stem.glm.out$rpm[i] = mean(y$rpm)

  # fitting
  fit.nb <- tryCatch( glm.nb(round(count) ~ replicate + pulldown + seqEffect + stemEffect + pulldown:seqEffect + pulldown:stemEffect, data=cur.data), error = function(err) { print("fit.nb"); print(err); return("err")})

  # test
  if (is.list(fit.nb)) {
    residual.nb <- round(cur.data$count)-fit.nb$fitted.values
    residual.z.nb <- residual.nb / sqrt(fit.nb$fitted.values)
    temp <- summary(fit.nb)
    stem.glm.out$nb.stem.beta[i] <- temp[["coefficients"]]["pulldownHA:stemEffect1", 1] # "Estimate"
    stem.glm.out$nb.stem.sd[i] <- temp[["coefficients"]]["pulldownHA:stemEffect1", 2] # "Std. Error"
    stem.glm.out$nb.stem.zval[i] <- temp[["coefficients"]]["pulldownHA:stemEffect1", 3] # "z value"
    stem.glm.out$nb.stem.pval[i] <- temp[["coefficients"]]["pulldownHA:stemEffect1", 4] # "Pr(>|z|)"
    
    stem.glm.out$nb.seq.beta[i] <- temp[["coefficients"]]["pulldownHA:seqEffect1", 1] # "Estimate"
    stem.glm.out$nb.seq.sd[i] <- temp[["coefficients"]]["pulldownHA:seqEffect1", 2] # "Std. Error"
    stem.glm.out$nb.seq.zval[i] <- temp[["coefficients"]]["pulldownHA:seqEffect1", 3] # "z value"
    stem.glm.out$nb.seq.pval[i] <- temp[["coefficients"]]["pulldownHA:seqEffect1", 4] # "Pr(>|z|)"

    stem.glm.out$nb.alpha[i] <- 1/temp$theta
  }
} 
```

Multiple test correction and select signficant stem
mutations

``` r
stem.glm.out$nb.padj <- p.adjust(stem.glm.out$nb.stem.pval, method="fdr")

stem.table$n <- stem.glm.out$n
stem.table$rpm <- stem.glm.out$rpm
stem.table$stem.log2FC <- log2(exp(stem.glm.out$nb.stem.beta))
stem.table$seq.log2FC <- log2(exp(stem.glm.out$nb.seq.beta))
stem.table$stem.pval <- stem.glm.out$nb.stem.pval
stem.table$seq.pval <- stem.glm.out$nb.seq.pval
stem.table$padj <- stem.glm.out$nb.padj
stem.table$sig <- factor(stem.table$padj<0.05, levels=c("TRUE", "FALSE"), labels=c("Sig", "Insig"))
```

Stem score plot (Figure 3C)

``` r
temp <- stem.table
temp$stemEffect <- 2^(temp$stem.log2FC)

p <- ggdotchart(temp, x = "stemID", y = "stemEffect", dot.size=3, add = "segments", palette=c("red", "black"))
p <- p + geom_hline(yintercept = c(1,2), linetype="dashed")
p <- p + xlab("Stem ID") + ylab("Stem effect")
p <- p + ggtitle("Stem score")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Fold change plot of “stem:8” (Figure 3C inset)

``` r
temp <- stem_log2FC_table("stem:8")
p <- ggboxplot(temp, x = "tileID", y = "log2FC", palette = "jco", add="jitter", add.params=list(size=1), outlier.shape=NA, legend = "none")
p <- p + stat_compare_means(comparisons = list( c("WT", "5' mut"), c("WT", "3' mut"), c("WT", "Comp. mut") ))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p <- p + xlab("Mutation") + ylab("Log2 (HA/Input)")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### Single mutation

Initialize a main
object.

``` r
single.table <- subset(tile.table, group=="Single")[,c("tileID", "mutType", "mutLoc")]
temp <- strsplit(single.table$mutType, split="to")
single.table$seqWT <- sapply(temp, "[[", 1)
single.table$seqMT <- sapply(temp, "[[", 2)
```

GLM
fitting

``` r
single.glm.out <- data.frame(tileID=single.table$tileID, n=NA, rpm=NA, poisson.beta=NA, poisson.sd=NA, poisson.zval=NA, poisson.pval=NA, poisson.overdisp=NA, poisson.overdisp.pval=NA, nb.beta=NA, nb.sd=NA, nb.zval=NA, nb.pval=NA, nb.alpha=NA, resid.diff=NA)
```

``` r
# WT object
X.wt <- data.frame(mutType = rep("WT", each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(subset(tile.table, group=="WT")$tileID, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(paste0("Barcode:", 1:BARCODE_NUM), each=REPLICATE_NUM*2),
                  replicate = rep(rep(c("Rep1", "Rep2", "Rep3"), each=2), BARCODE_NUM),
                  pulldown = rep(rep(c("Input", "HA"), REPLICATE_NUM), BARCODE_NUM),
                  betaMut = 0, 
                  stringsAsFactors = F)

for (i in 1:nrow(single.glm.out)) {
  
  # MT object
  X.mut <- data.frame(mutType = rep(single.table$mutType[i], each=BARCODE_NUM*REPLICATE_NUM*2), 
                      tileID = rep(single.table$tileID[i], each=BARCODE_NUM*REPLICATE_NUM*2),
                      barcodeID = rep(paste0("Barcode:", 1:BARCODE_NUM), each=REPLICATE_NUM*2),
                      replicate = rep(rep(c("Rep1", "Rep2", "Rep3"), each=2), BARCODE_NUM),
                      pulldown = rep(rep(c("Input", "HA"), REPLICATE_NUM), BARCODE_NUM),
                      betaMut = 1, 
                      stringsAsFactors = F)
  X <- rbind(X.wt, X.mut)
  X$oligoID <- paste(X$tileID, X$barcodeID, sep="_")
  
  y <- subset(count.table, tileID %in% unique(X$tileID) & pulldown %in% c("Input", "HA"))
  
  # filter barcodes by input rpm
  rpm.th <- 1
  temp <- dcast(subset(y, pulldown=="Input"), barcodeID~., value.var = "rpm", fun.aggregate = mean)
  y <- subset(y, barcodeID %in% temp[temp[, 2]>=rpm.th, 1])
  if (length(unique(y$barcodeID))<3) {next; }

  # input of glm
  idx <- match(paste(y$oligoID, y$replicate, y$pulldown, sep="_"), paste(X$oligoID, X$replicate, X$pulldown, sep="_"))
  cur.data <- data.frame(X[idx, ], count=y$normalized, stringsAsFactors = F)
  cur.data$pulldown <- factor(cur.data$pulldown, levels=c("Input", "HA"))
  cur.data$rep <- factor(cur.data$rep)
  cur.data$betaMut <- factor(cur.data$betaMut)
  #if (min(table(cur.data$rep, cur.data$pulldown)) < 2) {next;}

  single.glm.out$n[i] = nrow(y)
  single.glm.out$rpm[i] = mean(y$rpm)

  # fitting
  fit.nb <- tryCatch( glm.nb(round(count) ~ replicate + pulldown + betaMut + betaMut:pulldown, data=cur.data), error = function(err) { print("fit.nb"); print(err); return("err")})

  # test
  if (is.list(fit.nb)) {
    residual.nb <- round(cur.data$count)-fit.nb$fitted.values
    residual.z.nb <- residual.nb / sqrt(fit.nb$fitted.values)
    temp <- summary(fit.nb)
    single.glm.out$nb.beta[i] <- temp[["coefficients"]]["pulldownHA:betaMut1", 1] # "Estimate"
    single.glm.out$nb.sd[i] <- temp[["coefficients"]]["pulldownHA:betaMut1", 2] # "Std. Error"
    single.glm.out$nb.zval[i] <- temp[["coefficients"]]["pulldownHA:betaMut1", 3] # "z value"
    single.glm.out$nb.pval[i] <- temp[["coefficients"]]["pulldownHA:betaMut1", 4] # "Pr(>|z|)"
    single.glm.out$nb.alpha[i] <- 1/temp$theta
  }
}

single.glm.out$nb.padj <- p.adjust(single.glm.out$nb.pval, method="fdr")
```

Multiple test correction and select signficant mutations

``` r
single.table$n <- single.glm.out$n
single.table$rpm <- single.glm.out$rpm
single.table$log2FC <- log2(exp(single.glm.out$nb.beta))
single.table$pval <- single.glm.out$nb.pval
single.table$padj <- single.glm.out$nb.padj
single.table$sig <- factor(single.table$padj<0.05, levels=c("TRUE", "FALSE"), labels=c("Sig", "Insig"))
```

Summary plot (Suppl. Figure 5C)

``` r
temp <- single.table
temp$xlab <- paste(temp$mutLoc, temp$seqWT, sep=":")
temp$xlab <- factor(temp$xlab, levels=unique(temp$xlab))
temp$Sequence <- temp$seqMT
temp$P_value <- cut(temp$padj, breaks = c(0,0.1,0.2,0.3,1), labels=c("<=0.1", "<=0.2", "<=0.3", "<=1"))

p <- ggplot(temp, aes(x=xlab, y=2^log2FC, col=P_value, shape=Sequence))
p <- p + geom_segment(size=0.5, aes(xend=xlab), yend=0, color='grey80')
p <- p + geom_point(size=3)
p <- p + scale_color_manual(values = c("red", "orange", "black", "blue"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p <- p + xlab("Location on MS2") + ylab(expression(Delta* " (HA / Input)"))
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

# 4\. hTR MPRNA-IP

Load raw count data

``` r
load("~/Research/MPRNA/Manuscript/final_hTR_pipeline.Rdata")
```

### Enrichment barplot (Figure 4B)

``` r
temp <- subset(count.table, type %in% c("Random", "Control"))
temp <- dcast(temp, type+coord+tileID+barcodeID+replicate~pulldown, value.var = "rpm")
temp$log2fc <- log2(temp$IP / temp$Input)
temp <- subset(temp, Input>=1) # RPM>=1
temp$fc <- 2^(temp$log2fc)

temp1 <- dcast(temp, coord~. , value.var="fc", mean)
temp1$se <- dcast(temp, coord~. , value.var="fc", function(x) {sd(x)/sqrt(length(x))})[,"."]
colnames(temp1) <- c("Tile", "Mean", "SE")
temp1$Tile <- factor(temp1$Tile, levels=c("Random", "hTR:34-190", "hTR:210-366", "hTR:295-451"), labels = c("Random", "Tile 1", "Tile 2", "Tile 3"))

p <- ggplot(temp1, aes(x=Tile, y=Mean, fill=Tile)) 
p <- p + geom_bar(stat="Identity")
p <- p + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), size=0.8, width=0.5)
p <- p + scale_fill_manual(values=alpha(c("gray50", "orange", "red", "gray10"), 0.8))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

``` r
wilcox.test(subset(temp, coord=="Random")$fc, subset(temp, coord=="hTR:34-190")$fc)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  subset(temp, coord == "Random")$fc and subset(temp, coord == "hTR:34-190")$fc
    ## W = 64262, p-value = 0.004813
    ## alternative hypothesis: true location shift is not equal to 0

``` r
wilcox.test(subset(temp, coord=="Random")$fc, subset(temp, coord=="hTR:210-366")$fc)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  subset(temp, coord == "Random")$fc and subset(temp, coord == "hTR:210-366")$fc
    ## W = 15928, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0

``` r
wilcox.test(subset(temp, coord=="Random")$fc, subset(temp, coord=="hTR:295-451")$fc)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  subset(temp, coord == "Random")$fc and subset(temp, coord == "hTR:295-451")$fc
    ## W = 193492, p-value = 0.7561
    ## alternative hypothesis: true location shift is not equal to 0

### Deletion and replacement tiles

Initialize a main object.

``` r
temp <- subset(tile.table, type =="Deletion")
temp$bpLocPair <- paste(temp$bpLoc5prime, temp$bpLoc3prime, sep=":")
deletion.table <- dcast(temp, coord+region+bpLocPair~mutation, value.var = "tileID")
colnames(deletion.table)[(ncol(deletion.table)-1):ncol(deletion.table)] <- c("tileID.del", "tileID.rep")

temp <- subset(tile.table, type=="Control")
deletion.table$tileID.wt <- temp$tileID[match(deletion.table$coord, temp$coord)]
```

GLM fitting for deletion

``` r
for (i in 1:nrow(deletion.table)) {
  temp <- unlist(deletion.table[i, c("tileID.wt", "tileID.del"), drop=T])
  X <- data.frame(mutation = rep(names(temp), each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(temp, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(rep(BARCODE_LIST, each=REPLICATE_NUM*2), length(temp)),
                  replicate = rep(rep(rep(REPLICATE_LIST, each=2), BARCODE_NUM), length(temp)),
                  pulldown = rep(rep(rep(c("Input", "IP"), REPLICATE_NUM), BARCODE_NUM), length(temp)),
                  stringsAsFactors = F)
  X$mutation <- sapply(strsplit(X$mutation, split="[.]"), "[[", 2)
  X$oligoID <- paste(X$tileID, X$barcodeID, sep="_")
  
  X$deletion <- 0  
  X$deletion[which(X$mutation %in% c("del"))] <- 1

  y <- subset(count.table, tileID %in% unique(X$tileID))
  y$id <- paste(y$replicate, y$oligoID, sep="_")
  #nrow(y)/2 == length(unique(y$id))
  
  # filter
  temp <- dcast(y, replicate+oligoID~pulldown, value.var = "rpm")
  temp <- temp[apply(temp[,c("Input", "IP")], 1, max) >= 1, ]
  y <- subset(y, id %in% paste(temp$replicate, temp$oligoID, sep="_"))
  
  idx <- match(paste(y$id, y$pulldown, sep="_"), paste(X$replicate, X$oligoID, X$pulldown, sep="_"))
  X <- X[idx, ]
  
  X$pulldown <- factor(X$pulldown, levels=c("Input", "IP"))
  X$deletion <- factor(X$deletion, levels=c(0,1))
  #fit <- glm.nb(round(y$rpm) ~ X$pulldown + X$betaMut + X$betaMut : X$pulldown + X$betaStem + X$betaStem : X$pulldown)
  fit <- glm.nb(round(y$normalized) ~ X$replicate + X$pulldown + X$deletion + X$pulldown : X$deletion)
  #summary(fit)
  if (i==1) {
    out <- summary(fit)[["coefficients"]]["X$pulldownIP:X$deletion1", ]
  } else {
    out <- rbind(out, summary(fit)[["coefficients"]]["X$pulldownIP:X$deletion1", ])
  }
}
colnames(out) <- paste0("del.", c("beta", "sd", "statistic", "pval"))
deletion.table <- cbind(deletion.table, out)
```

GLM fitting for Replacement

``` r
for (i in 1:nrow(deletion.table)) {
  temp <- unlist(deletion.table[i, c("tileID.wt", "tileID.rep"), drop=T])
  X <- data.frame(mutation = rep(names(temp), each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(temp, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(rep(BARCODE_LIST, each=REPLICATE_NUM*2), length(temp)),
                  replicate = rep(rep(rep(REPLICATE_LIST, each=2), BARCODE_NUM), length(temp)),
                  pulldown = rep(rep(rep(c("Input", "IP"), REPLICATE_NUM), BARCODE_NUM), length(temp)),
                  stringsAsFactors = F)
  X$mutation <- sapply(strsplit(X$mutation, split="[.]"), "[[", 2)
  X$oligoID <- paste(X$tileID, X$barcodeID, sep="_")
  
  X$deletion <- 0  
  X$deletion[which(X$mutation %in% c("rep"))] <- 1

  y <- subset(count.table, tileID %in% unique(X$tileID))
  y$id <- paste(y$replicate, y$oligoID, sep="_")
  #nrow(y)/2 == length(unique(y$id))
  
  # filter
  temp <- dcast(y, replicate+oligoID~pulldown, value.var = "rpm")
  temp <- temp[apply(temp[,c("Input", "IP")], 1, max) >= 1, ]
  y <- subset(y, id %in% paste(temp$replicate, temp$oligoID, sep="_"))
  
  idx <- match(paste(y$id, y$pulldown, sep="_"), paste(X$replicate, X$oligoID, X$pulldown, sep="_"))
  X <- X[idx, ]
  
  X$pulldown <- factor(X$pulldown, levels=c("Input", "IP"))
  X$deletion <- factor(X$deletion, levels=c(0,1))
  #fit <- glm.nb(round(y$rpm) ~ X$pulldown + X$betaMut + X$betaMut : X$pulldown + X$betaStem + X$betaStem : X$pulldown)
  fit <- glm.nb(round(y$normalized) ~ X$replicate + X$pulldown + X$deletion + X$pulldown : X$deletion)
  #summary(fit)
  if (i==1) {
    out <- summary(fit)[["coefficients"]]["X$pulldownIP:X$deletion1", ]
  } else {
    out <- rbind(out, summary(fit)[["coefficients"]]["X$pulldownIP:X$deletion1", ])
  }
}
colnames(out) <- paste0("rep.", c("beta", "sd", "statistic", "pval"))
deletion.table <- cbind(deletion.table, out)
```

Multiple test
correction

``` r
deletion.table$del.padj <- p.adjust(deletion.table$del.pval, method="fdr")
deletion.table$rep.padj <- p.adjust(deletion.table$rep.pval, method="fdr")
```

A summary figure (Figure 4C)

``` r
temp <- subset(deletion.table, coord!="hTR:295-451")
temp$region <- factor(temp$region, levels=temp$region)
temp$Tile <- factor(temp$coord, levels=c("hTR:34-190", "hTR:210-366"), labels=c("Pesudoknot", "CR4/CR5"))
temp$log2fc1 <- log2(exp(temp$del.beta))
temp$Significance1 <- factor(temp$del.padj<0.05)
table(temp$Significance1)
```

    ## 
    ## FALSE  TRUE 
    ##    10     2

``` r
temp$log2fc2 <- log2(exp(temp$rep.beta))
temp$Significance2 <- factor(temp$rep.padj<0.05)
table(temp$Significance2)
```

    ## 
    ## FALSE 
    ##    12

``` r
p <- ggplot(temp, aes(x=region, y=log2fc1, col=Significance1), shape="circle")
p <- p + geom_point(size=3)
p <- p + geom_point(aes(y=log2fc2, col=Significance2), size=6, shape="triangle")
p <- p + scale_color_manual(values=c("black", "red"))
p <- p + geom_segment(size=0.7, aes(xend=region, yend=log2fc2), color='blue')
p <- p + geom_hline(yintercept = 0, linetype="dashed")
p <- p + xlab("Tile region") + ylab("Log2FC")
p <- p + background_grid() 
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

Effect size

``` r
temp <- subset(deletion.table, coord!="hTR:295-451")
temp$fc <- exp(temp$del.beta)
dcast(temp, coord~., value.var="fc", fun.aggregate=mean)
```

    ##         coord        .
    ## 1  hTR:34-190 1.036575
    ## 2 hTR:210-366 0.676701

### Compensatory mutation analysis

Initialize a main object.

``` r
temp <- subset(tile.table, type %in% c("Mut1bp", "Mut4bp"))
temp$bpLocPair <- paste(temp$bpLoc5prime, temp$bpLoc3prime, sep=":")
stem.table <- dcast(temp, type+coord+region+bpLocPair~mutation, value.var = "tileID")
colnames(stem.table)[(ncol(stem.table)-2):ncol(stem.table)] <- c("tileID.mut3p", "tileID.mut5p", "tileID.comp")

temp <- subset(tile.table, type=="Control")
stem.table$tileID.wt <- temp$tileID[match(stem.table$coord, temp$coord)]
```

GLM fitting

``` r
for (i in 1:nrow(stem.table)) {
  temp <- unlist(stem.table[i, c("tileID.wt", "tileID.mut5p", "tileID.mut3p", "tileID.comp"), drop=T])
  X <- data.frame(mutation = rep(names(temp), each=BARCODE_NUM*REPLICATE_NUM*2), 
                  tileID = rep(temp, each=BARCODE_NUM*REPLICATE_NUM*2),
                  barcodeID = rep(rep(BARCODE_LIST, each=REPLICATE_NUM*2), length(temp)),
                  replicate = rep(rep(rep(REPLICATE_LIST, each=2), BARCODE_NUM), length(temp)),
                  pulldown = rep(rep(rep(c("Input", "IP"), REPLICATE_NUM), BARCODE_NUM), length(temp)),
                  stringsAsFactors = F)
  X$mutation <- sapply(strsplit(X$mutation, split="[.]"), "[[", 2)
  X$oligoID <- paste(X$tileID, X$barcodeID, sep="_")

  X$seqEffect <- 0
  X$seqEffect[which(X$mutation %in% c("wt"))] <- 1

  X$stemEffect <- 0
  X$stemEffect[which(X$mutation %in% c("wt", "comp"))] <- 1

  y <- subset(count.table, tileID %in% unique(X$tileID))
  y$id <- paste(y$replicate, y$oligoID, sep="_")
  #nrow(y)/2 == length(unique(y$id))
  
  # filter
  temp <- dcast(y, replicate+oligoID~pulldown, value.var = "rpm")
  temp <- temp[apply(temp[,c("Input", "IP")], 1, max) >= 1, ]
  y <- subset(y, id %in% paste(temp$replicate, temp$oligoID, sep="_"))
  
  idx <- match(paste(y$id, y$pulldown, sep="_"), paste(X$replicate, X$oligoID, X$pulldown, sep="_"))
  X <- X[idx, ]
  
  X$pulldown <- factor(X$pulldown, levels=c("Input", "IP"))
  X$seqEffect <- factor(X$seqEffect, levels=c(0,1))
  X$stemEffect <- factor(X$stemEffect, levels=c(0,1))
  fit <- glm.nb(round(y$normalized) ~ X$replicate + X$pulldown + X$seqEffect + X$stemEffect + X$seqEffect : X$pulldown + X$stemEffect : X$pulldown)
  temp <- summary(fit)[["coefficients"]]
  if (i==1) {
    out <- c(mean(y$rpm), temp["X$pulldownIP:X$stemEffect1", ],temp["X$pulldownIP:X$seqEffect1", ])
  } else {
    out <- rbind(out, c(mean(y$rpm), temp["X$pulldownIP:X$stemEffect1", ],temp["X$pulldownIP:X$seqEffect1", ]))
  }
}
colnames(out) <- c("rpm.mean", paste0("bp.", c("beta", "sd", "statistic", "pval")), paste0("seq.", c("beta", "sd", "statistic", "pval")))
stem.table <- cbind(stem.table, out)
```

Multiple test correction and select signficant stem mutation

``` r
stem.table$bp.padj <- p.adjust(stem.table$bp.pval, method="fdr")

stem.table$bp.sig <- factor(stem.table$bp.padj<0.05, levels=c(TRUE, FALSE), labels=c("Sig", "Insig"))
table(stem.table$bp.sig)
```

    ## 
    ##   Sig Insig 
    ##    10   149

Dotplot of 4bp mutation beta (Suppl. Figure 7A)

``` r
temp <- subset(stem.table, coord!="hTR:295-451" & type=="Mut4bp")
temp$effect <- exp(temp$bp.beta)
temp$stemID <- paste(temp$region, temp$bpLocPair, sep="_")

p <- ggdotchart(temp, x = "stemID", y = "effect", col="bp.sig", dot.size=3, add = "segments", palette=c("red", "black"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p <- p + xlab("Stem ID") + ylab("Stem score")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

Dotplot of 1bp mutation beta (Suppl. Figure 8A)

``` r
temp <- subset(stem.table, coord!="hTR:295-451" & type=="Mut1bp")
temp$effect <- exp(temp$bp.beta)
temp$stemID <- paste(temp$region, temp$bpLocPair, sep="_")

p <- ggdotchart(temp, x = "stemID", y = "effect", col="bp.sig", dot.size=4, add = "segments", palette=c("red", "black"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p <- p + theme(axis.text.x=element_text(angle=70, hjust=1, vjust=1))
p <- p + xlab("Stem ID") + ylab("Stem score")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

Log2FC plot (Figure 4E)

``` r
temp <- hTR_stem_log2FC_table("Mut4bp", "hTR:210-366", "258:298")
temp <- subset(temp, Input>=1)
temp$log2FC <- log2(temp$IP+1) - log2(temp$Input+1)
p <- ggboxplot(temp, x = "tileID", y = "log2FC", outlier.shape=NA, legend = "none")
p <- p + stat_compare_means(comparisons = list( c("WT", "5' mut"), c("WT", "3' mut"), c("WT", "Comp. mut") ))
p <- p + geom_jitter(aes(col=replicate), size=2, width = 0.2, height=0)
p <- p + scale_color_manual(values=c("darkcyan", "black", "brown"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Log2FC plot (Suppl. Figure 7B)

``` r
temp <- hTR_stem_log2FC_table("Mut4bp", "hTR:210-366", "303:313")
temp <- subset(temp, Input>=1)
temp$log2FC <- log2(temp$IP+1) - log2(temp$Input+1)
p <- ggboxplot(temp, x = "tileID", y = "log2FC", outlier.shape=NA, legend = "none")
p <- p + stat_compare_means(comparisons = list( c("WT", "5' mut"), c("WT", "3' mut"), c("WT", "Comp. mut") ))
p <- p + geom_jitter(aes(col=replicate), size=2, width = 0.2, height=0)
p <- p + scale_color_manual(values=c("darkcyan", "black", "brown"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

Log2FC plot (Suppl. Figure 7C)

``` r
temp <- hTR_stem_log2FC_table("Mut4bp", "hTR:210-366", "243:326")
temp <- subset(temp, Input>=1)
temp$log2FC <- log2(temp$IP+1) - log2(temp$Input+1)
p <- ggboxplot(temp, x = "tileID", y = "log2FC", outlier.shape=NA, legend = "none")
p <- p + stat_compare_means(comparisons = list( c("WT", "5' mut"), c("WT", "3' mut"), c("WT", "Comp. mut") ))
p <- p + geom_jitter(aes(col=replicate), size=2, width = 0.2, height=0)
p <- p + scale_color_manual(values=c("darkcyan", "black", "brown"))
p <- p + geom_hline(yintercept = 1, linetype="dashed")
p
```

![](github_manuscript_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

# R session information

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggcorrplot_0.1.3     fdrtool_1.2.15       MASS_7.3-51.5       
    ##  [4] GenomicRanges_1.38.0 GenomeInfoDb_1.22.1  IRanges_2.20.2      
    ##  [7] S4Vectors_0.24.4     BiocGenerics_0.32.0  reshape2_1.4.3      
    ## [10] gridExtra_2.3        GGally_2.0.0         ggpubr_0.3.0        
    ## [13] ggrepel_0.9.1        ggplot2_3.3.3        cowplot_1.0.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyr_1.1.0            splines_3.6.3          carData_3.0-4         
    ##  [4] assertthat_0.2.1       GenomeInfoDbData_1.2.2 cellranger_1.1.0      
    ##  [7] yaml_2.2.1             pillar_1.6.0           backports_1.1.5       
    ## [10] lattice_0.20-40        glue_1.4.2             digest_0.6.27         
    ## [13] RColorBrewer_1.1-2     ggsignif_0.6.0         XVector_0.26.0        
    ## [16] colorspace_2.0-0       htmltools_0.5.1.1      Matrix_1.2-18         
    ## [19] plyr_1.8.6             pkgconfig_2.0.3        broom_0.5.6           
    ## [22] haven_2.2.0            zlibbioc_1.32.0        purrr_0.3.4           
    ## [25] scales_1.1.1           openxlsx_4.1.5         rio_0.5.16            
    ## [28] tibble_3.1.1           mgcv_1.8-31            generics_0.1.0        
    ## [31] farver_2.1.0           car_3.0-8              ellipsis_0.3.1        
    ## [34] withr_2.4.2            magrittr_2.0.1         crayon_1.4.1          
    ## [37] readxl_1.3.1           evaluate_0.14          fansi_0.4.2           
    ## [40] nlme_3.1-145           rstatix_0.5.0          forcats_0.5.0         
    ## [43] foreign_0.8-76         tools_3.6.3            data.table_1.13.6     
    ## [46] hms_1.0.0              lifecycle_1.0.0        stringr_1.4.0         
    ## [49] munsell_0.5.0          ggsci_2.9              zip_2.1.1             
    ## [52] compiler_3.6.3         rlang_0.4.10           grid_3.6.3            
    ## [55] RCurl_1.98-1.3         bitops_1.0-6           labeling_0.4.2        
    ## [58] rmarkdown_2.1          gtable_0.3.0           abind_1.4-5           
    ## [61] DBI_1.1.1              reshape_0.8.8          curl_4.3              
    ## [64] R6_2.5.0               knitr_1.28             dplyr_1.0.5           
    ## [67] utf8_1.2.1             stringi_1.5.3          Rcpp_1.0.6            
    ## [70] vctrs_0.3.7            tidyselect_1.1.0       xfun_0.12
