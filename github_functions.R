
# 1. Pipeline

### missing oligo count
mOligoCount <- function(cur.count.table, cur.sample, threshold) {
  sel.table <- subset(cur.count.table, grepl(cur.sample, sample))
  
  for (i in 1:length(threshold)) {
    temp <- subset(sel.table, raw<=threshold[i])
    temp <- table(temp$sample)
    temp <- as.data.frame.table(temp, stringsAsFactors = F)
    colnames(temp) <- c("sample", "count")
    if (i==1) {
      cur.res <- temp
    } else {
      if (all(cur.res$sample!=temp$sample)) {cat("Error"); return;}
      cur.res <- cbind(cur.res, temp[,2])
    }
  }
  
  colnames(cur.res)[-1] <- threshold
  return(cur.res)
}

### Correlation coefficients between two biological replicates 
### for variable barcode numbers
### perform n simulations for a given number of barcode numbers
simCorrLog2FC <- function(repList, bNumList, n) {
  # repList : c("Rep1", "Rep2")
  # bNumList: a vector of the number of barcodes
  # n: simulation/random sampling number
  out <- matrix(NA, nrow=length(bNumList), ncol=n, dimnames = list(bNumList, 1:n))
  
  for (i in 1:length(bNumList)) {
    #cat("Barcode number= ", i, "\n")
    barcode.random <- combn(1:BARCODE_NUM, bNumList[i])
    combIdx <- sample(1:ncol(barcode.random), n)
    for (j in 1:length(combIdx)) {
      cur.barcode <- BARCODE_LIST[barcode.random[, combIdx[j]]]
      simul1 <- subset(count.table, replicate==repList[1] & barcodeID %in% cur.barcode)
      simul1 <- dcast(simul1, tileID+barcodeID~pulldown, value.var = "rpm")
      simul1$log2fc <- log2(simul1$IP+1) - log2(simul1$Input+1)
      simul1 <- dcast(simul1, tileID~., value.var = "log2fc", fun.aggregate = mean)
      
      simul2 <- subset(count.table, replicate==repList[2] & barcodeID %in% cur.barcode)
      simul2 <- dcast(simul2, tileID+barcodeID~pulldown, value.var = "rpm")
      simul2$log2fc <- log2(simul2$IP+1) - log2(simul2$Input+1)
      simul2 <- dcast(simul2, tileID~., value.var = "log2fc", fun.aggregate = mean)
      
      if (any(simul1[,1]!=simul2[,1])) {cat("Error"); return;}
      out[i, j] <- cor(simul1[,2], simul2[,2])
    }
  }
  return(out)
}

### QC
### Correlation coefficients between two biological replicates of rpm

calMeanLog2fc <- function(r, f) {
  rpm.input <- subset(count.table, pulldown=="Input" & replicate==r)
  rpm.input <- dcast(rpm.input, tileID~barcodeID, value.var = "rpm")
  rpm.ip <- subset(count.table, pulldown==f & replicate==r)
  rpm.ip <- dcast(rpm.ip, tileID~barcodeID, value.var = "rpm")
  log2fc <- log2(rpm.ip[,-1]+1) - log2(rpm.input[,-1]+1)
  return(rowMeans(log2fc))
}

# 2. Enrichment

oligoFilter <- function(cur.oligo.count, inputFilterBy="rpm", inputTh=0, ipFilterBy="rpm", ipTh=0) {
  cur.oligo.count$oligoID <- paste(cur.oligo.count$tileID, oligo.count$barcodeID, sep="_")
  
  # Input filter
  temp <- subset(cur.oligo.count, pulldown=="Input" & get(inputFilterBy)>=inputTh)
  temp <- table(temp$oligoID)
  temp <- names(temp)[temp==REPLICATE_NUM] # oligoIDs should be expressed in all the replicates.
  #temp <- names(temp)[temp>=2] # oligoIDs should be expressed in at least 2 samples
  cur.oligo.count <- subset(cur.oligo.count, oligoID %in% temp)
  
  # IP filter
  temp <- subset(cur.oligo.count, pulldown=="IP" & get(ipFilterBy)>=ipTh)
  temp <- table(temp$oligoID)
  temp <- names(temp)[temp==REPLICATE_NUM] # oligoIDs should be expressed in all the replicates.
  #temp <- names(temp)[temp>=2] # oligoIDs should be expressed in at least 2 samples
  cur.oligo.count <- subset(cur.oligo.count, oligoID %in% temp)
  
  cat("Proportion of oligos after filtering: ", length(unique(cur.oligo.count$oligoID))/(nrow(tile.annot)*BARCODE_NUM), "\n")
  cur.oligo.count$oligoID <- NULL
  return(cur.oligo.count)
}


countPooler <- function(cur.oligo.count) {
  cur.oligo.count$sampleID <- paste(cur.oligo.count$repl, cur.oligo.count$pulldown, sep="_")
  
  # Pooling
  oligo.pooled <- dcast(cur.oligo.count, tileID~sampleID, value.var="raw", fun.aggregate = sum)
  rownames(oligo.pooled) <- oligo.pooled$tileID
  pooled.count <- data.matrix(oligo.pooled[,-1]) # Output as a matrix 
  
  cat("Proportion of tiles after filtering and pooling: ", nrow(pooled.count)/nrow(tile.annot), "\n")
  return(pooled.count)
}



### GLM
runGLMtest <- function(BARCODE_LIST, REPLICATE_LIST, IP, INPUT_RPM_TH, IP_COUNT_TH) {
  require(MASS)
  glm.out <- data.frame(tileID=tile.table$tileID, n=NA, rpm=NA, poisson.beta=NA, poisson.sd=NA, poisson.zval=NA, poisson.pval=NA, poisson.devPval=NA, poisson.overdisp=NA, poisson.overdisp.pval=NA, nb.beta=NA, nb.sd=NA, nb.zval=NA, nb.pval=NA, nb.devPval=NA, nb.alpha=NA, resid.diff=NA)
  
  for (i in 1:nrow(tile.table)) {
    #print(i)
    
    cur.tileID <- glm.out$tileID[i]
    
    # prepare X
    X <- data.frame(barcodeID = rep(BARCODE_LIST, each=length(REPLICATE_LIST)*2),
                    replicate = rep(rep(REPLICATE_LIST, each=2), length(BARCODE_LIST)),
                    pulldown = rep(rep(c("Input", IP), length(REPLICATE_LIST)), length(BARCODE_LIST)),
                    stringsAsFactors = F)
    X$oligoID <- paste(cur.tileID, X$barcodeID, sep="_")
    
    # prepare y
    y <- subset(count.table, tileID==cur.tileID & barcodeID %in% BARCODE_LIST & replicate %in% REPLICATE_LIST & pulldown %in% c("Input", IP))
    
    # filter barcodes by mean input rpm across replicates.
    # then at least 3 barcodes should exist in order to avoid bias by barcode.
    temp <- dcast(subset(y, pulldown=="Input"), barcodeID~., value.var = "rpm", fun.aggregate = mean)
    y <- subset(y, barcodeID %in% temp[temp[, 2]>=INPUT_RPM_TH, 1])
    if (length(unique(y$barcodeID))<3) {next; }
    
    # filter barcodes by merged IP raw count.
    # then at least 3 barcodes should exist in order to avoid bias by barcode.
    #temp <- dcast(subset(y, pulldown==IP), barcodeID~., value.var = "raw", fun.aggregate = sum)
    #y <- subset(y, barcodeID %in% temp[temp[, 2]>=IP_COUNT_TH, 1])
    #if (length(unique(y$barcodeID))<3) {next; }
    
    # filter by individual sample: every sample (repID*pulldown) should have at least 2 counted barcodes.
    #temp <- subset(y, raw>0)
    #temp <- table(temp$repID, temp$pulldown)
    #if (min(temp) < 2) {next;}
    
    glm.out$n[i] = nrow(y)
    glm.out$rpm[i] = mean(y$rpm)
    
    # prepare glm input by matching X and y
    idx <- match(paste(y$oligoID, y$replicate, y$pulldown, sep="_"), paste(X$oligoID, X$replicate, X$pulldown, sep="_"))
    cur.data <- data.frame(X[idx, ], count=y$normalized, stringsAsFactors = F) # use normalized count as y
    cur.data$pulldown <- factor(cur.data$pulldown, levels=c("Input", IP))
    cur.data$barcodeID <- factor(cur.data$barcodeID)
    cur.data$replicate <- factor(cur.data$replicate)
    
    # fitting
    fit.poisson <- tryCatch( glm(round(count) ~ replicate + barcodeID + pulldown, family=poisson, data=cur.data), error = function(err) { print("fit.poisson"); print(err); return("err")})
    
    # test
    if (is.list(fit.poisson)) {
      residual.poisson <- round(cur.data$count)-fit.poisson$fitted.values  
      residual.z.poisson <- residual.poisson / sqrt(fit.poisson$fitted.values)
      glm.out$poisson.overdisp[i] <- sum(residual.z.poisson^2) / fit.poisson$df.residual
      glm.out$poisson.overdisp.pval[i] <- 1-pchisq(sum(residual.z.poisson^2), fit.poisson$df.residual)
      
      temp <- summary(fit.poisson)
      glm.out$poisson.devPval[i] <- pchisq(temp$deviance, df=temp$df.residual, lower.tail=FALSE)
      
      temp <- temp[["coefficients"]][paste0("pulldown", IP),]
      glm.out$poisson.beta[i] <- temp[1] # "Estimate"
      glm.out$poisson.sd[i] <- temp[2] # "Std. Error"
      glm.out$poisson.zval[i] <- temp[3] # "z value"
      glm.out$poisson.pval[i] <- temp[4] # "Pr(>|z|)"
    }
    
  }
  return(glm.out)
}

# 3. Mutation

### A table of log2fc for a given stem
stem_log2FC_table <- function(sID) { # stemID
  idx <- which(stem.table$stemID==sID)
  cur.count <- subset(count.table, tileID %in% stem.table[idx, grep("tileID", colnames(stem.table))])
  cur.count <- rbind(cur.count, subset(count.table, tileID == subset(tile.table, group=="WT")$tileID))
  
  # filter barcodes by input rpm
  #rpm.th <- 1
  #temp <- dcast(subset(cur.count, pulldown=="Input"), barcodeID~., value.var = "rpm", fun.aggregate = mean)
  #cur.count <- subset(cur.count, barcodeID %in% temp[temp[, 2]>=rpm.th, 1])
  #length(unique(y$barcodeID))<3
  
  out <- dcast(cur.count, tileID+barcodeID+replicate~pulldown, value.var="rpm")
  out <- subset(out, Input>=1 | HA>=1)
  out$log2FC <- log2(out$HA+1)-log2(out$Input+1)
  #out$log2FC <- log2(out$HA/out$Input)
  
  temp <- subset(stem.table, stemID==sID)
  temp <- c(as.character(subset(tile.table, group=="WT")$tileID), unlist(c(temp[grep("tileID.mut5p", colnames(temp))], temp[grep("tileID.mut3p", colnames(temp))], temp[grep("tileID.comp", colnames(temp))])))
  
  out$tileID <- factor(out$tileID, levels=temp, labels=c("WT", "5' mut", "3' mut", "Comp. mut"))
  return(out)
}

### A table of log2fc for a given position
single_log2FC_table <- function(loc) { # location
  idx <- which(single.table$mutLoc==loc)
  cur.count <- subset(count.table, tileID %in% single.table[idx, grep("tileID", colnames(single.table))])
  cur.count <- rbind(cur.count, subset(count.table, tileID == subset(tile.table, group=="WT")$tileID))
  
  # filter barcodes by input rpm
  #rpm.th <- 1
  #temp <- dcast(subset(cur.count, pulldown=="Input"), barcodeID~., value.var = "rpm", fun.aggregate = mean)
  #cur.count <- subset(cur.count, barcodeID %in% temp[temp[, 2]>=rpm.th, 1])
  #length(unique(y$barcodeID))<3
  
  out <- dcast(cur.count, tileID+barcodeID+replicate~pulldown, value.var="rpm")
  out <- subset(out, Input>=1 | HA>=1)
  out$log2FC <- log2(out$HA+1)-log2(out$Input+1)
  #out$log2FC <- log2(out$HA/out$Input)
  
  temp <- c(as.character(subset(tile.table, group=="WT")$tileID), as.character(single.table$tileID[idx]))
  
  out$tileID <- factor(out$tileID, levels=temp, labels=c(paste("WT", unique(single.table$seqWT[idx]), sep=":"), single.table$seqMT[idx]))
  return(out)
}

### hTR log2FC table for a given stem
hTR_stem_log2FC_table <- function(cur.type, cur.coord, cur.loc) {
  idx <- which(stem.table$type==cur.type & stem.table$coord==cur.coord & stem.table$bpLocPair==cur.loc)
  temp <- stem.table[idx, grep("tileID", colnames(stem.table))]
  names(temp) <- sapply(strsplit(colnames(temp), split="[.]"), "[[", 2)
  
  cur.count <- subset(count.table, tileID %in% temp)
  cur.count$tileID <- factor(as.character(cur.count$tileID), levels=temp, labels=names(temp))
  
  out <- dcast(cur.count, tileID+barcodeID+replicate~pulldown, value.var="rpm")
  out <- out[apply(out[, c("Input", "IP")], 1, max)>=1, ]
  out$log2FC <- log2(out$IP)-log2(out$Input)
  
  out$tileID <- factor(as.character(out$tileID), levels=c("wt", "mut5p", "mut3p", "comp"), labels=c("WT", "5' mut", "3' mut", "Comp. mut"))
  return(out)
}



