# =================================
# 1. Enrichment
# =================================

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

# =================================
# 2. Mutation analyses
# =================================

runLimma.singleMut <- function(cur.single.table, cur.repl, rpm.th) {
  # the required columns in cur.single.table:
  # "mutID": mutation ID,
  # "tileID.mut": tileID for mutation,
  # "tileID.wt": tileID for WT.
  
  # oligo.count should be accessible globally.
  # the requited columns in oligo.count:
  # repl, pulldown, tileID, barcodeID, normed
  
  # 0. Check pulldown level
  if (length(table(oligo.count$pulldown))!=2) {stop("pulldown should have two levels of Input and IP.")}
  
  # 1. construct y and X
  for (i in 1:nrow(cur.single.table)) {
    cur.tileID <- setNames(as.character(cur.single.table[i, c("tileID.mut", "tileID.wt")]), c("mut", "wt"))
    
    temp <- lapply(cur.tileID, function(tID, sID=cur.repl) {subset(oligo.count, tileID==tID & repl==sID)[, c("normed", "barcodeID", "pulldown", "rpm")]})
    names(temp) <- names(cur.tileID)
    
    # (1) filter barcodes.
    temp <- lapply(temp, function(x) {temp0 <- table(x$barcodeID[which(x$rpm<rpm.th)]); ex.barcodes <- names(temp0)[which(temp0>=1)]; if (length(ex.barcodes)>0) {x$normed[which(x$barcodeID %in% ex.barcodes)]<-NA}; return(x)})
    
    # (2) y: y is normalized values, not RPM.
    cur.y <- as.vector(sapply(temp, function(x) {x$normed}))
    
    # (3) X
    cur.X <- data.frame(barcodeID=as.vector(sapply(temp, function(x) {x$barcodeID})), pulldown=as.vector(sapply(temp, function(x) {x$pulldown})))
    cur.X$mutation <- rep(names(temp), times=sapply(temp, nrow))
    
    if (i==1) {
      y <- cur.y
      X <- cur.X
    } else {
      y <- rbind(y, cur.y)
      if (!all(X==cur.X)) stop("X is different between mutations")
    }
  }  
  rownames(y) <- cur.single.table$mutID
  
  cat("Check the number of pseudo (barcode) samples: \n")
  print(table(X$mutation, X$pulldown))
  cat(" ==> Should be equal to:", BARCODE_NUM,"\n")
  
  # 2. Construct the linear model.
  X$mutation <- factor(X$mutation, levels=c("wt", "mut"), labels=c("WT", "MT"))
  
  limma.model <- model.matrix(data= X, ~ 0 + pulldown*mutation)
  # colnames(limma.model)
  
  # 3. Limma
  limma.fit <- limma::lmFit(log(y+2), limma.model) # limma cpm default is 2.
  limma.fit <- limma::eBayes(limma.fit, trend=TRUE)
  
  # Results
  limma.table <- limma::topTable(limma.fit, coef="pulldownIP:mutationMT", sort.by = "none", n = Inf)
  
  if(!(all(rownames(limma.table) == cur.single.table$mutID))) {stop("limma output is different from input in order.")}
  return(list(X=X, y=y, limma=limma.table))
}

runLimma.stemMut <- function(cur.stem.table, cur.repl, rpm.th) {
  # the required columns in cur.stem.table:
  # "stemID",
  # "tileID.wt",
  # "tileID.mut5p",
  # "tileID.mut3p",
  # "tileID.comp".
  
  # oligo.count should be accessible globally.
  # the requited columns in oligo.count:
  # repl, pulldown, tileID, barcodeID, normed
  
  # 0. Check pulldown level
  if (length(table(oligo.count$pulldown))!=2) {stop("pulldown should have two levels of Input and IP.")}
  
  # 1. construct y and X
  for (i in 1:nrow(cur.stem.table)) {
    cur.tileID <- setNames(as.character(cur.stem.table[i, c("tileID.mut5p", "tileID.mut3p", "tileID.comp", "tileID.wt")]), c("mut5p", "mut3p", "comp", "wt"))
    
    temp <- lapply(cur.tileID, function(tID, sID=cur.repl) {subset(oligo.count, tileID==tID & repl==sID)[, c("rpm", "normed", "barcodeID", "pulldown")]})
    names(temp) <- names(cur.tileID)
    
    # (1) filter barcodes.
    temp <- lapply(temp, function(x) {temp0 <- table(x$barcodeID[which(x$rpm<rpm.th)]); ex.barcodes <- names(temp0)[which(temp0>=1)]; if (length(ex.barcodes)>0) {x$normed[which(x$barcodeID %in% ex.barcodes)]<-NA}; return(x)})
    
    # (2) y: y is normalized values, not RPM.
    cur.y <- as.vector(sapply(temp, function(x) {x$normed}))
    
    # (3) X
    cur.X <- data.frame(barcodeID=as.vector(sapply(temp, function(x) {x$barcodeID})), pulldown=as.vector(sapply(temp, function(x) {x$pulldown})))
    cur.X$mutation <- rep(names(temp), times=sapply(temp, nrow))
    
    if (i==1) {
      y <- cur.y
      X <- cur.X
    } else {
      y <- rbind(y, cur.y)
      if (!all(X==cur.X)) stop("X is different between stems")
    }
  }  
  rownames(y) <- cur.stem.table$stemID
  cat("Check the number of pseudo (barcode) samples: \n")
  print(table(X$mutation, X$pulldown))
  cat(" ==> Should be equal to:", BARCODE_NUM,"\n")
  
  # 2. construct a linear model.
  
  # coding of stem effects
  X$seq5pEffect <- 0
  X$seq5pEffect[which(X$mutation %in% c("wt", "mut3p"))] <- 1
  X$seq3pEffect <- 0
  X$seq3pEffect[which(X$mutation %in% c("wt", "mut5p"))] <- 1
  X$stemEffect <- 0
  X$stemEffect[which(X$mutation %in% c("wt", "comp"))] <- 1
  X$seqEffect <- 0
  X$seqEffect[which(X$mutation %in% c("wt"))] <- 1
  
  X$seq5pEffect <- factor(X$seq5pEffect, levels=c(0,1))
  X$seq3pEffect <- factor(X$seq3pEffect, levels=c(0,1))
  X$stemEffect <- factor(X$stemEffect, levels=c(0,1))
  X$seqEffect <- factor(X$seqEffect, levels=c(0,1))
  
  # linear model and limma
  limma.model <- model.matrix(data= X, ~ 0 + pulldown * seq5pEffect + pulldown * seq3pEffect + pulldown * stemEffect)

  # 3. limma
  limma.fit <- limma::lmFit(log(y+2), limma.model) # limma cpm default is 2.
  limma.fit <- limma::eBayes(limma.fit, trend=TRUE)
  
  # Results
  limma.table <- limma::topTable(limma.fit, coef="pulldownIP:stemEffect1", sort.by = "none", n = Inf)
  if(!(all(rownames(limma.table) == cur.stem.table$stemID))) {stop("limma output is different from input in order.")}
  return(list(X=X, y=y, limma=limma.table))
}


mergeLimmaResults <- function(input.list) {
  # input.list: list of limma outputs
  # example: list(stem.rep1$limma, stem.rep2$limma, stem.rep3$limma)
  
  for (i in 2:length(input.list)) {
    if (!(all(rownames(input.list[[1]]) == rownames(input.list[[i]])))) {
      stop("ids are incompatible between input elements")
    }
  }
  
  out.logFC <- sapply(1:length(input.list), function(i) {input.list[[i]]$logFC})
  colnames(out.logFC) <- paste0("rep", 1:length(input.list), ".logFC")
  
  out.t <- sapply(1:length(input.list), function(i) {input.list[[i]]$t})
  colnames(out.t) <- paste0("rep", 1:length(input.list), ".stat")
  
  cat("Note that rank-based p-value assumption: ordinal scores follow an approximate normal distribution.")
  par(mfrow=c(1,ncol(out.t)))
  for (i in 1:ncol(out.t)) {
    qqnorm(out.t[,i], main=paste("Rep", i)); qqline(out.t[,i], col="red")  
  }
  
  out.test <- WGCNA::rankPvalue(out.t)
  
  out <- data.frame(out.logFC, out.t, row.names = rownames(input.list[[i]]))
  out$mean.logFC <- rowMeans(out.logFC)
  out$pval <- out.test$pValueExtremeScale
  out$qval <- out.test$qValueExtremeScale
  return(out)
}

stem_log2FC_tabler <- function(sID) {
  idx <- which(stem.table$stemID==sID)
  cur.count <- subset(oligo.count, tileID %in% as.character(stem.table[idx, grep("tileID", colnames(stem.table))]))
  
  out <- dcast(cur.count, tileID+barcodeID+repl~pulldown, value.var="rpm")
  out <- subset(out, Input>=1 & IP>=1) # Note !!
  out$log2FC <- log2(out$IP+1)-log2(out$Input+1)
  
  temp <- subset(stem.table, stemID==sID)
  temp <- unlist(c(temp[grep("tileID.wt", colnames(temp))], c(temp[grep("tileID.mut5p", colnames(temp))], temp[grep("tileID.mut3p", colnames(temp))], temp[grep("tileID.comp", colnames(temp))])))
  
  out$tileID <- factor(out$tileID, levels=temp, labels=c("WT", "5' mut", "3' mut", "Comp. mut"))
  return(out)
}

single_log2FC_tabler <- function(cur.table, tID) {
  idx <- which(cur.table$tileID==tID)
  cur.count <- subset(count.table, tileID %in% as.character(cur.table[idx, grep("tileID", colnames(cur.table))]))
  
  out <- dcast(cur.count, tileID+barcodeID+replicate~pulldown, value.var="rpm")
  out <- subset(out, Input>=1 & IP>=1) # Note !!
  out$log2FC <- log2(out$IP+1)-log2(out$Input+1)
  
  temp <- subset(cur.table, tileID==tID)
  temp <- unlist(c(temp[grep("tileID.wt", colnames(temp))], c(temp[grep("tileID.mut", colnames(temp))])))
  
  out$tileID <- factor(out$tileID, levels=temp, labels=c("WT", "Mut"))
  return(out)
}
