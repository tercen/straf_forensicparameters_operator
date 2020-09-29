library(tercen)
library(dplyr)
library(adegenet)
library(pegas)
library(hierfstat)

getIndicesFromGenind <- function(data, hw = FALSE, hwperm = 1000, ploidy = "Diploid") {
  
  freq <- apply(data@tab, 2, sum, na.rm = TRUE)
  nam <- strsplit(names(freq), split="[.]")
  loc <- as.factor(unlist(lapply(nam, function(x) x[1])))
  alle <- as.numeric(unlist(lapply(nam, function(x) sub("-", ".", x[2]))))
  DAT <- data.frame(freq, loc, alle)
  N <- tapply(DAT$freq, DAT$loc, sum)
  DAT$frequency <- DAT$freq / N[DAT$loc]
  
  PIC <- NULL
  for(i in unique(loc)) {
    
    FR <- c(DAT$frequency[names(DAT$frequency) == i])
    xu <- outer(FR, FR, "fu")
    som <- sum(xu[lower.tri(xu)])
    PIC[i] <-  1 - sum(FR ^ 2) - som
    
  } 
  
  Nall <- tapply(DAT[DAT$freq>0, ]$freq, DAT[DAT$freq>0, ]$loc, length)
  GD <- tapply(DAT$frequency, DAT$loc, function(x) 1 - sum(x ^ 2))
  GD <- GD * N / (N - 1) 
  PIC <- PIC[names(GD)]
  D2 <- genind2loci(data)
  
  sumloc <- summary(D2)[names(GD)]
  
  PM1 <- lapply(sumloc, function(x) {
    sum((x$genotype / sum(x$genotype)) ^ 2)
  })
  
  PM <- unlist(PM1)
  
  DF <- data.frame(
    locus = names(GD),
    N = N,
    Nall = Nall,
    GD = GD,
    PIC = PIC,
    PM = PM,
    PD = 1 - PM
  )
  
  if(ploidy == "Diploid") {
    
    DF$Hobs <- summary(data)$Hobs[names(GD)]
    DF$PE <- (DF$Hobs ^ 2) * (1 - 2 * (DF$Hobs) * ((1 - DF$Hobs) ^ 2))
    DF$TPI <- 1 / (2 * (1 - DF$Hobs))
    
  }
  
  
  if(length(unique(data@pop)) > 1 & length(locNames(data)) > 1) {
    
    basicstat <- basic.stats(data, diploid = switch(ploidy, Diploid = TRUE, Haploid = FALSE), digits = 4)
    Fst <- wc(data, diploid = switch(ploidy, Diploid = TRUE, Haploid = FALSE))$per.loc$FST
    
    names(Fst) <- as.character(unique(data@loc.fac))
    DF$Fst <- Fst[names(GD)]
    
    DF$Ht <- basicstat$perloc[names(GD), "Ht"]
    DF$Fis <- basicstat$perloc[names(GD), "Fis"]
    
  }
  
  if(ploidy == "Diploid" & hw) {
    withProgress(message = 'Performing HW test...', value = 0, {
      DF$pHW <- hw.test(data, B = hwperm)[names(GD), 4]
    })
  } 
  
  DF$population = as.character(data$pop[1])
  
  return(DF)
}
getIndicesAllPop <- function(data,
                             hw = FALSE,
                             hwperm = 1000,
                             ploidy = "Diploid") {
  
  ind <- list()
  # ind$all <-- getIndicesFromGenind(data, hw, hwperm, ploidy)
  
  for(popu in unique(data$pop)) {
    
    ind <- c(ind, x = NA)
    mat <- getIndicesFromGenind(data[data@pop == popu, ], hw, hwperm, ploidy)
    ind$x <- mat
    
    names(ind)[length(ind)] <- popu
    
  }
  
  return(ind)
}

fu <- function(a, b) { 2 * (a ^ 2) * (b ^ 2) }


ctx = tercenCtx()

freqTAB <- ctx %>% as.matrix %>% t
which.pop <- grep("population", names(ctx$cnames))
pop <- ctx$cselect(ctx$cnames[[which.pop]])[[1]]
which.samp <- grep("sample", names(ctx$cnames))
samp <- ctx$cselect(ctx$cnames[[which.samp]])[[1]]

loc <- ctx$rselect(ctx$rnames[[1]])[[1]]
all <- ctx$rselect(ctx$rnames[[2]])[[1]]
rownames(freqTAB) <- samp
colnames(freqTAB) <- paste(loc, all, sep =  ".")

dat2 <- as.genind(tab = freqTAB)
pop(dat2) <- pop

ind <- do.call(rbind, getIndicesAllPop(dat2))

#join with loc and pop
table <- tercen::dataframe.as.table(ctx$addNamespace(ind))
table$properties$name <- 'forensic_parameters'

for(i in 1:ncol(ind)) {
  tp <- ifelse(class(ind[, i]) %in% c('factor', 'integer'), 'character', 'double')
  table$columns[[i]]$type = tp
} 

relation <- SimpleRelation$new()
relation$id <- table$properties$name

join <- JoinOperator$new()
join$rightRelation <- relation

result <- OperatorResult$new()
result$tables <- list(table)
result$joinOperators <- list(join)

ctx$save(result) 

