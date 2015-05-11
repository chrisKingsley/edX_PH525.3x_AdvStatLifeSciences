# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 1

library(downloader)
library(genefilter)
library(devtools)
library(qvalue)

# load sample data
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset) 


# Question 1.1.1
sum(sampleInfo$date=="2005-06-27")

# Question 1.1.2
sum(geneAnnotation$CHR=="chrY", na.rm=T)

# Question 1.1.3
sample.idx <- which(sampleInfo$date=="2005-06-10")
gene.idx <- which(geneAnnotation$SYMBOL=="ARPC1A")
geneExpression[gene.idx, sample.idx]


set.seed(1)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population <- read.csv(filename)
pvals <- replicate(1000,{
    control <- sample(population[,1],12)
    treatment <- sample(population[,1],12)
    t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)


# Question 1.2.1
sum(pvals<0.05)/length(pvals)

# Question 1.2.2
sum(pvals<0.01)/length(pvals)


# Question 1.2.3
set.seed(100)
pvals <- replicate(20, {
    cases <- rnorm(10,30,2)
    controls <- rnorm(10,30,2)
    t.test(cases,controls)$p.val
})

sum(pvals<0.05)


# Question 1.2.4
set.seed(100)
sumSig_pvals <- replicate(1000, {
    pvals <- replicate(20, {
        cases <- rnorm(10,30,2)
        controls <- rnorm(10,30,2)
        t.test(cases,controls)$p.val
    })
    sum(pvals<0.05)
})
mean(sumSig_pvals)

# Question 1.2.5
sum(sumSig_pvals>0)/length(sumSig_pvals)


# Question 1.3.2
1-(1-0.05)^8793

# Question 1.3.3
1-(1-0.05)^(1/8793)


# Question 1.4.1
m <- 1000
alpha <- seq(0.01,0.25,0.01)
plot(alpha/m, 1-(1-alpha)^(1/m), pch=20, xlab="Bonferroni", ylab="Sidak")
abline(0, 1, col="red")

# Question 1.4.2
set.seed(1)
sumSigPvals <- replicate(10000, {
    pvals <- runif(8793,0,1)
    sum(pvals < 0.05/8753)
})
mean(sumSigPvals)

# Question 1.4.2
set.seed(1)
threshold <- 1-(1-0.05)^(1/8793)
sumSigPvals <- replicate(10000, {
    pvals <- runif(8793,0,1)
    sum(pvals < threshold)
})
mean(sumSigPvals)


# Question 1.5.1
tests <- rowttests(geneExpression, factor(sampleInfo$group))
sum(tests$p.value < 0.05)

# Question 1.5.2
sum(tests$p.value < 0.05/nrow(geneExpression))

# Question 1.5.3
qVals <- p.adjust(tests$p.value, method="fdr")
sum(qVals < 0.05)

# Question 1.5.4
qVals <- qvalue(tests$p.value)
sum(qVals$qvalues < 0.05)

# Question 1.5.5
qVals$pi0

# Question 1.5.7
n <- 24
m <- 8793
delta <- 2
positives <- 1:500
negatives <- 501:m
outcome <- factor(c(rep(0, n/2), rep(1, n/2)))

set.seed(1)
simData <- t(replicate(1000, {
    mat <- matrix(rnorm(n*m),m,n)
    mat[positives, 1:(n/2)] <- mat[positives, 1:(n/2)]+delta
    
    tests <- rowttests(mat, outcome)
    
    bonf_sig <- tests$p.value < 0.05/nrow(mat)
    bonf_fp <- intersect(which(bonf_sig), negatives)
    bonf_fn <- intersect(which(!bonf_sig), positives)
    
    bh_sig <- p.adjust(tests$p.value, method="fdr") < 0.05
    bh_fp <- intersect(which(bh_sig), negatives)
    bh_fn <- intersect(which(!bh_sig), positives)
      
    storey_sig <- qvalue(tests$p.value)$qvalues < 0.05
    storey_fp <- intersect(which(storey_sig), negatives)
    storey_fn <- intersect(which(!storey_sig), positives)
    
    c(length(bonf_fp)/length(negatives), length(bonf_fn)/length(positives),
      length(bh_fp)/length(negatives), length(bh_fn)/length(positives),
      length(storey_fp)/length(negatives), length(storey_fn)/length(positives)) 
}))

colnames(simData) <- c("bonf_fpr","bonf_fnr","bh_fpr","bh_fnr",
                       "storey_fpr","storey_fnr")
apply(simData, 2, mean)


# Question 1.6.1
library(SpikeInSubset)
data(mas133)

e <- exprs(mas133)
plot(e[,1], e[,2] ,main=paste0("corr=",signif(cor(e[,1],e[,2]),3)), cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

sum(e[,1]<k & e[,2]<k & e[,1]> -b & e[,1]> -b)/nrow(e)

# Question 1.6.2
plot(log2(e[,1]), log2(e[,2]),
     main=paste0("corr=", signif(cor(log2(e[,1]),log2(e[,2])),2)), cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

# Question 1.6.3
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2, e[,2]-e[,1], cex=0.5)

sd(e[,2]-e[,1])
sum(abs(e[,2]-e[,1]) > 1)
