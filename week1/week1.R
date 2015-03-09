# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 1

library(downloader)
library(genefilter)
library(devtools)

# load sample data
#install_github("genomicsclass/GSE5859Subset")
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
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
    control = sample(population[,1],12)
    treatment = sample(population[,1],12)
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
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.val
})

sum(pvals<0.05)


# Question 1.2.4
set.seed(100)
sumSig_pvals <- replicate(1000, {
    pvals <- replicate(20, {
        cases = rnorm(10,30,2)
        controls = rnorm(10,30,2)
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
