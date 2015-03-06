# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 1

library(downloader)

# load sample data
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




