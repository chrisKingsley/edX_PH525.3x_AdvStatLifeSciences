# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 3

library(dplyr)
library(Biobase)
library(qvalue)
library(genefilter)
library(Biobase)
library(sva)

library(GSE5859Subset)
data(GSE5859Subset)

# load admissions table
load('admissions.rda')

# load gene expression data
load("GSE5859.rda")
geneExpression <- exprs(e)
sampleInfo <- pData(e)


# Question 3.1.1
index = which(admissions$Gender==1)
maleAccepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
maleApplied = sum(admissions$Number[index])
maleAccepted/maleApplied

index = which(admissions$Gender==0)
femaleAccepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
femaleApplied = sum(admissions$Number[index])
femaleAccepted/femaleApplied


# Question 3.1.2
admTable <- matrix(c(maleAccepted, maleApplied, femaleAccepted, femaleApplied),
                   nrow=2)
chisq.test(admTable)


index = admissions$Gender==1
men = admissions[index,]
women = admissions[!index,]
print( data.frame( major=admissions[1:6,1],men=men[,3], women=women[,3]) )


# Question 3.1.3 & 3.1.4
adminTable <- data.frame(major=admissions$Major)
for(major in admissions$Major) {
    idx <- which(adminTable$major==major)
    adminTable[idx, "applied"] 
}

# calculate fraction of accepted applicants for each major
admByMajor <- admissions %>%
              group_by(Major) %>%
              summarise(applied=sum(Number),
                        accepted=sum(Number*Percent/100),
                        fracAccepted=accepted/applied)


# Question 3.1.4
# correlation of applicants per gender with fraction admitted for each major
cor(admissions$Number[admissions$Gender==1], admByMajor$fracAccepted)
cor(admissions$Number[admissions$Gender==0], admByMajor$fracAccepted)


# Question 3.2.1
sampleInfo$year <- format(sampleInfo$date,"%y")
ethPerYear <- sampleInfo %>%
    group_by(year) %>%
    summarise(numGroups=length(unique(ethnicity)))
sum(ethPerYear$numGroups > 1)


# Question 3.2.2
sampleInfo$month.year = format(sampleInfo$date,"%m%y")
ethPerMonthYear <- sampleInfo %>%
    group_by(month.year) %>%
    summarise(numGroups=length(unique(ethnicity)))
mean(ethPerMonthYear$numGroups > 1)


# Question 3.2.3
idx1 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="02")
idx2 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="03")
outcome <- rep(NA, nrow(sampleInfo))
outcome[idx1] <- 1
outcome[idx2] <- 2
pvals <- rowttests(geneExpression, as.factor(outcome))$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)
qvals$pi0


# Question 3.2.4
idx1 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="03")
idx2 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="04")
outcome <- rep(NA, nrow(sampleInfo))
outcome[idx1] <- 1
outcome[idx2] <- 2
pvals <- rowttests(geneExpression, as.factor(outcome))$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)


# Question 3.2.5
idx1 <- which(as.character(sampleInfo$ethnicity)=="CEU")
idx2 <- which(as.character(sampleInfo$ethnicity)=="ASN")
outcome <- rep(NA, nrow(sampleInfo))
outcome[idx1] <- 1
outcome[idx2] <- 2
pvals <- rowttests(geneExpression, as.factor(outcome))$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)


# Question 3.2.6
idx1 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="05")
idx2 <- which(as.character(sampleInfo$ethnicity)=="ASN" & sampleInfo$year=="05")
outcome <- rep(NA, nrow(sampleInfo))
outcome[idx1] <- 1
outcome[idx2] <- 2
pvals <- rowttests(geneExpression, as.factor(outcome))$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)


# Question 3.2.7
set.seed(3)
idx1 <- which(as.character(sampleInfo$ethnicity)=="CEU" & sampleInfo$year=="02")
idx1 <- sample(idx1, 3)
idx2 <- which(as.character(sampleInfo$ethnicity)=="ASN" & sampleInfo$year=="05")
outcome <- rep(NA, nrow(sampleInfo))
outcome[idx1] <- 1
outcome[idx2] <- 2
pvals <- rowttests(geneExpression, as.factor(outcome))$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.05)


# Question 3.3.1
sex <- as.factor(sampleInfo$group)
month <- factor(format(sampleInfo$date,"%m"))
table(sampleInfo$group, month)

pvals <- rowttests(geneExpression, sex)$p.value
qvals <- qvalue(pvals)
sum(qvals$qvalues < 0.1)


# Question 3.3.2
idx <- which(qvals$qvalues < 0.1)
sum(geneAnnotation$CHR[idx] %in% c("chrX","chrY"))/length(idx)


# Question 3.3.3
idx <- which(qvals$qvalues < 0.1 & !geneAnnotation$CHR %in% c("chrX","chrY"))
pvals <- rowttests(geneExpression[idx,], month)$p.value
sum(pvals < 0.05)/length(idx)


# Question 3.3.5
X <- model.matrix(~ sex + month)
pvals <- t( sapply(1:nrow(geneExpression), function(j) {
    y <- geneExpression[j,]
    fit <- lm(y ~ X - 1)
    summary(fit)$coef[2:3,4]
}) )
qvals <- qvalue(pvals[,1])
sum(qvals$qvalues < 0.1)


# Question 3.3.6
idx <- which(qvals$qvalues < 0.1)
sum(geneAnnotation$CHR[idx] %in% c("chrX","chrY"))/length(idx)

qvals <- qvalue(pvals[,2])
sum(qvals$qvalues < 0.1)


# Question 3.4.1
y <- geneExpression - rowMeans(geneExpression)
# samples ordered by gender
image(cor(y))
# samples ordered by date
idx <- order(sampleInfo$date)
image(cor(y[, idx]))


# Question 3.4.2 - 3.4.3
sex <- as.factor(sampleInfo$group)
month <- factor(format(sampleInfo$date,"%m"))
pcs <- svd(y)
plot(pcs$v[idx,1], col=month[idx])


# Question 3.4.4
d2 <- pcs$d^2
sum(d2/sum(d2) > 0.1)


# Question 3.4.5
pcMonthCor <- cor(pcs$v, as.numeric(month))
which(pcMonthCor==max(abs(pcMonthCor)))
max(abs(pcMonthCor))


# Question 3.4.6
pcGenderCor <- cor(pcs$v, sampleInfo$group)
which(pcGenderCor==max(abs(pcGenderCor)))
max(abs(pcGenderCor))


# Question 3.4.7
X <- model.matrix(~sex + pcs$v[,1:2])
pvals <- sapply(1:nrow(geneExpression), function(j) {
    y <- geneExpression[j,]
    fit <- lm(y ~ X - 1)
    summary(fit)$coef[2,4]
})

qvals <- qvalue(pvals)
idx <- which(qvals$qvalues < 0.1)
length(idx)
table(geneAnnotation$CHR[idx])


# Question 3.5.1
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group, s$v[,1])

sex <- sampleInfo$group
month <- factor(format(sampleInfo$date,"%m"))
mod <- model.matrix(~sex)
svafit <- sva(geneExpression, mod)
head(svafit$sv)

for(i in 1:ncol(svafit$sv)){
    print( cor(s$v[,i],svafit$sv[,i]) )
}

X <- model.matrix(~sex + svafit$sv)
pvals <- sapply(1:nrow(geneExpression), function(j) {
    y <- geneExpression[j,]
    fit <- lm(y ~ X - 1)
    summary(fit)$coef[2,4]
})

qvals <- qvalue(pvals)
idx <- which(qvals$qvalues < 0.1)
length(idx)
table(geneAnnotation$CHR[idx])
