# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 3

library(dplyr)
library(Biobase)
library(qvalue)
library(genefilter)


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
