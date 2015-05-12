# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 3

library(dplyr)

library(dagdata)
data(admissions)

# load admissions table
load('admissions.rda')


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
