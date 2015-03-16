# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 1


library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)

# Question 2.1.1
table(tissue)

# Question 2.1.2
d <- as.matrix(dist(t(e)))
d[3,45]

# Question 2.1.3
a <- e["210486_at",]
b <- e["200805_at",]
sqrt( crossprod(a-b) )

# Question 2.1.4
nrow(e)^2

# Question 2.1.5
length(dist(t(e)))


# Question 2.2.1

