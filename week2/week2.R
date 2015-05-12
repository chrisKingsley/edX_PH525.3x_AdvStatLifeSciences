# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 2

library(svd)
library(gplots)
library(genefilter)
library(devtools)
library(caret)
library(matrixStats)
library(class)

# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)

library(GSE5859Subset)
data(GSE5859Subset)


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
s <- svd(e)
m <- rowMeans(e)
cor(m, s$u[,1])


# Question 2.2.2
y <- e - rowMeans(e)
s2 <- propack.svd(y)


# Question 2.2.3
z <- s2$d * t(s2$v)
sqrt(crossprod(e[,3]-e[,45]))- sqrt(crossprod(z[1:2,3]-z[1:2,45]))


# Question 2.2.4
for(i in 2:10) {
    diff <- abs(sqrt(crossprod(e[,3]-e[,45])) - sqrt(crossprod(z[1:i,3]-z[1:i,45])))
    print(sprintf("dim:%d  diff:%f  error:%f", i, diff, 
                  diff/sqrt(crossprod(e[,3]-e[,45]))))
}


# Question 2.2.5
distances1 <- sqrt(apply(e[,-3]-e[,3],2,crossprod))
distances2 <- sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
cor(distances1, distances2, method = "spearman")


# Question 2.3.1
ftissue <- factor(tissue)
mypar2(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)

d <- dist(t(e))
mds <- cmdscale(d)
cor(z[1,], mds[,1])
cor(z[2,], mds[,2])


# Question 2.3.2
ftissue <- factor(tissue)
par(mfrow=c(1,3))
plot(-z[1,], -z[2,],col=as.numeric(ftissue))
plot(z[1,], z[2,],col=as.numeric(ftissue))
#legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
plot(mds[,1], mds[,2],col=as.numeric(ftissue))


# Question 2.3.3
s <- svd(geneExpression-rowMeans(geneExpression))
z <- s$d * t(s$v)
corrMat <- data.frame(dim=1:nrow(z),
                      cor=apply(z, 1, function(x) cor(x, sampleInfo$group)))
corrMat <- corrMat[order(corrMat$cor, decreasing=T),]


# Question 2.3.4
month <- format(sampleInfo$date, "%m")
month <- factor(month)
corrMat <- data.frame(dim=1:nrow(z),
                      cor=apply(z, 1, function(x) cor(x, as.numeric(month))))
corrMat <- corrMat[order(corrMat$cor, decreasing=T),]


# Question 2.3.5
boxplot(s$u[,6] ~ geneAnnotation$CHR)
boxplot(s$u[,6] ~ geneAnnotation$CHR, outline=F)


# Question 2.4.1
set.seed(1)
m <- 10000
n <- 24
x <- matrix(rnorm(m*n),m,n)
colnames(x)<-1:n

d <- dist(t(x))
plot(hclust(d))


# Question 2.4.2
set.seed(1)
m <- 10000
n <- 24

# get SE of number of clusters produced at cutheight of 143
numClusts <- replicate(100, {
    x <- matrix(rnorm(m*n),m,n)
    colnames(x) <- 1:n
    d <- dist(t(x))
    max(cutree(hclust(d), h=143))
})
sd(numClusts)


# Question 2.4.3
set.seed(10)
kclust <- kmeans(t(geneExpression), centers=5)

# examine relationship between cluster and different sample variables
table(kclust$cluster, sampleInfo$group)
table(kclust$cluster, sampleInfo$date)
table(kclust$cluster, sampleInfo$ethnicity)


# Question 2.5.1
rowIdx <- order(rowMads(geneExpression), decreasing=T)[1:25]
hmap <- heatmap.2(geneExpression[rowIdx,], trace="none", scale="row", key=F,
                  labCol=sampleInfo$date, labRow=geneAnnotation$CHR[rowIdx],
                  ColSideColors=as.character(sampleInfo$group))


# Question 2.5.2
set.seed(17)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
x <- matrix(rnorm(m*n),m,n)
g <- factor(sampleInfo$g )

row.sdIdx <- order(rowSds(x), decreasing=T)[1:50]
heatmap.2(x[row.sdIdx,], trace="none", scale="row", key=F, main="SD Filter",
          ColSideColors=as.character(g))

row.ttestIdx <- order(rowttests(x, g)$p.value)[1:50]
heatmap.2(x[row.ttestIdx,], trace="none", scale="row", key=F, 
          main="T-test Filter", ColSideColors=as.character(g))


# Question 2.6.1
n <- 10000
set.seed(1)
men <- rnorm(n,176,7) #height in centimeters
women <- rnorm(n,162,7) #height in centimeters
y <- c(rep(0,n),rep(1,n))
x <- round(c(men,women))
##mix it up
ind <- sample(seq(along=y))
y <- y[ind]
x <- x[ind]

mean(y[which(x==176)])


# Question 2.6.2
probs <- sapply(160:178, function(z) mean(y[which(x==z)]))
plot(160:178, probs)
abline(h=0.5, col="red")
max(seq(160,178)[probs >= 0.5])


# Question 2.7.1
set.seed(5)
N <- 250
ind <- sample(length(y),N)
Y <- y[ind]
X <- x[ind]

loessModel <- loess(formula = Y ~ X)
predict(loessModel, newdata=168)


# Question 2.7.2
set.seed(5)
predictions <- replicate(1000, {
    N <- 250
    ind <- sample(length(y),N)
    Y <- y[ind]
    X <- x[ind]
    
    loessModel <- loess(formula = Y ~ X)
    predict(loessModel, newdata=168)
})

sd(predictions)


# Question 2.8.1
y <- factor(sampleInfo$group)
X <- t(geneExpression)
out <- which(geneAnnotation$CHR %in% c("chrX","chrY"))
X <- X[,-out]

set.seed(1)
folds <- createFolds(y, k=10)


# Question 2.8.2
m <- 8
gene.t.tests <- colttests(X[ -folds[[2]], ], y[ -folds[[2]] ])
col.idx <- order(gene.t.tests$p.value)[1:m]

pred <- knn(train=X[ -folds[[2]], col.idx], test=X[ folds[[2]], col.idx],
            cl=y[ -folds[[2]] ], k=5)
sum(pred != y[ folds[[2]] ])


# Question 2.8.3
errors.knn <- sapply(seq_along(folds), function(i) {
    gene.t.tests <- colttests(X[ -folds[[i]], ], y[ -folds[[i]] ])
    col.idx <- order(gene.t.tests$p.value)[1:m]
    
    pred <- knn(train=X[ -folds[[i]], col.idx], test=X[ folds[[i]], col.idx],
                cl=y[ -folds[[i]] ], k=5)
    
    sum(pred != y[ folds[[i]] ])
})
sum(errors.knn)/length(y)
