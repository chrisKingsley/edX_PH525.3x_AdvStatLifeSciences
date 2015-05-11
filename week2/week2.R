# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 1

library(svd)

library(devtools)
install_github("genomicsclass/tissuesGeneExpression")

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

