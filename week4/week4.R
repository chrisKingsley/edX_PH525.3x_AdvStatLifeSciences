# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 4

library(genefilter)
library(limma)
library(R.utils)

library(tissuesGeneExpression)
data("tissuesGeneExpression")

library(dagdata)
data(hcmv)

library(rafalib)
mypar2()

# Question 4.1.1
dbinom(2, 4, 0.49)

# Question 4.1.2
dbinom(4, 10, 0.49)

sum(dbinom(11:20, 20, 0.4))
pbinom(10, 20, 0.4, lower.tail=F)

1-(1-1/175223510)^189000000
1-dbinom(0, 189000000, 1/175223510)

1-dbinom(0, 189000000, 1/175223510)-dbinom(1, 189000000, 1/175223510)
pbinom(1, 189000000, 1/175223510, lower.tail=F)


# Question 4.1.3
sum(dbinom(8:9, 20, 0.4))
pbinom(9, 20, 0.4)-pbinom(7, 20, 0.4)

sd <- sqrt(20*0.4*0.6)
pnorm(9, 8, sd=sd) - pnorm(7, 8, sd=sd)


# Question 4.1.4
binomProb <- pbinom(450, 1000, 0.4)-pbinom(350, 1000, 0.4)

sd <- sqrt(1000*0.4*0.6)
normProb <- pnorm(450, 400, sd=sd) - pnorm(350, 400, sd=sd)

abs(normProb - binomProb)


# Question 4.1.5
pdf("binomNormApprox.pdf")
par(mfrow=c(4, 5), mar=c(2,2,2,2))

for(N in c(5,10,30,100)) {
    for(p in c(0.01,0.10,0.5,0.9,0.99)) {
        exact <- dbinom(1:N, N, p)
        
        a <- (1:N + 0.5 - N*p)/sqrt(N*p*(1-p))
        b <- (1:N -0.5 - N*p)/sqrt(N*p*(1-p))
        approx = pnorm(a) - pnorm(b)
        
        plot(exact, col="blue", ylim=c(0, max(c(exact, approx))),
             main=sprintf("N:%s p:%s", N, p), pch=20)
        points(approx, col="red", pch=20)
        legend("topright", legend=c("binom","norm"), fill=c("blue","red"),
               cex=0.7, bty="n")
    }
}
dev.off()


# Question 4.1.6
N <- 189000000
p <- 1/175223510
dbinom(2,N,p)
dpois(2, N*p)

1 - dpois(0, N*p) - dpois(1, N*p)


# Question 4.2.1
breaks <- seq(0, 4000*round(max(locations)/4000), 4000)
tmp <- cut(locations, breaks)
counts <- as.numeric(table(tmp))
hist(counts)

logprobs <- dpois(counts, 4, log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

logLikPois <- function(lambda, counts) {
    logProbs <- dpois(counts, lambda, log=TRUE)
    sum(logProbs)
}
lambdas <- seq(0,15,len=300)
logLikelihoods <- sapply(lambdas, logLikPois, counts)
plot(lambdas, logLikelihoods)
idx <- which(logLikelihoods==max(logLikelihoods))
printf("maxLikelihood:%f   lambda:%f", logLikelihoods[idx], lambdas[idx])


# Question 4.2.2
binLocation <- (breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)
idx <- which(counts==max(counts))
printf("binLocation:%s counts:%s", binLocation[idx], counts[idx])


# Question 4.2.3
lambda <- mean(counts[ - which.max(counts) ])
ppois(13, lambda, lower.tail=F)


#  Question 4.2.4 
0.05/length(binLocation)


#  Question 4.2.5
ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean( counts[ -which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0, 1)


# Question 4.3.1
y <- e[,which(tissue=="endometrium")]
vars <- apply(y, 1, var, na.rm=T)
qqnorm(vars)
abline(0, 1)
qqnorm(sqrt(vars))
abline(0, 1)


# Question 4.3.2
f_fit <- fitFDist(vars, 14)
ps <- (seq(along=vars) - 0.5)/length(vars)
f_quantiles <- qf(ps, f_fit$scale, f_fit$df2)
qqplot(sqrt(vars), f_quantiles)
abline(0, 1)
qqplot(sqrt(vars[vars<quantile(vars, 0.95)]), f_quantiles, ylim=c(0,10))
abline(0, 1)
