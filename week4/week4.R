# code for HarvardX: PH525.3x Advanced Statistics for the Life Sciences - week 4

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


