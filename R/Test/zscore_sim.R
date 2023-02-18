num.subj <- 1000
pheno <- rnorm(num.subj)

beta.vec <- rep(0,num.subj)
std.err.vec <- rep(0,num.subj)
zscore.vec <- rep(0,num.subj)
af1.vec <- rep(0,num.subj)

for(i in 1:(num.subj-1)){
  pheno <- rnorm(num.subj)
  geno <- c(rep(2, i), rep(0, num.subj-i))
  mod <- lm(pheno~geno)
  out <- summary(mod)
  beta.vec[i] <- out$coefficients[2,1]
  std.err.vec[i] <- out$coefficients[2,2]
  zscore.vec[i] <- out$coefficients[2,3]
  af1.vec[i] <- sum(geno)/(2*length(geno))
}
plot(af1.vec, std.err.vec)
plot(af1.vec, zscore.vec)
plot(af1.vec, beta.vec)

mod <- lm(pheno~geno)
summary(mod)

mean(geno)
sd(geno)

norm.geno <- (geno-mean(geno))/sqrt(var(geno))
head(norm.geno)

p <- sum(geno)/(2*length(geno)) # allele frequency
p

2*p # mean
sqrt(2*p*(1-p)) # variance

