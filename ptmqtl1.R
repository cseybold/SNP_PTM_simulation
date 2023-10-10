#make github repo
library(purrr)
library(MASS)
set.seed(24601)


p <- 0.5 # prob of SNP
q <- 0.3 # prob of getting PTM w/o SNP
z <- 0.2 # effect of SNP on PTM
z <- min(z, 1-q)
z <- max(z, -q) # ensure 0<q+z<1
lambda_ <- 5
patient_count <- 200
id_count <- 1000

n <- rpois(patient_count, lambda = lambda_) # getting count of peptide
n # make random for each patient/ptmSNP pair

x <- matrix(0, nrow = id_count, ncol = patient_count)
y1 <- matrix(0, nrow = id_count, ncol = patient_count)
y2 <- matrix(0, nrow = id_count, ncol = patient_count)
y <- matrix(0, nrow = id_count, ncol = patient_count)

for (id_it in 1:id_count) {
  x[id_it,] <- rbinom(patient_count, n, p) # getting SNP count
  y1[id_it,] <- rbinom(patient_count, x[id_it,], q+z) # PTM w/ SNP
  y2[id_it,] <- rbinom(patient_count, n-x[id_it,], q) # PTM w/o SNP
}
y <- y1 + y2

x
y1
y2
y

dat <- matrix(0, nrow = sum(n) * id_count, ncol = 6)
colnames(dat) <- c("Patient", "SNP-PTM id", "Peptide", "SNP Exist?", "PTM Exist?", "Effect Size")

pos <- 1 #starting position iterator for each patient
id_pos <- 1 #starting SNP-PTM id position

for (id_iter in 1:id_count) {
  dat[id_pos:(id_pos + sum(n) - 1), 2] <- id_iter
  dat[id_pos:(id_pos + sum(n) - 1), 6] <- sample(c(-0.3, 0, 0.3), size=1, prob=c(0.1, 0.8, 0.1))
  id_pos <- id_pos + sum(n)
  for (i in 1:patient_count) {
    dat[(pos + (id_iter-1) * sum(n)):((n[i] + pos - 1) + (id_iter-1) * sum(n)), 1] <- i #patient number
    dat[(pos + (id_iter-1) * sum(n)):((x[id_iter, i] + pos - 1) + (id_iter-1) * sum(n)), 4] <- 1 #SNP exist (associate randomness in assignment order?)
    if (y1[id_iter, i] > 0) {
      dat[(pos + (id_iter-1) * sum(n)):((y1[id_iter, i] + pos - 1) + (id_iter-1) * sum(n)), 5] <- 1 #PTM exist? (PTM w/ SNP)
    }
    if (y2[id_iter, i] > 0) {
      dat[((pos + x[id_iter, i]) + (id_iter-1) * sum(n)):((y2[id_iter, i] + x[id_iter, i] + pos - 1) + (id_iter-1) * sum(n)), 5] <- 1 #PTM exist? (PTM w/o SNP)
    }
    pos2 <- pos #starting position iterator for each patient's peptide number
    for (j in 1:n[i]) {
      dat[pos2 + (id_iter-1) * sum(n), 3] <- j #peptide number
      pos2 <- pos2 + 1
    }
    pos <- n[i] + pos #next starting position
  }
  pos <- 1
  id_iter <- id_iter + 1
}
dat

#group by patient and SNP-PTM pairs, sum up peptide, sum up PTM exist, determine if 
# SNP exist is 0(all are 0) or 1(half 0 half 1) or 2(all 1). 
ptm_snp <- sum(y1)
nptm_snp <- sum(x) - sum(y1)
ptm_nsnp <- sum(y2)
nptm_nsnp <- id_count * sum(n) - sum(x) - sum(y2)

pdat <- matrix(c(ptm_snp, nptm_snp, ptm_nsnp, nptm_nsnp), nrow = 2, ncol = 2)
colnames(pdat) <- c("w/ SNP", "w/o SNP")
rownames(pdat) <- c("w/ PTM", "w/o PTM")
pdat


x_new <- dat[, 4]
y_new <- dat[, 5]
b <- dat[, 6]

neg_binom <- glm.nb(y_new ~ x_new * b, data = as.data.frame(dat)) #we don't know b
summary(neg_binom)


chi_test <- chisq.test(pdat) # do this for every SNP-PTM pair, then graph would be nice
chi_test
