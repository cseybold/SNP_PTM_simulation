options(max.print=1000000)
library(purrr)
library(MASS)
set.seed(24601)
library(tidyverse)


p <- 0.5 # prob of SNP
q <- 0.3 # prob of getting PTM w/o SNP
z <- 0.2 # effect of SNP on PTM
z <- min(z, 1-q)
z <- max(z, -q) # ensure 0<q+z<1
lambda_ <- 5 
patient_count <- 20
id_count <- 100


# n <- matrix(0, nrow = id_count, ncol = patient_count)
# x <- matrix(0, nrow = id_count, ncol = patient_count)
# y1 <- matrix(0, nrow = id_count, ncol = patient_count)
# y2 <- matrix(0, nrow = id_count, ncol = patient_count)
# y <- matrix(0, nrow = id_count, ncol = patient_count)
# 
# for (id_it in 1:id_count) {
#   n[id_it,] <- rpois(patient_count, lambda = lambda_) # getting count of peptide
#   x[id_it,] <- rbinom(patient_count, n[id_it,], p) # getting SNP count
#   y1[id_it,] <- rbinom(patient_count, x[id_it,], q+z) # PTM w/ SNP
#   y2[id_it,] <- rbinom(patient_count, n[id_it,]-x[id_it,], q) # PTM w/o SNP
# }
# y <- y1 + y2

n = rpois(patient_count*id_count, lambda = lambda_)
x = rbinom(patient_count*id_count, n, p)
y1 = rbinom(patient_count*id_count, x, q+z)
y2 = rbinom(patient_count*id_count, n-x, q)

# n = matrix(n, nrow = id_count, ncol=patient_count)
# x = matrix(x, nrow = id_count, ncol=patient_count)
# y1 = matrix(y1, nrow = id_count, ncol=patient_count)
# y2 = matrix(y2, nrow = id_count, ncol=patient_count)

# n[1:5, 1:5]
# x[1:5, 1:5]
# y1[1:5, 1:5]
# y2[1:5, 1:5]
# 
# n
# x
# y1
# y2
# y

# Convert the data to matrix form, this is the original form
dat = data.frame(patient=rep(1:patient_count,id_count), 
                 snp_ptm_id=rep(1:id_count,each=patient_count), peptide=n,
                 snp=x,non_snp=n-x,snp_ptm=y1,non_snp_ptm=y2,effect_size=z)
dat = as.tibble(dat)

# Convert to a coarser form
dat_coarse = dat
dat_coarse$snp_type = dat_coarse$snp/dat_coarse$peptide
dat_coarse$snp_type[dat_coarse$snp_type<0.33] = 0
dat_coarse$snp_type[dat_coarse$snp_type>0.66] = 1
dat_coarse$snp_type[dat_coarse$snp_type!=0 && dat_coarse$snp_type!=0] = 0.5
dat_coarse$ptm_count = dat_coarse$snp_ptm + dat_coarse$non_snp_ptm
dat_coarse = dat_coarse %>% select(patient, snp_ptm_id, peptide, snp_type,
                                   snp_ptm)
##################################################################
# Task 1, compare the performance using dat and dat_coarse
# Task 2, make the parameters sampled from distribution instead of using a fixed number
# for probability of getting a snp, we should make it a bit more realistic
# We first sample p from a defined distribution (prefer to have most of them as 
# small values). then for each patient, we do a binomial with sample size of 2
# mimicing the behavior of two chromosome copies. Note, p is associated with
# snp_ptm_id, so for different patient with same snp_ptm_id, p should be the same.
# Similarly, q should be associated with snp_ptm_id as well. no preference on distribution
# similarly, z should be associated with snp_ptm_id. make sure 0<=p+z<=1
# lambda_ associated with snp_ptm_id.
# compare the performance under this setting.





pos <- 1 #starting position iterator for each patient
id_pos <- 1 #starting SNP-PTM id position

# for (tr in 1:id_count) {
#   print(sum(n[tr,]))
# }


#db <- function() {
for (id_iter in 1:id_count) {
  dat[id_pos:(id_pos + sum(n[id_iter,]) - 1), 2] <- id_iter
  dat$snp_ptm_id[id_pos:(id_pos + sum(n[id_iter,]) - 1)] <- id_iter
  dat[id_pos:(id_pos + sum(n[id_iter,]) - 1), 6] <- sample(c(-0.3, 0, 0.3), size=1, prob=c(0.1, 0.8, 0.1))
  id_pos <- id_pos + sum(n[id_iter,])
  for (i in 1:patient_count) {
    dat[(pos + sum(n[0:(id_iter-1),])):((n[id_iter, i] + pos - 1) + sum(n[0:(id_iter-1),])), 1] <- i #patient number
    if (x[id_iter, i] > 0) {
      dat[(pos + sum(n[0:(id_iter-1),])):((x[id_iter, i] + pos - 1) + sum(n[0:(id_iter-1),])), 4] <- 1 #SNP exist
    }
    if (y1[id_iter, i] > 0) {
      dat[(pos + sum(n[0:(id_iter-1),])):((y1[id_iter, i] + pos - 1) + sum(n[0:(id_iter-1),])), 5] <- 1 #PTM exist? (PTM w/ SNP)
    }
    if (y2[id_iter, i] > 0) {
      dat[((pos + x[id_iter, i]) + sum(n[0:(id_iter-1),])):((y2[id_iter, i] + x[id_iter, i] + pos - 1) + sum(n[0:(id_iter-1),])), 5] <- 1 #PTM exist? (PTM w/o SNP)
    }
    pos2 <- pos #starting position iterator for each patient's peptide number
    for (j in 1:n[id_iter, i]) {
      dat[pos2 + sum(n[0:(id_iter-1),]), 3] <- j #peptide number
      pos2 <- pos2 + 1
    }
    pos <- n[id_iter, i] + pos #next starting position
  }
  pos <- 1
  id_iter <- id_iter + 1
}
dat #ascending order by SNP-PTM id
#dat[order(dat[,1],decreasing=FALSE),] #ascending order by patient number
#}
#debug(db)
#db()

#group by patient and SNP-PTM pairs, sum up peptide, sum up PTM exist, determine if 
# SNP exist is 0(all are 0) or 1(half 0 half 1) or 2(all 1). 

pdat <- matrix(0, nrow = patient_count*id_count*2, ncol = 4)
colnames(pdat) <- c("w/ SNP", "w/o SNP", "SNP_PTM id", "Patient #")
rownames(pdat) <- rep(c("w/ PTM", "w/o PTM"), patient_count * id_count)

snppat_label <- list()
label_iter <- 1
row_iter <- 1
for (ic in 1:id_count) {
  for (pc in 1:patient_count) {
    pdat[row_iter, 1] <- y1[ic, pc]
    pdat[row_iter + 1, 1] <- x[ic, pc] - y1[ic, pc]
    pdat[row_iter, 2] <- y2[ic, pc]
    pdat[row_iter + 1, 2] <- n[ic, pc] - x[ic, pc] - y2[ic, pc]
    pdat[row_iter:(row_iter + 1), 3] <- ic
    pdat[row_iter:(row_iter + 1), 4] <- pc
    snppat_label[[label_iter]] <- c(ic, pc)
    row_iter <- row_iter + 2
    label_iter <- label_iter + 1
  }
}
pdat


# snpchidat <- matrix(0, nrow = id_count*2, ncol = 3)
# colnames(snpchidat) <- c("w/ SNP", "w/o SNP", "SNP_PTM id")
# rownames(snpchidat) <- rep(c("w/ PTM", "w/o PTM"), id_count)
# 
# row_iter <- 1
# for (ic in 1:id_count) {
#   snpchidat[row_iter, 1] <- sum(y1[ic,])
#   snpchidat[row_iter + 1, 1] <- sum(x[ic,]) - sum(y1[ic,])
#   snpchidat[row_iter, 2] <- sum(y2[ic,])
#   snpchidat[row_iter + 1, 2] <- sum(n[ic,]) - sum(x[ic,]) - sum(y2[ic,])
#   snpchidat[row_iter:(row_iter + 1), 3] <- ic
#   row_iter <- row_iter + 2
# }
# snpchidat
# 
# 
# patchidat <- matrix(0, nrow = patient_count*2, ncol = 3)
# colnames(patchidat) <- c("w/ SNP", "w/o SNP", "Patient #")
# rownames(patchidat) <- rep(c("w/ PTM", "w/o PTM"), patient_count)
# 
# row_iter <- 1
# for (pc in 1:patient_count) {
#   patchidat[row_iter, 1] <- sum(y1[, pc])
#   patchidat[row_iter + 1, 1] <- sum(x[, pc]) - sum(y1[, pc])
#   patchidat[row_iter, 2] <- sum(y2[, pc])
#   patchidat[row_iter + 1, 2] <- sum(n[, pc]) - sum(x[, pc]) - sum(y2[, pc])
#   patchidat[row_iter:(row_iter + 1), 3] <- pc
#   row_iter <- row_iter + 2
# }
# patchidat


chipvals <- matrix(0, nrow = patient_count * id_count, ncol = 3)
colnames(chipvals) <- c("p-value", "SNP_PTM id", "Patient #")

row_it <- 1
for (rownum in seq(1, 2 * patient_count * id_count, 2)) {
  if (any(pdat[rownum:(rownum + 1), 1:2]) == 0) {
    chipvals[row_it, 1] <- 1 #best way to handle?
  }
  else {
    chipvals[row_it, 1] <- chisq.test(pdat[rownum:(rownum + 1), 1:2])$p.value
  }
  chipvals[row_it, 2] <- (row_it - 1) %/% patient_count + 1
  chipvals[row_it, 3] <- rep(seq(1, patient_count), id_count)[row_it]
  row_it <- row_it + 1
}

chipvals[is.nan(chipvals)] <- 1
chipvals
barplot(chipvals[, 1], xlab = "(SNP-PTM id, Patient #)", ylab = "p-value", main = "Significance of SNP-PTM associations", names.arg = snppat_label)


# x_new <- dat[, 4]
# y_new <- dat[, 5]
# b <- dat[, 6]
# 
# neg_binom <- glm.nb(y_new ~ x_new * b, data = as.data.frame(dat)) #we don't know b
# summary(neg_binom)
