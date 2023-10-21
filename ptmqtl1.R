library(purrr)
library(MASS)
set.seed(24601)
library(tidyverse)


patient_count <- 50
id_count <- 100

lambda_ <- rep(sample(2:8, size = id_count, replace = TRUE), patient_count)
p <- rep(rnorm(id_count, 0.4, 0.2), patient_count) # prob of SNP
p <- pmax(pmin(1, p), 0)
z <- rep(sample(c(-0.2, 0, 0.2), size = id_count, replace = TRUE, prob = c(0.15,0.7,0.15)), patient_count) # prob of getting PTM w/o SNP
q <- rep(rnorm(id_count, 0.3, 0.15), patient_count) # prob of getting PTM w/o SNP
q <- pmax(pmin(1, q), 0.2)


n = rpois(patient_count*id_count, lambda = lambda_)
n = replace(n, n==0, 1)
x = rbinom(patient_count*id_count, n, p)
y1 = rbinom(patient_count*id_count, x, q+z)
y2 = rbinom(patient_count*id_count, n-x, q)


# n = matrix(n, nrow = id_count, ncol=patient_count)
# x = matrix(x, nrow = id_count, ncol=patient_count)
# y1 = matrix(y1, nrow = id_count, ncol=patient_count)
# y2 = matrix(y2, nrow = id_count, ncol=patient_count)
# 
# n[1:5, 1:5]
# x[1:5, 1:5]
# y1[1:5, 1:5]
# y2[1:5, 1:5]
# 
# n[19]
# x[19]
# y1[19]
# y2[19]
# y[19]

# Convert the data to matrix form, this is the original form
dat = data.frame(patient=rep(1:patient_count,id_count), 
                 snp_ptm_id=rep(1:id_count,each=patient_count), peptide=n,
                 snp=x,non_snp=n-x,snp_ptm=y1,non_snp_ptm=y2,effect_size=z)
dat = as_tibble(dat)
print(dat, n = 100)

# Convert to a coarser form
dat_coarse = dat
dat_coarse$snp_type = dat_coarse$snp/dat_coarse$peptide
dat_coarse$snp_type[dat_coarse$snp_type<0.33] = 0
dat_coarse$snp_type[dat_coarse$snp_type>0.66] = 1
dat_coarse$snp_type[(dat_coarse$snp_type != 0) & (dat_coarse$snp_type != 1)] = 0.5
dat_coarse$ptm_count = dat_coarse$snp_ptm + dat_coarse$non_snp_ptm #not using?
dat_coarse = dat_coarse %>% select(patient, snp_ptm_id, peptide, snp_type, snp_ptm)
print(dat_coarse, n = 100)


coarse_lm <- lm(snp_ptm ~ snp_type, data = dat_coarse)
summary(coarse_lm)
plot(cooks.distance(coarse_lm), pch = 16, col = "blue")

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
