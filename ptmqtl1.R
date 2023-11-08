library(purrr)
library(MASS)
set.seed(24601)
library(tidyverse)


patient_count <- 50
id_count <- 100

lambda_ <- rep(sample(2:8, size = id_count, replace = TRUE), patient_count)

#probability of SNP
p <- rep(rnorm(id_count, 0.4, 0.2), patient_count)
p <- pmax(pmin(1, p), 0)

#effect size (q+z = prob of PTM with SNP)
z <- numeric(patient_count*id_count)
for (i in seq_len(patient_count)) {
  z[((i-1)*id_count + 1):(i*id_count)] <- rep(sample(c(-0.2, 0, 0.2),
    size = id_count, replace = TRUE, prob = c(0.15, 0.7, 0.15)), each = id_count)
}

#probability of PTM without SNP
q <- rep(rnorm(id_count, 0.3, 0.15), patient_count) # prob of getting PTM w/o SNP
q <- pmax(pmin(1, q), 0.2)

#peptide count
n = rpois(patient_count*id_count, lambda = lambda_)
n = replace(n, n==0, 1)

#peptide with SNP count
x = numeric(patient_count*id_count)
for (sampcnt in 1:(patient_count*id_count)) {
  x[sampcnt] = min(n[sampcnt], sample(c(0, 0.5, 1), size = 1, replace = TRUE,
    prob = c((1-p[sampcnt])^2, p[sampcnt]^2, max(1-((1-p[sampcnt])^2 + p[sampcnt]^2), 0))))
}
x[x == 1] <- n[x == 1]
x[x == 0.5] <- rbinom(length(x[x == 0.5]), n[x == 0.5], 0.5)

#peptide with PTM count
y1 = rbinom(patient_count*id_count, x, q+z) #with SNP
y2 = rbinom(patient_count*id_count, n-x, q) #without SNP

#x expanded as series of 1s (snp) and 0s (no snp) for new data table
snp_var <- numeric(0)
for (i in 1:(patient_count*id_count)) {
  snp_var <- c(snp_var, rep(1, x[i]), rep(0, n[i] - x[i]))
}

#expanded ptm count (y1 + y2) as 1s (ptm) and 0s (no ptm) for new data table
ptm_var <- numeric(0)
for (i in 1:(patient_count*id_count)) {
  ptm_var <- c(ptm_var, rep(1, y1[i] + y2[i]), rep(0, n[i] - (y1[i] + y2[i])))
}


# Convert the data to matrix form, this is the original form
dat = data.frame(patient=rep(1:patient_count,id_count),
                 snp_ptm_id=rep(1:id_count,each=patient_count), peptide=n,
                 snp=x,non_snp=n-x,snp_ptm=y1,non_snp_ptm=y2,effect_size=z)
dat = as_tibble(dat)
print(dat, n = 100)

#original form expanded out by number of peptides
new_dat = data.frame(patient=rep(rep(1:patient_count, id_count), n), 
                 snp_ptm_id=rep(rep(1:id_count,each=patient_count), n),
                 snp=snp_var, ptm=ptm_var,effect_size=rep(z, n))
new_dat = as_tibble(new_dat)
print(new_dat, n = 100)

# Convert to a coarser form
dat_coarse = dat
dat_coarse$snp_type = dat_coarse$snp/dat_coarse$peptide
dat_coarse$snp_type[dat_coarse$snp_type<0.33] = 0
dat_coarse$snp_type[dat_coarse$snp_type>0.66] = 1
dat_coarse$snp_type[(dat_coarse$snp_type != 0) & (dat_coarse$snp_type != 1)] = 0.5
dat_coarse$effect_size = dat_coarse$effect_size
dat_coarse = dat_coarse %>% select(patient, snp_ptm_id, peptide, snp_type, snp_ptm, effect_size)

dat_coarse$snp_ptm <- dat_coarse$snp_ptm / dat_coarse$peptide

print(dat_coarse, n = 100)

#coarse data linear regression plot
coarse_lm_pval <- numeric(id_count)
for (idc in 1:id_count) {
  coarse_lm <- lm(snp_ptm[(1+(patient_count*(idc-1))):(patient_count*idc)]
    ~ snp_type[(1+(patient_count*(idc-1))):(patient_count*idc)], data = dat_coarse)
  coarse_lm_pval[idc] <- summary(coarse_lm)$coefficients[2,4]
}

sig_col <- ifelse(coarse_lm_pval < 0.05, "red", "black")
plot(seq_len(id_count), coarse_lm_pval, main = "Significance of snp-ptm correlation (new)",
     xlab = "snp-ptm id", ylab = "p-value", pch = 20, col = sig_col)
abline(h = 0.05, col = "red4")


#old data linear regression plot
old_lm_pval <- numeric(id_count)
for (idc in 1:id_count) {
  old_lm <- lm(ptm[sum(n[0:(patient_count*(idc-1))]):sum(n[0:(patient_count*idc)])]
    ~ snp[sum(n[0:(patient_count*(idc-1))]):sum(n[0:(patient_count*idc)])], data = new_dat)
  old_lm_pval[idc] <- summary(old_lm)$coefficients[2,4]
}

sig_col_2 <- ifelse(old_lm_pval < 0.05, "blue", "black")
plot(seq_len(id_count), old_lm_pval, main = "Significance of snp-ptm correlation (old)",
     xlab = "snp-ptm id", ylab = "p-value", pch = 20, col = sig_col_2)
abline(h = 0.05, col = "blue4")


#coarse data analysis (id vs bvalue and snp-type vs normalized ptm count)
coarse_lm_b <- numeric(id_count)
for (idc in 1:id_count) {
  coarse_lm <- lm(snp_ptm[(1+(patient_count*(idc-1))):(patient_count*idc)]
                  ~ snp_type[(1+(patient_count*(idc-1))):(patient_count*idc)], data = dat_coarse)
  coarse_lm_b[idc] <- coef(coarse_lm)[2]
}
plot(dat_coarse$snp_type, dat_coarse$snp_ptm)
plot(1:100, coarse_lm_b, xlab = "snp-ptm_id", ylab = "B-value")
