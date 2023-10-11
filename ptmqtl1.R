options(max.print=1000000)
library(purrr)
library(MASS)
set.seed(24601)


p <- 0.5 # prob of SNP
q <- 0.3 # prob of getting PTM w/o SNP
z <- 0.2 # effect of SNP on PTM
z <- min(z, 1-q)
z <- max(z, -q) # ensure 0<q+z<1
lambda_ <- 5
patient_count <- 5
id_count <- 4


n <- matrix(0, nrow = id_count, ncol = patient_count)
x <- matrix(0, nrow = id_count, ncol = patient_count)
y1 <- matrix(0, nrow = id_count, ncol = patient_count)
y2 <- matrix(0, nrow = id_count, ncol = patient_count)
y <- matrix(0, nrow = id_count, ncol = patient_count)

for (id_it in 1:id_count) {
  n[id_it,] <- rpois(patient_count, lambda = lambda_) # getting count of peptide
  x[id_it,] <- rbinom(patient_count, n[id_it,], p) # getting SNP count
  y1[id_it,] <- rbinom(patient_count, x[id_it,], q+z) # PTM w/ SNP
  y2[id_it,] <- rbinom(patient_count, n[id_it,]-x[id_it,], q) # PTM w/o SNP
}
y <- y1 + y2

n
x
y1
y2
y


dat <- matrix(0, nrow = sum(n), ncol = 6)
colnames(dat) <- c("Patient", "SNP-PTM id", "Peptide", "SNP Exist?", "PTM Exist?", "Effect Size")

pos <- 1 #starting position iterator for each patient
id_pos <- 1 #starting SNP-PTM id position

for (tr in 1:id_count) {
  print(sum(n[tr,]))
}
#db <- function() {
for (id_iter in 1:id_count) {
  dat[id_pos:(id_pos + sum(n[id_iter,]) - 1), 2] <- id_iter
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
sum(x)
sum(dat[,4])
dat #ascending order by SNP-PTM id
dat[order(dat[,1],decreasing=FALSE),] #ascending order by patient number
#}
#debug(db)
#db()

#group by patient and SNP-PTM pairs, sum up peptide, sum up PTM exist, determine if 
# SNP exist is 0(all are 0) or 1(half 0 half 1) or 2(all 1). 

pdat <- matrix(0, nrow = patient_count*id_count*2, ncol = 4)
colnames(pdat) <- c("w/ SNP", "w/o SNP", "SNP_PTM id", "Patient #")
rownames(pdat) <- rep(c("w/ PTM", "w/o PTM"), patient_count * id_count)

row_iter <- 1
for (ic in 1:id_count) {
  for (pc in 1:patient_count) {
    pdat[row_iter, 1] <- y1[ic, pc]
    pdat[row_iter + 1, 1] <- x[ic, pc] - y1[ic, pc]
    pdat[row_iter, 2] <- y2[ic, pc]
    pdat[row_iter + 1, 2] <- n[ic, pc] - x[ic, pc] - y2[ic, pc]
    pdat[row_iter:(row_iter + 1), 3] <- ic
    pdat[row_iter:(row_iter + 1), 4] <- pc
    row_iter <- row_iter + 2
  }
}
pdat


chidat <- matrix(0, nrow = id_count*2, ncol = 3)
colnames(chidat) <- c("w/ SNP", "w/o SNP", "SNP_PTM id")
rownames(chidat) <- rep(c("w/ PTM", "w/o PTM"), id_count)

row_iter <- 1
for (ic in 1:id_count) {
  chidat[row_iter, 1] <- sum(y1[ic,])
  chidat[row_iter + 1, 1] <- sum(x[ic,]) - sum(y1[ic,])
  chidat[row_iter, 2] <- sum(y2[ic,])
  chidat[row_iter + 1, 2] <- sum(n[ic,]) - sum(x[ic,]) - sum(y2[ic,])
  chidat[row_iter:(row_iter + 1), 3] <- ic
  row_iter <- row_iter + 2
}
chidat


# chisq.test(pdat)$p.value # do this for every SNP-PTM pair, then graph would be nice
# barplot(H,xlab,ylab,main, names.arg,col)


x_new <- dat[, 4]
y_new <- dat[, 5]
b <- dat[, 6]

neg_binom <- glm.nb(y_new ~ x_new * b, data = as.data.frame(dat)) #we don't know b
summary(neg_binom)
