coef <- 0.5
p1 <- 0.2; N1 <- 200
p2 <- 0.9; N2 <- 500
hist(x1 <- rbeta(1e6, shape1.1 <- 2 * N1 * p1, shape2.1 <- 2 * N1 * (1 - p1)))
hist(x2 <- rbeta(1e6, shape1.2 <- 2 * N2 * p2, shape2.2 <- 2 * N2 * (1 - p2)))
hist(x <- coef * x1 + (1 - coef) * x2, "FD", freq = FALSE)
(param <- EnvStats::ebeta(x)$param)
curve(dbeta(x, param[1], param[2]), add = TRUE, col = "blue", lwd = 2)

a <- coef * shape1.1 / (shape1.1 + shape2.1) + (1 - coef) * shape1.2 / (shape1.2 + shape2.2)
c <- (1 - a) / a
b <- coef^2 * (shape1.1 * shape2.1) / (shape1.1 + shape2.1)^2 / (shape1.1 + shape2.1 + 1) +
  (1 - coef)^2 * (shape1.2 * shape2.2) / (shape1.2 + shape2.2)^2 / (shape1.2 + shape2.2 + 1)

param1 <- (c / b / (1 + c)^2 - 1) / (1 + c)
param2 <- c * param1
curve(dbeta(x, param1, param2), add = TRUE, col = "green", lwd = 2)


hist(rnorm(1e5))
hist(0.5 * rnorm(1e5) + 0.5 * rnorm(1e5))

###

library(dplyr)

all_freq <- bigreadr::fread2(
  runonce::download_file("https://figshare.com/ndownloader/files/31218445",
                         dir = "../bigsnpr/tmp-data", fname = "ref_freqs.csv.gz")
)

all_N <- bigreadr::fread2("all_prop.csv")$N
all_N[c(8, 11, 12, 15)] <- all_N[c(8, 11, 12, 15)] + c(99, 84, 104, 86)


# T2D from FinnGen
(N <- 29193 + 182573)
sumstats <- bigreadr::fread2("../freq-ancestry//tmp-data/summary_stats_finngen_R5_T2D.gz") %>%
  transmute(rsid = rsids, chr = `#chrom`, pos, a0 = ref, a1 = alt,
            freq = (2 * (n_hom_cases + n_hom_controls) +
                      (n_het_cases + n_het_controls)) / (2 * N)) %>%
  filter(chr %in% 1:22)


# Epilepsy
gz <- runonce::download_file(
  "http://www.epigad.org/gwas_ilae2018_16loci/all_epilepsy_METAL.gz",
  dir = "../bigsnpr/tmp-data")
readLines(gz, n = 3)
sumstats <- bigreadr::fread2(
  gz, select = c("CHR", "BP", "Allele2", "Allele1", "Freq1"),
  col.names = c("chr", "pos", "a0", "a1", "freq")
) %>%
  mutate_at(3:4, toupper)

# PAGE
gz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008053/WojcikG_PMID_invn_rheight_alls.gz",
  dir = "../bigsnpr/tmp-data")
sumstats <- bigreadr::fread2(
  gz, select = c("Chr", "Position_hg19", "Other-allele", "Effect-allele", "Effect-allele-frequency"),
  col.names = c("chr", "pos", "a0", "a1", "freq"))

pop <- c("African-American" = 17286, "Hispanic-Latino" = 22192, "Asian" = 4680,
         "Native Hawaiian" = 3939, "Native American" = 647,	Other	= 1052)
round(sort(pop / sum(pop), decreasing = TRUE), 3)
# Hispanic-Latino African-American    Asian  Native Hawaiian    Other  Native American
#           0.446            0.347    0.094            0.079    0.021            0.013


## CAD from Nikpay et al, NG 2015

txt <- runonce::download_file(
  "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
  dir = "../bigsnpr/tmp-data/", fname = "sumstats_CAD.txt")
sumstats <- bigreadr::fread2(
  txt, select = c("chr", "bp_hg19", "noneffect_allele", "effect_allele",
                  "beta", "se_dgc", "p_dgc", "effect_allele_freq",
                  "median_info", "n_studies", "het_pvalue"),
  col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se",
                "p", "freq", "info", "n_studies", "het_pvalue"))


## UKBB LD ref
sumstats <- mutate(readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  dir = "../bigsnpr/tmp-data", fname = "map_hm3_ldpred2.rds")),
  freq = af_UKBB)


## T2D in Africa

txt <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008114/ChenJ_31049640",
  dir = "tmp-data")
readLines(txt, n = 10)
sumstats <- bigreadr::fread2(txt, select = c("MarkerName", "Allele1", "Allele2", "Freq1"),
                   col.names = c("snpid", "a1", "a0", "freq")) %>%
  mutate(chr = as.integer(sub("^(.+):.+:.+_.+$", "\\1", snpid)),
         pos = as.integer(sub("^.+:(.+):.+_.+$", "\\1", snpid)),
         snpid = NULL, a1 = toupper(a1), a0 = toupper(a0))


matched <- bigsnpr::snp_match(
  mutate(sumstats, chr = as.integer(chr), beta = 1),
  all_freq[1:5], return_flip_and_rev = TRUE, #join_by_pos = FALSE
) %>%
  mutate(freq = ifelse(`_REV_`, 1 - freq, freq))

X <- as.matrix(all_freq[matched$`_NUM_ID_`, -(1:5)])
y <- matched$freq


SHAPE1 <- 1 + sweep(X, 2, 2 * N, '*')
SHAPE2 <- 1 + sweep(1 - X, 2, 2 * N, '*')
MEAN <- (SHAPE1 / (SHAPE1 + SHAPE2))
VAR <- ((SHAPE1 * SHAPE2) / (SHAPE1 + SHAPE2)^2 / (SHAPE1 + SHAPE2 + 1))


logLik <- function(coef) {
  A <- MEAN %*% coef
  C <- (1 - A) / A
  B <- VAR %*% (coef^2)
  PARAM1 <- (C / B / (1 + C)^2 - 1) / (1 + C)
  PARAM2 <- C * PARAM1
  mean(dbeta(y, PARAM1, PARAM2, log = TRUE), na.rm = TRUE)
}

(coef <- bigsnpr::snp_ancestry_summary(y, X)$coef)
logLik(coef)
coef2 <- coef * (coef > 0.2); coef2 <- coef2 / sum(coef2)
logLik(coef2)
# -632.1561 -> -519.7371
# -195.968 -> -310.2099


D <- crossprod(X)
eig_min <- min(eigen(D, symmetric = TRUE, only.values = TRUE)$values)

get_res <- function(coef) {
  diag(D) <- diag(D) - coef * eig_min
  res <- quadprog::solve.QP(
    Dmat = D,
    dvec = crossprod(y, X),
    Amat = cbind(1, diag(ncol(X))),
    bvec = c(1, rep(0, ncol(X))),
    meq  = 1
  )$solution
  round(res, 7)
}

coef <- get_res(0)
round(100 * setNames(coef, colnames(X)), 1)

SEQ <- seq(0, 0.999, length.out = 20)
all_ll <- sapply(SEQ, function(coeff) {
  res <- get_res(coeff)
  print(round(100 * setNames(res, colnames(X)), 1))
  logLik(res)
})
plot(SEQ, print(all_ll))

y_hat0 <- X %*% get_res(0)
mean(abs(y - y_hat0))
coef <- get_res(SEQ[which.max(all_ll)])
round(100 * setNames(coef, colnames(X)), 1)
y_hat <- X %*% coef
mean(abs(y - y_hat))
