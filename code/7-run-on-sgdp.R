library(bigsnpr)
library(dplyr)
library(ggplot2)

NCORES <- nb_cores()

prefix <- "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc"
bedfile <- runonce::download_file(paste0(prefix, ".bed"), dir = "tmp-data")
runonce::download_file(paste0(prefix, ".fam"), dir = "tmp-data")
bim_zip <- runonce::download_file(paste0(prefix, ".bim.zip"), dir = "tmp-data")
unzip(bim_zip, exdir = "tmp-data", overwrite = FALSE)

obj.bed <- bed(bedfile)
dim(obj.bed)  # 345 x 34,418,131

nbna_var <- bed_counts(obj.bed, ncores = NCORES)[4, ]
mean(keep <- (nbna_var == 0))  # 88.6%

all_freq <- bigreadr::fread2("ref_freqs.csv.gz")
projection <- bigreadr::fread2("projection.csv.gz")
correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

matched <- obj.bed$map %>%
  transmute(chr = chromosome, pos = physical.pos, a0 = allele2, a1 = allele1,
            beta = 1, keep) %>%
  snp_match(all_freq[1:5]) %>%
  filter(keep) %>%
  select(-keep) %>%
  print()
# 3,545,628 variants left

# project new individuals onto the PCA space
all_proj <- sapply(1:16, function(k) {
  bed_prodVec(
    obj.bed,
    projection[matched$`_NUM_ID_`, k + 5] / 2 * correction[k],
    ind.col = matched$`_NUM_ID_.ss`,
    # scaling to get G if beta = 1 and (2 - G) if beta = -1
    center = 1 - matched$beta,
    scale = matched$beta,
    ncores = NCORES
  )
})

# projected reference allele frequencies
X <- crossprod(as.matrix(projection[matched$`_NUM_ID_`, -(1:5)]),
               as.matrix(all_freq[matched$`_NUM_ID_`, -(1:5)]))
cp_X_pd <- Matrix::nearPD(crossprod(X), base.matrix = TRUE)
Amat <- cbind(1, diag(ncol(X)))
bvec <- c(1, rep(0, ncol(X)))

# solve a QP for each projected individual
all_res <- round(apply(all_proj, 1, function(y) {
  quadprog::solve.QP(
    Dmat = cp_X_pd$mat,
    dvec = crossprod(y, X),
    Amat = Amat,
    bvec = bvec,
    meq  = 1
  )$sol
}), 7)

rownames(all_res) <- colnames(X)

# get (population) information on individuals
fam2 <- bigreadr::fread2("https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt")
fam <- left_join(obj.bed$fam, fam2, by = c("sample.ID" = "Sample_ID(Aliases)"))

round(100 * all_res[, 1:5], 1)
info <- fam[c("Region", "Country", "Population_ID")]
Q <- t(all_res)
res <- cbind.data.frame(info, round(100 * Q, 1))


#### Also try ADMIXTURE ####

# prepare necessary files (bed and frequencies)
snp_readBed2(obj.bed$bedfile, ind.col = matched$`_NUM_ID_.ss`, ncores = NCORES)
obj.bigsnp <- snp_attach("tmp-data/cteam_extended.v4.maf0.1perc.rds")
G <- obj.bigsnp$genotypes
need_reverse <- matched$beta < 0
G[, need_reverse] <- 2 - G[, need_reverse]
sum(big_counts(G)[4, ])  # 0 missing values

ind_keep <- snp_clumping(G, infos.chr = obj.bigsnp$map$chromosome, thr.r2 = 0.2,
                         infos.pos = obj.bigsnp$map$physical.pos, ncores = NCORES)
# 258,925 variants kept

freq_kept <- all_freq[matched$`_NUM_ID_`[ind_keep], -(1:5)]

# verification
af_china <- big_scale()(G, ind.row = which(info$Country == "China"), ind.col = ind_keep)$center / 2
plot(af_china, freq_kept$`Asia (East)`)  # ok

snp_writeBed(obj.bigsnp, "tmp-data/sgdp_for_admixture.bed", ind.col = ind_keep)

data.table::fwrite(1 - freq_kept,  # needed to revert freqs (probably not using the same ref allele)
                   file = "tmp-data/sgdp_for_admixture.21.P.in",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

tgz <- runonce::download_file(
  "http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz",
  dir = "tmp-data")
untar(tgz, exdir = "tmp-data")

setwd("tmp-data")
admixture <- "dist/admixture_linux-1.3.0/admixture"
system(glue::glue(
  "{admixture} -P sgdp_for_admixture.bed 21 -j{nb_cores()}"))
# Converged in 4 iterations (179.67 sec)

Q2 <- read.table("sgdp_for_admixture.21.Q")
colnames(Q2) <- names(all_freq)[-(1:5)]
res2 <- cbind.data.frame(info, round(100 * Q2, 1))

ind <- as.vector(matrix(seq_len(2 * nrow(res)), nrow = 2, byrow = TRUE))
final_res <- rbind(
  cbind(Method = "snp_ancestry_summary", id = rows_along(res), res),
  cbind(Method = "admixture", id = rows_along(res2), res2)
)[ind, ]

bigreadr::fwrite2(final_res, "../ancestry_sgdp.csv")
