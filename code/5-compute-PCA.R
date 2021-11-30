all_eid <- readRDS("tmp-data/info-UKBB.rds")$eid
eid_sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_1[-1]
all_ind <- readRDS("data/list_ind_pop.rds")
all_ind2 <- match(all_eid[unlist(all_ind)], eid_sample)

all_freq <- readRDS("data/all_freq.rds")
all_snp_id <- with(all_freq, split(paste(chr, pos, a0, a1, sep = "_"), as.factor(1:22)[chr]))

library(bigsnpr)
snp_readBGEN(
  bgenfiles   = paste0("UKBB/bgen/ukb_imp_chr", 1:22, "_v3.bgen"),
  list_snp_id = all_snp_id,
  backingfile = "tmp-data/dosage_all_ind",
  ind_row     = all_ind2,
  ncores      = nb_cores()
)

ukbb <- snp_attach("tmp-data/dosage_all_ind.rds")
G <- ukbb$genotypes
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos
ind_chr <- with(ukbb$map, split(seq_along(chromosome), chromosome))
res_files <- paste0("tmp-data/ld_chr", 1:22, ".rds")

if (all(file.exists(res_files))) {
  future::plan("sequential")
} else {
  library(future.batchtools)
  NCORES <- 15
  plan(batchtools_slurm(workers = print(nrow(grid)), resources = list(
    t = "12:00:00", c = NCORES, mem = "120g",
    name = basename(rstudioapi::getSourceEditorContext()$path))))
}

ld <- unlist(furrr::future_map(1:22, function(chr) {
  runonce::save_run({
    ind <- ind_chr[[chr]]
    bigsnpr::snp_ld_scores(G, ind.col = ind, size = 1000, infos.pos = POS[ind],
                           ncores = NCORES)
  }, file = res_files[chr])
})) # 30 min for chr1 / 6 min for chr22


## 1000G

obj.bed <- bed("tmp-data/all_phase3.bed")
fam <- dplyr::left_join(
  obj.bed$fam,
  bigreadr::fread2("https://figshare.com/ndownloader/files/31080292"))
ind_row <- which(fam$Population %in% c("FIN", "BEB", "JPT", "PEL"))
map <- dplyr::select(obj.bed$map, chr = chromosome, pos = physical.pos,
                     a0 = allele2, a1 = allele1)
ind_col <- vctrs::vec_match(all_freq[1:4], map)


## PCA for 26 groups in UKBB and 4 pops in 1000G
NCORES <- nb_cores()
sqrt_ld <- sqrt(ld)

A <- function(x, args) {
  cat("_")
  c(
    big_prodVec(G, x, scale = sqrt_ld, ncores = NCORES),
    bed_prodVec(obj.bed, x, ind.row = ind_row, ind.col = ind_col,
                scale = sqrt_ld, ncores = NCORES)
  )
}
Atrans <- function(x, args) {
  cat("|")
  ind_G <- rows_along(G)
  big_cprodVec(G, x[ind_G], scale = sqrt_ld, ncores = NCORES) +
    bed_cprodVec(obj.bed, x[-ind_G], ind.row = ind_row, ind.col = ind_col,
                 scale = sqrt_ld, ncores = NCORES)
}

obj.svd <- runonce::save_run(
  structure(RSpectra::svds(A = A, Atrans = Atrans, k = 20,
                           dim = c(nrow(G) + length(ind_row), ncol(G))), class = "big_SVD"),
  file = "tmp-data/big_svd_all2.rds"
)

library(ggplot2)
plot(obj.svd)
# plot(obj.svd, type = "scores", scores = 1:20, coef = 0.4)
plot(obj.svd, type = "loadings", loadings = 1:20, coef = 0.3, ncol = 4)
# ggsave("figures/all_loadings.png", width = 8, height = 9)

all_freq <- readRDS("data/all_freq.rds")
project <- 2 * sweep(obj.svd$v, 1, sqrt_ld, '/')
colnames(project) <- paste0("PC", 1:20)
pred <- crossprod(as.matrix(all_freq[-(1:5)]), project)

pop <- c(rep(names(all_ind), lengths(all_ind)),
         unname(c("FIN" = "Finland", "BEB" = "Bangladesh", "JPT" = "Japan",
                  "PEL" = "South America")[fam$Population[ind_row]]))
pop <- factor(pop, unique(pop))

colors <- c("#b7ae38", "#426de0", "#78b436", "#ac56ce", "#4cbd64", "#d556b7",
            "#457836", "#8d7ce9", "#da9334", "#7652b2", "#97ae65", "#d9417f",
            "#54b38e", "#ce3a46", "#3abbcc", "#d95a2f", "#659ad7", "#a5572e",
            "#5466ad", "#806f29", "#b68fdd", "#d79e6a", "#975195", "#d67473",
            "#dc86b6", "#9f4566")
PC <- predict(obj.svd)

all_plots <- lapply(seq(0, 18, by = 2), function(offset) {

  k <- c(2, 1) + offset
  plot(obj.svd, type = "scores", scores = k, coef = 0.6) +
    coord_cartesian() + ggtitle(NULL) +
    aes(color = pop) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    geom_point(aes(x = pred[, k[1]], y = pred[, k[2]]), color = "black", size = 2) +
    # ggrepel::geom_text_repel(
    #   aes(x = c(PC[, k[1]], pred[, k[1]]), y = c(PC[, k[2]], pred[, k[2]]),
    #       label = c(sample(c(NA, ""), size = nrow(PC), replace = TRUE, prob = c(50, 1)),
    #                 rownames(pred))), na.rm = TRUE, size = 2.5, force = 2,
    #   color = "black", segment.linetype = 3, min.segment.length = 0, max.overlaps = Inf) +
    labs(color = "Ancestry group")
})

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
plot_grid2(plotlist = all_plots, ncol = 2, title_ratio = 0, legend_ratio = 0.3)
# ggsave("figures/all_scores.png", width = 8, height = 11)

# Use 16 PCs
# bigreadr::fwrite2(cbind.data.frame(all_freq[1:5], project[, 1:16]), "projection.csv.gz")


## Estimate shrinkage
data_1kg <- readRDS("tmp-data/info_af_1kg.rds")
ind <- vctrs::vec_match(all_freq[1:5], data_1kg[1:5])
pred_1kg <- crossprod(as.matrix(data_1kg[ind, -(1:5)]), project)

pop_ukbb <- c("Africa (West)", "Europe (South West)", "Sri Lanka", "Asia (East)", "United Kingdom", "Italy")
pop_1kg <- c("YRI", "IBS", "STU", "CHS", "GBR", "TSI")
coef <- sapply(cols_along(w), function(j)
  MASS::rlm(y = pred[pop_ukbb, j], x = pred_1kg[pop_1kg, j, drop = FALSE])$coef)
plot(coef[1:18]); points(coef_smooth <- pmax(lowess(coef[1:18])$y[1:16], 1), pch = 20, col = "red")
dput(round(coef_smooth, 3))
# c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099, 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
