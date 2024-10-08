library(bigreadr)
library(bigsnpr)
library(dplyr)
bigassertr::assert_dir("res_ancestry")

all_freq <- readRDS("data/all_freq.rds")
# fwrite2(all_freq, "ref_freqs.csv.gz")
projection <- as.matrix(fread2("projection.csv.gz", select = paste0("PC", 1:16)))

snp_match_ancestry <- function(sumstats, join_by_pos = TRUE) {

  matched <- bigsnpr::snp_match(
    transform(sumstats, chr = as.integer(chr), beta = 1),
    all_freq[1:5],
    join_by_pos = join_by_pos,
    match.min.prop = 0.1
  ) %>%
    mutate(freq = ifelse(beta > 0, freq, 1 - freq))


  freq_ref <- all_freq[matched$`_NUM_ID_`, -(1:5)]   #[c(1, 12, 15, 16, 19)]

  res <- snp_ancestry_summary(
    freq = matched$freq,
    info_freq_ref = freq_ref,
    projection = projection[matched$`_NUM_ID_`, ],
    correction = c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                   1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
  )


  # system.time(
  #   res2 <- unlist(
  #     Summix::summix(
  #       data = cbind(freq = matched$freq, freq_ref),
  #       reference = names(freq_ref),
  #       observed = "freq"
  #     )[names(freq_ref)]
  #   )
  # ) # > 3 min for 5 pops
  # Summix v1.99.2 does not converge when using too many reference populations


  # then this is equivalent to what Summix does (no PCA projection):
  X <- as.matrix(freq_ref)
  y <- matched$freq
  # solve QP problem using https://stats.stackexchange.com/a/21566/135793
  res3 <- quadprog::solve.QP(
    Dmat = crossprod(X),
    dvec = crossprod(y, X),
    Amat = cbind(1, diag(ncol(X))),
    bvec = c(1, rep(0, ncol(X))),
    meq  = 1
  )$solution  # < 1 sec
  # all.equal(res3, res2) # same


  # also try with a subset of 100K variants
  ind <- sample(nrow(freq_ref), 100e3)
  res4 <- snp_ancestry_summary(
    freq = matched$freq[ind],
    info_freq_ref = freq_ref[ind, ],
    projection = projection[matched$`_NUM_ID_`[ind], ],
    correction = c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                   1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)
  )


  list(res, setNames(round(res3, 7), names(res)), res4)
}


## T2D from Biobank Japan

# zip <- runonce::download_file(
#   "https://humandbs.biosciencedbc.jp/files/hum0197/hum0197.v3.BBJ.T2D.v1.zip",
#   dir = "tmp-data", fname = "sumstats_T2D_BBJ2.zip")
# unzip(zip, exdir = "tmp-data")
sumstats <- fread2(
  "../misspec/tmp-data/hum0197.v3.BBJ.T2D.v1/GWASsummary_T2D_Japanese_SakaueKanai2020.auto.txt.gz",
  select = c("CHR", "POS", "Allele2", "Allele1", "AF_Allele2", "imputationInfo",
             "BETA", "SE", "N"),
  col.names = c("chr", "pos", "a1", "a0", "freq", "info", "beta", "beta_se", "N")
) %>%
  filter(chr %in% 1:22)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/BBJ.rds")


## T2D from FinnGen

(N <- 29193 + 182573)
sumstats <- fread2("../atw-pred/tmp-data/summary_stats_finngen_R5_T2D.gz") %>%
  transmute(rsid = rsids, chr = `#chrom`, pos, a0 = ref, a1 = alt,
            freq = (2 * (n_hom_cases + n_hom_controls) +
                      (n_het_cases + n_het_controls)) / (2 * N)) %>%
  filter(chr %in% 1:22)

(res <- snp_match_ancestry(sumstats, join_by_pos = FALSE))
saveRDS(res, "res_ancestry/FinnGen.rds")


## GERA

tgz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007143/GERA-TC.tsv.gz",
  dir = "tmp-data")
sumstats <- fread2(
  tgz, select = c("chromosome", "position", "Allele 1", "Allele 2",
                  "Effect allele frequency (EAF)", "Sample size"),
  col.names = c("chr", "pos", "a1", "a0", "freq", "N")
) %>%
  filter(N > 85e3)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/GERA.rds")


## Height in PAGE study

gz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008053/WojcikG_PMID_invn_rheight_alls.gz",
  dir = "tmp-data")
sumstats <- fread2(
  gz, select = c("Chr", "Position_hg19", "Other-allele", "Effect-allele", "Effect-allele-frequency"),
  col.names = c("chr", "pos", "a0", "a1", "freq"))

pop <- c("African-American" = 17286, "Hispanic-Latino" = 22192, "Asian" = 4680,
         "Native Hawaiian" = 3939, "Native American" = 647,	Other	= 1052)
round(sort(pop / sum(pop), decreasing = TRUE), 3)
# Hispanic-Latino African-American    Asian  Native Hawaiian    Other  Native American
#           0.446            0.347    0.094            0.079    0.021            0.013

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/PAGE.rds")


## Urate

tgz <- runonce::download_file(
  "https://ckdgen.imbi.uni-freiburg.de/files/Tin2019/urate_chr1_22_LQ_IQ06_mac10_all_741_nstud37_summac400_rsid.txt.gz",
  dir = "tmp-data")
sumstats <- fread2(
  tgz, select = c("Chr", "Pos_b37", "Allele2", "Allele1", "Freq1", "n_total_sum"),
  col.names = c("chr", "pos", "a0", "a1", "freq", "N"), fill = TRUE
) %>%
  mutate_at(3:4, toupper)
hist(sumstats$N)
sumstats <- filter(sumstats, N > 4e5)

pop <- c("European" = 288649, "East Asian" = 125725, "African American" = 33671,
         "South Asian" = 9037, "Hispanic" = 608)
round(sort(pop / sum(pop), decreasing = TRUE), 3)
# European       East Asian African American      South Asian         Hispanic
#    0.631            0.275            0.074            0.020            0.001

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/urate.rds")


## Height in Arabic countries

tsv <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013304/GCST90013304_buildGRCh37.tsv",
  dir = "tmp-data")
readLines(tsv, n = 3)
sumstats <- fread2(
  tsv, select = c("chromosome", "base_pair_location", "other_allele", "effect_allele", "effect_allele_frequency"),
  col.names = c("chr", "pos", "a0", "a1", "freq"),
)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/arabic.rds")


## T2D in Africa

txt <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008114/ChenJ_31049640",
  dir = "tmp-data")
readLines(txt, n = 10)
sumstats <- fread2(txt, select = c("MarkerName", "Allele1", "Allele2", "Freq1"),
                   col.names = c("snpid", "a1", "a0", "freq")) %>%
  mutate(chr = as.integer(sub("^(.+):.+:.+_.+$", "\\1", snpid)),
         pos = as.integer(sub("^.+:(.+):.+_.+$", "\\1", snpid)),
         snpid = NULL, a1 = toupper(a1), a0 = toupper(a0))

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/t2d_africa.rds")


## COVID

tgz <- runonce::download_file(
  "https://storage.googleapis.com/covid19-hg-public/20210415/results/20210607/COVID19_HGI_C2_ALL_leave_23andme_20210607.b37.txt.gz",
  dir = "tmp-data")
sumstats <- fread2(
  tgz, select = c("#CHR", "POS", "REF", "ALT", "all_meta_AF", "all_meta_N"),
  col.names = c("chr", "pos", "a0", "a1", "freq", "N")
)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/covid.rds")


# Epilepsy (supposedly in EUR+EAS+AFR)

gz <- runonce::download_file(
  "http://www.epigad.org/gwas_ilae2018_16loci/all_epilepsy_METAL.gz",
  dir = "tmp-data")
sumstats <- fread2(
  gz, select = c("CHR", "BP", "Allele2", "Allele1", "Freq1"),
  col.names = c("chr", "pos", "a0", "a1", "freq")
) %>%
  mutate_at(3:4, toupper)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/epilepsy.rds")


## Body fat percentage in multi-ancestry

txt <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003435/body_fat_percentage_GWAS_PLUS_MC_ALL_ancestry_se_Sex_combined_for_locus_zoom_plot.TBL.txt",
  dir = "tmp-data")
sumstats <- fread2(
  txt, select = c("MarkerName", "Allele2", "Allele1", "Freq1", "SNPID", "N"),
  col.names = c("chr:pos", "a0", "a1", "freq", "rsid", "N")) %>%
  tidyr::separate(col = 1, into = c("chr", "pos"), sep = ":") %>%
  mutate(chr = readr::parse_number(chr), pos = as.integer(pos)) %>%
  mutate_at(3:4, toupper)
hist(sumstats$N)
sumstats <- filter(sumstats, 7e4 < N, N < 8e4)  # make sure same samples

(res <- snp_match_ancestry(sumstats, join_by_pos = FALSE))
saveRDS(res, "res_ancestry/body-fat.rds")


## Height in Peru

tgz <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002224/GCST90002224_buildGRCh37.tsv.gz",
  dir = "tmp-data")

sumstats <- fread2(
  tgz, select = c("chromosome", "base_pair_location", "other_allele", "effect_allele", "effect_allele_frequency"),
  col.names = c("chr", "pos", "a0", "a1", "freq"),
)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/height-peru.rds")


#### Eczema EAGLE ####

txt <- runonce::download_file(
  "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003184/EAGLE_AD_no23andme_results_29072015.txt",
  dir = "tmp-data")
sumstats <- fread2(
  txt, select = c("chromosome", "position", "reference_allele", "other_allele", "eaf"),
  col.names = c("chr", "pos", "a1", "a0", "freq"))

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/eczema.rds")


## PrCa from Schumacher et al., NG 2018

# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "tmp-data/sumstats_PRCA.zip")
# unzip("tmp-data/sumstats_PRCA.zip", exdir = "tmp-data")
sumstats <- fread2(
  "tmp-data/meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect", "StdErr",
             "Pvalue", "Freq1", "OncoArray_imputation_r2"),
  col.names = c("chr", "pos", "a1", "a0", "beta", "beta_se", "p", "freq", "info")
) %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1)) %>%
  filter(info > 0.4, beta_se > 0)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/PrCa.rds")


## BrCa iCOGs

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "tmp-data/sumstats_BRCA.txt.gz")
# R.utils::gunzip("tmp-data/sumstats_BRCA.txt.gz")
sumstats <- fread2("tmp-data/sumstats_BRCA.txt", na.strings = "NULL",
                   select = c("chr", "position_b37", "a0", "a1",
                              "bcac_icogs2_eaf_controls",
                              "bcac_icogs2_beta",
                              "bcac_icogs2_se",
                              "bcac_icogs2_r2"),
                   col.names = c("chr", "pos", "a0", "a1", "freq",
                                 "beta", "beta_se", "info")) %>%
  filter(info > 0.4, beta_se > 0)

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/BrCa.rds")


## CAD from Nikpay et al, NG 2015

# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "tmp-data/sumstats_CAD.txt")
sumstats <- fread2("../misspec/tmp-data/sumstats_CAD.txt",
                   select = c("chr", "bp_hg19", "noneffect_allele",
                              "effect_allele", "beta", "se_dgc", "p_dgc",
                              "effect_allele_freq", "median_info", "n_studies", "het_pvalue"),
                   col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se",
                                 "p", "freq", "info", "n_studies", "het_pvalue"))

(res <- snp_match_ancestry(sumstats))
saveRDS(res, "res_ancestry/CAD.rds")


## All of them

group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"
pop_names <- unique(group)

library(dplyr)
pop_N0 <- lengths(readRDS("data/list_ind_pop.rds"))
pop_N <- setNames(pop_N0[pop_names], pop_names)
pop_N["Africa (East)"] <- pop_N0["Africa (East 1)"]
pop_N["Europe (East)"] <- sum(pop_N0[c("Europe (South East)", "Europe (North East)")])
pop_N["Europe (North West)"] <- sum(pop_N0[c("Scandinavia", "United Kingdom", "Ireland")])
# TODO: add 1000G N later (in the table)

files <- list.files("res_ancestry", full.names = TRUE)
all_res <- setNames(lapply(files, readRDS), sub("\\.rds$", "", basename(files)))


library(dplyr)
purrr::iwalk(setNames(1:3, c("", "_summix", "_100K")), function(.x, .y) {
  all_res %>%
    purrr::map_dfc(function(three_res) tapply(three_res[[.x]], factor(group, pop_names), sum)) %>%
    bind_cols(Population = pop_names, N = pop_N, .) %>%
    select(1:2, BBJ, FinnGen, `height-peru`, arabic, t2d_africa, GERA, PAGE,
           BrCa, PrCa, CAD, everything()) %>%
    print(n = Inf, width = Inf) %>%
    fwrite2(paste0("all_prop", .y, ".csv")) %>%
    print()
})


# Model for exporting to LaTeX
df <- fread2("all_prop.csv") %>%
  as_tibble() %>%
  mutate_at(-(1:2), ~ round(100 * ., 1)) %>%
  print(n = Inf, width = Inf)

df %>%
  mutate_at(-(1:2), ~ ifelse(. == 0, "", as.character(.))) %>%
  xtable::xtable() %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 4, 5, 11, 12, 15, 18),
        include.rownames = FALSE, file = "ancestry-proportions.tex")
