#### Prepare variants to use ####

snp_qc <- bigreadr::fread2(runonce::download_file(
  "https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt",
  dir = "tmp-data"
))

af_qc <- bigreadr::fread2(runonce::download_file(
  "http://kunertgraf.com/data/files/snp_maf_comparison.csv",
  dir = "tmp-data"
))


library(dplyr)
snp_qc$qc_all <- rowMeans(select(snp_qc, ends_with("_qc")))
snp_qc2 <- snp_qc %>%
  select(-ends_with("_qc"), -ends_with("_loading")) %>%
  left_join(af_qc, by = c("rs_id" = "rsid"))

plot(as.factor(snp_qc2$oob), snp_qc2$qc_all)

snp_qc_final <- snp_qc %>%
  select(-2, -3, -ends_with("_qc"), -ends_with("_loading")) %>%
  as_tibble() %>%
  filter(array == 2, qc_all == 1) %>%
  anti_join(filter(af_qc, oob == 1), by = c("rs_id" = "rsid")) %>%
  select(2:5)

library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {
    paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
      bigreadr::fread2(showProgress = FALSE) %>%
      filter(V6 > 0.01, V8 == 1) %>%   ## MAF > 1% & INFO = 1
      dplyr::semi_join(
        filter(snp_qc_final, chromosome == !!chr),
        by = c("V3" = "position", "V4" = "allele1_ref", "V5" = "allele2_alt")
      ) %>%
      with(paste(chr, V3, V4, V5, sep = "_"))
  }
}) # 40 sec
stopCluster(cl)

sum(lengths(list_snp_id))  # 586,534

obj.bed <- bed(download_1000G("tmp-data"))
list_snp_id_1kg <- with(obj.bed$map, paste(chromosome, physical.pos, allele2, allele1, sep = "_"))

list_snp_id2 <- lapply(list_snp_id, intersect, y = list_snp_id_1kg)
sum(lengths(list_snp_id2))  # 470,878
lengths(list_snp_id2)
# 36992 37176 31338 29663 27834 33759 25869 24112 21171 23708 23058
# 21881 16323 15268 15206 16697 15554 14302 13234 12546  7228  7959

#### Prepare individual genotypes ####

all_eid <- readRDS("tmp-data/info-UKBB2.rds")$eid
eid_sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_1[-1]
all_ind <- readRDS("data/list_ind_pop.rds")

bigsnpr::snp_readBGEN(
  bgenfiles = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
  list_snp_id = list_snp_id2,
  backingfile = "tmp-data/UKBB_geno",
  ind_row = match(all_eid[unlist(all_ind)], eid_sample),
  read_as = "random",
  ncores = nb_cores()
)

ukbb <- snp_attach("tmp-data/UKBB_geno.rds")
G <- ukbb$genotypes
ind_keep <- snp_clumping(G, as.integer(ukbb$map$chromosome), ncores = nb_cores())
ind_keep2 <- match(unlist(list_snp_id2)[ind_keep], list_snp_id_1kg)


to_write <- big_apply(G, function(X, ind, ind_keep, obj.bed, ind_keep2) {
  get_1kg <- bigsnpr:::read_bed_scaled(obj.bed, rows_along(obj.bed), ind_keep2[ind],
                                       rep(0, length(ind)), rep(1, length(ind)))
  apply(rbind(X[, ind_keep[ind]], get_1kg), 2, paste, collapse = "")
}, ind = seq_along(ind_keep), ind_keep = ind_keep, obj.bed = obj.bed, ind_keep2 = ind_keep2,
a.combine = 'c', block.size = 100, ncores = nb_cores())

writeLines(to_write, "tmp-data/for_snmf.geno")

obj.snmf <- LEA::snmf("tmp-data/for_snmf.geno", K = 15,
                      CPU = nb_cores(), iterations = 1000)

fam2 <- bigreadr::fread2(paste0(obj.bed$prefix, ".fam2"))
group <- unlist(list(rep(names(all_ind), lengths(all_ind)), fam2$Population))
group <- ordered(group, levels = unique(group))
Q <- read.table("tmp-data/for_snmf.snmf/K15/run1/for_snmf_r1.15.Q")
main_comp_by_group <- by(Q, group, function(x) which.max(colSums(x)))
val_main_comp_group <- Q[cbind(seq_along(group), main_comp_by_group[group])]

library(dplyr)
tbl <-
  data.frame(grp = group, val = val_main_comp_group) %>%
  bind_cols(Q) %>%
  group_by(grp) %>%
  mutate(id = rank(-val, ties.method = "first"), val = NULL) %>%
  tidyr::pivot_longer(-c(grp, id))

library(ggplot2)
cbPalette <- c("#64cb64", "#e7be28", "#b94a7b", "#6687d2", "#add042", "#cc4eb7",
               "#e6831f", "#8d5edf", "#5ccda0", "#613d9a", "#c79335", "#46d0e5",
               "#d14f29", "#60924e", "#a45441")

ggplot(tbl) +
  geom_col(aes(id, value, color = name, fill = name)) +
  theme_bw(13) +
  scale_x_continuous(breaks = NULL) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  theme(legend.position = "none") +
  facet_wrap(~ grp, nrow = 9, scales = "free_x") +
  labs(x = "Individual # (ordered by main component of group)",
       y = "Ancestry proportion", color = "Ancestry", fill = "Ancestry")
# ggsave("figures/admixture.pdf", width = 13, height = 16)
