# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
pgen_zst <- runonce::download_file(
  "https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1",
  dir = "tmp-data")
pvar_zst <- runonce::download_file(
  "https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1",
  dir = "tmp-data")
psam <- runonce::download_file(
  "https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1",
  dir = "tmp-data", fname = "all_phase3.psam")

library(bigsnpr)
plink2 <- download_plink2("tmp-data")

sysglue <- function(...) system(glue::glue(...))
pgen <- sub('\\.zst', '', pgen_zst)
if (!file.exists(pgen)) sysglue("{plink2} --zst-decompress {pgen_zst} > {pgen}")
pvar <- sub('\\.zst', '', pvar_zst)
if (!file.exists(pvar)) sysglue("{plink2} --zst-decompress {pvar_zst} > {pvar}")

# https://dx.doi.org/10.1016%2Fj.ajhg.2017.03.004
outliers <- c("NA20314", "HG00731", "HG00732", "HG01880", "HG01882",
              "HG01944", "HG02497", "NA20320", "NA20321")
bigsnpr:::write.table2(cbind.data.frame(0, outliers), (tmp <- tempfile()))

sysglue("{plink2} --pfile tmp-data/all_phase3",
        " --remove {tmp} --keep-founders",
        " --max-alleles 2 --mac 10 --autosome",
        " --make-bed --out tmp-data/all_phase3",
        " --threads {nb_cores()[[1]]} --memory 50e3")

obj.bed <- bed("tmp-data/all_phase3.bed")

library(doParallel)
cl <- makeCluster(22)
snp_info <- do.call("rbind", parLapply(cl, 1:22, function(chr) {
  library(dplyr)
  paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
    bigreadr::fread2(showProgress = FALSE) %>%
    filter(V6 > 0.001, V8 > 0.3) %>%
    { .[!vctrs::vec_duplicate_detect(.$V2), ] } %>%
    transmute(chr, pos = V3, a0 = V4, a1 = V5)
}))
stopCluster(cl)

library(dplyr)
matched <- snp_match(
  transmute(obj.bed$map, rsid = marker.ID,
            chr = as.integer(chromosome), pos = physical.pos,
            a0 = allele2, a1 = allele1, beta = 1),
  snp_info
) %>%
  filter(beta > 0) %>%
  select(-beta)
# 23,675,459 variants to be matched.
# 1,999,870 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 12,373,833 variants have been matched; 0 were flipped and 167 were reversed.


fam <- dplyr::left_join(
  obj.bed$fam,
  bigreadr::fread2("https://figshare.com/ndownloader/files/31080292"))
table(fam$Population)
# ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD
#  92  53  86  93  99 103 105  94  99  99  91 101 113
# IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
# 107 102 104  99  97  85  64  84  96 102  99 107 108
# saveRDS(table(fam$Population), "tmp-data/pop_1kg_size.rds")

all_freq_1000G <- do.call(cbind, tapply(rows_along(obj.bed), fam$Population, function(ind) {
  bed_MAF(obj.bed, ind.row = ind, ind.col = matched$`_NUM_ID_.ss`, ncores = nb_cores())$af
}))
# saveRDS(cbind.data.frame(matched[1:5], all_freq_1000G), "tmp-data/info_af_1kg.rds")
