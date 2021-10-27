data_1kg <- readRDS("tmp-data/info_af_1kg.rds")
## Get all SNP IDs and split them in blocks
all_snp_id <- with(data_1kg, split(paste(chr, pos, a0, a1, sep = "_"), as.factor(1:22)[chr]))
lengths(all_snp_id)
#      1       2       3       4       5       6       7       8       9      10      11
# 981656 1036649  878387  899012  782944  819453  724662  674884  530864  627816  610113
#     12      13      14      15      16      17      18      19      20      21      22
# 598814  448686  404642  358726  379725  346916  352132  293684  275652  172731  175518

library(dplyr)
grid <- purrr::map_dfr(1:22, function(chr) {
  bgenfile <- paste0("UKBB/bgen/ukb_imp_chr", chr, "_v3.bgen")
  all_splits <- bigparallelr::split_vec(all_snp_id[[chr]], block_len = 50e3)
  tibble::tibble(bgenfile, snp_id = all_splits)
})
nrow(grid) # 257 parts

all_eid <- readRDS("tmp-data/info-UKBB2.rds")$eid
eid_sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_1[-1]
all_ind <- lapply(readRDS("data/list_ind_pop.rds"), function(ind) {
  match(all_eid[ind], eid_sample)
})

library(future.batchtools)
NCORES <- 15
plan(batchtools_slurm(workers = print(nrow(grid)), resources = list(
  t = "12:00:00", c = NCORES, mem = "120g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("tmp-data/bgen")
bigassertr::assert_dir("all_freq_ukbb")
bigassertr::assert_dir("log")

furrr::future_pwalk(mutate(grid, id = row_number()), function(bgenfile, snp_id, id) {

  for (ipop in seq_along(all_ind)) {

    res_file <- paste0("all_freq_ukbb/pop", ipop, "_part", id, ".rds")

    runonce::save_run({

      tmp <- tempfile(tmpdir = "tmp-data/bgen")
      on.exit(file.remove(paste0(tmp, c(".bk", ".rds"))), add = TRUE)

      library(bigsnpr)
      snp_attach(snp_readBGEN(
        bgenfiles   = bgenfile,
        list_snp_id = list(snp_id),
        backingfile = tmp,
        ind_row     = all_ind[[ipop]],
        ncores      = NCORES
      ))$map[c("freq", "info")]

    }, file = print(res_file))
  }

})

all_map <- lapply(seq_along(all_ind), function(ipop) {
  print(ipop)
  res_files <- paste0("all_freq_ukbb/pop", ipop, "_part", rows_along(grid), ".rds")
  purrr::map_dfr(res_files, readRDS)
})
# saveRDS(setNames(all_map, names(all_ind)), "tmp-data/info_af_ukbb.rds")
