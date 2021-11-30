data_1kg <- readRDS("tmp-data/info_af_1kg.rds")
data_ukbb <- readRDS("tmp-data/info_af_ukbb.rds")

INFO <- do.call(cbind, lapply(data_ukbb, function(.) .$info))
round(sort(colMeans(INFO, na.rm = TRUE)), 3)
# Japan         Asia (East)         Philippines          Bangladesh           Sri Lanka
# 0.628               0.659               0.680               0.748               0.752
# Pakistan    Africa (Central)      Africa (South)       Africa (West)     Africa (East 1)
#    0.768               0.785               0.792               0.799               0.811
# Africa (East 2) Europe (South East) Europe (North East)               India               Italy
#           0.813               0.815               0.820               0.820               0.821
# Europe (Central)             Ireland         Scandinavia         Middle East           Ashkenazi
#            0.822               0.825               0.835               0.835               0.846
# United Kingdom           Caribbean      Africa (North) Europe (South West)             Finland
#          0.856               0.869               0.874               0.876               0.878
# South America
#         0.909

mean_info <- rowMeans(INFO)
min_info <- matrixStats::rowMins(INFO)
hist(mean_info)
hist(min_info)
ind_keep <- which(min_info > 0.4)
length(ind_keep) # 5,840,630

X_ukbb <- do.call(cbind, lapply(data_ukbb, function(.) .$freq[ind_keep]))
summary(rowMeans(X_ukbb))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.004086 0.136399 0.274173 0.335182 0.497979 0.992190
X_1kg <- as.matrix(data_1kg[ind_keep, -(1:5)])

pop_1kg_size  <- readRDS("tmp-data/pop_1kg_size.rds")
pop_ukbb_size <- lengths(readRDS("data/list_ind_pop.rds"))

library(bigsnpr)

X_FBM <- as_FBM(cbind(X_ukbb, X_1kg), backingfile = "tmp-data/all_af")$save()
all_N <- c(pop_ukbb_size, pop_1kg_size)

npop <- length(all_N)
all_fst <- FBM(npop, npop, init = 0)
library(doFuture)
plan("multisession", workers = 15)
registerDoFuture()

foreach(i = 2:npop, .export = "all_fst") %:%
  foreach(j = seq_len(i - 1)) %dopar% {
    fst <- snp_fst(list(
      data.frame(af = X_FBM[, i], N = all_N[i]),
      data.frame(af = X_FBM[, j], N = all_N[j])
    ), overall = TRUE)
    all_fst[i, j] <- all_fst[j, i] <- print(fst)
    NULL
  }

all_fst <- all_fst[]
colnames(all_fst) <- rownames(all_fst) <- names(all_N)
# write.table(all_fst, "all_fst.csv", sep = ",", quote = FALSE)
# all_fst <- as.matrix(read.csv("all_fst.csv", check.names = FALSE))

print_fst <- function(pop, nb = 7) {
  smallest_fst <- sort(all_fst[, pop])[1:nb + 1L]
  paste0(" - [", pop, "] ",paste(purrr::imap_chr(
    smallest_fst, ~ glue::glue("{.y}: {signif(.x, 2)}")), collapse = " || "))
}
purrr::walk(names(pop_ukbb_size), ~ cat(print_fst(.), "\n"))
# - [Japan] JPT: 0.00033 || CHB: 0.007 || CHS: 0.0089 || Asia (East): 0.01 || KHV: 0.014 || CDX: 0.017 || Philippines: 0.021
# - [Asia (East)] CHS: 0.0009 || KHV: 0.0021 || CHB: 0.0029 || CDX: 0.0034 || Philippines: 0.01 || Japan: 0.01 || JPT: 0.011
# - [Philippines] Asia (East): 0.01 || KHV: 0.01 || CDX: 0.011 || CHS: 0.012 || CHB: 0.015 || Japan: 0.021 || JPT: 0.021
# - [Africa (South)] Africa (Central): 0.0017 || LWK: 0.0035 || Caribbean: 0.0043 || Africa (West): 0.0044 || Africa (East 2): 0.005 || YRI: 0.005 || ESN: 0.0055
# - [Africa (North)] PUR: 0.01 || Middle East: 0.01 || Italy: 0.01 || Europe (South West): 0.011 || TSI: 0.012 || IBS: 0.012 || Ashkenazi: 0.012
# - [Africa (East 1)] Africa (North): 0.018 || ASW: 0.027 || Africa (East 2): 0.028 || PUR: 0.029 || Middle East: 0.036 || ACB: 0.037 || CLM: 0.039
# - [Africa (East 2)] LWK: 0.0033 || Caribbean: 0.0045 || ASW: 0.0048 || ACB: 0.0048 || Africa (South): 0.005 || Africa (Central): 0.0056 || Africa (West): 0.0074
# - [Caribbean] ACB: 0.00048 || Africa (West): 0.0011 || YRI: 0.0022 || Africa (Central): 0.0023 || ESN: 0.0028 || ASW: 0.003 || MSL: 0.0042
# - [Africa (West)] YRI: 0.00061 || Caribbean: 0.0011 || ESN: 0.0015 || Africa (Central): 0.0019 || ACB: 0.0023 || MSL: 0.0028 || Africa (South): 0.0044
# - [Africa (Central)] Africa (South): 0.0017 || Africa (West): 0.0019 || Caribbean: 0.0023 || YRI: 0.0024 || ESN: 0.0028 || ACB: 0.0036 || LWK: 0.0045
# - [Middle East] Italy: 0.0043 || TSI: 0.0059 || Europe (South East): 0.007 || Ashkenazi: 0.007 || Europe (South West): 0.008 || India: 0.0088 || IBS: 0.0089
# - [United Kingdom] Scandinavia: 0.00043 || Ireland: 0.00069 || CEU: 0.001 || GBR: 0.0011 || Europe (Central): 0.0012 || Europe (South East): 0.0021 || Europe (South West): 0.0023
# - [Ireland] United Kingdom: 0.00069 || Scandinavia: 0.0015 || GBR: 0.0017 || CEU: 0.0018 || Europe (Central): 0.0023 || Europe (South West): 0.0033 || Europe (South East): 0.0033
# - [Finland] FIN: 0.00041 || Europe (North East): 0.0052 || Scandinavia: 0.0054 || Europe (Central): 0.0056 || United Kingdom: 0.0066 || CEU: 0.0067 || GBR: 0.0072
# - [Scandinavia] United Kingdom: 0.00043 || Europe (Central): 0.0012 || CEU: 0.0013 || Ireland: 0.0015 || GBR: 0.0015 || Europe (South East): 0.0026 || Europe (North East): 0.0027
# - [Europe (South West)] IBS: 0.00079 || Italy: 0.0014 || Europe (South East): 0.002 || TSI: 0.0021 || United Kingdom: 0.0023 || Europe (Central): 0.003 || CEU: 0.003
# - [Italy] TSI: 0.0011 || Europe (South West): 0.0014 || Europe (South East): 0.0017 || IBS: 0.0023 || Europe (Central): 0.004 || United Kingdom: 0.0041 || Middle East: 0.0043
# - [Europe (South East)] Europe (Central): 0.00075 || Italy: 0.0017 || Europe (North East): 0.0017 || Europe (South West): 0.002 || United Kingdom: 0.0021 || Scandinavia: 0.0026 || TSI: 0.0027
# - [Europe (Central)] Europe (North East): 0.00053 || Europe (South East): 0.00075 || Scandinavia: 0.0012 || United Kingdom: 0.0012 || CEU: 0.0019 || GBR: 0.0023 || Ireland: 0.0023
# - [Europe (North East)] Europe (Central): 0.00053 || Europe (South East): 0.0017 || Scandinavia: 0.0027 || United Kingdom: 0.003 || CEU: 0.0036 || Ireland: 0.0039 || GBR: 0.0039
# - [Pakistan] PJL: 0.0026 || India: 0.0041 || GIH: 0.0053 || Bangladesh: 0.0059 || Sri Lanka: 0.0065 || ITU: 0.0068 || STU: 0.0074
# - [Sri Lanka] STU: 0.00034 || ITU: 0.0015 || Bangladesh: 0.0023 || BEB: 0.0025 || PJL: 0.0037 || GIH: 0.0043 || Pakistan: 0.0065
# - [Bangladesh] BEB: 0.00068 || Sri Lanka: 0.0023 || STU: 0.0028 || ITU: 0.003 || PJL: 0.0037 || GIH: 0.0044 || Pakistan: 0.0059
# - [South America] MXL: 0.0023 || CLM: 0.0053 || PUR: 0.014 || PEL: 0.019 || India: 0.022 || Pakistan: 0.027 || Europe (South East): 0.029
# - [Ashkenazi] Italy: 0.0051 || TSI: 0.0064 || Europe (South West): 0.0069 || Middle East: 0.007 || Europe (South East): 0.0072 || IBS: 0.008 || Europe (Central): 0.0097
# - [India] Pakistan: 0.0041 || PJL: 0.007 || Middle East: 0.0088 || Europe (South East): 0.0091 || Bangladesh: 0.0095 || Europe (Central): 0.0096 || United Kingdom: 0.0097


all_fst_1kg <- all_fst[1:26, 27:52]
min_fst <- tibble::tibble(
  pop_1kg = colnames(all_fst_1kg),
  pop_ukbb = rownames(all_fst_1kg)[apply(all_fst_1kg, 2, which.min)],
  fst = apply(all_fst_1kg, 2, min)
)
print(dplyr::arrange(min_fst, fst), n = Inf)
#    pop_1kg pop_ukbb                 fst
#  1 JPT     Japan               0.000327
#  2 STU     Sri Lanka           0.000340
#  3 FIN     Finland             0.000411
#  4 ACB     Caribbean           0.000483
#  5 YRI     Africa (West)       0.000615
#  6 BEB     Bangladesh          0.000675
#  7 IBS     Europe (South West) 0.000794
#  8 CHS     Asia (East)         0.000897
#  9 CEU     United Kingdom      0.00101
# 10 TSI     Italy               0.00108
# 11 GBR     United Kingdom      0.00112
# 12 ITU     Sri Lanka           0.00148
# 13 ESN     Africa (West)       0.00149
# 14 KHV     Asia (East)         0.00210
# 15 MXL     South America       0.00235
# 16 PJL     Pakistan            0.00260
# 17 MSL     Africa (West)       0.00280
# 18 CHB     Asia (East)         0.00288
# 19 ASW     Caribbean           0.00295
# 20 LWK     Africa (East 2)     0.00330
# 21 CDX     Asia (East)         0.00337
# 22 GIH     Sri Lanka           0.00426
# 23 GWD     Africa (West)       0.00495
# 24 CLM     South America       0.00527
# 25 PUR     Europe (South West) 0.00941
# 26 PEL     South America       0.0191

all_scaled_fst <- purrr::pmap(min_fst, function(pop_1kg, pop_ukbb, fst) {
  var_fst <- snp_fst(list(
    data.frame(af = X_1kg[, pop_1kg],   N = pop_1kg_size[[pop_1kg]]),
    data.frame(af = X_ukbb[, pop_ukbb], N = pop_ukbb_size[[pop_ukbb]])
  ))
  var_fst / fst
})
diff_score <- unname(Reduce(`+`, all_scaled_fst))
hist(diff_score, "FD", xlim = c(-100, 400))
(thr <- bigutilsr::tukey_mc_up(diff_score)) # 280
sum(diff_score >= thr) # 24,040

ind <- which(diff_score < thr)

X <- as.data.frame(X_ukbb[ind, ])
# Increase sample size by merging with 1000G
X[, "Finland"] <- (143 * X[, "Finland"] + 99 * X_1kg[ind, "FIN"]) / (143 + 99)
X[, "Bangladesh"] <- (223 * X[, "Bangladesh"] + 86 * X_1kg[ind, "BEB"]) / (223 + 86)
X[, "Japan"] <- (240 * X[, "Japan"] + 104 * X_1kg[ind, "JPT"]) / (240 + 104)
X[, "South America"] <- (473 * X[, "South America"] + 84 * X_1kg[ind, "PEL"]) / (473 + 84)

# Remove some groups
X2 <- dplyr::select(X, -c(`Africa (East 2)`, Caribbean, `Africa (Central)`, `Europe (Central)`, India),
                    `Africa (East)` = `Africa (East 1)`)

ord <- c("Africa (West)", "Africa (South)", "Africa (East)", "Africa (North)", "Middle East",
         "Ashkenazi", "Italy", "Europe (South East)", "Europe (North East)", "Finland",
         "Scandinavia", "United Kingdom", "Ireland", "Europe (South West)", "South America",
         "Sri Lanka", "Pakistan", "Bangladesh", "Asia (East)", "Japan", "Philippines")
stopifnot(ncol(X2) == length(ord))
dim(X2[ord])

# saveRDS(cbind.data.frame(data_1kg[ind_keep[ind], 1:5], X2[ord]), "data/all_freq.rds")
