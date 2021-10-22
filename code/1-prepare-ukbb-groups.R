library(dplyr)
library(bigreadr)

#### Get PC / ancestry info from UKBB ####

# file.symlink("~/NCRR-PRS/faststorage/UKBB/", ".")

df_UKBB <- runonce::save_run({

  # Self-reported ancestry
  code_ancestry <- fread2("UKBB/coding1001.tsv")
  # Country of birth
  code_country   <- filter(fread2("UKBB/coding89.tsv"), selectable == "Y")
  code_continent <- filter(fread2("UKBB/coding89.tsv"), selectable == "N")

  eid_imputed <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")$ID_1[-1]

  fread2(
    "UKBB/ukb41181.csv",
    select = c("eid", "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16)),
    col.names = c("eid", "ancestry", "country", paste0("PC", 1:16))
  ) %>%
    filter(eid %in% eid_imputed) %>%
    mutate(
      ancestry = factor(ancestry, levels = code_ancestry$coding,
                        labels = code_ancestry$meaning) %>%
        forcats::fct_recode("Other White" = "Any other white background",
                            "Other Asian" = "Any other Asian background",
                            "Other Asian" = "Asian or Asian British",
                            "Other Black" = "Any other Black background",
                            "Other Black" = "Black or Black British") %>%
        forcats::fct_other(drop = c("Prefer not to answer", "Do not know", "Mixed",
                                    "Other ethnic group", "Any other mixed background",
                                    "White and Asian", "White and Black Caribbean",
                                    "White and Black African"),
                           other_level = NA) %>%
        droplevels(exclude = NA),

      continent = factor(country, levels = code_country$coding,
                         labels = code_country$parent_id) %>%
        factor(labels = code_continent$meaning),

      country = factor(country, levels = code_country$coding,
                       labels = code_country$meaning)
    )
}, file = "tmp-data/info-UKBB2.rds")
df_UKBB$country[df_UKBB$ancestry == "British"] <- "United Kingdom"
df_UKBB$country[df_UKBB$ancestry == "Irish"]   <- "Ireland"
df_UKBB$continent[df_UKBB$country %in% c("United Kingdom", "Ireland")] <- "Europe"

print_table <- function(x, min_count = 5) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  tab <- tab[tab >= min_count]
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}

PC_UKBB <- as.matrix(select(df_UKBB, PC1:PC16))
continent <- df_UKBB$continent
country <- df_UKBB$country
ancestry <- df_UKBB$ancestry

## A bit of QC on the country of birth so that it better match the actual ancestry
table(ancestry, continent, exclude = NULL)

# Many from Africa have Asian ancestry
asian_in_africa <- which(continent == "Africa" &
                           ancestry %in% c("Indian", "Other Asian", "Pakistani", "Chinese"))
country[asian_in_africa] <- NA

# Check for European or UK ancestry
center_UK <- bigutilsr::geometric_median(PC_UKBB[which(country == "United Kingdom"), ])
dist_to_UK <- sqrt(bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center_UK, '-')))

library(ggplot2)
POP <- c("Morocco", "Tunisia", "Algeria", "Brazil", "Turkey", "Australia", "South Africa")
POP <- c("United Kingdom", "France", "Belgium", "Denmark", "Sweden", "Ireland", "Greece", "Russia", "Portugal")
qplot(log(dist_to_UK[country %in% POP]), geom = "density",
      fill = country[country %in% POP], alpha = I(0.4)) +
  theme_bw() +
  labs(fill = "Country of birth", x = "log-distance to UK") +
  geom_vline(xintercept = 4, linetype = 3) +
  # geom_vline(xintercept = c(2.5, 5), linetype = 3) +
  theme(legend.position = c(0.2, 0.7))
# ggsave("figures/noeur-dist-uk.pdf", width = 8, height = 6, scale = 0.8)
# ggsave("figures/eur-dist-uk.pdf", width = 8, height = 6, scale = 0.8)

# Probably European:
country[log(dist_to_UK) < 4 & continent != "Europe"] <- NA
# Maybe UK:
country[log(dist_to_UK) < 2.5 & continent == "Europe" & country != "United Kingdom"] <- NA
# Probably not European:
country[log(dist_to_UK) > 5 & continent == "Europe"] <- NA


# Keep only African ancestry for South Africa
qplot(log(dist_to_UK[country == "South Africa"]), geom = "density",
      fill = ancestry[country == "South Africa"], alpha = I(0.4)) +
  theme_bw() +
  labs(fill = "Self-reported ancestry", x = "log-distance to UK") +
  theme(legend.position = c(0.25, 0.7))
country[country == "South Africa" & log(dist_to_UK) < 6] <- NA

ind_nona <- which(!is.na(country))
table(ancestry[ind_nona], continent[ind_nona], exclude = NULL)


country2 <- forcats::fct_lump_min(country, 50, other_level = NA) %>%
  droplevels(exclude = NA)

print_table(country2)
# United Kingdom: 430432 || NA: 21594 || Ireland: 11738 || India: 3352 || Caribbean: 2490 ||
# Pakistan: 1368 || Germany: 1100 || Nigeria: 1036 || Ghana: 883 || Italy: 746 || France: 734 ||
# Sri Lanka: 671 || Poland: 595 || Iran: 504 || Hong Kong: 490 || Malaysia: 447 || Netherlands: 411 ||
# China: 380 || Spain: 339 || Philippines: 324 || Zimbabwe: 323 || Iraq: 309 || The Guianas: 298 ||
# Portugal: 290 || Barbados: 262 || USA: 260 || Japan: 249 || Bangladesh: 240 || Colombia: 220 ||
# Denmark: 219 || Sierra Leone: 217 || Sweden: 199 || Mauritius: 194 || Kenya: 179 || Singapore: 173 ||
# Turkey: 166 || Nepal: 161 || Brazil: 160 || Russia: 155 || Uganda: 153 || Finland: 153 || Congo: 151 ||
# Egypt: 147 || Cyprus: 142 || Switzerland: 141 || Austria: 136 || Greece: 114 || Czech Republic: 113 ||
# Somalia: 111 || Afghanistan: 111 || Norway: 111 || Hungary: 100 || Malta: 100 || Belgium: 96 ||
# Sudan: 90 || Myanmar (Burma): 90 || Thailand: 90 || Morocco: 78 || Vietnam: 73 || Algeria: 72 ||
# Lithuania: 71 || Zambia: 67 || Bulgaria: 67 || Chile: 67 || Mexico: 66 || Ethiopia: 65 || Romania: 63 ||
# Ukraine: 62 || Israel: 60 || Yemen: 59 || South Africa: 56 || Lebanon: 55 || Peru: 54 || Latvia: 52 ||
# Angola: 51 || Eritrea: 50 || Serbia/Montenegro: 50

sqdist_to_pop <- function(ind_pop, thr = 7) {
  X <- PC_UKBB[ind_pop, ]
  center_pop <- bigutilsr::geometric_median(X)
  sq_dist_to_pop <- bigutilsr:::rowSumsSq(sweep(X, 2, center_pop, '-'))
  X2 <- X[sq_dist_to_pop < exp(thr), ]
  center_pop2 <- bigutilsr::geometric_median(X2)
  bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center_pop2, '-'))
}

pop <- "Afghanistan"; thr <- 5.5
length(ind_pop <- which(country == pop))
hist(log(dist_to_UK[ind_pop]))
sq_dist_to_pop <- sqdist_to_pop(ind_pop)
qplot(log(sq_dist_to_pop[ind_pop]), geom = "density",
      fill = ancestry[ind_pop], alpha = I(0.4)) +
  theme_bw() + labs(fill = "Self-reported ancestry")
length(ind <- setdiff(which(sq_dist_to_pop < exp(thr)), which(country == "United Kingdom")))
print_table(country[ind], 2)
#

center_ashkenazi <- c(20.9739546693973, -11.4341656011295, 28.1040767857325, -85.6976204805128,
                      16.3211859722377, -8.62497185585874, -28.3339726063373, -18.400229766216,
                      0.127012311278296, 3.00741482315339, 0.384465030376849, 3.03032525906422,
                      0.654442583800136, -0.225132565608239, 0.0971148832455303, -0.061088005330329)
sq_dist_to_ash <- bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center_ashkenazi, '-'))
print_table(country[which(sq_dist_to_ash < 100)])
# United Kingdom: 1182 || NA: 572 || USA: 115 || Israel: 27 || Hungary: 12 || Ireland: 8 || France: 7 || Canada: 6 || Russia: 5

print_table(country[PC_UKBB[, 6] > 50 & PC_UKBB[, 8] < -30])
# Colombia: 173 || Chile: 57 || Mexico: 51 || Peru: 50 || Ecuador: 33 || Bolivia: 21 || NA: 20 || Venezuela: 15 || Brazil: 10 || United Kingdom: 8 || USA: 8 || Argentina: 5

grid <- tibble::tribble(
  ~pop, ~thr, ~name,
  "Japan", 5.5, "Japan",
  "Hong Kong", 5.5, "Asia (East)",
  "Philippines", 7, "Philippines",
  "South Africa", 6.5, "Africa (South)",
  "Morocco", 7.5, "Africa (North)",
  "Ethiopia", 8.5, "Africa (East 1)",
  "Kenya", 8.5, "Africa (East 2)",
  "Barbados", 5, "Caribbean",
  "Sierra Leone", 4.5, "Africa (West)",
  "Cameroon", 6, "Africa (Central)",
  "Iraq", 6, "Middle East",
  "United Kingdom", 4, "United Kingdom",
  "Ireland", 4, "Ireland",
  "Finland", 6, "Finland",
  "Norway", 4.5, "Scandinavia",
  "Spain", 5, "Europe (South West)",
  "Italy", 5.5, "Italy",
  "Serbia/Montenegro", 5.5, "Europe (South East)",
  "Czech Republic", 5, "Europe (Central)",
  "Russia", 5, "Europe (North East)",
  "Pakistan", 4.5, "Pakistan",
  "Sri Lanka", 4.5, "Sri Lanka",
  "Bangladesh", 6, "Bangladesh"
)

set.seed(1)
all_ind <- c(
  setNames(purrr::pmap(grid, function(pop, thr, name) {
    ind_pop <- which(country == pop)
    sq_dist_to_pop <- sqdist_to_pop(ind_pop)
    ind <- which(sq_dist_to_pop < exp(thr))
    if (pop == name) ind <- intersect(ind, ind_pop)
    if (pop != "United Kingdom") ind <- setdiff(ind, which(country == "United Kingdom"))
    `if`(length(ind) > 2000, sort(sample(ind, 2000)), ind)
  }), grid$name),
  list("South America" = unname(which(PC_UKBB[, 6] > 50 & PC_UKBB[, 8] < -30)),
       "Ashkenazi" = unname(which(sq_dist_to_ash < 100)))
)
c(length(Reduce(union, all_ind)), sum(lengths(all_ind)))  # 15224 // 15227


all_dist <- sapply(all_ind, function(ind) sqdist_to_pop(ind))
all_min_dist <- matrixStats::rowMins(all_dist)
hist(log(all_min_dist), "FD"); abline(v = 7.8, col = "red")
ind_no_match <- setdiff(which(all_min_dist > exp(7.8)), unlist(all_ind))
print_table(country[ind_no_match], 10)
# NA: 2632 || United Kingdom: 831 || Caribbean: 573 || India: 193 || Brazil: 97 || Nepal: 88 ||
# The Guianas: 73 || Mauritius: 68 || Zimbabwe: 56 || Myanmar (Burma): 54 || Malaysia: 45 ||
# Colombia: 42 || Barbados: 38 || Hong Kong: 30 || Singapore: 30 || Seychelles: 26 || Pakistan: 25 ||
# Ghana: 24 || Thailand: 24 || Kenya: 23 || USA: 23 || Afghanistan: 21 || Angola: 20 || Sri Lanka: 20 ||
# Somalia: 18 || Nigeria: 16 || Yemen: 16 || Ireland: 15 || Malawi: 14 || Sierra Leone: 13 ||
# New Zealand: 13 || Mozambique: 12 || Sudan: 12 || Uganda: 12 || Venezuela: 12 || Argentina: 11 ||
# Tanzania: 10 || Netherlands: 10 || Mexico: 10

thr <- 6.5
length(ind_pop <- ind_no_match[which(country[ind_no_match] == "India")])
sq_dist_to_pop <- sqdist_to_pop(ind_pop)
qplot(log(sq_dist_to_pop[ind_pop]), geom = "density",
      fill = ancestry[ind_pop], alpha = I(0.4)) +
  theme_bw() + labs(fill = "Self-reported ancestry")
length(ind <- setdiff(which(sq_dist_to_pop < exp(thr)), which(country == "United Kingdom")))
print_table(country[ind], 2)
all_ind[["India"]] <- ind

plot(eulerr::euler(all_ind[c(21:23, 26)]))  # disjoint
# saveRDS(all_ind, "data/list_ind_pop.rds")

purrr::iwalk(all_ind, ~ {cat("\n -", .y, ": "); print_table(country[.x])})
# - Japan : Japan: 240
# - Asia (East) : Hong Kong: 390 || Malaysia: 179 || NA: 122 || China: 114 || Singapore: 58 || Vietnam: 33 || Indonesia: 14 || Taiwan: 10 || Thailand: 8 || Brunei: 6 || Macau (Macao): 6
# - Philippines : Philippines: 295
# - Africa (South) : Zimbabwe: 251 || Uganda: 57 || South Africa: 47 || Zambia: 41 || Malawi: 10 || NA: 10 || Tanzania: 6 || Congo: 5
# - Africa (North) : Algeria: 69 || Morocco: 62 || Egypt: 54 || Libya: 36 || NA: 18 || Tunisia: 13
# - Africa (East 1) : Somalia: 81 || Ethiopia: 58 || Sudan: 54 || Eritrea: 46 || NA: 27
# - Africa (East 2) : Kenya: 104 || Uganda: 60 || Burundi: 18 || NA: 16 || Sudan: 14 || Tanzania: 14 || Rwanda: 13 || Nigeria: 8 || Somalia: 7 || Angola: 5 || Zimbabwe: 5
# - Caribbean : Caribbean: 327 || NA: 310 || Barbados: 73 || The Guianas: 23 || Antigua and Barbuda: 7
# - Africa (West) : Nigeria: 217 || Ghana: 175 || NA: 122 || Sierra Leone: 98 || Caribbean: 80 || Ivory Coast: 7 || Liberia: 7 || The Guianas: 7 || Togo: 6
# - Africa (Central) : Congo: 117 || Cameroon: 42 || Nigeria: 33 || NA: 33 || Angola: 21 || Caribbean: 8 || The Guianas: 7
# - Middle East : Iraq: 240 || Iran: 172 || Turkey: 55 || NA: 31 || Syria: 11
# - United Kingdom : United Kingdom: 2000
# - Ireland : Ireland: 2000
# - Finland : Finland: 143
# - Scandinavia : Denmark: 142 || Norway: 82 || Sweden: 75 || NA: 68 || Germany: 25 || Netherlands: 17 || Iceland: 5
# - Europe (South West) : Spain: 279 || Portugal: 198 || NA: 91 || France: 25 || Gibraltar: 7
# - Italy : Italy: 345
# - Europe (South East) : NA: 115 || Romania: 47 || Serbia/Montenegro: 41 || Bosnia and Herzegovina: 34 || Bulgaria: 33 || Croatia: 28 || Hungary: 17 || Italy: 11 || Macedonia: 8
# - Europe (Central) : NA: 196 || Germany: 111 || Czech Republic: 79 || Poland: 50 || Austria: 33 || Hungary: 27 || Slovakia: 25 || Croatia: 10 || Slovenia: 8 || France: 7
# - Europe (North East) : Russia: 88 || NA: 87 || Poland: 72 || Ukraine: 34 || Latvia: 10 || Germany: 7 || Kazakhstan: 6 || Lithuania: 5
# - Pakistan : Pakistan: 400
# - Sri Lanka : Sri Lanka: 372
# - Bangladesh : Bangladesh: 223
# - South America : Colombia: 173 || Chile: 57 || Mexico: 51 || Peru: 50 || Ecuador: 33 || Bolivia: 21 || NA: 20 || Venezuela: 15 || Brazil: 10 || United Kingdom: 8 || USA: 8 || Argentina: 5
# - Ashkenazi : United Kingdom: 1182 || NA: 572 || USA: 115 || Israel: 27 || Hungary: 12 || Ireland: 8 || France: 7 || Canada: 6 || Russia: 5
# - India : NA: 258 || India: 108 || Sri Lanka: 16 || Pakistan: 11

sort(lengths(all_ind))
#  Finland          Bangladesh               Japan      Africa (North)     Africa (East 1)
#      143                 223                 240                 268                 276
# Africa (Central)     Africa (East 2)         Philippines Europe (North East)               Italy
#              276                 279                 295                 320                 345
# Europe (South East)           Sri Lanka            Pakistan         Scandinavia               India
#                 347                 372                 400                 416                 418
# Africa (South)       South America         Middle East    Europe (Central) Europe (South West)
#            449                 473                 523                 556                 603
# Africa (West)           Caribbean         Asia (East)           Ashkenazi      United Kingdom
#           735                 752                 961                1975                2000
# Ireland
#    2000
