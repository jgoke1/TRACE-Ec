snp <- read.delim("snp_dist_matrix.tsv")

snpps <- Full_amr %>% select(ERR, Source, ST)

unique_snpps <- unique(snpps$ERR)

df_unique <- snpps %>%
  distinct(ERR, .keep_all = TRUE)


write.csv(latest_mlst, "mlst.csv", row.names = FALSE)


ACTHMAN_MLST <- read_csv("ACTHMAN_MLST.csv")
View(ACTHMAN_MLST)

ACTHMAN_MLST <- ACTHMAN_MLST %>% select(ERR, ST)

df_unique <- df_unique %>% select(ERR, Source)

latest_mlst <- left_join(ACTHMAN_MLST, df_unique, by = "ERR")
