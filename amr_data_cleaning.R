library(dplyr)
library(tidyverse)



africa <- Africa_esbl %>% select(ERR, `Organism Name`, YEAR, Country, Source, )
vivli <- vivli_africa_metadata %>% select(IDS, ST)

vivli <- rename(vivli, ERR = IDS)
viv_africa <- left_join(africa, vivli, by = "ERR")
viv_africa <- left_join(viv_africa, longread, by = "ERR")
viv_africa <- viv_africa %>% 
  filter(!grepl("longread", TYPE))
viv_africa <- distinct(viv_africa)


amr <- read.delim("combined_amrfinder.tsv")
amr_type <- amr %>% filter(Type == "AMR")
amr_beta <- amr_type %>% filter(Class == "BETA-LACTAM")
amr_cleaned <- amr_beta %>% 
  filter(grepl("CEPHALOSPORIN|Carbapenem", Subclass, ignore.case = TRUE))


amr_cleaned <- rename(amr_cleaned, Gene = Element.symbol)

amr_cleaned <- amr_cleaned %>% select(ERR, Gene, Class, Subclass)

amr_wide <- pivot_wider(amr_cleaned, names_from = Gene, values_from = Subclass)

Full_amr <- left_join(viv_africa, amr_cleaned, by = "ERR")

Full_amr <- Full_amr %>% select(-TYPE)

write.csv(Full_amr, "cleaned_amr.csv", row.names = FALSE)
