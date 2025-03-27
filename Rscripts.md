**R code for data anylses** ##update w Tara's code 

R script for Variant Density Plot:
```
snpdens <- read_tsv("snpdens.tsv")
str(snpdens)
unique(snpdens$CHROM)
snpdens$CHROM <- as.factor(snpdens$CHROM)
any(is.na(snpdens$CHROM))
snpdens <- as.data.frame(snpdens)

snpdenchr <- snpdens %>%
  filter(!str_starts(CHROM, "NW")) %>%
  mutate(CHROM = case_when(
    CHROM == "NC_014690.1" ~ "MT",
    CHROM == "NC_055957.1" ~ "01",
    CHROM == "NC_055958.1" ~ "02",
    CHROM == "NC_055959.1" ~ "03",
    CHROM == "NC_055960.1" ~ "04",
    CHROM == "NC_055961.1" ~ "05",
    CHROM == "NC_055962.1" ~ "06",
    CHROM == "NC_055963.1" ~ "07",
    CHROM == "NC_055964.1" ~ "08",
    CHROM == "NC_055965.1" ~ "09",
    CHROM == "NC_055966.1" ~ "10",
    CHROM == "NC_055967.1" ~ "11",
    CHROM == "NC_055968.1" ~ "12",
    CHROM == "NC_055969.1" ~ "13",
    CHROM == "NC_055970.1" ~ "14",
    CHROM == "NC_055971.1" ~ "15",
    CHROM == "NC_055972.1" ~ "16",
    CHROM == "NC_055973.1" ~ "17",
    CHROM == "NC_055974.1" ~ "18",
    CHROM == "NC_055975.1" ~ "19",
    CHROM == "NC_055976.1" ~ "20",
    CHROM == "NC_055977.1" ~ "21",
    CHROM == "NC_055978.1" ~ "22",
    CHROM == "NC_055979.1" ~ "23",
    CHROM == "NC_055980.1" ~ "24",
    TRUE ~ CHROM 
  )) %>%  filter(CHROM != "MT")
snpplot <- ggplot(snpdenchr, aes(x = BIN_START, y = `VARIANTS/KB`, color = CHROM)) +
  geom_col() +
  scale_color_viridis_d(name = "CHROM") + 
  labs(
      x = "Position", 
      y = "Variant Density",
      title = "Variant Density Across Genome") +
  theme_bw () +
  facet_wrap(~CHROM, scales = "free_x")

print(snpplot)
```

#Tara insert data 





