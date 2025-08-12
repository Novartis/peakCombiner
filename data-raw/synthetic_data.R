# Load library
library("tidyverse")
library("GenomicRanges")

# Define column names
names <- c("chrom", "start", "end", "name", "score", "strand", "center", "sample_name")

# Create the entire synthetic data
synthetic_data <- tibble("chr1", 200, 900, NA, 100, ".", 500, "treatment_rep1") |>
  rename_all(.fun = ~names) |>
  rbind(
    tibble("chr1", 0, 900, NA, 97, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 100, 300, NA, 94, ".", 200, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 300, 900, NA, 94, ".", 500, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 200, 900, NA, 100, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 300, 900, NA, 98, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 300, 1000, NA, 96, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 300, 1100, NA, 93, ".", 500, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 1300, 1600, NA, 97, ".", 1400, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 1900, 2200, NA, 98, ".", 2000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 2500, 3100, NA, 97, ".", 2800, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 2500, 3400, NA, 98, ".", 3000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 2600, 3200, NA, 99, ".", 2800, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 3500, 4200, NA, 44, ".", 3800, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 3500, 4400, NA, 95, ".", 3800, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 3600, 4400, NA, 43, ".", 3900, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4500, 5000, NA, 97, ".", 4800, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4500, 5200, NA, 60, ".", 4700, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4500, 5200, NA, 59, ".", 5000, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4500, 5300, NA, 98, ".", 4800, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4500, 5300, NA, 98, ".", 5100, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4600, 5100, NA, 93, ".", 4900, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4600, 5200, NA, 94, ".", 4800, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4700, 5300, NA, 46, ".", 4900, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4700, 5300, NA, 45, ".", 5100, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 5600, 6100, NA, 26, ".", 5700, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 5700, 6400, NA, 98, ".", 6200, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 5800, 6300, NA, 30, ".", 6100, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 25, ".", 7000, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 44, ".", 7000, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 43, ".", 7000, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 29, ".", 7000, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 98, ".", 7000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 6700, 7400, NA, 97, ".", 7000, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr10", 100, 800, NA, 95, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr10", 100, 900, NA, 80, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr10", 200, 900, NA, 95, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr10", 300, 1000, NA, 75, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr10", 300, 1000, NA, 90, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr10", 300, 1000, NA, 90, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 100, 800, NA, 30, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr2", 100, 900, NA, 10, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 200, 900, NA, 50, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr2", 300, 1000, NA, 50, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr2", 300, 1000, NA, 10, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 300, 1000, NA, 30, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("Chr2", 100, 800, NA, 80, ".", 700, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4 2", 300, 1000, NA, 30, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4 2", 100, 800, NA, 25, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr4 2", 300, 1000, NA, 35, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr4-2", 400, 1100, NA, 20, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4-2", 200, 900, NA, 30, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr4?2", 100, 900, NA, 25, ".", 400, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr4|2", 100, 800, NA, 80, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr42", 300, 1000, NA, 90, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names)
  )
write_tsv(synthetic_data, "data-raw/synthetic_data.bed")

#syn_data_tibble <- synthetic_data |> mutate(summit = center - start +1) |> select(-name, -center) |> arrange(sample_name, chrom, start)
#usethis::use_data(syn_data_tibble, overwrite = TRUE)
###Save as GRanges
#syn_data_granges <- as(syn_data_tibble, "GRanges") |> data.frame() |> select(-name) #|>
#write_tsv("lists/synthetic_data_as_granges.bed")
#usethis::use_data(syn_data_granges)


###Save as BED
#syn_data_bed <- synthetic_data |> select(chrom, start, end, sample_name) #|>
#write_tsv("lists/synthetic_data_as_bed.bed", col_names = FALSE)
#usethis::use_data(syn_data_bed)


#syn_data_control01 <- synthetic_data |> filter(sample_name == "control_rep1") |> select(-name, -sample_name) #|>
#readr::write_tsv("lists/synthetic_data_C1.bed", col_names = FALSE)
#usethis::use_data(syn_data_control01, overwrite = TRUE)


#synthetic_data |> filter(sample_name == "control_rep2") |> select(-name, -sample_name) |>
#  readr::write_tsv("lists/synthetic_data_C2.bed", col_names = FALSE)
#synthetic_data |> filter(sample_name == "control_rep3") |> select(-name, -sample_name) |>
#  readr::write_tsv("lists/synthetic_data_C3.bed", col_names = FALSE)


#syn_data_treatment01 <- synthetic_data |> filter(sample_name == "treatment_rep1") |> select(-name, -sample_name) #|>
#readr::write_tsv("lists/synthetic_data_T1.bed", col_names = FALSE)
#usethis::use_data(syn_data_treatment01, overwrite = TRUE)

#synthetic_data |> filter(sample_name == "treatment_rep2") |> select(-name, -sample_name) |>
#  readr::write_tsv("lists/synthetic_data_T2.bed", col_names = FALSE)
#synthetic_data |> filter(sample_name == "treatment_rep3") |> select(-name, -sample_name) |>
#  readr::write_tsv("lists/synthetic_data_T3.bed", col_names = FALSE)

