# Load library
library("tidyverse")
library("GenomicRanges")

# Define column names
names <- c("chrom", "start", "end", "name", "score", "strand", "center", "sample_name")

# Create the entire synthetic data
synthetic_data <- tibble("chr1", 200, 900, NA, 100, ".", 500, "treatment_rep1") |>
  rename_all(.fun = ~names) |>
  rbind(
    tibble("chr1", 1, 900, NA, 97, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 101, 300, NA, 94, ".", 200, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 301, 900, NA, 94, ".", 500, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 201, 900, NA, 100, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 301, 900, NA, 98, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 301, 1000, NA, 96, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 301, 1100, NA, 93, ".", 500, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 1301, 1600, NA, 97, ".", 1400, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 1901, 2200, NA, 98, ".", 2000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 2501, 3100, NA, 97, ".", 2800, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 2501, 3400, NA, 98, ".", 3000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 2601, 3200, NA, 99, ".", 2800, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 3501, 4200, NA, 44, ".", 3800, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 3501, 4400, NA, 95, ".", 3800, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 3601, 4400, NA, 43, ".", 3900, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4501, 5000, NA, 97, ".", 4800, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4501, 5200, NA, 60, ".", 4700, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4501, 5200, NA, 59, ".", 5000, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4501, 5300, NA, 98, ".", 4800, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4501, 5300, NA, 98, ".", 5100, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4601, 5100, NA, 93, ".", 4900, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 4601, 5200, NA, 94, ".", 4800, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 4701, 5300, NA, 46, ".", 4900, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 4701, 5300, NA, 45, ".", 5100, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 5601, 6100, NA, 26, ".", 5700, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 5701, 6400, NA, 98, ".", 6200, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 5801, 6300, NA, 30, ".", 6100, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 25, ".", 7000, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 44, ".", 7000, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 43, ".", 7000, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 29, ".", 7000, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 98, ".", 7000, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr1", 6701, 7400, NA, 97, ".", 7000, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr10", 101, 800, NA, 95, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr10", 101, 900, NA, 80, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr10", 201, 900, NA, 95, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr10", 301, 1000, NA, 75, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr10", 301, 1000, NA, 90, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("chr10", 301, 1000, NA, 90, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 101, 800, NA, 30, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr2", 101, 900, NA, 10, ".", 500, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 201, 900, NA, 50, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr2", 301, 1000, NA, 50, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr2", 301, 1000, NA, 10, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr2", 301, 1000, NA, 30, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names),
    tibble("Chr2", 101, 800, NA, 80, ".", 700, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4 2", 301, 1000, NA, 30, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4 2", 101, 800, NA, 25, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr4 2", 301, 1000, NA, 35, ".", 600, "control_rep3") |> rename_all(.fun = ~names),
    tibble("chr4-2", 401, 1100, NA, 20, ".", 600, "control_rep1") |> rename_all(.fun = ~names),
    tibble("chr4-2", 201, 900, NA, 30, ".", 500, "treatment_rep1") |> rename_all(.fun = ~names),
    tibble("chr4?2", 101, 900, NA, 25, ".", 400, "treatment_rep3") |> rename_all(.fun = ~names),
    tibble("chr4|2", 101, 800, NA, 80, ".", 400, "control_rep2") |> rename_all(.fun = ~names),
    tibble("chr42", 301, 1000, NA, 90, ".", 600, "treatment_rep2") |> rename_all(.fun = ~names)
  )
write_tsv(synthetic_data, "data-raw/synthetic_data.bed")
