# Load the entire synthetic data
synthetic_data <- read_tsv("data-raw/synthetic_data.bed", show_col_types = FALSE)

# Filter synthetic data for treatment rep 3
synthetic_data |>
  filter(sample_name == "treatment_rep3") |>
  select(-sample_name) |>
  mutate(
    name = ".",
    signalValue = score,
    pValue = -1,
    qValue = c(log10(score * score)),
    peak = center - start,
    score = 0
  ) |>
  select(-center) |>
  write_tsv("data-raw/syn_treatment_rep3.narrowPeak", col_names = FALSE)
