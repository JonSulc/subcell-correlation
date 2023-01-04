common_features <- Reduce(
  intersect,
  sample_list |>
    lapply(rownames)
)

do.call(
  rbind,
  lapply(
    names(sample_list),
    \(technology) {
      expression <- SummarizedExperiment::assay(sample_list[[technology]][common_features,]) |>
        t() |>
        data.table::as.data.table()
      expression$Technology <- technology
      expression
    }
  )
) |>
  saveRDS("data/all_expression_data.rds")
