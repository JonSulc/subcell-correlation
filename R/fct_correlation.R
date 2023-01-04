library(ggplot2)
sum_counts_per_bin <- function(spe,
                               nbins_x = 200,
                               nbins_y = 200,
                               normalize = TRUE) {
  spe <- shinySTRegister:::add_bins(spe, nbins_x = nbins_x, nbins_y = nbins_y)

  counts <- data.table::data.table(t(assay(spe))) |>
    cbind(data.table::as.data.table(colData(spe)[, c("binx", "biny", "sample_id")]))
  counts <- counts[
    ,
    lapply(.SD, sum),
    .SDcols = -c("binx", "biny", "sample_id"),
    by = .(binx, biny, sample_id)
  ][
    CJ(biny = seq_len(nbins_y),
       binx = seq_len(nbins_x),
       sample_id = unique(colData(spe)$sample_id)),
    on = .(binx, biny, sample_id),
    .SDcols = -c("binx", "biny", "sample_id")
  ]
  if (!normalize)
    return(counts)

  counts[
    ,
    log1p(.SD) |>
      c(.(binx = binx,
          biny = biny,
          sample_id = sample_id)),
    .SDcols = -c("binx", "biny", "sample_id")
  ]
}

split_counts_by_sample <- function(counts) {
  split(counts[, -c("binx", "biny", "sample_id")],
        counts$sample_id)
}
get_genewise_correlation_from_counts <- function(counts) {
  split_counts_by_sample(counts) |>
    get_genewise_correlation_from_split_counts()
}
get_genewise_correlation_from_split_counts <- function(split_counts) {
  mapply(cor,
         as.data.frame(split_counts[[1]]),
         as.data.frame(split_counts[[2]]),
         use = "pair")
}
get_binwise_correlation_from_counts <- function(counts) {
  split_counts_by_sample(counts) |>
    get_binwise_correlation_from_split_counts()
}

get_binwise_correlation_from_split_counts <- function(split_counts) {
  mapply(cor,
         as.data.frame(split_counts[[1]] |> t()),
         as.data.frame(split_counts[[2]] |> t()),
         use = "pair")
}

get_binwise_correlation_from_spe <- function(
    spe,
    nbins_x = 100,
    nbins_y = nbins_x
) {
  sum_counts_per_bin(spe, nbins_x = nbins_x, nbins_y = nbins_y) |>
    get_binwise_correlation_from_counts()
}

get_shuffled_correlation_from_counts <- function(counts) {
  split_counts <- split_counts_by_sample(counts)
  shuffled_index <- sum(!is.na(split_counts[[2]])) |>
    seq_len() |>
    sample()
  split_counts[[2]] <- `[<-.data.frame`(
    split_counts[[2]],
    !is.na(split_counts[[2]]),
    value = `[.data.frame`(split_counts[[2]], !is.na(split_counts[[2]]))[shuffled_index])

  get_binwise_correlation_from_split_counts(split_counts)
}

plot_counts_correlations <- function(counts, plot_shuffled = TRUE) {
  true_correlation <- get_binwise_correlation_from_counts(counts)
  data <- data.frame(correlation = unname(unlist(true_correlation)),
                     type = "true")
  if (plot_shuffled)
    data <- rbind(data,
                  data.frame(correlation = counts |>
                               get_shuffled_correlation_from_counts() |>
                               unlist() |>
                               unname(),
                             type = "shuffled"))
  ggplot(data, aes(x = correlation, color = type)) +
    geom_density()
}

get_correlation_data_frame <- function(
    sample_list,
    dimension = "spot",
    correlation_function = ifelse(dimension == "spot",
                                  get_binwise_correlation_from_counts,
                                  get_genewise_correlation_from_counts),
    nbins = 200,
    only_common_features = FALSE
  ) {
  if (only_common_features)
    sample_list <- subset_common_features(sample_list)
  Reduce(
    rbind,
    lapply(
      names(sample_list),
      \(sample_name) {
        data.table(
          Feature = rownames(sample_list[[sample_name]]),
          Correlation = sum_counts_per_bin(sample_list[[sample_name]], nbins) |>
            correlation_function(),
          Technology = sample_name
        )
      }
    )
  )
}

subset_common_features <- function(sample_list) {
  common_features <- Reduce(
    intersect,
    lapply(
      sample_list,
      rownames
    )
  )

  lapply(
    sample_list,
    \(sample) sample[common_features, ]
  )
}

plot_counts_correlation_2d <- function(counts,
                                       minimal = FALSE) {
  counts_plot <- data.frame(
    Correlation = get_binwise_correlation_from_counts(counts),
    x = counts$binx[counts$sample_id == counts$sample_id[1]],
    y = counts$biny[counts$sample_id == counts$sample_id[1]]
  ) |>
    ggplot(aes(x = x, y = y, fill = Correlation)) +
    geom_raster() +
    # viridis::scale_fill_viridis(trans = ifelse(log1p, "log1p", "identity")) +
    viridis::scale_fill_viridis() +
    theme_minimal()

  if (minimal)
    counts_plot <- counts_plot +
      theme(axis.text = element_blank(), legend.position = "none")

  counts_plot
}

plot_spe_correlation_2d <- function(spe, nbins = 200, ...) {
  spe |>
    sum_counts_per_bin(nbins, nbins) |>
    plot_counts_correlation_2d()
}
