## GRAPHICAL ABSTRACT

## Install rayshader (if needed)
if (!requireNamespace("rayshader", quietly = TRUE))
  install.packages("rayshader")

## Read XRD data
## (see analysis/supplementary-materials/01-data.R)
xrd_counts <- here::here("analysis/data/derived_data/xrd_clean.csv") |>
  read.table(sep = ",", dec = ".", header = TRUE)

## Remove the first column (2 theta positions)
theta_scale  <- xrd_counts[, 1]
xrd_counts <- xrd_counts[, -1]

## Transpose (one sample per line)
xrd_counts <- xrd_counts |> t() |> as.data.frame()
colnames(xrd_counts) <- paste0("T", seq_len(ncol(xrd_counts)))

## Read archaeological data
corpus <- here::here("analysis/data/raw_data/corpus.csv") |>
  read.table(sep = ",", dec = ".", header = TRUE, encoding = "UTF-8")

## Mas de Moreno: ceramics + unfired sherds
# all(corpus$sample == rownames(xrd_counts))
index_moreno <- corpus$site == "Mas de Moreno"

## Correspondance analysis
ca_moreno <- dimensio::ca(xrd_counts[index_moreno, ])

## Get CA order
ca_coords <- dimensio::get_coordinates(ca_moreno)
ca_index <- order(ca_coords[, 1], decreasing = TRUE)

## Reorder matrix
xrd_counts2 <- xrd_counts[ca_index, ]
colnames(xrd_counts2) <- as.character(theta_scale)
xrd_counts2$sample <- factor(rownames(xrd_counts2), levels = rownames(xrd_counts2))

gg_ray <- xrd_counts2 |>
  tidyr::gather(key = "theta", value = "value", -sample) |>
  dplyr::filter(theta >= 25 & theta <= 32) |>
  dplyr::mutate(
    theta = as.numeric(theta),
    sample2 = as.numeric(sample),
    # value = sqrt(value),
    value = value / 1000
  ) |>
  ggplot2::ggplot() +
  ggplot2::geom_tile(ggplot2::aes(x = theta, y = sample2, fill = value)) +
  ggplot2::scale_x_continuous(breaks = seq(25, 32, 1), expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 70, 5), expand = c(0, 0)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = expression(2*theta),
    fill = expression("count ("*10^3*")")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colourbar(
      title.vjust = 0.8,
      frame.colour = "black"
    )
  ) +
  khroma::scale_fill_YlOrBr()
gg_ray

rayshader::plot_gg(
  gg_ray +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "bottom")),
  scale = 300,
  shadow_intensity = 0.7,
  sunangle = 135,
  theta = 45,
  phi = 45,
  zoom = 0.9,
  windowsize = 600
)

here::here("analysis/figures/abstract.png") |>
  rayshader::render_snapshot(clear = TRUE)
