conditions_palette <- c(
  `control_mock` = "#A6CEE3",
  `inactivated_control` = "#B2DF8A",
  `infected_control` = "#1F78B4",
  `heparin` = "#33A02C",
  `inactivated_heparin` = "#B15928",
  `infected_heparin` = "#FB9A99",
  `plasma_healthy` = "#E31A1C",
  `inactivated_healthy` = "#FF7F00",
  `infected_healthy` = "#FDBF6F",
  `plasma_covid19` = "#CAB2D6",
  `inactivated_covid19` = "#FFFF99",
  `infected_covid19` = "#6A3D9A"
)

plasma_palette <- c(
  `no_plasma` = '#239B56',
  `plasma_healthy` = '#3498DB',
  `plasma_covid` = '#E74C3C')

heparin_palette <- c(
  `no_heparin` = '#5499C7',
  `heparin` = '#9B59B6')

sars_palette <- c(
  `no_sars` = '#239B56',
  `inactivated_sars` = '#3498DB',
  `infectious_sars` = '#E74C3C')

batch_palette <- c(
  `b1` = '#239B56',
  `b2` = '#3498DB',
  `b3` = '#E74C3C')

# Palette used for each grouping variable in PCA plots.
pca_palettes <- list(
  condition = conditions_palette,
  batch     = batch_palette,
  plasma    = plasma_palette,
  heparin   = heparin_palette,
  sars      = sars_palette)