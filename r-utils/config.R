# R/config.R

# set the path to the project.
project_path <- '~/3d4covid-bulkRNAseq'

# aliner name.
mapper <- 'STAR'

# Canonical sample IDs and their intended order.
samples_id <- c(
  '1_CM', '2_CS', '3_CH', '4_hCM', '5_hCS',
  '6_21', '7_21S', '8_21H', '9_14', '10_14S',
  '11_14H', '12_28', '13_28S', '14_28H', '15_48',
  '16_48S', '17_48H', '18_CM', '19_CS', '20_CH',
  '21_hCM', '22_hCS', '23_hCH', '24_08', '25_08S',
  '26_08H', '27_10', '28_10S', '29_10H', '30_76',
  '31_76S', '32_76H', '33_27', '34_27S', '35_27H',
  '36_CM', '37_CS', '38_CH', '39_hCM', '40_hCS',
  '41_hCH', '42_29', '43_29S', '44_29H', '45_11',
  '46_11S', '47_11H', '48_22', '49_22S', '50_22H',
  '51_58', '52_58S', '53_58H')
