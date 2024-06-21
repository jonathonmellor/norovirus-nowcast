# script to install relevant dependencies locally
install.packages(setdiff("pacman", rownames(installed.packages())))
# note: some packages have been experimentally commented out to reduce
# the amount of unnecessary dependencies.
# If your code breaks, this may be why
pacman::p_load(
  janitor,
  aws.ec2metadata,
  aws.s3,
  tidyverse,
  patchwork,
  cmdstanr,
  data.table,
  yaml,
  zoo,
  # readxl,
  # kableExtra,
  ggthemes,
  patchwork,
  glue,
  scales,
  # actuar,
  # cplm,
  # Boom,
  bsts,
  # fUnitRoots,
  str2str,
  gt,
  webshot
)

remotes::install_github("epinowcast/epinowcast", dependencies = TRUE, ref = "v0.2.2")
library(epinowcast)
