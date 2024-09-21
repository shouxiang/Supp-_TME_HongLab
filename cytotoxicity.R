
# package loading ---------------------------------------------------------
{
  library(tidyverse)
}

cytotoxic <- read.csv(
  "D:/Publication/TME/sourceData/biocompatibility/alarmaBlue_TME.csv",
  header = TRUE,
  stringsAsFactors = FALSE)

maxValue <- mean(cytotoxic$neg)
minValue <- mean(cytotoxic$dead)

viability <- ((cytotoxic - minValue) / (maxValue - minValue)) %>%
  dplyr::select(-neg, -dead)

s0 <- viability %>%
  pivot_longer(
    colnames(viability),
    names_to = "sample",
    values_to = "intensity") %>%
  separate(sample, 
           into = c("conc", "stain", "time"),
           sep = "_",
           convert = TRUE)

point_30 <- s0 %>%
  filter(time == 30)

point_60 <- s0 %>%
  filter(time == 60)

s1 <- s0 %>%
  group_by(conc, stain, time) %>%
  summarise(
    weight = mean(intensity),
    sd = sd(intensity),
    n = n()
  )

s1$stain <- factor(s1$stain, levels = c("DMSO", "TME"))

