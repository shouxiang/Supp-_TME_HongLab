
# package loading ---------------------------------------------------------
{
  library(tidyverse)
  rm(list = ls())
}

# Load Data ---------------------------------------------------------------

{
  # BLG kinetics
  blg_kinetics_0 <- read.csv("E:/Publication/TME/sourceData/in_vitro/BLG/blg_kinetics.csv",
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
    pivot_longer(paste0(rep(LETTERS[1:12], each = 4), 1:4),
                 names_to = "sample",
                 values_to = "intensity"
    ) %>%
    separate(sample, 
             into = c("experiment", "rep"),
             sep = 1,
             convert = TRUE) %>%
    mutate(name = case_when(experiment == "A"  ~ "BLG only",
                            experiment == "B" ~ "BLG only + U",
                            experiment == "C" ~ "GSH only",
                            experiment == "D" ~ "GSH only + U",
                            experiment == "E" ~ "TME only",
                            experiment == "F" ~ "TME only + U",
                            experiment == "G" ~ "Folded",
                            experiment == "H" ~ "Unfolded",
                            experiment == "I" ~ "GSH",
                            experiment == "J" ~ "GSH + U",
                            experiment == "K" ~ "Folded + NMM",
                            experiment == "L" ~ "NMM + Unfolded")) %>%
    filter(name %in% c("BLG only", "TME only",
                       "Folded", "Unfolded", "GSH", "NMM + Unfolded")) %>%
    mutate(time = time / 3600) %>%
    select(time, name, rep, intensity)
  
  blg_kinetics <- blg_kinetics_0 %>%
    group_by(name, time) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  blg_kinetics$name <- factor(blg_kinetics$name,
                              levels = c("Unfolded", "Folded", "GSH",
                                         "NMM + Unfolded", "BLG only", "TME only"))
  rm(blg_kinetics_0)
  
  #BLG fluorometer
  
  blg_fluorometer <- read.csv("E:/Publication/TME/sourceData/in_vitro/BLG/blg_fluorometer.csv",
                              header = TRUE,
                              stringsAsFactors = FALSE) %>%
    pivot_longer(LETTERS[1:12],
                 names_to = "experiment",
                 values_to = "intensity"
    ) %>%
    mutate(name = case_when(experiment == "A"  ~ "BLG only",
                            experiment == "B" ~ "BLG only + U",
                            experiment == "C" ~ "GSH only",
                            experiment == "D" ~ "GSH only + U",
                            experiment == "E" ~ "TME only",
                            experiment == "F" ~ "TME only + U",
                            experiment == "G" ~ "Folded",
                            experiment == "H" ~ "Unfolded",
                            experiment == "I" ~ "GSH",
                            experiment == "J" ~ "GSH + U",
                            experiment == "K" ~ "Folded + NMM",
                            experiment == "L" ~ "NMM + Unfolded")) %>%
    filter(name %in% c("BLG only", "TME only",
                       "Folded", "Unfolded", "GSH", "NMM + Unfolded")) %>%
    select(Wavelength, name, intensity)
  
  blg_fluorometer$name <- factor(blg_fluorometer$name,
                                 levels = c("Unfolded", "Folded", "GSH",
                                            "NMM + Unfolded", "BLG only", "TME only"))
  
  #BLG GSH damping
  blg_gsh0 <- read.csv("E:/Publication/TME/sourceData/in_vitro/BLG/blg_gsh.csv",
                       header = TRUE,
                       stringsAsFactors = FALSE) %>%
    pivot_longer(paste0(rep(LETTERS[1:3], each = 4), 1:4),
                 names_to = "sample",
                 values_to = "intensity") %>%
    separate(sample, 
             into = c("experiment", "rep"),
             sep = 1,
             convert = TRUE) %>%
    mutate(name = case_when(experiment == "A"  ~ "BLG only",
                            experiment == "B" ~ "BLG + TME",
                            experiment == "C" ~ "TME only"))
  
  blg_gsh <- blg_gsh0 %>%
    group_by(name, GSH) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  rm(blg_gsh0)
  
  # Dsb kinetics
  
  dsb_kinetics0 <- read.csv("E:/Publication/TME/sourceData/in_vitro/Dsb/dsb_kinetics.csv",
                            header = TRUE,
                            stringsAsFactors = FALSE) %>%
    pivot_longer(paste0(rep(LETTERS[1:6], each = 4), 1:4),
                 names_to = "sample",
                 values_to = "intensity"
    ) %>%
    separate(sample, 
             into = c("experiment", "rep"),
             sep = 1,
             convert = TRUE) %>%
    mutate(name = case_when(experiment == "A" ~ "TME only",
                            experiment == "B" ~ "TME only + U",
                            experiment == "C" ~ "Folded DsbA",
                            experiment == "D" ~ "Folded C33A",
                            experiment == "E" ~ "Unfolded DsbA",
                            experiment == "F" ~ "Unfolded C33A")) %>%
    filter(name %in% c("TME only",
                       "Folded DsbA", "Folded C33A", "Unfolded DsbA")) %>%
    mutate(time = time / 3600)
  
  dsb_kinetics <- dsb_kinetics0 %>%
    group_by(name, time) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  dsb_kinetics$name <- factor(
    dsb_kinetics$name,
    levels = c("Unfolded DsbA",
               "Folded DsbA",
               "Folded C33A",
               "TME only"))
  
  rm(dsb_kinetics0)
  
  #BLG PAGE fluorescence
  
  blg_page0 <- read.csv("E:/Publication/TME/sourceData/in_vitro/BLG/blg_page.csv",
                        header = TRUE,
                        stringsAsFactors = FALSE) %>%
    pivot_longer(paste0(rep(c("F", "U"), each = 4), 1:4),
                 names_to = "sample",
                 values_to = "intensity"
    ) %>%
    separate(sample, 
             into = c("experiment", "rep"),
             sep = 1,
             convert = TRUE)
  
  blg_page <- blg_page0 %>%
    group_by(experiment, time) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  rm(blg_page0)
  
  # Dsb PAGE fluorescence
  
  dsb_page0 <- read.csv("E:/Publication/TME/sourceData/in_vitro/Dsb/dsb_page.csv",
                        header = TRUE,
                        stringsAsFactors = FALSE) %>%
    pivot_longer(paste0(rep(LETTERS[1:4], each = 4), 1:4),
                 names_to = "sample",
                 values_to = "intensity") %>%
    separate(sample, 
             into = c("experiment", "rep"),
             sep = 1,
             convert = TRUE) %>%
    mutate(name = case_when(
      experiment == "A" ~ "Unfolded C33A",
      experiment == "B" ~ "Folded C33A",
      experiment == "C" ~ "Unfolded DsbA",
      experiment == "D" ~ "Folded DsbA"))
  
  dsb_page <- dsb_page0 %>%
    group_by(name, time) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  rm(dsb_page0)
}

# fitG --------------------------------------------------------------------
{
  {
    blg_deltaG0 <- read.csv(
      "E:/Publication/TME/sourceData/in_vitro/BLG/blg_deltaG.csv",
      header = TRUE,
      stringsAsFactors = FALSE) %>%
      pivot_longer(c("trp_1","trp_2","trp_3",
                     "tme_1","tme_2","tme_3"),
                   names_to = "stain",
                   values_to = "intensity") %>%
      separate(stain, 
               into = c("grp", "rep"),
               sep = "_",
               convert = TRUE)
    
    blg_deltaG0$grp <- factor(blg_deltaG0$grp,
                              levels = c("trp", "tme"))
    
    fitG <- function(df, tLower = 3.5, tUpper = 6.75){
      
      left <- df %>%
        filter(urea >= 0 & urea <= tLower)
      
      middle <- df %>%
        filter(urea > tLower & urea < tUpper)
      
      right <- df %>%
        filter(urea >= tUpper & urea <= 9)
      
      leftMod <- lm(intensity ~ urea, data = left)
      left <- left %>%
        modelr::add_predictions(leftMod)
      k_f <- as.double(left[left$urea == tLower,]$pred)
      
      rightMod <- lm(intensity ~ urea, data = right)
      right <- right %>%
        modelr::add_predictions(rightMod)
      k_u <- as.double(right[right$urea == tUpper,]$pred)
      
      R <- 1.987
      T <- 298
      
      middle <- middle %>%
        filter(intensity >= k_f & intensity <= k_u)
      
      middle <- middle %>%
        mutate (K = (intensity - k_f)/(k_u - intensity)) %>%
        mutate(appG = -1 * R * T * log( K, exp(1)))
      
      lmod <- lm(appG ~ urea, data = middle)
      corr <- cor(middle$urea, middle$appG) %>% abs()
      
      middle <- middle %>%
        modelr::add_predictions(lmod) %>%
        mutate(delG = lmod$coefficients[1]/1000,
               corr = corr)
    }
    
    blg_deltaG <- blg_deltaG0 %>%
      split(list(.$rep, .$grp)) %>%
      map(fitG, 3.5, 6.75) %>%
      reduce(rbind)
    
    rm("fitG")
    
    plotG <- function(df){
      
      stopifnot("Input data does not contain correct column names" =
                  c("urea", "appG", "pred", "corr") %in%
                  colnames(df))
      
      ggplot(df) +
        geom_point(aes(urea, appG/1000),
                   size = .4) +
        geom_line(aes(urea, pred/1000),
                  colour = "red",
                  size = .5,
                  alpha = .5) +
        scale_y_continuous(breaks = seq(-2, 2, 1),
                           labels = seq(-2, 2, 1),
                           limits = c(-2.3, 2.5)) +
        geom_text(
          label = paste0("| r | = ", format(unique(df$corr), 
                                            digits = 2,
                                            nsmall = 2)),
          x = 4.7,
          y = -0.5,
          size = 6/.pt
        ) +
        xlab("[Urea] (M)") +
        ylab(expression(paste("Apparent ", Delta, "G (Cal/mol)"))) +
        theme0
    }
    
    fitTrp <- list()
    for (i in 1 : 3) {
      
      df <- blg_deltaG %>%
        filter(grp == "trp" & rep == i)
      
      name <- paste0("p", i)
      fitTrp[[name]] <- plotG(df) + 
        ggtitle(paste0("Trp - Replicate " , i))
    }
    
    fitTME <- list()
    for (i in 1 : 3) {
      
      df <- blg_deltaG %>%
        filter(grp == "tme" & rep == i)
      
      name <- paste0("p", i)
      fitTME[[name]] <- plotG(df) + 
        ggtitle(paste0("TME - Replicate " , i))
    }
  }
  
  # deltaG deviation ------------------------------------------------------
  
  blg_deltaG_1 <- 
    blg_deltaG %>%
    select(grp, rep, delG) %>%
    distinct()
  
  blg_deltaG_2 <-
    blg_deltaG_1%>%
    group_by(grp) %>%
    summarise(
      weight = mean(delG),
      sd = sd(delG),
      n = n()
    )
  
  blg_deltaG_dev <- 
    ggplot(blg_deltaG_1,
           aes(grp, delG)) +
    geom_jitter(
      width = 0.05,
      size = .5) +
    geom_errorbar(
      data = blg_deltaG_2,
      aes(grp, weight,
          ymin = weight - sd,
          ymax = weight + sd),
      colour = "black",
      width = .15,
      size = .2) +
    scale_x_discrete(
      labels=c("trp" = "Trp", "tme" = "TME")) +
    scale_y_continuous(
      breaks = c(7:9),
      labels = c(7:9),
      limits = c(6.5, 9),
      guide = "prism_offset"
    ) +
    ylab(expression(paste(Delta, "G (Cal/mol)"))) +
    theme0 +
    theme(axis.title.x = element_blank())
  
  # deltaG sensitivity ---------------------------------------------------------
  
  blg_sensitivity0 <-
    blg_deltaG0 %>%
    group_by(grp, rep) %>%
    summarise(
      sensitivity = max(intensity)/min(intensity)
    )
  
  blg_sensitivity_1 <-
    blg_sensitivity0%>%
    group_by(grp) %>%
    summarise(
      weight = mean(sensitivity),
      sd = sd(sensitivity),
      n = n()
    )
  
  blg_deltaG_sensitivity <- 
    ggplot(blg_sensitivity0,
           aes(grp, sensitivity)) +
    geom_jitter(
      width = 0.15,
      size = .5) +
    geom_errorbar(
      data = blg_sensitivity_1,
      aes(grp, weight,
          ymin = weight - sd,
          ymax = weight + sd),
      colour = "black",
      width = .15,
      size = .2) +
    scale_x_discrete(
      labels=c("trp" = "Trp", "tme" = "TME")) +
    scale_y_continuous(
      breaks = seq(0,8,2),
      labels = seq(0,8,2),
      limits = c(0, 8),
      guide = "prism_offset"
    ) +
    ylab("Max / Min") +
    theme0 +
    theme(axis.title.x = element_blank())
  
  # Denaturation curve ------------------------------------------------------
  
  s <- blg_deltaG0 %>%
    group_by(grp, urea) %>%
    summarise(
      weight = mean(intensity),
      sd = sd(intensity),
      n = n()
    )
  
  tLower = 3.5
  tUpper = 6.75
  
  fitTail <- function(df, tLower = 3.5, tUpper = 6.75){
    
    left <- df %>%
      filter(urea >= 0 & urea <= tLower)
    
    middle <- df %>%
      filter(urea > tLower & urea < tUpper)
    
    right <- df %>%
      filter(urea >= tUpper & urea <= 9)
    
    leftMod <- lm(weight ~ urea, data = left)
    left <- left %>%
      modelr::add_predictions(leftMod)
    
    rightMod <- lm(weight ~ urea, data = right)
    right <- right %>%
      modelr::add_predictions(rightMod)
    
    rbind(left, middle) %>%
      rbind(right)
  }
  
  s1 <- rbind(fitTail(s[s$grp == "trp",]), fitTail(s[s$grp == "tme",]))
}

