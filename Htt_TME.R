
# Package loading ---------------------------------------------------------
{
  library(ggprism)
  library(conicfit)
  library(sf)
  library(tidyverse)
  rm(list = ls())
}

# Color setting -----------------------------------------------------------

{
  colour_3stage <- c(
    stage1 = "#FF9595", 
    stage2 = "#6AA3FF", 
    stage3 = "#FFA3F8")
  
  color_mutant <- c(
    mCherry = "#BDC3CB",
    `25Q` = "#bcff14",
    `46Q` = "#ade6f4",
    `97Q` = "#cb99c9"
  )
  
  color_H <- c(
    head = "#806da5", 
    tail = "#f4a016")
  
  color_W <- c(
    head = "#E9738D", 
    tail = "#71AACD")
  
}

# Gating parameters -------------------------------------------------------

{
  incluPoints <- rbind(
    c(50 * 1e+03, 200),
    c(80 * 1e+03, 1000),
    c(125 * 1e+03, 100 * 1e+03),
    c(50 * 1e+03, 100 * 1e+03),
    c(50 * 1e+03, 200)) %>% as.data.frame()
  
  arrow_1stage <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 35000, 
                        y1 = 10^2.5, 
                        x2 = 2*50000, 
                        y2 = 10^4.8),
      arrow = arrow(length = unit(2, "mm"), angle = 20),
      size = .5,
      curvature = .2,
      colour = colour_3stage[["stage1"]]
    )
  
  arrow_3stage1.1 <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 35000,
                        y1 = 10^2, 
                        x2 = 1.3*50000, 
                        y2 = 10^2.5),
      arrow = arrow(length = unit(1.5, "mm"), angle = 15),
      size = .5,
      curvature = .2,
      colour = colour_3stage[["stage1"]]
    )
  
  arrow_3stage1.2 <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 2.2*50000, 
                        y1 = 10, 
                        x2 = 3.1*50000, 
                        y2 = 10^2),
      arrow = arrow(length = unit(1.5, "mm"), angle = 15),
      size = .5,
      curvature = .2,
      colour = colour_3stage[["stage1"]]
    )
  
  arrow_3stage2.1 <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 1.3*50000,
                        y1 = 10^2.6, 
                        x2 = 1.15*50000, 
                        y2 = 10^3.5),
      arrow = arrow(length = unit(1.5, "mm"), angle = 15),
      size = .5,
      curvature = .5,
      colour = colour_3stage[["stage2"]]
    )
  
  arrow_3stage2.2 <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 3.1*50000,
                        y1 = 10^2.2, 
                        x2 = 2.7*50000, 
                        y2 = 10^3.5),
      arrow = arrow(length = unit(1.5, "mm"), angle = 15),
      size = .5,
      curvature = .5,
      colour = colour_3stage[["stage2"]]
    )
  
  arrow_3stage3.1 <-
    geom_curve(
      aes(x = x1, y = y1, xend = x2, yend = y2),
      data = data.frame(x1 = 1.15*50000,
                        y1 = 10^3.6, 
                        x2 = 2*50000, 
                        y2 = 10^5.25),
      arrow = arrow(length = unit(1.5, "mm"), angle = 15),
      size = .5,
      curvature = -.5,
      colour = colour_3stage[["stage3"]]
    )
  
  gate_46q_1 <-
    geom_path(
      aes(x = x, y = y),
      data = data.frame(
        x = 50000*c(0.7, 1.5, 3.4, 3.3, 1.3, 0.7),
        y = 10 ^  c(1.9, 1.5, 1.3, 1.8, 2.1, 1.9)
      ),
      alpha = .8,
      size = .2
    )
  
  gate_46q_2b <-
    geom_path(
      aes(x = x, y = y),
      data = data.frame(
        x = 50000*c(1.3, 1.9, 2.2, 2.0, 1.3),
        y = 10 ^  c(2.2, 2.1, 3.3, 3.7, 2.2)
      ),
      alpha = .8,
      size = .2
    )
  
  gate_46q_2a <-
    geom_path(
      aes(x = x, y = y),
      data = data.frame(
        x = 50000*c(2.0, 2.8, 2.6, 2.3, 2.0),
        y = 10 ^  c(2.1, 2.0, 4.2, 3.5, 2.1)
      ),
      alpha = .8,
      size = .2
    )
  
  gate_46q_3 <-
    geom_path(
      aes(x = V1, y = V2),
      data = incluPoints,
      alpha = .8,
      size = .2
    )
  
  # text_46q <-
  #   draw_text(
  #     text = c("\U2160",paste0("\U2161","a"),paste0("\U2161","b"),"\U2162"),
  #     x = 50000*c(3.8, 3.0, 0.7, 0.7),
  #     y = 10 ^  c(2.0, 3.0, 2.2, 4.5),
  #     size = 6
  #   )
  
  gate_97q_2a <-
    geom_path(
      aes(x = x, y = y),
      data = data.frame(
        x = 50000*c(2.0, 2.8, 2.6, 2.3, 2.0),
        y = 10 ^  c(2.1, 2.0, 3.5, 3.5, 2.1)
      ),
      alpha = .8,
      size = .2
    )

}

# helper ------------------------------------------------------------------

{
  # bin a parameter p within group, number = N, name = pName
  binP <- function(df, p, group, N){
    
    df %>%
      split(.[[group]]) %>%
      map(function(df) df %>%
            mutate(!!quo_name(paste(p, ".l", sep = "")) := cut_number(
              x = .[[p]], n = N, labels = c(1:N)))) %>%
      reduce(rbind)
  }
  
  genBreaks <- function(df) {
    
    summary <- df %>%
      mutate(x = cut_number(
        mCherry.A, n = N, labels = c(1:N))) %>%
      group_by(x) %>% 
      summarise(
        `#` = n(),
        min = min(mCherry.A),
        max = max(mCherry.A)
      )
    
    c(min(df$mCherry.A) - 1,
      summary$max[1: N-1],
      summary$max[N] +1)
  }
  
  
  tTest <- function(d0, d1){ # Statistical analysis for two inputs
    
    stopifnot("dataframe does not contain correct column names" =
                c("mutant",
                  "mCherry.A.l",
                  "medianV450") %in%
                colnames(d0))
    
    stopifnot("dataframe does not contain correct column names" =
                c("mutant",
                  "mCherry.A.l",
                  "medianV450") %in%
                colnames(d1))
    
    s0 <- d0 %>%
      group_by(mutant, mCherry.A.l) %>%
      summarise(
        weight = mean(medianV450),
        sd = sd(medianV450),
        n = n(),
        se = sd / sqrt(n)
      ) %>%
      mutate(p = c(1.0))
    
    s1 <- d1 %>%
      group_by(mutant, mCherry.A.l) %>%
      summarise(
        weight = mean(medianV450),
        sd = sd(medianV450),
        n = n(),
        se = sd / sqrt(n)
      ) %>%
      mutate(p = c(1.0))
    
    v <- var.test(d0$medianV450, d1$medianV450)
    
    if(v$p.value <= .05) {
      p <- t.test(d0$medianV450, d1$medianV450, var.equal = FALSE)$p.value
    } else {
      p <- t.test(d0$medianV450, d1$medianV450, var.equal = TRUE)$p.value
    }
    
    if(s0$weight > s1$weight){
      s0$p <- p
    }else{
      s1$p <- p
    }
    
    rbind(s0, s1) %>%
      mutate(sig = cut(p, 
                       breaks = c(0, .0001, .001, .01, .05, Inf),
                       labels = c("****",
                                  "***",
                                  "**",
                                  "*",
                                  ""),
                       right = FALSE)
      )
  }
  
  plot_expression <- function(df){
    
    stopifnot("dataframe does not contain correct column names" =
                c("mCherry.W",
                  "mCherry.H",
                  "mCherry.A") %in%
                colnames(df))
    
    ggplot(df, aes(mCherry.W, mCherry.H)) +
      geom_point(aes(colour = mCherry.A), alpha = 1/15, size = .5)  +
      scale_x_continuous(
        name = "mCherry.W x 50000 (a.u.)",
        limits = c(30000, 200000),
        breaks = 50000 * c(1:4),
        labels =  c(1:4)) +
      scale_y_continuous(
        name = expression(log[10]*"(mCherry.H) (a.u.)"),
        trans = "log10",
        limits = c(-Inf, 10^5.3),
        breaks = c(10^c(2:5)),
        labels =  c(2:5)) + 
      scale_colour_gradientn(
        trans = "log10", colours=rainbow(10),
        name = expression(log[10]*"(mCherry.A) (a.u.)"),
        limits = c(-Inf, 10^5),
        breaks = c(10^c(2:5)),
        labels =  c(2:5),
        guide = guide_colourbar(title.position = "left"),) +
      geom_path(data = incluPoints,
                aes(V1, V2),
                alpha = .8,
                size = .2)
  }
  
  plot_tme <- function(df){
    
    stopifnot("dataframe does not contain correct column names" =
                c("mCherry.W",
                  "mCherry.H",
                  "V450.A") %in%
                colnames(df))
    
    ggplot(df, aes(mCherry.W, mCherry.H)) +
      geom_point(aes(color = V450.A), alpha = 1/20, size = .5) +
      scale_x_continuous(
        name = "Pulse width x 50000 (a.u.)",
        limits = c(30000, 200000),
        breaks = 50000 * c(1:4),
        labels =  c(1:4)) +
      scale_y_continuous(
        name = expression(log[10]*"(Pulse height) (a.u.)"),
        trans = "log10",
        breaks = c(10^c(2:5)),
        labels =  c(2:5)) + 
      scale_colour_gradient2(
        name = "TME signal x 5000 (a.u.)",
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 19000,
        limits = c(0, 35000),
        breaks = c(5000*c(0:7)),
        labels =  c(0:7),
        guide = guide_colourbar(title.position = "left"),) +
      geom_path(data = incluPoints,
                aes(V1, V2),
                color = "black",
                alpha = .8,
                size = .2)
  }
  
  plot_x_tme <- function(df) {
    
    ggplot() +
      geom_rect(
        aes(xmin = 40, xmax = 10^2.4,
            ymin = 35000, ymax = 40000),
        fill = colour_3stage[["stage1"]],
        alpha = .7) +
      geom_rect(
        aes(xmin = 10^2.4, xmax = 10^4,
            ymin = 35000, ymax = 40000),
        fill = colour_3stage[["stage2"]],
        alpha = .7) +
      geom_rect(
        aes(xmin = 10^4, xmax = 10^4.8,
            ymin = 35000, ymax = 40000),
        fill = colour_3stage[["stage3"]],
        alpha = .7) +
      geom_point(
        data = df, 
        aes(medianmCherry, weight, colour = mutant),
        alpha = .5,
        size = .2) +
      geom_errorbar(
        data = df,
        aes(medianmCherry, weight,
            ymin = weight - se,
            ymax = weight + se), 
        width =.08,
        alpha = .5,
        size = .25,
        colour = "black") +
      
      scale_x_continuous(trans = "log10",
                         name = expression(log[10]*"(mCherry) (a.u.)"),
                         breaks = 10^c(2:4),
                         labels = c(2:4)) +
      scale_y_continuous(name = "TME signal x 10000 (a.u.)",
                         limits = c(0, 40000),
                         breaks = 10000*c(1:3),
                         labels = c(1:3)) +
      geom_line(
        data = df,
        aes(medianmCherry, weight, colour = mutant),
        alpha = 0.8,
        show.legend = FALSE,
        size = .5) +
      geom_text(
        data = df,
        aes(x = medianmCherry,
            y = weight + se + 4000,
            label = sig),
        colour = "black",  size = 2) +
      scale_color_manual(
        name = NULL,
        labels = c("mCherry", "25Q",
                   "46Q", "97Q"),
        values = c(
          mcherry = color_mutant[["mCherry"]], 
          `25q` = color_mutant[["25Q"]], 
          `46q` = color_mutant[["46Q"]], 
          `97q` = color_mutant[["97Q"]])) +
      theme(
        strip.text = element_blank(),
        legend.title = element_blank()
      )
  }
  
  sepUL <- function(df) {
    df <- df %>%
      mutate(pos = c("middle"))
    
    df <- arrange(df, mCherry.H)
    df[1: (nrow(df)/4), ]$pos <- "tail"
    
    df <- arrange(df, desc(mCherry.H))
    df[1: (nrow(df)/4), ]$pos <- "head"
    df
  }
  
  shape <- function(df){# p value of Upper and Lower
    
    df$mCherry.A.l <- as.integer(df$mCherry.A.l)
    df$mutant <- as.character(df$mutant)
    
    s <- df %>%
      group_by(mutant, mCherry.A.l, pos) %>%
      summarise(
        weight = mean(medianV450),
        sd = sd(medianV450),
        n = n(),
        se = sd / sqrt(n)
      ) %>%
      mutate(p = c(1.0), sig = c(" "))
    
    d1 <- df[df$pos == "head", ]$medianV450
    d2 <- df[df$pos == "tail", ]$medianV450
    
    v <- var.test(d1, d2)
    
    if(v$p.value <= .05) {
      p <- t.test(d1, d2, var.equal = FALSE)$p.value
    } else {
      p <- t.test(d1, d2, var.equal = TRUE)$p.value
    }
    
    if(p >= .05) {
      sig <- " "
    }else if(p >= .01) {
      sig <- "*"
    }else if(p >= .001) {
      sig <- "**"
    }else if(p >= .0001) {
      sig <- "***"
    }else{
      sig <- "****"
    }
    if(s[s$pos == "head",][["weight"]] > s[s$pos == "tail",][["weight"]]){
      s[s$pos == "head",][["p"]] <- p
      s[s$pos == "head",][["sig"]] <- sig
    }else{
      s[s$pos == "tail",][["p"]] <- p
      s[s$pos == "tail",][["sig"]] <- sig
    }
    s
  }
  
  plot_shape <- function(df) {
    
    stopifnot("dataframe does not contain correct column names" =
                c("medianmCherry",
                  "weight",
                  "se",
                  "pos",
                  "sig") %in%
                colnames(df))
    
    ggplot() +
      geom_rect(
        aes(xmin = 40, xmax = 10^2.4,
            ymin = 47000, ymax = 52000),
        fill = colour_3stage[["stage1"]],
        alpha = .7) +
      geom_rect(
        aes(xmin = 10^2.4, xmax = 10^4,
            ymin = 47000, ymax = 52000),
        fill = colour_3stage[["stage2"]],
        alpha = .7) +
      geom_rect(
        aes(xmin = 10^4, xmax = 10^5,
            ymin = 47000, ymax = 52000),
        fill = colour_3stage[["stage3"]],
        alpha = .7) +
      geom_point(
        data = df,
        aes(medianmCherry, weight),
        alpha = .5, size = .2) +
      geom_errorbar(
        data = df,
        aes(medianmCherry, weight,
            ymin = weight - se,
            ymax = weight + se), 
        width =.08,
        alpha = .5,
        size = .25,
        colour = "black") +
      scale_x_continuous(trans = "log10",
                         name = expression(log[10]*"(mCherry) (a.u.)"),
                         breaks = 10^c(2:5),
                         labels = c(2:5)) +
      scale_y_continuous(name = "TME signal x 10000 (a.u.)",
                         breaks = 10000*c(1:5),
                         labels = c(1:5)) +
      geom_line(
        data = df,
        aes(medianmCherry, weight, colour = pos),
        size = .5,
        alpha = 0.8) +
      geom_text(
        data = df,
        aes(medianmCherry,
            y = weight + se + 5000,
            label = sig),
        colour = "black",
        size = 2,
        alpha = .6)
  }
  
  sepLR <- function(df) {
    df <- df %>%
      mutate(pos = c("middle"))
    
    df <- arrange(df, mCherry.W)
    df[1: (nrow(df)/4), ]$pos <- "tail"
    
    df <- arrange(df, desc(mCherry.W))
    df[1: (nrow(df)/4), ]$pos <- "head"
    df
  }
  
  processFlowCsv <- function(df) {
    
    isInGate <- function(df, gatePoints, x, y) {
      
      gatePoints[nrow(gatePoints), ] <- gatePoints[1, ]
      
      gateSfc <- gatePoints %>% 
        list() %>% 
        st_polygon() %>% 
        st_sfc()
      
      xySfc <- cbind(df[[x]], df[[y]]) %>% 
        st_multipoint() %>%
        st_sfc() %>%
        st_cast("POINT")
      
      st_intersects(xySfc, gateSfc, sparse =  FALSE)[, 1]
    }
    
    #p1 the main population
    p1Points <-calculateEllipse(x = 140 * 1e+03, y = 110 * 1e+03,
                                a = 130 * 1e+03, b= 90 * 1e+03,
                                angle = 60, steps = 50)
    p1 <- isInGate(df, p1Points,
                   "FSC.A", "SSC.A")
    
    #Singlet from FSC
    sFPoints <- rbind(c(75 * 1e+3, 0),
                      c(75 * 1e+03, 100 * 1e+03),
                      c(225 * 1e+03, 100 * 1e+03),
                      c(225 * 1e+03, 0 * 1e+03),
                      c(75 * 1e+03, 0 * 1e+03)) 
    sF <- isInGate(df, sFPoints, 
                   "FSC.W", "FSC.H")
    
    #Singlet from SSC
    sSPoints <- rbind(c(75 * 1e+3, 0),
                      c(75 * 1e+03, 100 * 1e+03),
                      c(200 * 1e+03, 100 * 1e+03),
                      c(200 * 1e+03, 0 * 1e+03),
                      c(75 * 1e+03, 0 * 1e+03))
    sS <- isInGate(df, sSPoints,
                   "SSC.W", "SSC.H")
    
    #highV450
    stainPoints <- rbind(c(3, 10^3.5),
                         c(3, 10^5),
                         c(10^5.5, 10^5),
                         c(10^5.5, 10^3.5),
                         c(3, 10^3.5))
    stain <- isInGate(df, stainPoints,
                      "mCherry.A", "V450.A")
    
    df %>% mutate(p1 = p1,
                  sF = sF,
                  sS = sS,
                  stain = stain)
  }
  
  readFlowCsv <- function(
    fileName,
    inputDir = "D:/Publication/TME/sourceData/flowCytometry/Htt/csv/singleStain/") {
    
    file <- read.csv(paste0(inputDir, fileName),
                     header = TRUE, stringsAsFactors = FALSE,
                     nrows = 5000)
    file %>%
      mutate("mutant" = str_remove_all(fileName, "(_001.exported.FCS3.csv)|(Specimen_001_)"))
  }
  
  plot_singleStain <- function(df, percentText){
    
    ggplot() +
      geom_path(
        data = stainPoints,
        aes(V1, V2),
        alpha = .8,
        size = .2) +
      geom_rect(
        aes(xmin = 40, xmax = 10^2.4,
            ymin = 10^1, ymax = 10^5.2),
        fill = "red",
        alpha = .3) +
      geom_rect(
        aes(xmin = 10^2.4, xmax = 10^4,
            ymin = 10^1, ymax = 10^5.2),
        fill = "blue",
        alpha = .3) +
      geom_rect(
        aes(xmin = 10^4, xmax = 10^5.7,
            ymin = 10^1, ymax = 10^5.2),
        fill = "purple",
        alpha = .3) +
      geom_point(
        data = df,
        aes(mCherry.A, V450.A),
        alpha = 1/10, 
        size = .5) +
      geom_text(
        aes(x = 2*10^3, y = 2*10^4, label = percentText),
        size = 7/.pt) +
      scale_x_continuous(
        trans = "log10",
        limits = c(2, 10^5.8),
        breaks = 10^c(1:5),
        labels =  c(1:5)) +
      scale_y_continuous(
        trans = "log10",
        breaks = 10^c(1:5),
        labels =  c(1:5)) +
      xlab(expression(log[10]*"(mCherry) (a.u.)")) + 
      ylab(expression(log[10]*"(TME) (a.u.)"))
  }
}

# Load data ---------------------------------------------------------------
{
  a0 <- read.csv(
    "D:/Publication/TME/sourceData/flowCytometry/Htt/Htt_TME_gated.csv",
    header = TRUE,
    stringsAsFactors = FALSE) %>%
    filter(mutant %in% c("mcherry", "25q", "46q", "97q")) %>%
    select(mCherry.A, mCherry.W, mCherry.H,
           inclu,
           mutant, 
           V450.A, V500.A,
           replicate)
  
  a0$mutant <- factor(a0$mutant,
                      levels = c("mcherry", "25q", "46q", "97q")
  )
  
  #3 stages
  a_3stage <- a0 %>%
    filter(replicate == 1) %>%
    split(.$mutant) %>%
    map(function(df) df %>%
          head(5000)) %>%
    reduce(rbind)
  
  #4 phases
  a_4phase <- a0 %>%
    filter(replicate == 1) %>%
    split(.$mutant) %>%
    map(function(df) df %>%
          head(10000)) %>%
    reduce(rbind)
  a_4phase$V450.A <- ifelse(
    a_4phase$V450.A < 35000, a_4phase$V450.A, 35000)
  
  {## Upper&Lower -------------------------------------------------------------

    a_UL <- a0 %>%
      binP(p = "mCherry.A",
           group = "mutant", N = 25)
    
    a_UL <- a_UL %>%
      split(list(.$mutant, .$mCherry.A.l, .$replicate)) %>%
      map(sepUL) %>%
      reduce(rbind)
    
    mCherryList_UL <- a_UL %>%
      group_by(mutant, mCherry.A.l) %>% 
      summarise(
        medianmCherry = median(mCherry.A)
      )
    
    TPEList_UL <- a_UL %>% 
      filter(pos %in% c("head", "tail")) %>%
      group_by(mutant, mCherry.A.l, replicate, pos) %>% 
      summarise(
        medianV450 = median(V450.A)
      )
    
    s_UL <- TPEList_UL %>%
      split(list(.$mutant, .$mCherry.A.l)) %>%
      map(shape) %>%
      reduce(rbind)
    
    s_UL$mCherry.A.l <- factor(s_UL$mCherry.A.l)
    s_UL$mutant <- factor(
      s_UL$mutant, levels = c("mcherry", "25q", "46q", "97q"))
    f_UL <- s_UL %>% 
      left_join(mCherryList_UL, by = c("mutant", "mCherry.A.l"))
  }
  
  {## Left&Right --------------------------------------------------------------
    a_LR <- a0 %>%
      binP(p = "mCherry.A",
           group = "mutant", N = 25)
    
    a_LR <- a_LR %>%
      split(list(.$mutant, .$mCherry.A.l, .$replicate)) %>%
      map(sepLR) %>%
      reduce(rbind)
    
    mCherryList_LR <- a_LR %>%
      group_by(mutant, mCherry.A.l) %>% 
      summarise(
        medianmCherry = median(mCherry.A)
      )
    
    TPEList_LR <- a_LR %>% 
      filter(pos %in% c("head", "tail")) %>%
      group_by(mutant, mCherry.A.l, replicate, pos) %>% 
      summarise(
        medianV450 = median(V450.A)
      )
    
    s_LR <- TPEList_LR %>%
      split(list(.$mutant, .$mCherry.A.l)) %>%
      map(shape) %>%
      reduce(rbind)
    
    s_LR$mCherry.A.l <- factor(s_LR$mCherry.A.l)
    s_LR$mutant <- factor(
      s_LR$mutant, levels = c("mcherry", "25q", "46q", "97q"))
    f_LR <- s_LR %>% 
      left_join(mCherryList_LR, by = c("mutant", "mCherry.A.l"))
  }
  
  #Expression ~ inclusion%
  a_x_inclusion <- a0 %>%
    filter(replicate == 1) %>%
    binP(p = "mCherry.A",
         group = "mutant", N = 100)
  
  mCherry_Inclu <- a_x_inclusion %>%
    group_by(mutant, mCherry.A.l) %>% 
    summarise(
      count = n(),
      medianmCherry.A = median(mCherry.A),
      inclusion = mean(inclu)
    )
  
  {## Expression~TME ----------------------------------------------------------
    
    # expression histogram of all
    # ggplot(a0, aes(mCherry.A)) +
    #   geom_freqpoly(aes(y= ..density.., colour = mutant)) +
    #   scale_x_continuous(trans = "log10")
    
    a_x_tme <- a0 %>%
      select(-inclu, -mCherry.H, -mCherry.W)
    
    N = 25
    a_x_tme <- a_x_tme %>% # bin expression based on the breaks, discard outliers
      mutate(mCherry.A.l = cut(mCherry.A,
                               breaks = genBreaks(a_x_tme[a_x_tme$mutant == "97q",]),
                               labels = c(1:N)))
    a_x_tme <- a_x_tme[!is.na(a_x_tme$mCherry.A.l),]
    
    mCherry_l <- a_x_tme %>%
      group_by(mCherry.A.l) %>% 
      summarise(
        medianmCherry = median(mCherry.A)
      )
    
    TPE_l <- a_x_tme %>%
      group_by(mutant, replicate, mCherry.A.l) %>% 
      summarise(
        medianV450 = median(V450.A)
      )
    
    {
      q2_m <- TPE_l %>% #25q vs mcherry
        filter(mutant %in% c("25q", "mcherry")) %>%
        split(.$mCherry.A.l)
      q2m <- vector(mode = "list", length(q2_m))
      for(i in seq_along(q2_m)){
        q2m[[i]] <- tTest(
          q2_m[[i]] %>% filter(mutant == "25q"),
          q2_m[[i]] %>% filter(mutant == "mcherry"))
      }
      q2m <- reduce(q2m, rbind) %>%
        mutate(grouping = "25-mcherry")
      
      
      q2_4 <- TPE_l %>% #25q vs 46q
        filter(mutant %in% c("25q", "46q")) %>%
        split(.$mCherry.A.l)
      q24 <- vector(mode = "list", length(q2_4))
      for(i in seq_along(q2_4)){
        q24[[i]] <- tTest(
          q2_4[[i]] %>% filter(mutant == "25q"),
          q2_4[[i]] %>% filter(mutant == "46q"))
      }
      q24 <- reduce(q24, rbind) %>%
        mutate(grouping = "25-46q")
      
      
      q2_9 <- TPE_l %>% #25q vs 97q
        filter(mutant %in% c("25q", "97q")) %>%
        split(.$mCherry.A.l)
      q29 <- vector(mode = "list", length(q2_9))
      for(i in seq_along(q2_9)){
        q29[[i]] <- tTest(
          q2_9[[i]] %>% filter(mutant == "25q"),
          q2_9[[i]] %>% filter(mutant == "97q"))
      }
      q29 <- reduce(q29, rbind) %>%
        mutate(grouping = "25-97q")
      
      #whole dataset
      f_x_tme <- rbind(q2m, q24) %>%
        rbind(q29) %>%
        left_join(mCherry_l, by = c("mCherry.A.l"))
      
      f_x_tme$mutant <- factor(
        f_x_tme$mutant,
        levels = c("97q","46q","mcherry","25q"))
      }
    }
}

# singleStainCheck --------------------------------------------------------

{
  
  fileName <- list.files(
    path = "D:/Publication/TME/sourceData/flowCytometry/Htt/csv/singleStain/",
    pattern = ".exported.FCS3.csv$")
  
  samples <- data.frame()
  
  for (i in 1:length(fileName)) {
    samples <- samples %>% rbind(readFlowCsv(fileName[i]))
  }
  samples <- processFlowCsv(samples)
  
  singlet <- samples[samples$p1 & 
                       samples$sF &
                       samples$sS,] %>%
    filter(mCherry.A > 3)
  
  singlet$mutant <- factor(singlet$mutant,
                           levels = c("untransfect_unstained", 
                                      "untransfect_emi", 
                                      "mcherry_unstained",
                                      "25q_unstained",
                                      "46q_unstained",
                                      "97q_unstained"))
  
  singleStain_summary <- singlet %>% 
    group_by(mutant) %>% 
    summarise(
      
      `singlet #` = sum(p1),
      `% of highV450` = mean(stain) * 100,
      `parent #` = n()
    )
  View(singleStain_summary)
  
  #############################     Plot
  stainPoints <- rbind(c(3, 10^3.5),
                       c(3, 10^5),
                       c(10^5.5, 10^5),
                       c(10^5.5, 10^3.5),
                       c(3, 10^3.5)) %>% as.data.frame()
  
  plot_single1 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "untransfect_unstained"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "untransfect_unstained"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("Blank")
  
  plot_single2 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "untransfect_emi"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "untransfect_emi"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("TME")
  
  plot_single3 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "mcherry_unstained"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "mcherry_unstained"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("mCherry")
  
  plot_single4 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "25q_unstained"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "25q_unstained"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("25Q")
  
  plot_single5 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "46q_unstained"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "46q_unstained"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("46Q")
  
  plot_single6 <- plot_singleStain(
    singlet %>% dplyr::filter(mutant == "97q_unstained"),
    paste0(
      format(
        round((singleStain_summary%>%filter(mutant == "97q_unstained"))$`% of highV450`,2), 
        nsmall = 2) ,
      "%")
  ) +
    ggtitle("97Q")
  
  
  plot_ext8_col1 <- plot_grid(
    plot_single1,
    plot_single2,
    plot_single3,
    plot_single4,
    plot_single5,
    plot_single6,
    labels = c("a", "","","","",""),
    label_size = 7,
    nrow = 2,
    ncol = 3,
    rel_widths = c(30,30,30),
    rel_heights = c(32,32)
  )

}

# Upper Lower Left Right points Check ------------------------------------------------
{
  plotHbinExample <- function(df,  binLevel, mutant = "97q"){
    
    p <- df[df$mutant == mutant &
              df$replicate == 1 &
              df$mCherry.A.l == binLevel,]
    
    pEg <- ggplot(p, aes(mCherry.W, mCherry.H)) +
      geom_path(data = incluPoints, 
                aes(V1, V2),
                alpha = .8,
                size = .2) +
      geom_point(alpha = 1/100, size = .1) +
      scale_x_continuous(name = "Pulse width x 50000 (a.u.)",
                         limits = c(30000, 200000),
                         breaks = 50000 * c(1:4),
                         labels =  c(1:4)) +
      scale_y_continuous(name = expression(log[10]*"(Pulse height) (a.u.)"),
                         trans = "log10",
                         limits = c(10, 100000),
                         breaks = c(10^c(2:5)),
                         labels =  c(2:5)) +
      
      geom_point(data=p[p$pos == "head",], colour=color_H[["head"]], size=.01) +
      geom_point(data=p[p$pos == "tail",], colour=color_H[["tail"]], size=.01)
    
    pEg
  }
  
  plotWbinExample <- function(df,  binLevel, mutant = "97q"){
    
    p <- df[df$mutant == mutant &
              df$replicate == 1 &
              df$mCherry.A.l == binLevel,]
    
    pEg <- ggplot(p, aes(mCherry.W, mCherry.H)) +
      geom_path(data = incluPoints, 
                aes(V1, V2),
                alpha = .8,
                size = .2) +
      geom_point(alpha = 1/20, size = .1) +
      scale_x_continuous(name = "Pulse width x 50000 (a.u.)",
                         limits = c(30000, 200000),
                         breaks = 50000 * c(1:4),
                         labels =  c(1:4)) +
      scale_y_continuous(name = expression(log[10]*"(Pulse height) (a.u.)"),
                         trans = "log10",
                         limits = c(10, 100000),
                         breaks = c(10^c(2:5)),
                         labels =  c(2:5)) +
      
      geom_point(data=p[p$pos == "head",], colour=color_W[["head"]], size=.01) +
      geom_point(data=p[p$pos == "tail",], colour=color_W[["tail"]], size=.01)
    
    pEg
  }
  
  
  plot_H_7 <- plotHbinExample(a_UL, 7) + ggtitle("Bin level = 7")
  plot_H_17 <- plotHbinExample(a_UL, 17) + ggtitle("Bin level = 17")
  plot_H_25 <- plotHbinExample(a_UL, 25) + ggtitle("Bin level = 25")
  
  plot_W_7 <- plotWbinExample(a_LR, 7) + ggtitle("Bin level = 7")
  plot_W_17 <- plotWbinExample(a_LR, 17) + ggtitle("Bin level = 17")
  plot_W_25 <- plotWbinExample(a_LR, 25) + ggtitle("Bin level = 25")
  
}