
# Package loading ---------------------------------------------------------
{
  library(tidygraph)
  library(conicfit)
  library(sf)
  library(tidyverse)
  rm(list = ls())
}

# helper ------------------------------------------------------------------
{
  
  read <- function(address){ #address inside ms folder
    
    inputDir <- paste0(
      "D:/Publication/TME/sourceData/ms/",
      address)
    
    read.csv(inputDir,
             header = TRUE,
             stringsAsFactors = FALSE)
  }
}

# Load data ---------------------------------------------------------------
{
  mouse <- OmicsXZ::mouse
  
  string_db_mouse <- 
    STRINGdb::STRINGdb$new(
      version = "11",
      species = 10090,
      score_threshold = 200)
  
  f_T_lysate <- read("stressor/processed/T_lysate_imputed.csv") %>%
    mutate(C_mean = (impute_LFQ.intensity.U1 +
                       impute_LFQ.intensity.U2 +
                       impute_LFQ.intensity.U3)/3,
           T_mean = (impute_LFQ.intensity.T1 +
                       impute_LFQ.intensity.T2 +
                       impute_LFQ.intensity.T3)/3) %>%
    arrange(C_mean) %>%
    mutate(C_Rank = row_number()/nrow(.)) %>%
    arrange(T_mean) %>%
    mutate(T_Rank = row_number()/nrow(.)) %>%
    select(proteinID, C_Rank, T_Rank)
  
  f_T_CD_CE <- 
    read("stressor/processed/T_CD_CE_imputed.csv") %>% 
    filter(Cys > 0) %>%
    left_join(
      f_T_lysate %>% select(proteinID, C_Rank),
      by = "proteinID") %>%
    rename(rankInLysate = C_Rank)
  
  f_T_TD_TE <- 
    read("stressor/processed/T_TD_TE_imputed.csv") %>% 
    filter(Cys > 0) %>%
    left_join(
      f_T_lysate %>% select(proteinID, T_Rank),
      by = "proteinID") %>%
    rename(rankInLysate = T_Rank)
  
  f_T_CE_TE <- read("stressor/processed/T_CE_TE_imputed.csv") %>% filter(Cys > 0)
  
  T_unfold <- 
    inner_join(f_T_TD_TE %>%
                 filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
                 select(proteinID, Cys),
               f_T_CE_TE %>%
                 filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
                 select(proteinID, Cys))
  
  T_ctrBinder <- 
    read("stressor/binder/T_ctrBinder.csv") %>%
    select(proteinID)
  
  T_tuniBinder <- 
    read("stressor/binder/T_tuniBinder.csv") %>%
    select(proteinID)
  
  # Enrichment
  
  T_ctrBinder_enrich <- read("stressor/enrichment/T_ctrBinder_enrich.csv")
  T_tuniBinder_enrich <- read("stressor/enrichment/T_tuniBinder_enrich.csv")
  
  #relation of term ~ functional categories. contains all such relations 

  ID_funCat <- list(
    ID = c("GO:0007100", "R-MMU-8854050", "R-MMU-69242", 
           "R-MMU-69206", "R-MMU-69202", "R-MMU-69656", "R-MMU-69580", "R-MMU-69275", 
           "R-MMU-69613", "R-MMU-453276", "GO:0045182", "GO:0006417", "GO:0006413", 
           "GO:0005844", "GO:0042788", "GO:0043024", "R-MMU-72702", "R-MMU-72649", 
           "R-MMU-72613", "R-MMU-72737", "R-MMU-73856", "R-MMU-72662", "GO:0006397", 
           "GO:0071013", "R-MMU-72203", "R-MMU-72163", "R-MMU-72172", "GO:0003729", 
           "GO:0003727", "GO:0003725", "GO:0090079", "GO:0022618", "GO:0008064", 
           "GO:0015629", "GO:0030863", "GO:0046785", "GO:0006457", "GO:0051087", 
           "GO:0031072", "GO:0016234", "GO:0016239", "GO:0000502", "mmu03050", 
           "GO:0033962", "GO:0036464"), 
    funCat = c("Cell cycle", "Cell cycle", 
               "Cell cycle", "Cell cycle", "Cell cycle", "Cell cycle", "Cell cycle", 
               "Cell cycle", "Cell cycle", "Cell cycle", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Translation initiation and regulation", "mRNA processing", "mRNA processing", 
               "mRNA processing", "mRNA processing", "mRNA processing", "Nucleotide binding", 
               "Nucleotide binding", "Nucleotide binding", "Nucleotide binding", 
               "Nucleotide binding", "Cytoskeleton organization", "Cytoskeleton organization", 
               "Cytoskeleton organization", "Cytoskeleton organization", "Quality control", 
               "Quality control", "Quality control", "Quality control", "Quality control", 
               "Quality control", "Quality control", "Biomolecular condensates", 
               "Biomolecular condensates")) %>% as.data.frame()
  
  ID_funCat$funCat <- 
    factor(ID_funCat$funCat,
           levels = c("Cell cycle",
                      "Translation initiation and regulation",
                      "mRNA processing",
                      "Nucleotide binding",
                      "Cytoskeleton organization",
                      "Quality control",
                      "Biomolecular condensates"))
  
  #flow cytometry
  flow_tunicamycin <- read.csv(
    "D:/Publication/TME/sourceData/flowCytometry/drugs/tunicamycin_gated.csv",
    header = TRUE,
    stringsAsFactors = FALSE)
  
  flow_mg132 <- read.csv(
    "D:/Publication/TME/sourceData/flowCytometry/drugs/MG132_gated.csv",
    header = TRUE,
    stringsAsFactors = FALSE)
  
  #noMatchBetweenRuns
  
  T_NoMatch <- read.csv(
    "D:/Publication/TME/sourceData/ms/stressor/raw/T_NoMatchForAll.csv",
    header = TRUE, stringsAsFactors = FALSE) %>%
    filter(Only.identified.by.site != "+" &
             Reverse != "+" &
             Potential.contaminant != "+") %>%
    rename(gene = `Gene.names`) %>%
    mutate(DD1 = log2(Intensity.A1),
           DD2 = log2(Intensity.A2),
           DD3 = log2(Intensity.A3),
           TD1 = log2(Intensity.B1),
           TD2 = log2(Intensity.B2),
           TD3 = log2(Intensity.B3),
           DE1 = log2(Intensity.C1),
           DE2 = log2(Intensity.C2),
           DE3 = log2(Intensity.C3),
           TE1 = log2(Intensity.D1),
           TE2 = log2(Intensity.D2),
           TE3 = log2(Intensity.D3)) %>%
    select(gene, id,
           DD1, DD2, DD3,
           TD1, TD2, TD3,
           DE1, DE2, DE3,
           TE1, TE2, TE3)
  is.na(T_NoMatch) <- sapply(T_NoMatch, is.infinite)
  
  T_NoMatch_long <- T_NoMatch %>%
    pivot_longer(c("DD1","DD2","DD3",
                   "TD1","TD2","TD3",
                   "DE1","DE2","DE3",
                   "TE1","TE2","TE3"),
                 names_to = "run",
                 values_to = "intensity") %>%
    separate(run, 
             into = c("experiment", "rep"),
             sep = 2,
             convert = TRUE) %>%
    separate(experiment, 
             into = c("treatment_tmp", "stain_tmp"),
             sep = 1,
             convert = TRUE) %>%
    mutate(treatment = fct_recode(treatment_tmp,
                                  dmso = "D",
                                  tunicamycin = "T"),
           stain = fct_recode(stain_tmp,
                              dmso = "D",
                              emi = "E")) %>%
    select(-treatment_tmp, -stain_tmp)
}

# flowCytometry -----------------------------------------------------------
{
  # tunicamycin dataset
  scatter_tunicamycin <- flow_tunicamycin %>%
    group_by(treatment, rep) %>%
    summarize(
      TPE = median(V500.A),
    ) %>%
    filter(treatment %in% c("Control", "Tunicamycin")) 
  
  scatter_tunicamycin$TPE <- scatter_tunicamycin$TPE / mean((scatter_tunicamycin%>%filter(treatment=="Control"))$TPE)
  
  compare_means(TPE~treatment, data = scatter_tunicamycin, method = 't.test')
  
  # mg132 dataset
  scatter_mg132 <- flow_mg132 %>%
    group_by(treatment, rep) %>%
    summarize(
      TPE = median(V500.A),
    ) %>%
    filter(treatment %in% c("Control", "MG132")) 
  
  scatter_mg132$TPE <- scatter_mg132$TPE / mean((scatter_mg132%>%filter(treatment=="Control"))$TPE)
  
  compare_means(TPE~treatment, data = scatter_mg132, method = 't.test')
  
  # two datasets
  scatter_2 <- scatter_mg132 %>%
    rbind(scatter_tunicamycin%>%filter(treatment == "Tunicamycin"))
  
  scatter_2$treatment <- factor(
    scatter_2$treatment,
    levels = c("Control","Tunicamycin","MG132"))
  
  compare_means(TPE~treatment, data = scatter_2, method = 't.test')
  
  error_tunicamycin <- scatter_tunicamycin %>%
    group_by(treatment) %>%
    summarise(
      weight = mean(TPE),
      sd = sd(TPE),
      n = n()
    )
  
  error_mg132 <- scatter_mg132 %>%
    group_by(treatment) %>%
    summarise(
      weight = mean(TPE),
      sd = sd(TPE),
      n = n()
    )
  
  error_2 <- error_mg132 %>%
    rbind(error_tunicamycin %>% filter(treatment == "Tunicamycin"))
  
  flow_2 <-
    ggplot(scatter_2, aes(treatment, TPE)) +
    geom_point(
      size = .4,
      position = position_jitter(w = 0.1, h = 0)) +
    geom_errorbar(
      data = error_2,
      aes(treatment, weight,
          ymin = weight - sd,
          ymax = weight + sd),
      size = .2,
      width = .15) +
    scale_y_continuous(
      limits = c(.5, 1.25),
      breaks = c(.5, .75, 1, 1.25)) +
    stat_compare_means(
      aes(label = ..p.signif..),
      comparisons = 
        list(c("Control", "Tunicamycin"), 
             c("Control", "MG132")),
      label.y = c(1.14, 1.22),
      tip.length = .01,
      size = 5/.pt,
      bracket.size = 0.1,
      vjust = .5,
      method = 't.test'
    ) +
    ylab("Fluorescence fold change") +
    ggtitle("TME signals - Stress vs Control") +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank())
}

# proteinCounts -----------------------------------------------------------

{
  
  pr_QTY_T <- T_NoMatch_long %>%
    group_by(id, treatment, stain) %>% 
    summarise(
      na_QTY = sum(is.na(intensity)),
    )
  
  pNum_T <- data.frame(
    treatment = c("DD", "TD", "DE", "TE"),
    proteinNumber = c(nrow(pr_QTY_T %>%
                             filter(treatment == "dmso" & 
                                      stain == "dmso" & 
                                      na_QTY == 0)),
                      nrow(pr_QTY_T %>%
                             filter(treatment == "tunicamycin" & 
                                      stain == "dmso" & 
                                      na_QTY == 0)),
                      nrow(pr_QTY_T %>%
                             filter(treatment == "dmso" & 
                                      stain == "emi" & 
                                      na_QTY == 0)),
                      nrow(pr_QTY_T %>%
                             filter(treatment == "tunicamycin" & 
                                      stain == "emi" & 
                                      na_QTY == 0)))
  )
  
  pNum_T$treatment <- factor(pNum_T$treatment,
                             c("DD", "TD", "DE", "TE"))
  plot_pNum_T <- 
    ggplot(pNum_T, aes(treatment, proteinNumber)) +
    geom_bar(stat = 'identity') +
    geom_text(aes(y = proteinNumber + 150,
                  label = proteinNumber),
              size = 6/.pt) +
    scale_y_continuous(breaks = c(100*c(5,10,15,20)), labels = c(5,10,15,20)) +
    ylab("Protein counts x 100") +
    ggtitle("Protein counts in tunicamycin dataset") +
    theme(
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank())
}

# Gating strategies -------------------------------------------------------

{
  process <- function(df) {
    
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
    p1Points <-conicfit::calculateEllipse(x = 100000, y = 80000,
                                          a = 60000, b= 45000,
                                          angle = 60, steps = 10)
    p1 <- isInGate(df, p1Points,
                   "FSC.A", "SSC.A")
    
    # to generate a rectangular polygon points 1e+03
    recRange <- function(xmin, xmax, ymin, ymax) {
      
      rbind(c(xmin * 1e+03,ymin * 1e+03),
            c(xmin * 1e+03,ymax * 1e+03),
            c(xmax * 1e+03,ymax * 1e+03),
            c(xmax * 1e+03,ymin * 1e+03),
            c(xmin * 1e+03,ymin * 1e+03))
    }
    
    #Singlet from FSC
    sFPoints <- recRange(75, 150, 10, 90)
    sFPoints <- 
      sF <- isInGate(df, sFPoints, 
                     "FSC.W", "FSC.H")
    
    #Singlet from SSC
    sSPoints <- recRange(95, 175, 1, 70)
    sS <- isInGate(df, sSPoints, "SSC.W", "SSC.H")
    
    #live from Singlet
    livePoints <- rbind(c(4e04, 5),
                        c(4e04, 4000),
                        c(18e04, 4000),
                        c(18e04, 5),
                        c(4e04, 5))
    live <- isInGate(df, livePoints, "FSC.A", "APC.A")
    
    df %>% mutate(p1 = p1,
                  sF = sF,
                  sS = sS,
                  live = live)
  }
  
  readFile <- function(fileName) {
    
    file <- read.csv(fileName,
                     header = TRUE, stringsAsFactors = FALSE,
                     nrow = 5000)
    file %>% dplyr::select(FSC.H, FSC.W, FSC.A, 
                    SSC.H, SSC.W, SSC.A,
                    V450.A, APC.A) %>%
      mutate("sampleName" = str_remove_all(fileName,
                                           "Specimen_001_|.exported.FCS3.csv")) %>%
      separate(sampleName, into = c("treatment","rep"),
               sep = "_", convert = TRUE)
  }
  
  setwd("D:/Publication/TME/sourceData/flowCytometry/drugs/Tunicamycin/csv")
  
  fileName <- list.files(pattern = "^Specimen_001_Blank")
  samples <- data.frame()
  
  for (i in 1:length(fileName)) {
    samples <- samples %>% rbind(readFile(fileName[i]))
  }
  
  samples <- process(samples)
  
  # P1 in all cells 

  df <- samples
  summary <- df %>% 
    group_by(treatment) %>% 
    summarise(
      
      `P1 #` = sum(p1),
      `% of the parent` = mean(p1) * 100,
      `parent #` = n()
    )
  
  p1Points <-calculateEllipse(x = 80000, y = 80000,
                              a = 60000, b= 45000,
                              angle = 60, steps = 10) %>% as.data.frame()
  
  plot_gating1 <- 
    ggplot(df, aes(FSC.A, SSC.A)) +
    geom_point(alpha = 1/20,
               size = .5) +
    scale_x_continuous(limits = c(0 , 20* 1e04),
                       breaks = c(0,100000, 200000),
                       labels = c(0, 10, 20),
                       name = "FSC.A x 10000") +
    scale_y_continuous(limits = c(0 , 20* 1e04),
                       breaks = c(0,100000, 200000),
                       labels = c(0, 10, 20),
                       name = "SSC.A x 10000") +
    geom_path(data = p1Points, aes(V1, V2), size = .2) +
    ggtitle("General gating strategies for flow cytometry")
  

  ## see FSC single in P1 
  df <- samples[samples$p1 ,]
  
  summary <- df %>% 
    group_by(treatment) %>% 
    summarise(
      
      `sF #` = sum(sF),
      `% of the parent` = mean(sF) * 100,
      `parent #` = n()
    )
  
  # to generate a rectangular polygon points 1e+03
  recRange <- function(xmin, xmax, ymin, ymax) {
    
    rbind(c(xmin * 1e+03,ymin * 1e+03),
          c(xmin * 1e+03,ymax * 1e+03),
          c(xmax * 1e+03,ymax * 1e+03),
          c(xmax * 1e+03,ymin * 1e+03),
          c(xmin * 1e+03,ymin * 1e+03))
  }
  
  sFPoints <- recRange(75, 150, 10, 90) %>% as.data.frame()
  
  plot_gating2 <- 
    ggplot(df, aes(FSC.W, FSC.H)) +
    geom_point(alpha = 1/20, size = .5) +
    scale_x_continuous(limits = c(5e+04 , 22e+04),
                       breaks = 10000*c(5,10,15,20),
                       labels = c(5,10,15,20),
                       name = "FSC.W x 10000") +
    scale_y_continuous(limits = c(0e+04 , 10e+04),
                       breaks = 25000*c(1,2,3,4),
                       labels = c(1,2,3,4),
                       name = "FSC.H x 25000") +
    geom_path(data = sFPoints, aes(V1, V2), size = .2) +
    ggtitle(" ")
  
  
  ################################### see SSC single in FSC
  df <- samples[samples$p1 & 
                  samples$sF,]
  
  summary <- df %>%
    group_by(treatment) %>% 
    summarise(
      
      `sS #` = sum(sS),
      `% of the parent` = mean(sS) * 100,
      `parent #` = n()
    )
  
  sSPoints <- recRange(95, 200, 1, 70) %>% as.data.frame()
  
  plot_gating3 <- 
    ggplot(df, aes(SSC.W, SSC.H)) +
    geom_point(alpha = 1/20, size = .5) +
    scale_x_continuous(limits = c(9e04 , 30e04),
                       breaks = 50000*c(2:6),
                       labels = c(2:6),
                       name = "SSC.W x 50000") +
    scale_y_continuous(limits = c(0e04 , 8e04),
                       breaks = 20000*c(1:4),
                       labels = c(1:4),
                       name = "SSC.H x 20000") +
    geom_path(data = sSPoints, aes(V1, V2), size = .2) +
    ggtitle(" ")
    
  
  ##############################   see live in singlets
  df <- samples[samples$p1 & 
                  samples$sF &
                  samples$sS,]
  
  summary <- df %>% 
    group_by(treatment) %>% 
    summarise(
      
      'live #' = sum(live),
      '% of the parent' = mean(live) * 100,
      'parent #' = n()
    )
  
  livePoints <- rbind(c(4e04, 5),
                      c(4e04, 3000),
                      c(16e04, 3000),
                      c(16e04, 5),
                      c(4e04, 5)) %>% as.data.frame()
  
  plot_gating4 <- 
    ggplot(df, aes(FSC.A, APC.A)) +
    geom_point(alpha = 1/20, size = .5) +
    scale_x_continuous(limits = c(3e04 , 20e04),
                       breaks = 50000*c(1:4),
                       labels = c(1:4),
                       name = "FSC.A x 50000") +
    scale_y_continuous(trans = "log10",
                       breaks = 10^c(1,3,5),
                       labels = c(1,3,5),
                       name = expression(log[10]*"(APC.A) (a.u.)")) +
    geom_path(data = livePoints, aes(V1, V2), size = .2) +
    ggtitle(" ")
  
  plot_gating <- plot_grid(
    plot_gating1,
    plot_gating2,
    plot_gating3,
    plot_gating4,
    align = 'vh',
    axis = 'bl',
    nrow = 2,
    ncol = 2,
    rel_widths = c(45, 45)
  )
}

# Load data PD ---------------------------------------------------------------
{
  human <- OmicsXZ::human
  
  string_db_human <-
    STRINGdb::STRINGdb$new(
      version = "11",
      species = 9606,
      score_threshold = 200)
  
  f_P_lysate <- read("pd/processed/P_lysate_imputed.csv") %>%
    mutate(H_mean = (impute_LFQ.intensity.C1 +
                       impute_LFQ.intensity.C2 +
                       impute_LFQ.intensity.C3 +
                       impute_LFQ.intensity.C4)/4,
           P_mean = (impute_LFQ.intensity.P1 +
                       impute_LFQ.intensity.P2 +
                       impute_LFQ.intensity.P3 +
                       impute_LFQ.intensity.P4)/4) %>%
    arrange(H_mean) %>%
    mutate(H_Rank = row_number()/nrow(.)) %>%
    arrange(P_mean) %>%
    mutate(P_Rank = row_number()/nrow(.))
  
  f_P_HD_HE <-
    read("pd/processed/P_HD_HE_imputed.csv") %>%
    filter(Cys > 0) %>%
    left_join(
      f_P_lysate %>% select(proteinID, H_Rank),
      by = "proteinID") %>%
    rename(rankInLysate = H_Rank)
  
  f_P_PD_PE <-
    read("pd/processed/P_PD_PE_imputed.csv") %>%
    filter(Cys > 0) %>%
    left_join(
      f_P_lysate %>% select(proteinID, P_Rank),
      by = "proteinID") %>%
    rename(rankInLysate = P_Rank)
  
  f_P_HE_PE <- read("pd/processed/P_HE_PE_imputed.csv")
  
  P_moreUnfold <-
    inner_join(f_P_PD_PE %>%
                 filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
                 select(proteinID, Cys),
               f_P_HE_PE %>%
                 filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
                 select(proteinID, Cys))
  
  P_lessUnfold <-
    inner_join(f_P_HD_HE %>%
                 filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
                 select(proteinID, Cys),
               f_P_HE_PE %>%
                 filter(p.adjust < 0.05 & log2fc < -log2(1.5)) %>%
                 select(proteinID, Cys))
  
  P_healBinder <-
    read("pd/binder/healBinder.csv") %>%
    select(proteinID)
  
  P_pdBinder <-
    read("pd/binder/pdBinder.csv") %>%
    select(proteinID)
  
  ID_funCat <- list(
    ID = c("GO:0070628", "GO:0005844", "R-HSA-156902", 
           "R-HSA-156842", "R-HSA-72764", "GO:0097526", "GO:0000502", "R-HSA-72613", 
           "GO:0034063", "R-HSA-72649", "GO:0006413", "hsa03010", "hsa03050", 
           "GO:0042255", "R-HSA-72187", "GO:0022618", "GO:0071826", "GO:0008135", 
           "R-HSA-72165", "R-HSA-72766", "GO:0019843", "R-HSA-72163", "R-HSA-72172", 
           "GO:0003725", "GO:0045182", "GO:0043021", "hsa03040", "R-HSA-72203", 
           "R-HSA-453276", "GO:0006457", "GO:0000377", "GO:0000398", "GO:0051082", 
           "GO:0003729", "GO:0003727", "R-HSA-69206", "GO:0010972", "R-HSA-453279"), 
    funCat = c("Quality control", "Quality control", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "mRNA processing", "Quality control", "Translation initiation and regulation", 
               "Biomolecular condensates", "Translation initiation and regulation", 
               "Translation initiation and regulation", "Translation initiation and regulation", 
               "Quality control", "Translation initiation and regulation", "mRNA processing", 
               "Biomolecular condensates", "Biomolecular condensates", "Nucleotide binding", 
               "mRNA processing", "Translation initiation and regulation", "Nucleotide binding", 
               "mRNA processing", "mRNA processing", "Nucleotide binding", "Translation initiation and regulation", 
               "Biomolecular condensates", "mRNA processing", "mRNA processing", 
               "Cell cycle", "Quality control", "mRNA processing", "mRNA processing", 
               "Quality control", "Nucleotide binding", "Nucleotide binding", 
               "Cell cycle", "Cell cycle", "Cell cycle")) %>% as.data.frame()
  
  ID_funCat$funCat <- 
    factor(ID_funCat$funCat,
           levels = c("Cell cycle",
                      "Translation initiation and regulation",
                      "mRNA processing",
                      "Nucleotide binding",
                      "Quality control",
                      "Biomolecular condensates"))
  
}

# PCA ---------------------------------------------------------------------

{
  m.ls <-
    f_P_lysate %>%
    select(contains("impute_LFQ.intensity."))
  rownames(m.ls) <- f_P_lysate$proteinID
  colnames(m.ls) <- sub("impute_LFQ\\.intensity\\.", "", colnames(m.ls))
  colnames(m.ls) <- sub("C", "H", colnames(m.ls))
  
  m.ap <-
    select(f_P_HE_PE, contains("impute_LFQ.intensity."))
  rownames(m.ap) <- f_P_HE_PE$proteinID
  colnames(m.ap) <- sub("impute_LFQ\\.intensity\\.", "", colnames(m.ap))
  colnames(m.ap) <- sub("E", "", colnames(m.ap))
  
  # PCA
  pca.ls <- prcomp(t(m.ls), scale. = TRUE)
  pca.ap <- prcomp(t(m.ap), scale. = TRUE)
  
  #Scree Plot
  
  var.ls <- pca.ls$sdev^2
  var.per.ls <- data.frame(per = round(var.ls/sum(var.ls)*100, 1),
                           com = 1:8)
  
  var.ap <- pca.ap$sdev^2
  var.per.ap <- data.frame(per = round(var.ap/sum(var.ap)*100, 1),
                           com = 1:8)
  
}

