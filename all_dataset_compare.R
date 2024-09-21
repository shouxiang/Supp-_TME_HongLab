# Package loading ---------------------------------------------------------
{
  library(ggprism)
  library(cowplot)
  library(ggraph)
  library(tidygraph)
  library(gtable)
  library(grid)
  library(gridExtra)
  library(tidyverse)
  library(ggpubr)
  library(ggrepel)
  library(idpr)
  rm(list = ls())
}

# utils---------------------------------------------------------------
{
  read_ms <- function(address){ 
    
    inputDir <- paste0(
      BASE_DIR, '/ms',
      address)
    read.csv(inputDir,
             header = TRUE,
             stringsAsFactors = FALSE)
  }
  EXPORT_ADDRESS <- "C:/Users/zhang/Desktop"
  BASE_DIR <- "E:/Publication/TME/sourceData" # dir address to sourceData
}

# Load data ---------------------------------------------------------------
{
  mouse <- OmicsXZ::mouse
  string_db_mouse <- 
    STRINGdb::STRINGdb$new(
      version = "11",
      species = 10090,
      score_threshold = 200)
  
  mouse_stringid <- string_db_mouse$map(
    data.frame(UniProt = OmicsXZ::mouse$proteinID),
    "UniProt",
    removeUnmappedRows = TRUE) %>%
    dplyr::rename(proteinID = UniProt)
  
  paxdb_mouse_whole <- read_ms('/paxdb/10090-WHOLE_ORGANISM-integrated.csv') %>%
    dplyr::rename(STRING_id = string_external_id) %>%
    dplyr::inner_join(mouse_stringid) %>%
    group_by(proteinID) %>%
    summarize(
      ppm = max(log(abundance, base = 10))
    )
  
  paxdb_mouse_brain <- read_ms('/paxdb/10090-BRAIN-integrated.csv') %>%
    dplyr::rename(STRING_id = string_external_id) %>%
    dplyr::inner_join(mouse_stringid) %>%
    group_by(proteinID) %>%
    summarize(
      ppm = max(log(abundance, base = 10))
    )
  
  # Tunicamycin dataset
  T_ctrBinder <- 
    read_ms("/stressor/binder/T_ctrBinder.csv")
  
  T_tuniBinder <- 
    read_ms("/stressor/binder/T_tuniBinder.csv")
  
  T_lysate <- read_ms('/stressor/raw/T_lysate_ibaq.csv') %>%
    dplyr::filter(Only.identified.by.site == '') %>%
    dplyr::filter(Reverse == '') %>%
    dplyr::filter(Potential.contaminant == '') %>%
    dplyr::select(Majority.protein.IDs, iBAQ) %>%
    dplyr::mutate(proteinID = str_split_fixed(Majority.protein.IDs, ";", 2)[, 1]) %>%
    dplyr::select(proteinID, iBAQ) %>%
    dplyr::filter(iBAQ !=0) %>%
    dplyr::inner_join(
      read_ms('/stressor/processed/T_lysate_imputed.csv') %>%   # The imputed data
        dplyr::select(proteinID)
    ) %>%
    dplyr::mutate(iBAQ = log10(iBAQ))
  
  T_ctrBinder_undetectable <- 
    dplyr::left_join(T_ctrBinder, T_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  T_tuniBinder_undetectable <- 
    dplyr::left_join(T_tuniBinder, T_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  # MG132 dataset
  M_ctrBinder <- 
    read_ms("/stressor/binder/M_ctrBinder.csv")
  
  M_mgBinder <- 
    read_ms("/stressor/binder/M_mgBinder.csv")
  
  M_lysate <- read_ms('/stressor/raw/M_lysate_ibaq.csv') %>%
    dplyr::filter(Only.identified.by.site == '') %>%
    dplyr::filter(Reverse == '') %>%
    dplyr::filter(Potential.contaminant == '') %>%
    dplyr::select(Majority.protein.IDs, iBAQ) %>%
    dplyr::mutate(proteinID = str_split_fixed(Majority.protein.IDs, ";", 2)[, 1]) %>%
    dplyr::select(proteinID, iBAQ) %>%
    dplyr::filter(iBAQ !=0) %>%
    dplyr::inner_join(
      read_ms('/stressor/processed/M_lysate_imputed.csv') %>%   # The imputed data
        dplyr::select(proteinID)
    ) %>%
    dplyr::mutate(iBAQ = log10(iBAQ))
  
  M_ctrBinder_undetectable <-
    dplyr::left_join(M_ctrBinder, M_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  M_mgBinder_undetectable <-
    dplyr::left_join(M_mgBinder, M_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  # human ref
  human <- OmicsXZ::human
  string_db_human <-
    STRINGdb::STRINGdb$new(
      version = "11",
      species = 9606,
      score_threshold = 200)
  
  human_stringid <- string_db_human$map(
    data.frame(UniProt = OmicsXZ::human$proteinID),
    "UniProt",
    removeUnmappedRows = TRUE) %>%
    dplyr::rename(proteinID = UniProt)
  
  paxdb_human_whole <- read_ms('/paxdb/9606-WHOLE_ORGANISM-integrated.csv') %>%
    dplyr::rename(STRING_id = string_external_id) %>%
    dplyr::inner_join(human_stringid) %>%
    group_by(proteinID) %>%
    summarize(
      ppm = max(log(abundance, base = 10))
    )
  
  paxdb_human_lymph <- read_ms('/paxdb/9606-LYMPH_NODE-integrated.csv') %>%
    dplyr::rename(STRING_id = string_external_id) %>%
    dplyr::inner_join(human_stringid) %>%
    group_by(proteinID) %>%
    summarize(
      ppm = max(log(abundance, base = 10))
    )
  
  # PD dataset
  P_healBinder <-
    read_ms("/pd/binder/healBinder.csv")
  
  P_pdBinder <-
    read_ms("/pd/binder/pdBinder.csv")
  
  P_lysate <- read_ms('/pd/raw/P_lysate_ibaq.csv') %>%
    dplyr::filter(Only.identified.by.site == '') %>%
    dplyr::filter(Reverse == '') %>%
    dplyr::filter(Potential.contaminant == '') %>%
    dplyr::select(Majority.protein.IDs, iBAQ) %>%
    dplyr::mutate(proteinID = str_split_fixed(Majority.protein.IDs, ";", 2)[, 1]) %>%
    dplyr::select(proteinID, iBAQ) %>%
    dplyr::filter(iBAQ !=0) %>%
    dplyr::inner_join(
      read_ms('/pd/processed/P_lysate_imputed.csv') %>%   # The imputed data
        dplyr::select(proteinID)
    ) %>%
    dplyr::mutate(iBAQ = log10(iBAQ))
  
  P_healBinder_undetectable <- 
    dplyr::left_join(P_healBinder, P_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  P_pdBinder_undetectable <-
    dplyr::left_join(P_pdBinder, P_lysate) %>%
    dplyr::filter(is.na(iBAQ))
  
  #n2a dataset
  n2a <- read_ms('/paxdb/n2a-abundance/n2a.csv') %>%
    dplyr::select(Majority.protein.IDs, iBAQ, rank_N2a) %>%
    dplyr::mutate(proteinID = str_split_fixed(Majority.protein.IDs, ";", 2)[, 1]) %>%
    dplyr::filter(rank_N2a != 'NaN') %>%
    dplyr::mutate(ppm = iBAQ) %>%
    dplyr::select(proteinID, ppm)
  
}

# Protein Abundance -----------------------------------------------------------
{
  
  ## Fun -----------------------------------------------------------
  {
    gen_compare_abundance <- function(ref, lysate, binder, undetected, title){
      
      color <- c(
        Lysate = "#AAACAF",
        Binder = "#71AACD",
        Undetectable = "#E9738D"
      )
      
      lysate <- dplyr::inner_join(ref, lysate%>%dplyr::select(proteinID)) %>%
        dplyr::mutate(group = 'Lysate')
      
      Undetectable <- dplyr::inner_join(ref, undetected%>%dplyr::select(proteinID)) %>%
        dplyr::mutate(group = 'Undetectable')
      
      binder_exclude <- binder %>% 
        dplyr::select(proteinID) %>% 
        dplyr::filter(!proteinID %in% (Undetectable %>% dplyr::select(proteinID) %>% pull(proteinID))) %>%
        dplyr::inner_join(ref) %>%
        dplyr::mutate(group = 'Binder')
      
      data <- 
        rbind(lysate, binder_exclude) %>%
        rbind(Undetectable)
      
      data$group <- factor(
        data$group, 
        levels = c('Lysate', 'Binder', 'Undetectable'))
      
      plot <- 
        ggplot(data, aes(group, ppm)) +
        geom_violin(
          aes(fill = group, color = group), 
          width = 1, 
          trim = FALSE, 
          size = .3) +
        scale_fill_manual(values = color) +
        scale_color_manual(values = color) +
        geom_boxplot(
          width = .15,
          lwd = .2,
          fill = "black", 
          color = "black", 
          outlier.colour = NA) +
        stat_summary(
          fun = median, 
          geom = "point", 
          fill = "white", 
          shape = 21,
          size = 1.5) +
        scale_y_continuous(
          limits = c(-2.5, 5.5),
          breaks = seq(-5, 5, 2.5),
          labels = seq(-5, 5, 2.5)
        ) +
        stat_compare_means(
          aes(label = ..p.signif..),
          comparisons =
            list(c("Lysate", "Undetectable"),
                 c("Binder", "Undetectable")),
          label.y = c(4.8, 4.2),
          tip.length = .01,
          size = 5/.pt,
          bracket.size = 0.1,
          vjust = .5
        ) +
        ylab(expression(paste(log[10],"(ppm)",sep=""))) +
        ggtitle(title) +
        THEME0 +
        theme(
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none'
        )
      plot
    }
    
    gen_compare_abundance_ibaq <- function(ref, lysate, binder, undetected, title){
      
      color <- c(
        Lysate = "#AAACAF",
        Binder = "#71AACD",
        Undetectable = "#E9738D"
      )
      
      lysate <- dplyr::inner_join(ref, lysate%>%dplyr::select(proteinID)) %>%
        dplyr::mutate(group = 'Lysate')
      
      Undetectable <- dplyr::inner_join(ref, undetected%>%dplyr::select(proteinID)) %>%
        dplyr::mutate(group = 'Undetectable')
      
      binder_exclude <- binder %>% 
        dplyr::select(proteinID) %>% 
        dplyr::filter(!proteinID %in% (Undetectable %>% dplyr::select(proteinID) %>% pull(proteinID))) %>%
        dplyr::inner_join(ref) %>%
        dplyr::mutate(group = 'Binder')
      
      data <- 
        rbind(lysate, binder_exclude) %>%
        rbind(Undetectable)
      
      data$group <- factor(
        data$group, 
        levels = c('Lysate', 'Binder', 'Undetectable'))
      
      plot <- 
        ggplot(data, aes(group, ppm)) +
        geom_violin(
          aes(fill = group, color = group), 
          width = 1, 
          trim = FALSE, 
          size = .3) +
        scale_fill_manual(values = color) +
        scale_color_manual(values = color) +
        geom_boxplot(
          width = .15,
          lwd = .2,
          fill = "black", 
          color = "black", 
          outlier.colour = NA) +
        stat_summary(
          fun = median, 
          geom = "point", 
          fill = "white", 
          shape = 21,
          size = 1.5) +
        # scale_y_continuous(
        #   limits = c(-2.5, 5.5),
        #   breaks = seq(-5, 5, 2.5),
        #   labels = seq(-5, 5, 2.5)
        # ) +
        stat_compare_means(
          aes(label = ..p.signif..),
          comparisons =
            list(c("Lysate", "Undetectable"),
                 c("Binder", "Undetectable")),
          label.y = c(40, 37),
          tip.length = .01,
          size = 5/.pt,
          bracket.size = 0.1,
          vjust = .5
        ) +
        ylab(expression(paste(log[2],"(iBAQ)",sep=""))) +
        ggtitle(title) +
        THEME0 +
        theme(
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none'
        )
      plot
    }
    
    
    gen_percentage_bar_plot <- function(ref, lysate, binder, undetectable, title, ymax, n=20){
      gen_distribution <- function(ref, df, n, group){
        
        complete_percentiles <- data.frame(
          percentile = as.character(c(1:n, "U"))
        )
        
        ref <- ref %>%
          mutate(percentile = ntile(ppm, n))
        
        df_percentiles <- dplyr::left_join(df, ref) %>%
          dplyr::mutate(
            percentile = ifelse(
              is.na(percentile),
              "U", 
              as.character(percentile)))
        
        df_distribution <- df_percentiles %>%
          count(percentile, name = "count") %>%
          mutate(percentage = count / nrow(df_percentiles))
        
        complete_distribution <- complete_percentiles %>%
          left_join(df_distribution, by = "percentile") %>%
          mutate(
            count = ifelse(is.na(count), 0, count),
            percentage = ifelse(is.na(percentage), 0, percentage)
          ) %>%
          mutate(group = group)
        complete_distribution
      }
      
      binder_exclude <- binder %>% 
        dplyr::filter(!proteinID %in% (undetectable %>% dplyr::select(proteinID) %>% pull(proteinID)))
      
      lysate_distribution <- gen_distribution(ref, lysate, n, 'Lysate')
      binder_exclude_distribution <- gen_distribution(ref, binder_exclude, n, 'Binder')
      undetectable_distribution <- gen_distribution(ref, undetectable, n, 'Undetectable')
      
      a <- rbind(lysate_distribution, binder_exclude_distribution) %>%
        rbind(undetectable_distribution) %>%
        dplyr::mutate(group = factor(group, levels = c("Undetectable", "Binder", 'Lysate')))
      
      p <- 
        ggplot(a, aes(x = factor(percentile, levels = as.character(c(1:20, "U"))), y = percentage, fill = group)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_x_discrete(labels = c(5 * (1:20), "U")) +
        
        labs(x = "Percentile", y = "Percent of Proteins", title = title) +
        scale_fill_manual(values = COLOR_ABUNDANCE) + 
        THEME0 + 
        theme(
          legend.position = c(0.3, 0.8),
          legend.key.size = unit(3, 'mm'),
          legend.title = element_blank(),
          legend.spacing.y = unit(5, 'mm'),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45,size = 3)
        )
      
      custom_labels <- function(x) {
        percent(x, accuracy = 1)
      }
      
      if(ymax>1){
        p <- p + scale_y_continuous(
          limits = c(0, ymax),
          labels = scales::percent,
          breaks = seq(0, 1, 0.25),
          expand = expansion(mult = c(0.02, 0.05)))
      } else{
        p <- p + scale_y_continuous(
          limits = c(0, ymax),
          labels = scales::percent,
          expand = expansion(mult = c(0.02, 0.05))) 
      }
      
      p
    }
  }
  
# Repeatability -----------------------------------------------------------
{
  ## Load data -----------------------------------------------------------
  {
    T_NoMatch <- read_ms("/stressor/raw/T_NoMatchForAll.csv") %>%
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

    M_NoMatch <- read_ms("/stressor/raw/M_NoMatchForAll.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(DD1 = log2(Intensity.A1),
             DD2 = log2(Intensity.A2),
             DD3 = log2(Intensity.A3),
             MD1 = log2(Intensity.B1),
             MD2 = log2(Intensity.B2),
             MD3 = log2(Intensity.B3),
             DE1 = log2(Intensity.C1),
             DE2 = log2(Intensity.C2),
             DE3 = log2(Intensity.C3),
             ME1 = log2(Intensity.D1),
             ME2 = log2(Intensity.D2),
             ME3 = log2(Intensity.D3)) %>%
      select(gene, id,
             DD1, DD2, DD3,
             MD1, MD2, MD3,
             DE1, DE2, DE3,
             ME1, ME2, ME3)
    is.na(M_NoMatch) <- sapply(M_NoMatch, is.infinite)

    PD_NoMatch <- read_ms("/pd/raw/NoMatchForAll_AP.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(
        HD1 = log2(Intensity.HD1),
        HD2 = log2(Intensity.HD2),
        HD3 = log2(Intensity.HD3),
        HD4 = log2(Intensity.HD4),

        HE1 = log2(Intensity.HE1),
        HE2 = log2(Intensity.HE2),
        HE3 = log2(Intensity.HE3),
        HE4 = log2(Intensity.HE4),

        PD1 = log2(Intensity.PD1),
        PD2 = log2(Intensity.PD2),
        PD3 = log2(Intensity.PD3),
        PD4 = log2(Intensity.PD4),

        PE1 = log2(Intensity.PE1),
        PE2 = log2(Intensity.PE2),
        PE3 = log2(Intensity.PE3),
        PE4 = log2(Intensity.PE4)) %>%
      dplyr::select(gene,
                    HD1, HD2, HD3, HD4,
                    HE1, HE2, HE3, HE4,
                    PD1, PD2, PD3, PD4,
                    PE1, PE2, PE3, PE4)
    is.na(PD_NoMatch) <- sapply(PD_NoMatch, is.infinite)

    htt_T_NoMatch <- read_ms("/htt/raw/total.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(
        P1_1 = log2(Intensity.P1_1),
        P1_2 = log2(Intensity.P1_2),
        P1_3 = log2(Intensity.P1_3),

        P2a_1 = log2(Intensity.P2a_1),
        P2a_2 = log2(Intensity.P2a_2),
        P2a_3 = log2(Intensity.P2a_3),

        P2b_1 = log2(Intensity.P2b_1),
        P2b_2 = log2(Intensity.P2b_2),
        P2b_3 = log2(Intensity.P2b_3),

        P3_1 = log2(Intensity.P3_1),
        P3_2 = log2(Intensity.P3_2),
        P3_3 = log2(Intensity.P3_3)
      ) %>%
      dplyr::select(
        P1_1,P1_2,P1_3,
        P2a_1,P2a_2, P2a_3,
        P2b_1, P2b_2, P2b_3,
        P3_1, P3_2 ,P3_3
      )
    is.na(htt_T_NoMatch) <- sapply(htt_T_NoMatch, is.infinite)

    htt_T_NoMatch <- read_ms("/htt/raw/total.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(
        P1_1 = log2(Intensity.P1_1),
        P1_2 = log2(Intensity.P1_2),
        P1_3 = log2(Intensity.P1_3),

        P2a_1 = log2(Intensity.P2a_1),
        P2a_2 = log2(Intensity.P2a_2),
        P2a_3 = log2(Intensity.P2a_3),

        P2b_1 = log2(Intensity.P2b_1),
        P2b_2 = log2(Intensity.P2b_2),
        P2b_3 = log2(Intensity.P2b_3),

        P3_1 = log2(Intensity.P3_1),
        P3_2 = log2(Intensity.P3_2),
        P3_3 = log2(Intensity.P3_3)
      ) %>%
      dplyr::select(
        P1_1,P1_2,P1_3,
        P2a_1,P2a_2, P2a_3,
        P2b_1, P2b_2, P2b_3,
        P3_1, P3_2 ,P3_3
      )
    is.na(htt_T_NoMatch) <- sapply(htt_T_NoMatch, is.infinite)

    htt_S_NoMatch <- read_ms("/htt/raw/super.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(
        P1_1 = log2(Intensity.P1_1),
        P1_2 = log2(Intensity.P1_2),
        P1_3 = log2(Intensity.P1_3),

        P2a_1 = log2(Intensity.P2a_1),
        P2a_2 = log2(Intensity.P2a_2),
        P2a_3 = log2(Intensity.P2a_3),

        P2b_1 = log2(Intensity.P2b_1),
        P2b_2 = log2(Intensity.P2b_2),
        P2b_3 = log2(Intensity.P2b_3),

        P3_1 = log2(Intensity.P3_1),
        P3_2 = log2(Intensity.P3_2),
        P3_3 = log2(Intensity.P3_3)
      ) %>%
      dplyr::select(
        P1_1,P1_2,P1_3,
        P2a_1,P2a_2, P2a_3,
        P2b_1, P2b_2, P2b_3,
        P3_1, P3_2 ,P3_3
      )
    is.na(htt_S_NoMatch) <- sapply(htt_S_NoMatch, is.infinite)

    htt_P_NoMatch <- read_ms("/htt/raw/pellet.csv") %>%
      filter(Only.identified.by.site != "+" &
               Reverse != "+" &
               Potential.contaminant != "+") %>%
      rename(gene = `Gene.names`) %>%
      mutate(
        P1_1 = log2(Intensity.P1_1),
        P1_2 = log2(Intensity.P1_2),
        P1_3 = log2(Intensity.P1_3),

        P2a_1 = log2(Intensity.P2a_1),
        P2a_2 = log2(Intensity.P2a_2),
        P2a_3 = log2(Intensity.P2a_3),

        P2b_1 = log2(Intensity.P2b_1),
        P2b_2 = log2(Intensity.P2b_2),
        P2b_3 = log2(Intensity.P2b_3),

        P3_1 = log2(Intensity.P3_1),
        P3_2 = log2(Intensity.P3_2),
        P3_3 = log2(Intensity.P3_3)
      ) %>%
      dplyr::select(
        P1_1,P1_2,P1_3,
        P2a_1,P2a_2, P2a_3,
        P2b_1, P2b_2, P2b_3,
        P3_1, P3_2 ,P3_3
      )
    is.na(htt_P_NoMatch) <- sapply(htt_P_NoMatch, is.infinite)
  }

  ## Fun -----------------------------------------------------------
  {
    len = 50

    genRepeat <- function(arr, den_range) {

      dat <- arr
      n <- dim(dat)[2]

      ## List of scatterplots
      scatter <- list()

      for (i in 2:n) {
        for (j in 1:(i-1)) {

          # Data frame
          df.point <- na.omit(data.frame(cbind(x = dat[ , j], y = dat[ , i])))

          # Plot
          p <- ggplot(df.point, aes(x, y)) +
            geom_point(colour = "grey",
                       alpha = .4,
                       size = .05) +
            theme_classic() +
            theme(
              strip.background = element_rect(colour = "white"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line=element_blank())

          name <- paste0("Item", j, i)
          scatter[[name]] <- p
        }
      }

      ## List of density plots
      den <- list()
      for(i in 1:n) {

        # Data frame
        den.df <- data.frame(dat[ , i])
        names(den.df) <- c("x")
        den.df <- den.df %>%
          filter(!is.na(x))

        # Plot
        p <- ggplot(den.df, aes(x)) +
          geom_density(linewidth = .2) +
          theme_classic() +
          theme(
            strip.background = element_rect(colour = "white"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line=element_blank()) +
          scale_x_continuous(limits = c(den_range[1],den_range[2]))

        name <- paste0("Item", i)
        den[[name]] <- p
      }

      ## List of tiles
      tile <- list()
      colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')

      for (i in 1:(n-1)) {
        for (j in (i+1):n) {

          # Data frame
          df.point <- na.omit(data.frame(cbind(x = dat[ , j], y = dat[ , i])))

          x = df.point[, 1]
          y = df.point[, 2]
          correlation = cor.test(x, y)
          cor <- data.frame(estimate = correlation$estimate,
                            statistic = correlation$statistic,
                            p.value = correlation$p.value)
          cor$cor = paste0("r = ", sprintf("%.2f", cor$estimate))

          fill <- colFn(100)[findInterval(cor$estimate, seq(0, 1, length=100))]

          # Plot
          p <- ggplot(cor, aes(x = 1, y = 1)) +
            geom_tile(fill = fill) +
            theme_classic() +
            theme(
              strip.background = element_rect(colour = "white"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line=element_blank())

          name <- paste0("Item", j, i)
          tile[[name]] <- p
        } }

      # Convert the ggplots to grobs,
      # and select only the plot panels
      denGrob <- lapply(den, ggplotGrob)
      denGrob <- lapply(denGrob, gtable_filter, "panel")

      scatterGrob <- lapply(scatter, ggplotGrob)
      scatterGrob <- lapply(scatterGrob, gtable_filter, "panel")

      tileGrob <- lapply(tile, ggplotGrob)
      tileGrob <- lapply(tileGrob, gtable_filter, "panel")


      ## Set up the gtable layout
      gt <- gtable(unit(rep(1, n), "null"), unit(rep(1, n), "null"))

      ## Add the plots to the layout
      # Bar plots along the diagonal
      for(i in 1:n) {
        gt <- gtable_add_grob(gt, denGrob[[i]], t=i, l=i)
      }

      # Scatterplots in the lower half
      k <- 1
      for (i in 2:n) {
        for (j in 1:(i-1)) {
          gt <- gtable_add_grob(gt, scatterGrob[[k]], t=i, l=j)
          k <- k+1
        }
      }

      # Tiles in the upper half
      k <- 1
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          gt <- gtable_add_grob(gt, tileGrob[[k]], t=i, l=j)
          k <- k+1
        }
      }

      gt

    }
  }
  ## Tunicamycin -----------------------------------------------------------
  {
    T_repeatability <- genRepeat(T_NoMatch[,3:14], c(13,32))
  }
  ## MG132 -----------------------------------------------------------
  {
    M_repeatability <- genRepeat(M_NoMatch[,3:14], c(13, 35))
  }

  ## PD -----------------------------------------------------------
  {
    PD_repeatability <- genRepeat(PD_NoMatch[,2:17], c(13, 32))
  }

  ## htt_T -----------------------------------------------------------
  {
    htt_T_repeatability <- genRepeat(htt_T_NoMatch[,1:12], c(16, 42))
  }

  ## htt_S -----------------------------------------------------------
  {
    htt_S_repeatability <- genRepeat(htt_S_NoMatch[,1:12], c(16, 42))
  }

  ## htt_P -----------------------------------------------------------
  {
    htt_P_repeatability <- genRepeat(htt_P_NoMatch[,1:12], c(17, 44))
  }
}