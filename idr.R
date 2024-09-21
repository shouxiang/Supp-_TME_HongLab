# Package loading ---------------------------------------------------------
{
  library(idpr)
  library(ggprism)
  library(cowplot)
  library(tidyverse)
  rm(list = ls())
}

# params ---------------------------------------------------------
{
  mouse_dir <- paste0(getwd(), '/idr_analysis/mouse')
  
  type_0 <- "no_cys"
  type_1 <- "structured"
  type_2 <- "cys_in_only_folded"
  type_3 <- "cys_in_only_disordered"
  type_4 <- "cys_in_both"
  type_5 <- "idp"
}

# utils-----------------------------------------------------------------
{
  identify_disordered_regions <- function(df, disordered_length) {
    fullLength <- paste0(df$AA,collapse = "")
    # Identify disordered positions
    disordered <- df$IUPred2 > 0.5
    
    # Initialize variables to track regions
    regions <- list()
    current_region <- NULL
    
    # Loop through the disordered vector to identify continuous regions
    for (i in seq_along(disordered)) {
      if (disordered[i]) {
        if (is.null(current_region)) {
          # Start a new region
          current_region <- list(df$Position[i], df$Position[i])
        } else {
          # Extend the current region
          current_region[[2]] <- df$Position[i]
        }
      } else {
        if (!is.null(current_region)) {
          if(current_region[[2]] - current_region[[1]] + 1 >= disordered_length){
            # End the current region and store it as a list of start and end positions
            regions <- c(regions, list(current_region))
          }
          current_region <- NULL
        }
      }
    }
    
    # Add the last region if it ended on the last position
    if (!is.null(current_region)) {
      if(current_region[[2]] - current_region[[1]] + 1 >= disordered_length){
        # End the current region and store it as a list of start and end positions
        regions <- c(regions, list(current_region))
      }
    }
    
    disordered_aa <- list()
    
    for (region in regions){
      aa_seq <- substr(fullLength, region[[1]], region[[2]])
      disordered_aa <- c(disordered_aa, list(aa_seq))
    }
    
    disordered_aa_string <- paste(disordered_aa, collapse = " ")
    
    return(disordered_aa_string)
  }
  
  analyse_proteome <- function(env_path, disordered_threshold, save_path,  save_name) {
    
    env <- readRDS(env_path)
    keys <- ls(envir = env)
    
    results <- list()
    
    for (key in keys){
      
      fullLength <- paste0(env[[key]]$AA,collapse = "")
      
      disordered_regions <- 
        identify_disordered_regions(env[[key]], disordered_threshold)
      
      cys_total <- stringr::str_count(fullLength, "C")
      
      cys_in_disorder <- 0
      
      for (region in disordered_regions) {
        cys_in_disorder <- cys_in_disorder + stringr::str_count(region, "C")
      }
      
      aa_if_disordered <- env[[key]]$IUPred2 > 0.5
      
      percent_disordered <- sum(aa_if_disordered) / length(aa_if_disordered)
      
      type <- ""
      if(cys_total == 0){
        type <- type_0
      } else{
        if(percent_disordered < 0.05){
          type <- type_1
        }
        else if(percent_disordered <=0.8){
          if(cys_in_disorder == cys_total){
            type <- type_3
          }
          else if(cys_in_disorder == 0){
            type <- type_2
          }
          else{
            type <- type_4
          }
        }
        else{
          type <- type_5
        }
      }
      
      results[[key]] <- list(
        proteinID = key,
        fullLength = fullLength,
        disordered_regions = disordered_regions,
        cys_total = cys_total,
        cys_in_disorder = cys_in_disorder,
        percent_disordered =percent_disordered,
        type = type
      )
    }
    
    results_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
    row.names(results_df) <- 1:nrow(results_df)
    
    saveRDS(
      results_df,
      file = file.path(
        save_path,
        save_name))
  }
  
  get_aa_composition <- function(vector){
    
    common_amino_acids <- 
      c("A", "R", "N", "D", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "C")
    
    
    all_sequences <- paste(vector, collapse = "")
    filtered_sequences <- strsplit(all_sequences, NULL)[[1]]
    filtered_sequences <- filtered_sequences[filtered_sequences %in% common_amino_acids]
    
    
    amino_acid_counts <- table(filtered_sequences)
    
    amino_acid_counts_df <- as.data.frame(amino_acid_counts, stringsAsFactors = FALSE)
    colnames(amino_acid_counts_df) <- c("aa", "count")
    
    total_amino_acid_count <- sum(amino_acid_counts_df$count)
    
    amino_acid_counts_df$composition <- amino_acid_counts_df$count / total_amino_acid_count
    
    return(amino_acid_counts_df)
  }
  
  fisher <- function(cat_in, cat_total, proteome_cat, proteome_total){
    
    structured_aa <- proteome_cat - cat_in
    disordered_not_aa <- cat_total - cat_in
    structured_not_aa <- proteome_total - cat_in - structured_aa - disordered_not_aa
    
    contingency_table <- matrix(
      c(cat_in, disordered_not_aa, structured_aa, structured_not_aa), 
      nrow = 2, byrow = TRUE)
    
    dimnames(contingency_table) <- list(
      disordered = c("in", "not-in"),
      aa = c("aa", "not-this-aa")
    )
    
    fisher_test_result <- fisher.test(contingency_table)
    
    return(fisher_test_result$p.value)
  }
  
  gen_disorder_type_plot <- function(base, df1, df2, df1_name, df2_name, color, title){
    
    df1 <- dplyr::left_join(
      df1, base)
    
    df2 <- dplyr::left_join(
      df2, base)
    
    
    y.position <- c(rep(0.1, 5), rep(0.2, 5))
    xmin <- rep((1:5) - 0.3, 2)
    xmax <- c((1:5), (1:5) + 0.3)
    
    base_count <- c(
      sum(base$type==type_1),
      sum(base$type==type_2),
      sum(base$type==type_3),
      sum(base$type==type_4),
      sum(base$type==type_5)
    )
    
    base_total <- sum(base_count)
    base_ratio <- base_count/base_total
    
    
    df1_count <- c(
      sum(df1$type==type_1),
      sum(df1$type==type_2),
      sum(df1$type==type_3),
      sum(df1$type==type_4),
      sum(df1$type==type_5)
    )
    df1_total <- sum(df1_count)
    df1_ratio <- df1_count/df1_total
    
    
    
    df2_count <- c(
      sum(df2$type==type_1),
      sum(df2$type==type_2),
      sum(df2$type==type_3),
      sum(df2$type==type_4),
      sum(df2$type==type_5)
    )
    df2_total <- sum(df2_count)
    df2_ratio <- df2_count/df2_total
    
    type <- paste0("Type ", 1:5)
    
    group1 <- rep("Proteome", 5)
    group2 <- rep(df1_name, 5)
    foldEnrichment <- df1_ratio/base_ratio
    pvalue <- c(
      fisher(df1_count[[1]], df1_total, base_count[[1]], base_total),
      fisher(df1_count[[2]], df1_total, base_count[[2]], base_total),
      fisher(df1_count[[3]], df1_total, base_count[[3]], base_total),
      fisher(df1_count[[4]], df1_total, base_count[[4]], base_total),
      fisher(df1_count[[5]], df1_total, base_count[[5]], base_total)
    )
    pstar <- ifelse(pvalue < 0.0001, "****",
                    ifelse(pvalue < 0.001, "***",
                           ifelse(pvalue < 0.01, "**",
                                  ifelse(pvalue < 0.05, "*", "ns"))))
    label <- sprintf("%.2f, %s", foldEnrichment, pstar)
    
    stat.pair <- data.frame(
      group1,
      group2,
      annotation = type,
      label
    )
    
    group1 <- rep("Proteome", 5)
    group2 <- rep(df2_name, 5)
    foldEnrichment <- df2_ratio/base_ratio
    pvalue <- c(
      fisher(df2_count[[1]], df2_total, base_count[[1]], base_total),
      fisher(df2_count[[2]], df2_total, base_count[[2]], base_total),
      fisher(df2_count[[3]], df2_total, base_count[[3]], base_total),
      fisher(df2_count[[4]], df2_total, base_count[[4]], base_total),
      fisher(df2_count[[5]], df2_total, base_count[[5]], base_total)
    )
    pstar <- ifelse(pvalue < 0.0001, "****",
                    ifelse(pvalue < 0.001, "***",
                           ifelse(pvalue < 0.01, "**",
                                  ifelse(pvalue < 0.05, "*", "ns"))))
    label <- sprintf("%.2f, %s", foldEnrichment, pstar)
    
    
    stat.pair <- rbind(stat.pair, data.frame(
      group1,
      group2,
      annotation = type,
      label
    ))
    
    stat.pair <- stat.pair %>%
      dplyr::mutate(stress = group2) %>%
      dplyr::mutate(y.position = y.position) %>%
      dplyr::mutate(xmin = xmin) %>%
      dplyr::mutate(xmax = xmax)
    
    r1 <- data.frame(
      stress = c(rep('Proteome', 5), rep(df1_name, 5), rep(df2_name, 5)),
      annotation = rep(type, 3),
      percentage = c(base_ratio, df1_ratio, df2_ratio)
    )
    
    r1$stress <- factor(r1$stress, levels = c("Proteome", df1_name, df2_name))
    
    p <- 
      ggplot(r1, aes(x = annotation, y = percentage, fill = stress)) +
      geom_bar(position = "dodge", stat = "identity",
               aes(fill = stress)) +
      scale_fill_manual(values = color) +
      scale_x_discrete() +
      scale_y_continuous(labels = scales::percent, limits = c(0,0.5)) +
      ylab("Percent of proteins") +
      scale_color_manual(labels = c("Proteome", df1_name, df2_name)) +
      stat_pvalue_manual(
        stat.pair, 
        label = "label",
        tip.length = 0.01,
        size = 5/.pt,
        bracket.size = 0.1,
        vjust = -.25
      ) +
      labs(title = title) +
      THEME0 +
      theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(.7, .8),
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(3, 'mm'),
        legend.title = element_blank())
    p
  }
}

# Load data -----------------------------------------------------------------
{
  T_ctrBinder <- 
    read_ms("/stressor/binder/T_ctrBinder.csv")
  
  T_tuniBinder <- 
    read_ms("/stressor/binder/T_tuniBinder.csv")
  
  M_ctrBinder <- 
    read_ms("/stressor/binder/M_ctrBinder.csv")
  
  M_mgBinder <- 
    read_ms("/stressor/binder/M_mgBinder.csv")
  
  mouse_result_15 <- readRDS(paste0(mouse_dir, '/mouse_result_15.rds'))
  mouse_result_30 <- readRDS(paste0(mouse_dir, '/mouse_result_30.rds'))
  
  T_unfold <- 
    dplyr::inner_join(
      read_ms("/stressor/processed/T_TD_TE_imputed.csv") %>% 
        dplyr::filter(Cys > 0) %>%
        dplyr::filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
        dplyr::select(proteinID, geneName),
      read_ms("/stressor/processed/T_CE_TE_imputed.csv") %>% 
        dplyr::filter(Cys > 0) %>%
        dplyr::filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
        dplyr::select(proteinID, geneName)) %>%
    dplyr::left_join(
      mouse_result_15 %>% dplyr::select(proteinID, type)) %>%
    dplyr::rename(type_15 = type) %>%
    dplyr::left_join(
      mouse_result_30 %>% dplyr::select(proteinID, type)) %>%
    dplyr::rename(type_30 = type)
  
  M_unfold <- 
    inner_join(
      read_ms("/stressor/processed/M_MD_ME_imputed.csv") %>% 
        filter(Cys > 0) %>%
        filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
        select(proteinID, geneName),
      read_ms("/stressor/processed/M_CE_ME_imputed.csv") %>% 
        filter(Cys > 0) %>%
        filter(p.adjust < 0.05 & log2fc > log2(1.5)) %>%
        select(proteinID, geneName)) %>%
    dplyr::left_join(
      mouse_result_15 %>% dplyr::select(proteinID, type)) %>%
    dplyr::rename(type_15 = type) %>%
    dplyr::left_join(
      mouse_result_30 %>% dplyr::select(proteinID, type)) %>%
    dplyr::rename(type_30 = type)
}


# mouse AA composition -----------------------------------------------------------------
{
  
  mouse_proteome_composition <- get_aa_composition(mouse_result_30$fullLength) %>%
    dplyr::mutate(group = "proteome")
  mouse_disordered_15_composition <- get_aa_composition(mouse_result_15$disordered_regions) %>%
    dplyr::mutate(group = "disordered_15")
  mouse_disordered_30_composition <- get_aa_composition(mouse_result_30$disordered_regions) %>%
    dplyr::mutate(group = "disordered_30")
  
  mouse_aa_composition <- rbind(mouse_proteome_composition, mouse_disordered_15_composition) %>%
    rbind(mouse_disordered_30_composition)
  
  group_order <- c("proteome", "disordered_15", "disordered_30")
  
  COLOR_AA_COMPOSITION <- c(
    proteome = "#AAACAF", 
    disordered_15 = "#6E9ECE", 
    disordered_30 = "#E6928F")
  
  mouse_aa_composition$group <- factor(
    mouse_aa_composition$group, levels = group_order)
  
  aa_order <- 
    c("S", "P", "E", "A", "G", "R", "K", "Q", "D", "T",
      "H", "L", "V", "I", "F", "N", "M", "Y", "W", "C")
  
  
  mouse_aa_composition$aa <- factor(mouse_aa_composition$aa, levels = aa_order)

}