
# Package loading ---------------------------------------------------------
{
  library(ggprism)
  library(circlize)
  library(ComplexHeatmap)
  library(ggraph)
  library(tidygraph)
  library(magick)
  library(tidyverse)
  library(ggpubr)
  rm(list = ls())
}

# Color setting -----------------------------------------------------------

{
  color_stress <- c(
    Proteome = "#AAACAF",
    Descending = "#6E9ECE",
    Ascending = "#E6928F"
  )
  
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
  
 
  
  tTest <- function(g0, g1){
    
    v <- var.test(g0, g1)
    
    if(v$p.value <= .05) {
      p <- t.test(g0, g1, var.equal = FALSE)$p.value
    } else {
      p <- t.test(g0, g1, var.equal = TRUE)$p.value
    }
    p
  }
  
  znorm <- function(df) {
    
    for (i in 1:nrow(df)) {
      
      tmp <- unlist(df[i,])
      
      result <- (tmp-mean(tmp))/sd(tmp)
      
      if(i == 1){
        a <- as.matrix(result)
      } else{
        a <- cbind(a, result)
      }
    }
    
    b <- t(a)
    rownames(b) <- rownames(df)
    b
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
  
  populationName <- c("P1", "P2a", "P2b", "P3")
  colNames <- paste0(rep(populationName, each = 3), rep(1:3))
  
  # heatmap
  
  {
    pellet_P2a_P2b <- read("htt/processed/pellet_P2a_P2b.csv")
    
    pellet_heatmap <- read("htt/processed/pellet_heatmap.csv")
    rownames(pellet_heatmap) <- pellet_heatmap$proteinID
    
    pellet.norm <- znorm(pellet_heatmap %>% dplyr::select(all_of(colNames)))
    
    pellet_row_group <- factor(
      pellet_heatmap$group, levels = c("high1", "high3"))
    names(pellet_row_group) <- pellet_heatmap$proteinID
    
    pellet_col_group <- factor(
      rep(populationName, each = 3), levels = c("P1", "P2a", "P2b", "P3"))
    names(pellet_col_group) <- colNames
  }
  
  {
    total_P2a_P2b <- read("htt/processed/total_P2a_P2b.csv")
    
    total_heatmap <- read("htt/processed/total_heatmap.csv")
    rownames(total_heatmap) <- total_heatmap$proteinID
    
    total.norm <- znorm(total_heatmap %>% dplyr::select(all_of(colNames)))
    
    total_row_group <- factor(
      total_heatmap$group, levels = c("high1", "high3"))
    names(total_row_group) <- total_heatmap$proteinID
    
    total_col_group <- factor(
      rep(populationName, each = 3), levels = c("P1", "P2a", "P2b", "P3"))
    names(total_col_group) <- colNames
  }
  
  {
    super_P2a_P2b <- read("htt/processed/super_P2a_P2b.csv")
    
    super_heatmap <- read("htt/processed/super_heatmap.csv")
    rownames(super_heatmap) <- super_heatmap$proteinID
    
    super.norm <- znorm(super_heatmap %>% dplyr::select(all_of(colNames)))
    
    super_row_group <- factor(
      super_heatmap$group, levels = c("high1", "high3"))
    names(super_row_group) <- super_heatmap$proteinID
    
    super_col_group <- factor(
      rep(populationName, each = 3), levels = c("P1", "P2a", "P2b", "P3"))
    names(super_col_group) <- colNames
  }
  
  # heatmap enrichment
  T_ctrBinder_enrich <- read("stressor/enrichment/T_ctrBinder_enrich.csv")
  pellet_high3_enrich <- read("htt/enrichment/pellet_high3_enrich.csv")
  pellet_high1_enrich <- read("htt/enrichment/pellet_high1_enrich.csv")
  
  total_high3_enrich <- read("htt/enrichment/total_high3_enrich.csv")
  total_high1_enrich <- read("htt/enrichment/total_high1_enrich.csv")
  
  super_high3_enrich <- read("htt/enrichment/super_high3_enrich.csv")
  super_high1_enrich <- read("htt/enrichment/super_high1_enrich.csv")
  
  {
    ID_funCat_pellet_high3 <- list( #use dput
      
      ID = c("R-MMU-8854050", "R-MMU-69242", "R-MMU-69206", 
             "R-MMU-69202", "R-MMU-69656", "R-MMU-69580", "R-MMU-69275", "R-MMU-69613", 
             "R-MMU-453276", "GO:0045182", "GO:0090079", "GO:0051087", "GO:0031072", 
             "GO:0016234", "GO:0000502", "mmu03050", "GO:0036464", "GO:0043248", 
             "GO:0090083", "GO:0051085", "GO:0070841", "GO:0061077", "GO:0051082", 
             "GO:0005776", "GO:0003697", "R-MMU-2262752", "R-MMU-8953897", 
             "GO:0003729", "GO:0005525"),
      
      funCat = c("Cell cycle", "Cell cycle", "Cell cycle", "Cell cycle", 
                 "Cell cycle", "Cell cycle", "Cell cycle", "Cell cycle", "Cell cycle", 
                 "Translation initiation and regulation", 
                 "Nucleotide binding", "Quality control", "Quality control", "Quality control", 
                 "Quality control", "Quality control", "Biomolecular condensates",
                 "Quality control",
                 "Quality control",
                 "Quality control",
                 "Quality control",
                 "Quality control",
                 "Quality control",
                 "Quality control",
                 "Nucleotide binding",
                 "Quality control",
                 "Quality control",
                 "mRNA processing",
                 "Nucleotide binding")
      
    ) %>% as.data.frame()
    
    ID_funCat_pellet_high3$funCat <- 
      factor(ID_funCat_pellet_high3$funCat,
             levels = c("Cell cycle",
                        "Translation initiation and regulation",
                        "mRNA processing",
                        "Nucleotide binding",
                        "Quality control",
                        "Biomolecular condensates"))
    
    cat_color_pellet_high3 = c(
      `Cell cycle` = "#ef7770", 
      `Translation initiation and regulation` = "#bf9a14", 
      `mRNA processing` = "#60b304",
      `Nucleotide binding` = "#40bf93",
      # `Stress response` = "#42b6ea",
      `Quality control` = "#a58bfe",
      `Biomolecular condensates` = "#f264d7"
    )
  }
  
  {
    ID_funCat_total_high3 <- list( #use dput
      
      ID = c("GO:0005838", "GO:0043248", "GO:0000502", 
             "mmu03050", "GO:0070628", "R-MMU-69610", "R-MMU-69613", "R-MMU-8932339", 
             "R-MMU-69541", "R-MMU-69563", "R-MMU-69580", "R-MMU-1234174", 
             "R-MMU-69206", "R-MMU-9711123", "GO:0061077", "GO:0006413", "GO:0010498", 
             "GO:0045182", "R-MMU-2262752", "R-MMU-8953897", "GO:0034976"),
      
      funCat = c("Quality control", "Quality control", "Quality control", 
                 "Quality control", "Quality control", "Stress response", 
                 "Stress response", "Stress response", "Stress response", 
                 "Stress response", "Stress response", "Stress response", 
                 "Cell cycle", "Stress response", "Quality control", "Translation initiation and regulation", 
                 "Quality control", "Translation initiation and regulation", 
                 "Stress response", "Stress response", "Stress response")
      
    ) %>% as.data.frame()
    
    ID_funCat_total_high3$funCat <- 
      factor(ID_funCat_total_high3$funCat,
             levels = c("Stress response",
                        "Quality control",
                        "Translation initiation and regulation",
                        "Cell cycle"))
    
    cat_color_total_high3 = c(
      
      `Stress response` = "#42b6ea",
      `Quality control` = "#a58bfe",
      `Translation initiation and regulation` = "#bf9a14", 
      `Cell cycle` = "#ef7770"
    )
  }
  
  {
    ID_funCat_super_high3 <- list( #use dput
      
      ID = c("GO:0005761", "GO:0003735", "R-MMU-5389840", 
             "GO:0044391", "R-MMU-5419276", "R-MMU-5368287", "R-MMU-72766", 
             "GO:0032543", "R-MMU-72649", "R-MMU-72702", "R-MMU-72613", "R-MMU-72737", 
             "GO:0005844", "R-MMU-72312", "R-MMU-8868773", "GO:0006364", "GO:0005681", 
             "GO:0008380"),
      
      funCat = c("Ribosome", "Ribosome", "Translation initiation and regulation", 
                 "Ribosome", "Translation initiation and regulation", "Translation initiation and regulation", 
                 "Translation initiation and regulation", "Translation initiation and regulation", 
                 "Translation initiation and regulation", "Translation initiation and regulation", 
                 "Translation initiation and regulation", "Translation initiation and regulation", 
                 "Translation initiation and regulation", "RNA processing", "RNA processing", 
                 "RNA processing", "RNA processing", "RNA processing")
      
    ) %>% as.data.frame()
    
    ID_funCat_super_high3$funCat <- 
      factor(ID_funCat_super_high3$funCat,
             levels = c("Translation initiation and regulation",
                        "RNA processing",
                        "Ribosome"))
    
    cat_color_super_high3 = c(
      
      `Translation initiation and regulatio` = "#42b6ea",
      `RNA processing` = "#a58bfe",
      `Ribosome` = "#bf9a14"
    )
  }
  
  #volcano enrichment
  pellet_P2a_P2b_right_enrich <- read("htt/enrichment/pellet_P2a_P2b_right_enrich.csv") %>% modifyFunctionalTerms()
  pellet_P2a_P2b_left_enrich <- read("htt/enrichment/pellet_P2a_P2b_left_enrich.csv") %>% modifyFunctionalTerms()
  
}



# Venn cluster D and A in P, in T and in S-----------------------------------------------------

{
  
  {## pellet ------------------------------------------------------------------
    lt_pellet_13 =
      list(
        TctrBinder = T_ctrBinder_enrich$ID,
        high1 = pellet_high1_enrich$ID,
        high3 = pellet_high3_enrich$ID
      )
    
    plot_venn0 <- 
      ggVennDiagram::ggVennDiagram(
        lt_pellet_13, 
        edge_size = .1, 
        set_size = 6/.pt,
        label = "count", label_alpha = 0, label_size = 6/.pt,
        
        category.names = c("TME Binder",
                           "P-cluster-D",
                           "P-cluster-A")) +
      scale_fill_gradient(
        low = "#F4FAFE", high = "#4981BF",
        limits = c(0, 100),
        breaks = c(0, 50, 100),
        name = "Count") +
      scale_color_manual(values = c("black","black","black"))
    
    legend_venn <- get_legend(
      plot_venn0 +
        theme(
          legend.key.width = unit(3, 'mm'),
          legend.key.height = unit(2.5, 'mm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5)
        )
    )
    
    plot_venn_pellet_13 <- 
      ggdraw(plot_venn0 + theme(legend.position = "none")) +
      draw_plot(
        legend_venn,
        x = .95,
        y = .95,
        width = .1,
        height = .3,
        hjust = 1,
        vjust = 1
      )
  }
  
  {## total -------------------------------------------------------------------
    lt_total_13 =
      list(
        TctrBinder = T_ctrBinder_enrich$ID,
        high1 = total_high1_enrich$ID,
        high3 = total_high3_enrich$ID
      )
    
    plot_venn0 <- 
      ggVennDiagram::ggVennDiagram(
        lt_total_13, 
        edge_size = .3, 
        set_size = 6/.pt,
        label = "count", label_alpha = 0, label_size = 6/.pt,
        
        category.names = c("TME Binder",
                           "T-cluster-D",
                           "T-cluster-A")) +
      scale_fill_gradient(
        low = "#F4FAFE", high = "#4981BF",
        limits = c(0, 100),
        breaks = c(0, 50, 100),
        name = "Count") +
      scale_color_manual(values = c("black","black","black"))
    
    legend_venn <- get_legend(
      plot_venn0 +
        theme(
          legend.key.width = unit(3, 'mm'),
          legend.key.height = unit(2.5, 'mm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5)
        )
    )
    
    plot_venn_total_13 <- 
      ggdraw(plot_venn0 + theme(legend.position = "none")) +
      draw_plot(
        legend_venn,
        x = .95,
        y = .95,
        width = .1,
        height = .3,
        hjust = 1,
        vjust = 1
      )
    
    # ggsave(filename = paste0("venn-",
    #                          format(Sys.time(), format = "%Y%m%d-%H%M%S") %>%
    #                            as.character(),
    #                          ".pdf"),
    #        plot = plot_venn_total_13, path = "C:/Users/cc/Desktop",
    #        width = 45, height = 45, units = c("mm"))
    
  }
  
  {## super -------------------------------------------------------------------

    lt_super_13 =
      list(
        TctrBinder = T_ctrBinder_enrich$ID,
        high1 = "",
        high3 = super_high3_enrich$ID
      )
    
    plot_venn0 <- 
      ggVennDiagram::ggVennDiagram(
        lt_super_13, 
        edge_size = .3, 
        set_size = 6/.pt,
        label = "count", label_alpha = 0, label_size = 6/.pt,
        
        category.names = c("TME Binder",
                           "S-cluster-D",
                           "S-cluster-A")) +
      scale_fill_gradient(
        low = "#F4FAFE", high = "#4981BF",
        limits = c(0, 200),
        breaks = c(0, 100, 200),
        name = "Count") +
      scale_color_manual(values = c("black","black","black"))
    
    legend_venn <- get_legend(
      plot_venn0 +
        theme(
          legend.key.width = unit(3, 'mm'),
          legend.key.height = unit(2.5, 'mm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5)
        )
    )
    
    plot_venn_super_13 <- 
      ggdraw(plot_venn0 + theme(legend.position = "none")) +
      draw_plot(
        legend_venn,
        x = .95,
        y = .95,
        width = .1,
        height = .3,
        hjust = 1,
        vjust = 1
      )
    
    # ggsave(filename = paste0("venn-",
    #                          format(Sys.time(), format = "%Y%m%d-%H%M%S") %>%
    #                            as.character(),
    #                          ".pdf"),
    #        plot = plot_venn_super_13, path = "C:/Users/cc/Desktop",
    #        width = 45, height = 45, units = c("mm"))
  }
  
} 


# Venn T, S and P in A ----------------------------------------------

{
  {## in cluster A ------------------------------------------------------------------
    
    lt_cluserD_TSP =
      list(
        Total = total_high3_enrich$ID,
        Supernatant = super_high3_enrich$ID,
        Pellet = pellet_high3_enrich$ID
      )
    
    plot_venn0 <- 
      ggVennDiagram::ggVennDiagram(
        lt_cluserD_TSP, 
        edge_size = .3, 
        set_size = 6/.pt,
        label = "count", label_alpha = 0, label_size = 6/.pt,
        category.names = c("T-cluster-A",
                           "S-cluster-A",
                           "P-cluster-A")) +
      scale_fill_gradient(
        low = "#F4FAFE", high = "#4981BF",
        limits = c(0, 150),
        breaks = c(0, 75, 150),
        name = "Count") +
      scale_color_manual(values = c("black","black","black"))
    
    legend_venn <- get_legend(
      plot_venn0 +
        theme(
          legend.key.width = unit(3, 'mm'),
          legend.key.height = unit(2.5, 'mm'),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5)
        )
    )
    
    plot_venn_clusterD_TSP <- 
      ggdraw(plot_venn0 + theme(legend.position = "none")) +
      draw_plot(
        legend_venn,
        x = .95,
        y = .95,
        width = .1,
        height = .3,
        hjust = 1,
        vjust = 1
      )
    
    # ggsave(filename = paste0("venn-",
    #                          format(Sys.time(), format = "%Y%m%d-%H%M%S") %>%
    #                            as.character(),
    #                          ".pdf"),
    #        plot = plot_venn_clusterD_TSP, path = "C:/Users/cc/Desktop",
    #        width = 45, height = 45, units = c("mm"))
  }
}

# Htt abundance in 4 populations ------------------------------------------
{
  {## pellet ------------------------------------------------------------------

    htt_pellet <- pellet_heatmap %>% 
      filter(proteinID == "ZZZZZ2") %>%
      dplyr::select(all_of(colNames), proteinID) %>%
      pivot_longer(
        all_of(colNames),
        names_to = "tmp",
        values_to = "abundance") %>%
      separate(
        tmp,
        into = c("population", "rep"),
        sep = -1)
    
    htt_pellet$population <- factor(
      htt_pellet$population,
      levels = c("P1", "P2a", "P2b", "P3"))
    
    plot_htt_pellet <-
      ggplot(htt_pellet, aes(population, abundance)) +
      geom_point(
        size = .5,
        position = position_jitter(w = 0.1, h = 0)) +
      scale_x_discrete(labels = c("P1" = "1", "P2a" = "n", "P2b" = "a", "P3" = "3")) +
      scale_y_continuous(
        limits = c(30, 40),
        breaks = c(30,35,40),
        guide = "prism_minor") +
      xlab("") +
      ylab(expression(log[2]*"(46Q abundance) (a.u.)")) +
      ggtitle("Pellet") +
      theme(axis.ticks.x = element_blank())
  }
  
  {## total ------------------------------------------------------------------
    
    htt_total <- total_heatmap %>% 
      filter(proteinID == "ZZZZZ2") %>%
      dplyr::select(all_of(colNames), proteinID) %>%
      pivot_longer(
        all_of(colNames),
        names_to = "tmp",
        values_to = "abundance") %>%
      separate(
        tmp,
        into = c("population", "rep"),
        sep = -1)
    
    htt_total$population <- factor(
      htt_total$population,
      levels = c("P1", "P2a", "P2b", "P3"))
    
    plot_htt_total <-
      ggplot(htt_total, aes(population, abundance)) +
      geom_point(
        size = .5,
        position = position_jitter(w = 0.1, h = 0)) +
      scale_x_discrete(labels = c("P1" = "1", "P2a" = "n", "P2b" = "a", "P3" = "3")) +
      scale_y_continuous(
        limits = c(28, 40),
        breaks = c(30,35, 40),
        guide = "prism_minor") +
      xlab("") +
      ylab(expression(log[2]*"(46Q abundance) (a.u.)")) +
      ggtitle("Total cell lysate") +
      theme(axis.ticks.x = element_blank())
  }
  
  {## super ------------------------------------------------------------------
    
    htt_super <- super_heatmap %>% 
      filter(proteinID == "ZZZZZ2") %>%
      dplyr::select(all_of(colNames), proteinID) %>%
      pivot_longer(
        all_of(colNames),
        names_to = "tmp",
        values_to = "abundance") %>%
      separate(
        tmp,
        into = c("population", "rep"),
        sep = -1)
    
    htt_super$population <- factor(
      htt_super$population,
      levels = c("P1", "P2a", "P2b", "P3"))
    
    plot_htt_super <-
      ggplot(htt_super, aes(population, abundance)) +
      geom_point(
        size = .5,
        position = position_jitter(w = 0.1, h = 0)) +
      scale_x_discrete(labels = c("P1" = "1", "P2a" = "n", "P2b" = "a", "P3" = "3")) +
      scale_y_continuous(
        limits = c(27, 35),
        breaks = c(30,35),
        guide = "prism_minor") +
      xlab("") +
      ylab(expression(log[2]*"(46Q abundance) (a.u.)")) +
      ggtitle("Supernatant") +
      theme(axis.ticks.x = element_blank())
  }
}
