## ---- garnish_plot
#' Graph summary results using ggplot2.
#'
#' Plots ADN, CDN, priority, frameshift, and fusion derived neoepitopes for Class I and Class II MHC by sample.
#' See `garnish_summary` for more information on neoepitope classification by antigen.garnish.
#'
#' @param dt Data table. Prediction output from `garnish_predictions`. Either dt or jdt must be provided as input.
#' @param jdt Data table. Prediction output from `garnish_jaffa`` %>% `garnish_predictions`. Either dt or jdt must be provided as input.
#' @param save Optional: Should the plot be saved to the working directory as antigen.garnish_summary.pdf? Default is TRUE.
#'
#' @examples
#'\dontrun{
#'library(magrittr)
#'library(antigen.garnish)
#'
#'  # download an example VCF
#'    g <- "antigen.garnish_example.vcf" %T>%
#'    utils::download.file("http://get.rech.io/antigen.garnish_example.vcf", .) %>%
#'
#'  # extract variants
#'    antigen.garnish::garnish_variants %>%
#'
#'  # add test MHC types
#'      .[, MHC := c("HLA-A*02:01 HLA-DRB1*14:67",
#'                   "H-2-Kb H-2-IAd",
#'                  "HLA-A*01:47 HLA-DRB1*03:07")] %>%
#'
#'  # predict neoepitopes
#'    antigen.garnish::garnish_predictions %>%
#'
#'  # plot it
#'    antigen.garnish::garnish_plot
#'   
#'  # here is our ggplot object that can be further customized 
#'   class(g)
#'  
#'  # save the plot with a custom ggplot theme 
#'  g + my_ggplot_custom_theme_object %>% cowplot::ggsave("antigen.garnish_custom_themed_plot.pdf")
#'   
#'  # show location of where default plot was saved
#'   list.files(pattern = "summary\\.pdf", full.names = TRUE)   
#'   
#'}
#'
#' @return
#'
#' A list of ggplot2 objects and, if save = TRUE, saves each of these objects pdfs in the working directory.
#' The first plot is a column graph indicating the number of peptides in each sample classified as ADN, CDN, priority, frameshift-derived, and fusion-derived neoepitopes
#' (if any) for each sample, faceted by Class I vs Class II MHC. The threshold for fusion and frameshift derived neoepitopes is Consensus_scores < 1000nM.
#' See `garnish_summary`` documentation for further explanation of neoepitope classification of ADN, CDN, and priority.
#' 
#' The second plot and third plots are returned only if frameshift and/or fusion mutants (respectively) are present in the input table. These plot
#' the number of peptides per sample by MHC class I or II and binned by binding affinity (1000 - 500nM, 500 - 50nM, and < 50 nM). 
#'
#' @export garnish_plot
#'
#' @md

garnish_plot <- function(dt = NULL, jdt = NULL, save = TRUE){
 
  ag_gg_theme <-
    ggplot2::theme(line = ggplot2::element_line(colour = "#000000")) +
    ggplot2::theme(axis.line = ggplot2::element_line(color="black")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(2*1.2),
                                                      colour = "#000000", lineheight = 0.9, face = "bold", vjust = 0)) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "#000000",
                                                        size = ggplot2::rel(2*1.425), face = "bold")) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(colour = "#000000",
                                                        size = ggplot2::rel(2*1.425), face = "bold", angle = 90, vjust = 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "#000000",
                                                       size = ggplot2::rel(2*1.425), face = "bold", angle = 30, hjust = 0.9, vjust = 0.92)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "#000000",
                                                       size = ggplot2::rel(2*1.425), face = "bold", hjust = 0.9, vjust = 0.92, angle = 30)) +
    ggplot2::theme(legend.text = ggplot2::element_text(colour = "#000000",
                                                       size = ggplot2::rel(2*1))) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#eeeeee")) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold")) +
    ggplot2::theme(legend.key = ggplot2::element_blank()) +
    ggplot2::theme(legend.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::theme(plot.margin=grid::unit(c(0, 0, 0, 0 ), "cm"))
  
  ag_colors <-   c("#ff80ab",
                   "#b388ff",
                   "#82b1ff",
                   "#a7ffeb",
                   "#b9f6ca",
                   "#f4ff81",
                   "#ffe57f",
                   "#ff9e80",
                   "#ff8a80",
                   "#ea80fc",
                   "#8c9eff",
                   "#80d8ff",
                   "#84ffff",
                   "#ccff90",
                   "#ffff8d",
                   "#ffd180")
  
  ##check function input
  if (missing(dt) & missing(jdt)) stop("at least one dt from garnish_predictions must be provided.")
  
  if (!missing(dt) & !missing(jdt)) dt <- cat_tables(dt, jdt)
  
  if(!missing(dt)){
  
  if (!(c("nmer",
          "MHC",
          "sample_id",
          "DAI",
          "Consensus_scores") %chin% names(dt)) %>% any)
    stop("'sample_id', 'nmer', 'MHC', 'frameshift', 'DAI', and 'Consensus_scores' columns are required in dt")
  
  dt <- dt[pep_type != "wt"] %>% unique(by = c("nmer",
                                               "MHC",
                                               "sample_id",
                                               "DAI",
                                               "Consensus_scores"))
  
  dt <- dt %>% .[Consensus_scores < 5000] %>%
    .[pep_type == "mutnfs" & Consensus_scores < 50, type := "CDN"] %>%
      .[pep_type == "mutnfs" & DAI > 10, type := "ADN"] %>%
        .[pep_type == "mutnfs" & Consensus_scores < 50 & DAI > 10, type := "priority"] %>%
          .[pep_type == "fus" & Consensus_scores < 1000, type := "fusion"] %>%
            .[pep_type == "mut_other" & Consensus_scores < 1000, type := "frameshift"]
  
  if (nrow(dt) < 1) stop("No neoeptiopes with Consensus_scores < 5000nM ")
  
  dt_pl <- dt[MHC %like% "(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])", MHC := "Class I"] %>%
    .[MHC %like% "(HLA-D[A-Z0-9]+\\*)|(H-2-[A-Z]{2}[a-z])", MHC := "Class II"] %>%
      .[!is.na(type)]
  
  if (nrow(dt_pl) < 1) stop("No neoepitopes meet classification criteria")
  
  gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "type")]
  
  gdt <- dt_pl %>% (function(dt){
    
    ns <- dt[, sample_id %>% unique %>% length]
    
    gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                      MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                      type = c(replicate(ns * 2, "ADN"), replicate(ns * 2, "CDN"), replicate(ns * 2, "priority"),
                                replicate(ns * 2, "frameshift"), replicate(ns * 2, "fusion")),
                      N = 0) %>% unique
    return(gdt)
  })
  
  gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>% 
    .[, N := max(N), by = c("sample_id", "type", "MHC")] %>% unique
  
  if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){
    
    message("Sample_id has >7 characters, shortening names for aesthetics.
    To circumvent this, change sample_ids to less than 7 characters in the input data.table.")
    
    for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){
      
    gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
          sample_id := sample_id %>%
               substr(1, 7) %>% paste0(., "_", i)]
                                                                  }
                                                          }
  
 g1 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
    ggplot2::geom_col(ggplot2::aes(fill = type), col = "black", position = "dodge") +
    ggplot2::facet_wrap(~MHC) +
    ag_gg_theme +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = ag_colors) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
    ggplot2::ylab("peptides") +
    ggplot2::xlab("") +
    ggplot2::ggtitle(paste0("antigen.garnish summary"))
 
 if(save == "TRUE") cowplot::ggsave(plot = g1, "antigen.garnish_summary.pdf", height = 6, width = 9)
  
 if (!(dt[, effect_type %>% unique] %like% "frameshift" %>% any)) g2 <- NA
 
 if (dt[, effect_type %>% unique] %like% "frameshift" %>% any){
 # frameshift plot
   
   dt_pl <- dt[MHC %like% "(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])", MHC := "Class I"] %>%
     .[MHC %like% "(HLA-D[A-Z0-9]+\\*)|(H-2-[A-Z]{2}[a-z])", MHC := "Class II"] %>%
        .[effect_type %like% "frameshift" & Consensus_scores < 1000]
   
   dt_pl[, binding := "<1000nM"] %>% 
     .[Consensus_scores < 500, binding := "<500nM"] %>%
        .[Consensus_scores < 50, binding := "<50nM"]
   
gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")] 

gdt <- dt_pl %>% (function(dt){
  
  ns <- dt[, sample_id %>% unique %>% length]
 
   gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                  MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                  binding = c(replicate(ns * 2, "<50nM"), replicate(ns * 2, "<500nM"), replicate(ns * 2, "<1000nM")),
                  N = 0) %>% unique
        })

gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>% 
  .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique
 
if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){
  
  message("Sample_id has >7 characters, shortening names for aesthetics.
          To circumvent this, change sample_ids to less than 7 characters in the input data.table.")
  
  for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){
    
    gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
          sample_id := sample_id %>%
            substr(1, 7) %>% paste0(., "_", i)]
  }
}

 g2 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
   ggplot2::geom_col(ggplot2::aes(fill = binding), col = "black", position = "dodge") +
   ggplot2::facet_wrap(~MHC) +
   ag_gg_theme +
   ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
   ggplot2::scale_fill_manual(values = ag_colors[1:3]) +
   ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
   ggplot2::ylab("peptides") +
   ggplot2::xlab("") +
   ggplot2::ggtitle(paste0("Frameshift neoepitopes"))
 
 if (save == "TRUE") cowplot::ggsave(plot = g2, "antigen.garnish_Frameshifts_summary.pdf", height = 6, width = 9)
 
 }
  }
  
 # fusion plot
 
  if(!missing(jdt) & missing(dt)) jdt <- jdt
  
  if(!missing(jdt) & !missing(dt)) jdt <- dt
  
 if (!"fusion genes" %chin% (jdt %>% names)) g3 <- NA
 
 if ("fusion genes" %chin% (jdt %>% names)){
  
   
   dt_pl <- jdt[MHC %like% "(HLA-[ABC]\\*)|(H-2-[A-Z][a-z])", MHC := "Class I"] %>%
     .[MHC %like% "(HLA-D[A-Z0-9]+\\*)|(H-2-[A-Z]{2}[a-z])", MHC := "Class II"] %>%
     .[!is.na(fusion_uuid) & Consensus_scores < 1000]
   
   dt_pl[, binding := "<1000nM"] %>% 
     .[Consensus_scores < 500, binding := "<500nM"] %>%
     .[Consensus_scores < 50, binding := "<50nM"]
   
   gg_dt <- dt_pl[, .N, by = c("sample_id", "MHC", "binding")] 
   
   gdt <- dt_pl %>% (function(dt){
     
     ns <- dt[, sample_id %>% unique %>% length]
     
     gdt <- data.table(sample_id = dt[, sample_id %>% unique],
                       MHC = c(replicate(ns, "Class I"), replicate(ns, "Class II")),
                       binding = c(replicate(ns * 2, "<50nM"), replicate(ns * 2, "<500nM"), replicate(ns * 2, "<1000nM")),
                       N = 0) %>% unique
   })
   
   gg_dt <- merge(gg_dt, gdt, by = intersect(names(gg_dt), names(gdt)), all = TRUE) %>% 
     .[, N := max(N), by = c("sample_id", "binding", "MHC")] %>% unique
   
   if (any(gg_dt[, sample_id %>% unique %>% nchar] > 7)){
     
     message("Sample_id has >7 characters, shortening names for aesthetics.
             To circumvent this, change sample_ids to less than 7 characters in the input data.table.")
     
     for (i in gg_dt[nchar(sample_id) > 7, sample_id %>% unique] %>% seq_along){
       
       gg_dt[sample_id == gg_dt[nchar(sample_id) > 7, sample_id %>% unique][i],
             sample_id := sample_id %>%
               substr(1, 7) %>% paste0(., "_", i)]
     }
   }
   
   g3 <- ggplot2::ggplot(gg_dt, ggplot2::aes(x = sample_id, y = N)) +
     ggplot2::geom_col(ggplot2::aes(fill = binding), col = "black", position = "dodge") +
     ggplot2::facet_wrap(~MHC) +
     ag_gg_theme +
     ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
     ggplot2::scale_fill_manual(values = ag_colors[1:3]) +
     ggplot2::theme(strip.text.x = ggplot2::element_text(size = ggplot2::rel(2))) +
     ggplot2::ylab("peptides") +
     ggplot2::xlab("") +
     ggplot2::ggtitle(paste0("Fusion neoepitopes"))
   
   if (save == "TRUE") cowplot::ggsave(plot = g3, "antigen.garnish_Fusions_summary.pdf", height = 6, width = 9)
   
 }
  
  if (missing(dt) & !missing(jdt)) return(g3)
  
  if (missing(jdt) & !missing(dt)) return(list(g1,g2))
 
 glist <- list(g1, g2, g3)
 
 return(glist)
 
}

 






## ---- cat_tables
#' Internal function to combine garnish_predictions outputs from garnish_jaffa
#' and garnish_variants or direct input to garnish_predictions.
#'
#' @export cat_tables
#' @md
#' 
cat_tables <- function(dt1, dt2){
  ##get the jaffa_dt
  if ("fus_tx" %chin% names(dt1) %>% any) jdt <- dt1
  if ("fus_tx" %chin% names(dt2) %>% any) jdt <- dt2
  
  ##get the vcf predictions output dt  
  if (grepl("mutnfs", (dt1[, .SD %>% unique, .SDcols = "pep_type"] %>%
                       unlist)) %>% any) vdt <- dt1
  if (grepl("mutnfs", (dt2[, .SD %>% unique, .SDcols = "pep_type"] %>%
                       unlist)) %>% any) vdt <- dt2
  
  ##subset tables for what we care about for downstream merge
  
  vdt <- vdt[, .SD, .SDcols = c("nmer", "sample_id", "MHC", "Consensus_scores", "DAI",
                                "Lower.CI", "Upper.CI", "ensembl_transcript_id",
                                "effect_type", "external_gene_name", "frameshift", "protein_change",
                                "cDNA_change", "mutant_index", "pep_mut", "pep_wt", "pep_type", "var_uuid", "dai_uuid", "nmer_uuid")]
  
  jdt <- jdt[, .SD, .SDcols = c("nmer", "sample_id", "MHC", "Consensus_scores", "Lower.CI", "Upper.CI",
                                "fusion genes", "frameshift", "mutant_index", "pep_mut", "pep_gene_1", "var_uuid",
                                "fusion_uuid", "nmer_uuid")]   
  
  ##add fus as pep_type for later analysis, add mutfs for frameshifts
  jdt[, pep_type := "fus"]
  
  vdt[is.na(pep_type) & effect_type %like% "frameshift",
      pep_type := "mutfs"]
  
  ##change names to overlap in the vcf dt
  strings <- jdt[, sample_id %>% unique]
  
  ###contruct our regex
  strings <- paste("(", strings, ")", sep = "", collapse = "|")
  
  ##rename stings in vcf dt from full BAM name
  vdt[, sample_id := sample_id %>% stringr::str_extract(pattern = strings)]
  
  ##final merge
  dto <- merge(vdt, jdt, all = TRUE, by = intersect(names(vdt), names(jdt)))
  return(dto)
}

  

