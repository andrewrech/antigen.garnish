library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("get_snpeff_annot", {

  # load test data

  dt <- data.table::data.table(INFO = c("DP=547;ECNT=5;IN_PON;NLOD=27.52;N_ART_LOD=29.20;POP_AF=1.000e-03;P_GERMLINE=-2.422e+01;TLOD=71.63;ANN=C|missense_variant|MODERATE|Phf13|ENSMUSG00000047777|transcript|ENSMUST00000055688.9|protein_coding|3/4|c.431A>G|p.Tyr144Cys|825/3205|431/891|144/296||,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759|transcript|ENSMUST00000036680.7|protein_coding||c.-3347A>G|||||3306|,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759|transcript|ENSMUST00000105665.2|protein_coding||c.-3347A>G|||||3333|,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759|transcript|ENSMUST00000146471.1|retained_intron||n.-3323A>G|||||3323|",
  "DP=349;ECNT=8;NLOD=41.54;N_ART_LOD=-2.143e+00;POP_AF=1.000e-03;P_GERMLINE=-3.824e+01;TLOD=92.73;ANN=A|missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474|transcript|ENSMUST00000074241.8|protein_coding|4/10|c.332C>T|p.Ala111Val|469/1821|332/1452|111/483||,A|missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474|transcript|ENSMUST00000160918.7|protein_coding|3/9|c.332C>T|p.Ala111Val|493/2077|332/1452|111/483||,A|non_coding_transcript_exon_variant|MODIFIER|Thnsl2|ENSMUSG00000054474|transcript|ENSMUST00000170455.1|processed_transcript|2/6|n.150C>T||||||")) %>%
  get_snpeff_annot


    dt[, .SD, .SDcols = c("effect_type", "ensembl_transcript_id",
    "ensembl_gene_id", "protein_change", "cDNA_change",
    "protein_coding")] %>%

    testthat::expect_equal(.,
      structure(list(effect_type = c("missense_variant", "missense_variant"),
                     ensembl_transcript_id = c("ENSMUST00000055688",
                                               "ENSMUST00000074241"),
                     ensembl_gene_id = c("ENSMUSG00000047777",
                                         "ENSMUSG00000054474"),
                     protein_change = c("p.Tyr144Cys",
                                        "p.Ala111Val"),
                     cDNA_change = c("c.431A>G",
                                     "c.332C>T"),
                     protein_coding = c(TRUE, TRUE)),
                     .Names = c("effect_type",
                                "ensembl_transcript_id",
                                "ensembl_gene_id",
                                "protein_change",
                                "cDNA_change",
                                "protein_coding"),
                     class = c("data.table", "data.frame"),
                     row.names = c(NA, -2L))
              )

    })