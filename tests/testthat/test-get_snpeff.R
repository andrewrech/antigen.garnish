testthat::test_that("get_vcf_snpeff_dt", {

  # load test data

  dt <- data.table::data.table(ANN = c(
    "C|missense_variant|MODERATE|Phf13|ENSMUSG00000047777.1|transcript|ENSMUST00000055688.9|protein_coding|3/4|c.431A>G|p.Tyr144Cys|825/3205|431/891|144/296||,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000036680.7|protein_coding||c.-3347A>G|||||3306|,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000105665.2|protein_coding||c.-3347A>G|||||3333|,C|upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000146471.1|retained_intron||n.-3323A>G|||||3323|",
    "A|missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474.1|transcript|ENSMUST00000074241.8|protein_coding|4/10|c.332C>T|p.Ala111Val|469/1821|332/1452|111/483||,A|missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474.1|transcript|ENSMUST00000160918.7|protein_coding|3/9|c.332C>T|p.Ala111Val|493/2077|332/1452|111/483||,A|non_coding_transcript_exon_variant|MODIFIER|Thnsl2|ENSMUSG00000054474|transcript|ENSMUST00000170455.1|processed_transcript|2/6|n.150C>T||||||"
  )) %>%
    # run test
    antigen.garnish:::get_vcf_snpeff_dt()

  dtNames <- dt %>% names()
  dtNames %<>% .[!dtNames %like% "uuid"]

  correctlyParsed <- structure(list(ANN = c(
    "missense_variant|MODERATE|Phf13|ENSMUSG00000047777.1|transcript|ENSMUST00000055688.9|protein_coding|3/4|c.431A>G|p.Tyr144Cys|825/3205|431/891|144/296||",
    "upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000036680.7|protein_coding||c.-3347A>G|||||3306|",
    "upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000105665.2|protein_coding||c.-3347A>G|||||3333|",
    "upstream_gene_variant|MODIFIER|Thap3|ENSMUSG00000039759.1|transcript|ENSMUST00000146471.1|retained_intron||n.-3323A>G|||||3323|",
    "missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474.1|transcript|ENSMUST00000074241.8|protein_coding|4/10|c.332C>T|p.Ala111Val|469/1821|332/1452|111/483||",
    "missense_variant|MODERATE|Thnsl2|ENSMUSG00000054474.1|transcript|ENSMUST00000160918.7|protein_coding|3/9|c.332C>T|p.Ala111Val|493/2077|332/1452|111/483||",
    "non_coding_transcript_exon_variant|MODIFIER|Thnsl2|ENSMUSG00000054474|transcript|ENSMUST00000170455.1|processed_transcript|2/6|n.150C>T||||||"
  ), transcript_id = c(
    "ENSMUST00000055688.9", "ENSMUST00000036680.7",
    "ENSMUST00000105665.2", "ENSMUST00000146471.1", "ENSMUST00000074241.8",
    "ENSMUST00000160918.7", "ENSMUST00000170455.1"
  ), effect_type = c(
    "missense_variant",
    "upstream_gene_variant", "upstream_gene_variant", "upstream_gene_variant",
    "missense_variant", "missense_variant", "non_coding_transcript_exon_variant"
  ), putative_impact = c(
    "MODERATE", "MODIFIER", "MODIFIER", "MODIFIER",
    "MODERATE", "MODERATE", "MODIFIER"
  ), gene = c(
    "Phf13", "Thap3",
    "Thap3", "Thap3", "Thnsl2", "Thnsl2", "Thnsl2"
  ), gene_id = c(
    "ENSMUSG00000047777.1",
    "ENSMUSG00000039759.1", "ENSMUSG00000039759.1", "ENSMUSG00000039759.1",
    "ENSMUSG00000054474.1", "ENSMUSG00000054474.1", "ENSMUSG00000054474"
  ), feature_type = c(
    "transcript", "transcript", "transcript",
    "transcript", "transcript", "transcript", "transcript"
  ), feature_id = c(
    "ENSMUST00000055688.9",
    "ENSMUST00000036680.7", "ENSMUST00000105665.2", "ENSMUST00000146471.1",
    "ENSMUST00000074241.8", "ENSMUST00000160918.7", "ENSMUST00000170455.1"
  ), transcript_bioptype = c(
    "protein_coding", "protein_coding",
    "protein_coding", "retained_intron", "protein_coding", "protein_coding",
    "processed_transcript"
  ), exon_intron_rank = c(
    "3/4", "", "",
    "", "4/10", "3/9", "2/6"
  ), cDNA_change = c(
    "c.431A>G", "c.-3347A>G",
    "c.-3347A>G", "n.-3323A>G", "c.332C>T", "c.332C>T", "n.150C>T"
  ), protein_change = c(
    "p.Tyr144Cys", "", "", "", "p.Ala111Val",
    "p.Ala111Val", ""
  ), cDNA_position_cDNA_len = c(
    "825/3205", "",
    "", "", "469/1821", "493/2077", ""
  ), CDS_position_CDS_len = c(
    "431/891",
    "", "", "", "332/1452", "332/1452", ""
  ), Protein_position_Protein_len = c(
    "144/296",
    "", "", "", "111/483", "111/483", ""
  ), Distance_to_feature = c(
    "",
    "3306", "3333", "3323", "", "", ""
  ), protein_coding = c(
    TRUE,
    FALSE, FALSE, FALSE, TRUE, TRUE, FALSE
  )), row.names = c(
    NA,
    -7L
  ), class = c("data.table", "data.frame"))

  dt %>%
    .[, .SD, .SDcols = dtNames] %>%
    testthat::expect_equal(
      .,
      correctlyParsed
    )
})
