library(testthat)
library(antigen.garnish)
library(data.table)
library(magrittr)
library(dt.inflix)

testthat::test_that("garnish_predictions_assemble_generate", {

  list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove
  on.exit(list.files(pattern = "(netMHC|nuggets|flurry).*-.*-.*\\.csv") %>% file.remove)

  # load test data
    dt <- data.table::data.table(
        CHROM = c("7", "11", "11"),
        POS = c("25348381", "60554341", "60731323"),
        ID = c(NA_character_, NA_character_, NA_character_),
        REF = c("C", "T", "T"),
        ALT = c("T", "G", "A"),
        QUAL = c(NA_character_, NA_character_, NA_character_),
        FILTER = c("PASS", "PASS", "PASS"),
        INFO = c("SOMATIC;QSS=124;TQSS=1;NT=ref;QSS_NT=124;TQSS_NT=1;SGT=CC->CT;DP=492;MQ=57.28;MQ0=9;ReadPosRankSum=-6.35;SNVSB=0.00;EVS=21.34;ANN=T|missense_variant|MODERATE|Megf8|ENSMUSG00000045039|transcript|ENSMUST00000128119.1|protein_coding|28/41|c.4988C>T|p.Ala1663Val|5331/10040|4988/8370|1663/2789||,T|non_coding_transcript_exon_variant|MODIFIER|Megf8|ENSMUSG00000045039|transcript|ENSMUST00000153077.1|retained_intron|5/18|n.963C>T||||||",
        "SOMATIC;QSS=2219;TQSS=1;NT=ref;QSS_NT=3070;TQSS_NT=1;SGT=TT->GT;DP=1165;MQ=50.88;MQ0=119;ReadPosRankSum=0.00;SNVSB=0.00;EVS=25.57;ANN=G|missense_variant|MODERATE|Alkbh5|ENSMUSG00000042650|transcript|ENSMUST00000044250.3|protein_coding|4/4|c.1114T>G|p.Ser372Ala|1559/5730|1114/1188|372/395||,G|non_coding_transcript_exon_variant|MODIFIER|Alkbh5|ENSMUSG00000042650|transcript|ENSMUST00000134770.1|processed_transcript|4/4|n.598T>G||||||",
        "SOMATIC;QSS=487;TQSS=1;NT=ref;QSS_NT=3070;TQSS_NT=1;SGT=TT->AT;DP=171;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;EVS=34.24;ANN=A|missense_variant|MODERATE|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000018743.4|protein_coding|4/4|c.718T>A|p.Trp240Arg|894/2522|718/1365|240/454||,A|upstream_gene_variant|MODIFIER|Flii|ENSMUSG00000002812|transcript|ENSMUST00000002889.4|protein_coding||c.-4134A>T|||||4060|,A|upstream_gene_variant|MODIFIER|Flii|ENSMUSG00000002812|transcript|ENSMUST00000137226.1|retained_intron||n.-4100A>T|||||4100|,A|upstream_gene_variant|MODIFIER|Gm23341|ENSMUSG00000089581|transcript|ENSMUST00000158956.1|snRNA||n.-4538T>A|||||4538|,A|downstream_gene_variant|MODIFIER|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000154890.1|retained_intron||n.*327T>A|||||327|,A|downstream_gene_variant|MODIFIER|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000146159.1|processed_transcript||n.*660T>A|||||660|,A|downstream_gene_variant|MODIFIER|Mir5100|ENSMUSG00000092734|transcript|ENSMUST00000174993.1|miRNA||n.*2597T>A|||||2597|"),
        FORMAT = c("DP:FDP:SDP:SUBDP:AU:CU:GU:TU", "DP:FDP:SDP:SUBDP:AU:CU:GU:TU", "DP:FDP:SDP:SUBDP:AU:CU:GU:TU"),
        NORMAL = c("149:0:0:0:0,0:149,149:0,0:0,6", "584:1:0:0:0,0:0,0:1,1:582,810", "116:0:0:0:0,0:0,0:0,0:116,116"),
        TUMOR = c("270:1:0:0:0,0:233,234:0,0:36,103", "252:4:0:0:0,0:0,1:248,350:0,3", "55:0:0:0:55,55:0,0:0,0:0,0"),
        sample_id = c("normal_tumor.bam", "normal_tumor.bam", "normal_tumor.bam"),
        vcf_type = c("other", "other", "other"),
        se = c("missense_variant|MODERATE|Megf8|ENSMUSG00000045039|transcript|ENSMUST00000128119.1|protein_coding|28/41|c.4988C>T|p.Ala1663Val|5331/10040|4988/8370|1663/2789||,T|non_coding_transcript_exon_variant|MODIFIER|Megf8|ENSMUSG00000045039|transcript|ENSMUST00000153077.1|retained_intron|5/18|n.963C>T||||||",
        "missense_variant|MODERATE|Alkbh5|ENSMUSG00000042650|transcript|ENSMUST00000044250.3|protein_coding|4/4|c.1114T>G|p.Ser372Ala|1559/5730|1114/1188|372/395||,G|non_coding_transcript_exon_variant|MODIFIER|Alkbh5|ENSMUSG00000042650|transcript|ENSMUST00000134770.1|processed_transcript|4/4|n.598T>G||||||",
        "missense_variant|MODERATE|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000018743.4|protein_coding|4/4|c.718T>A|p.Trp240Arg|894/2522|718/1365|240/454||,A|upstream_gene_variant|MODIFIER|Flii|ENSMUSG00000002812|transcript|ENSMUST00000002889.4|protein_coding||c.-4134A>T|||||4060|,A|upstream_gene_variant|MODIFIER|Flii|ENSMUSG00000002812|transcript|ENSMUST00000137226.1|retained_intron||n.-4100A>T|||||4100|,A|upstream_gene_variant|MODIFIER|Gm23341|ENSMUSG00000089581|transcript|ENSMUST00000158956.1|snRNA||n.-4538T>A|||||4538|,A|downstream_gene_variant|MODIFIER|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000154890.1|retained_intron||n.*327T>A|||||327|,A|downstream_gene_variant|MODIFIER|Mief2|ENSMUSG00000018599|transcript|ENSMUST00000146159.1|processed_transcript||n.*660T>A|||||660|,A|downstream_gene_variant|MODIFIER|Mir5100|ENSMUSG00000092734|transcript|ENSMUST00000174993.1|miRNA||n.*2597T>A|||||2597|"),
        snpeff_uuid = c("e74947ee-7e8b-11e7-9090-12e0492b8512", "dc4ffadd-fda9-4efc-9604-a6a31498c59d", "0db9e696-aeef-4272-b184-0e3eaad29242"),
        effect_type = c("missense_variant", "missense_variant", "missense_variant"),
        ensembl_transcript_id = c("ENSMUST00000128119", "ENSMUST00000044250", "ENSMUST00000018743"),
        ensembl_gene_id = c("ENSMUSG00000045039", "ENSMUSG00000042650", "ENSMUSG00000018599"),
        protein_change = c("p.Ala1663Val", "p.Ser372Ala", "p.Trp240Arg"),
        cDNA_change = c("c.4988C>T", "c.1114T>G", "c.718T>A"),
        protein_coding = c(TRUE, TRUE, TRUE),
        MHC = c("HLA-A*02:01 HLA-DRB1*14:67", "HLA-A*02:01 HLA-DRB1*14:67", "HLA-A*03:01 HLA-DRB1*03:01"))

  # run test
    dto <- garnish_predictions(dt, predict = FALSE)

  testthat::expect_equal(dto %>% nrow, 551)
  testthat::expect_equal(dto %>% length, 54)
  testthat::expect_true(dto$cDNA_locs %>% unique %in% c(4988, 718, 1114) %>% all)

})
