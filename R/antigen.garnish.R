## -------- garnish_variants
#' Intakes variants and returns an intersected data table for epitope prediction.
#'
#' Intakes variants from a \href{https://github.com/broadinstitute/gatk}{MuTect2}/\href{https://github.com/Illumina/strelka}{Strelka2} - \href{https://github.com/pcingola/SnpEff}{SnpEff} variant annotation pipeline and filters for neoepitope prediction. Hg38 (human) and GRCm38 (murine) variant calls are required. Mutect2 and Strelka variant threshold prior to intersection were empirically established to limit false positives.
#'
#' @param vcfs Character vector. Strelka2 or Mutect2 VFC files to import.
#'
#' @return A list of three data tables: 1) all passing variants and 2) intersected variants, 3) intersected missense mutations for neoepitope prediction.
#'
#' @export garnish_variants

garnish_variants <- function(vcfs) {

  ### package building
    invisible(dt.inflix::allduplicated(data.table::data.table(a="a")))
    invisible(testthat::compare(1, 1))
  ###

  # NGS data loaded in parallel

  ivfdtl <- lapply(vcfs %>% seq_along, function(ivf){

  # load dt
      vcf <-  vcfR::read.vcfR(vcfs[ivf], verbose = TRUE)

  # extract sample names from Mutect2 and Strelka command line for intersection

      sample_id <- c((vcf@meta %include% "[Cc]ommand" %>% stringr::str_extract_all("[^ ]+\\.bam") %>% unlist) %include% (vcf@meta %include% "[Cc]ommand" %>% stringr::str_extract("(?<=--tumorSampleName )[^ ]*")), vcf@meta %>% stringr::str_extract("(?<=tumorBam )[^ ]*") %>% basename %>% stringr::str_replace("\\.bam", "")) %>%
                                  stats::na.omit %>%
                                  basename %>%
                                  stringr::str_replace("\\.bam", "")

      # extract vcf type
      vcf_type <- vcf@meta %>% unlist %>% stringr::str_extract(stringr::regex("Strelka|Mutect", ignore_case = TRUE)) %>% stats::na.omit %>% unlist %>% data.table::first

  # return a data table of variants

  vdt <- vcf@fix %>%data.table::as.data.table

  if (vcf@gt %>% length > 0) vdt <- cbind(vdt, vcf@gt %>%data.table::as.data.table)

  if(vdt %>% nrow < 1) return(data.table::data.table(sample_id = sample_id))

  # filter passing Strelka2 variants
  if(vcf_type == "Strelka") vdt <- vdt[FILTER == "PASS"]
  if(vcf_type == "Mutect") vdt <- vdt[INFO %>%
                                      stringr::str_extract("(?<=TLOD=)[0-9\\.]") %>%
                                      as.numeric > 6.0]
  vdt[, sample_id := sample_id]
  vdt[, vcf_type := vcf_type]

  # extract full snpeff annotation

  vdt[, se := INFO %>%
    stringr::str_extract("ANN.*") %>%
    stringr::str_replace("ANN=[^\\|]+\\|", "")]

  # add a variant identifier
  suppressWarnings(vdt[, uuid :=
                lapply(1:nrow(vdt),
                uuid::UUIDgenerate) %>% unlist])

  # abort if no variants passed filtering
  if (vdt %>% nrow < 1) return(NULL)

  # spread SnpEff annotation over rows
  vdt %>% tidyr::separate_rows("se", sep = ",")
  vdt %<>% tidyr::separate_rows("se", sep = ",")

  # extract info from snpeff annotation
    vdt[, effect_type := se %>%
        stringr::str_extract("^[a-z0-9][^\\|]+")]
    vdt[, transcript_affected_v := se %>%
        stringr::str_extract("(?<=\\|)(ENSMUST|ENST)[0-9.]+(?=\\|)")]
    vdt[, transcript_affected := se %>%
        stringr::str_extract("(?<=\\|)(ENSMUST|ENST)[0-9]+")]
    vdt[, gene_affected := se %>%
        stringr::str_extract("(?<=\\|)(ENSMUSG|ENSG)[0-9.]+(?=\\|)")]
    vdt[, protein_change := se %>%
        stringr::str_extract("c\\.[^\\|]+")]
    vdt[, dna_change := se %>%
        stringr::str_extract("c\\.[^\\|]+")]
    vdt[, aa_mutation := se %>%
        stringr::str_extract("(?<=p\\.)[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}") %>%
        aa_convert]
    vdt[, protein_coding := se %>%
        stringr::str_detect("protein_coding")]

  return(vdt)
  })

  ivfdt <- ivfdtl %>% data.table::rbindlist

  merge_vcf <- function(dt, dt2){

    # a function to intersect annotated variants across VCFs using SnpEff

      sdt <- merge(dt, dt2[, .SD,
             .SDcols = c("CHROM",
              "POS",
              "REF",
              "dna_change")],
               by = c("CHROM",
              "POS",
              "REF",
              "dna_change")) %>% .[, vcf_type := "intersect"]
      return(sdt)

    }

  # return an intersected data table of variants

  sdt <- parallel::mclapply(ivfdt[, sample_id %>% unique], function(sn){

    # find data tables with matching sample names
    sdt <- lapply(ivfdtl, function(dt){

     dt[, sample_id %>% .[1]] == sn

    }) %>% unlist

  # merge all data tables with matching sample names
    if (ivfdtl[sdt] %>% length == 1) return(ivfdtl[[sdt]])
    if (ivfdtl[sdt] %>% length > 1) return(ivfdtl[sdt] %>% Reduce(merge_vcf, .))


  }) %>% data.table::rbindlist

if (ivfdt$transcript_affected %>%
    stats::na.omit %>%
    .[1] %like% "ENSMUST") bmds <- "mmusculus_gene_ensembl"

if (ivfdt$transcript_affected %>%
    stats::na.omit %>%
    .[1] %like% "ENST") bmds <- "hsapiens_gene_ensembl"


  # add metadata

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                             dataset = bmds,
                                             host = 'ensembl.org')

  trneff <- ivfdt$transcript_affected %>%
                              sort %>%
                              unique %>%
                              stats::na.omit

  if (trneff %>% length > 1){
    var_dt <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                  "external_gene_name", "ensembl_gene_id", "description", "chromosome_name",
                  "start_position", "end_position", "transcript_start", "transcript_end",
                  "transcript_length", "refseq_mrna"),
                       filters = c("ensembl_transcript_id"),
                       values = list(trneff),
                       mart = mart) %>%
          data.table::as.data.table %>%
          data.table::setnames("ensembl_transcript_id", "transcript_affected")

  ivfdt <- merge(ivfdt, var_dt, by = "transcript_affected", all.x = TRUE)
  sdt <- merge(sdt, var_dt, by = "transcript_affected", all.x = TRUE)
 }

agdt <- sdt[!transcript_affected %>% is.na &
            !aa_mutation %>% is.na &
            protein_coding == TRUE &
            effect_type == "missense_variant"]

return(
       list(
        all_variants = ivfdt,
        all_intersected_variants = sdt,
        antigen.garnish_input = agdt
            ))
}



## -------- garnish_predictions
#' Performs epitope prediction.
#'
#' Performs epitope prediction on a data table of missence mutations.
#'
#' @param mhc_dt Data table. Input data frame from garnish_variants. Must contain three columns: sample_id, transcript_affected (ensembl_transcript_id), aa_mutation (single letter abbreviations), and space-seperated MHC (MHC type, e.g. 'HLA-A*02:01 HLA-A*03:01', 'H-2-Kb H-2-Kb' or HLA-DRB1*11:07 [second type]')
#' @param assemble Logical. Assemble data table?
#' @param generate Logical. Generate peptides and commands?
#' @param predict Logical. Predict binding affinities?
#'
#' @return A data table of neoepitopes.
#' @export garnish_predictions

garnish_predictions <- function(mhc_dt,
                                   assemble = TRUE,
                                   generate = TRUE,
                                   predict = TRUE) {

      if (!c("sample_id", "transcript_affected", "MHC") %chin% (mhc_dt %>% names) %>% all) stop("mhc_dt must contain these three columns: sample_id, transcript_affected (ensembl_transcript_id), aa_mutation (single letter abbreviations), and space-seperated MHC (MHC type, e.g. 'HLA-A*02:01 HLA-A*03:01', 'H-2-Kb H-2-Kb' or HLA-DRB1*11:07 [second type]')")

if (assemble){

  ## --------  load peptide databse

  if (mhc_dt[, transcript_affected %>% unique %>% stringr::str_detect("ENSMUST") %>% any]){

      if (!file.exists("Mus_musculus.GRCm38"))utils::download.file(method = "wget", url = "http://get.rech.io/Mus_musculus.GRCm38", destfile = "Mus_musculus.GRCm38")

      if (!file.exists("Mus_musculus.GRCm38.pep"))utils::download.file(method = "wget", url = "http://get.rech.io/Mus_musculus.GRCm38.pep", destfile = "Mus_musculus.GRCm38.pep")

        qdt <- readRDS("Mus_musculus.GRCm38")
        db_pep <- readRDS("Mus_musculus.GRCm38.pep")

    } else {

      if (!file.exists("Homo_sapiens.GRCh38"))utils::download.file(method = "wget", url = "http://get.rech.io/Homo_sapiens.GRCh38", destfile = "Homo_sapiens.GRCh38")
      if (!file.exists("Homo_sapiens.GRCh38.pep"))utils::download.file(method = "wget", url = "http://get.rech.io/Homo_sapiens.GRCh38.pep", destfile = "Homo_sapiens.GRCh38.pep")

      qdt <- readRDS("Homo_sapiens.GRCh38")
      db_pep <- readRDS("Homo_sapiens.GRCh38.pep")

       }



  ## -------- peptide generation for epitope prediction

  mhc_dt[, mutant_loc := aa_mutation %>% stringr::str_extract("[0-9]+") %>% as.integer]
  mhc_dt[, mutate_from := aa_mutation %>% stringr::str_extract("^[A-Z]")]

  # create peptide lengths
  mhc_dt[, frag_begin := mutant_loc - 14]
  mhc_dt[, frag_end := mutant_loc + 14]
  mhc_dt[frag_begin >= 1, mut_pep_loc := 15 %>% as.integer]
  mhc_dt[frag_begin < 1, `:=` (frag_begin = 1 %>% as.integer, mut_pep_loc = mutant_loc)]

  # directly replace mutant amino acid

  mhc_dt[, mutate_to := aa_mutation %>% stringr::str_extract("[A-Z]$")]

  if ("protein_index" %chin% (mhc_dt %>% names)) mhc_dt[, "protein_index" := NULL]

  mhc_dt <- merge(mhc_dt, qdt, by = "transcript_affected", all.x = TRUE)

  # function to report mutant amino acid from peptide database

    get_prediction_aa <- function(db_pep, protein_index, begin, end){

    parallel::mclapply(1:length(protein_index), function(x){
    vec <- try(seqinr::getFrag(db_pep[protein_index[x]], begin = begin[x], end = end[x]) %>% unlist %>% toupper)

    if ("try-error" %chin% (vec %>% class)) vec <- NA
    return(vec)
    }) %>% unlist
    }

  # function to report total length of mutant protein for peptide length adjustment

    get_protein_length <- function(db_pep, protein_index){

    parallel::mclapply(1:length(protein_index), function(x){
      vec <- try(db_pep[[protein_index[x]]] %>% length)
      if ("try-error" %chin% (vec %>% class)) vec <- NA
      return(vec)
    }) %>% unlist
    }

  # function to return a peptide fragment

    get_pred_frag <- function(db_pep, protein_index, begin, end){

    parallel::mclapply(1:length(protein_index), function(x){
      vec <- try(seqinr::getFrag(db_pep[protein_index[x]],
                                 begin = begin[x],
                                end = end[x]) %>%
                                 unlist %>%
                                 paste(collapse = "") %>%
                                 toupper)

      if ("try-error" %chin% (vec %>% class)) vec <- NA
      return(vec)

    }) %>% unlist
    }

  # function to produce the mutant peptide fragment

    get_mut_frag <- function(pep_wt, mut_pep_loc, mutate_to){

    parallel::mclapply(1:length(pep_wt), function(x){

      if(pep_wt[x] %>% is.na) return(NA)

      y <- (pep_wt[x] %>% strsplit("") %>% unlist)
      y[mut_pep_loc[x]] <- mutate_to[x]
    return(y %>% paste(collapse = ""))
    }) %>% unlist
    }

  mhc_dt[!protein_index %>% is.na, pep_db_aa_from := get_prediction_aa(db_pep, protein_index, mutant_loc, mutant_loc)]
  mhc_dt[, pep_db_aa_length := get_protein_length(db_pep, protein_index)]


  mhc_dt[, frag_begin := frag_begin %>% as.integer]
  mhc_dt[, frag_end := frag_end %>% as.integer]
  mhc_dt[frag_end > pep_db_aa_length, frag_end := pep_db_aa_length]

  mhc_dt[, pep_wt := get_pred_frag(db_pep, protein_index, frag_begin, frag_end)]
  mhc_dt[, pep_mut := get_mut_frag(pep_wt, mut_pep_loc, mutate_to)]

}

if (generate) {


  # generation a uuid for each unique variant
  suppressWarnings(mhc_dt[, variant_uuid :=
                  parallel::mclapply(1:nrow(mhc_dt),
                  uuid::UUIDgenerate) %>% unlist])

  # generate a datatable of unique variants for peptide generation
  basepep_dt <- mhc_dt[pep_db_aa_from == mutate_from &
                !pep_mut %>% is.na &
                !pep_wt %>% is.na &
                !mut_pep_loc %>% is.na,
                .(pep_mut, pep_wt, mut_pep_loc, variant_uuid)] %>%
                unique

if (basepep_dt %>% nrow == 0) stop("No unique missense variants for peptide generation")

  peptide_dt <- garnish_predictions_worker(basepep_dt = basepep_dt)

  # generation a uuid for each unique peptide
  suppressWarnings(peptide_dt[, pep_uuid :=
                  parallel::mclapply(1:nrow(peptide_dt),
                  uuid::UUIDgenerate) %>% unlist])

  mhc_dt <- merge(mhc_dt,
                  peptide_dt,
                    by = "variant_uuid",
                    all.x = TRUE)


}

if (predict) {


# test that prediction tools are in PATH
  if (suppressWarnings(system('which mhcflurry-predict', intern = TRUE)) %>%
        length == 0) stop("mhcflurry-predict is not in PATH")
  if (suppressWarnings(system('which netMHC', intern = TRUE)) %>%
        length == 0) stop("netMHC is not in PATH")
  if (suppressWarnings(system('which netMHCpan', intern = TRUE)) %>%
        length == 0) stop("netMHCpan is not in PATH")
  if (suppressWarnings(system('which netMHCII', intern = TRUE)) %>%
        length == 0) stop("netMHCII is not in PATH")
  if (suppressWarnings(system('which netMHCIIpan', intern = TRUE)) %>%
        length == 0) stop("netMHCIIpan is not in PATH")

  if (mhc_dt[, MHC %>% unique] %>% stringr::str_detect(" ") %>% any) mhc_dt %<>% tidyr::separate_rows("MHC", sep = " ")

# function for writing input peptides for prediction to disk in segments

  write_peps <- function(dt, type){

    if (dt %>% nrow == 0) return(NULL)

    combs <- data.table::CJ(dt[, get(type)] %>% unique,
                            dt[, pep_length] %>% unique)

    dto <- parallel::mclapply(1:nrow(combs), function(i) {

      dts <- dt[get(type) == combs$V1[i] & pep_length == combs$V2[i]]

      chunks <- ((dts %>% nrow)/100) %>% ceiling

      dto <- parallel::mclapply(dts %>% split(1:chunks), function(dtw){

        filename <- paste0(type, "_",
                    uuid::UUIDgenerate(), ".csv")


        # write out unique peptides for MHC type, length
        data.table::fwrite(dtw[, .(peptide)] %>% unique,
                          filename,
                          col.names = FALSE)

        return(data.table::data.table(
               type = type,
               allele = combs$V1[i],
               pep_length = combs$V2[i],
               filename = filename))


        }) %>% data.table::rbindlist
      return(dto)

    }) %>% data.table::rbindlist

  return(dto)
  }


  # function to replacement HLA with matching type for netMHC

 detect_hla <- function(x){

  prog <- deparse(substitute(x))

  for (hla in (x %>% unique)) {

    allele <- alleles[type == prog, allele %>%
      .[stringi::stri_detect_fixed(., hla) %>% which]]

      if (allele %>% length == 0) allele <- NA

      x %<>% stringr::str_replace_all(stringr::fixed(hla), allele)
      }

    return(x)
  }

  # get available MHC alleles for predictions

  alleles <- data.table::rbindlist(
                         list(
      system.file("extdata",
          "netMHC_alleles.txt", package = "antigen.garnish") %>%
                          data.table::fread(header = FALSE) %>%
                          data.table::setnames("V1", "allele") %>%
                          .[, type := "netMHC"],
      system.file("extdata",
          "netMHCpan_alleles.txt", package = "antigen.garnish") %>%
                          data.table::fread(header = FALSE) %>%
                          data.table::setnames("V1", "allele") %>%
                          .[, type := "netMHCpan"],
      system.file("extdata",
          "mhcflurry_alleles.txt", package = "antigen.garnish") %>%
                          data.table::fread(header = FALSE) %>%
                          data.table::setnames("V1", "allele") %>%
                          .[, type := "mhcflurry"],
      system.file("extdata",
          "netMHCII_alleles.txt", package = "antigen.garnish") %>%
                          data.table::fread(header = FALSE) %>%
                          data.table::setnames("V1", "allele") %>%
                          .[, type := "netMHCII"],
      system.file("extdata",
          "netMHCIIpan_alleles.txt", package = "antigen.garnish") %>%
                          data.table::fread(header = FALSE) %>%
                          data.table::setnames("V1", "allele") %>%
                          .[, type := "netMHCIIpan"]
                          ))

  # generate csv for mhcflurry predictions

      mhc_dt[MHC %chin% alleles[type == "mhcflurry", allele] &
      pep_length < 15,
            .SD, .SDcols = c("MHC", "peptide")] %>%
      data.table::copy %>%
      data.table::setnames("MHC", "allele") %>%
      unique %>%
      data.table::fwrite("mhcflurry_input.csv")

  # generate matchable MHC substring for netMHC tools

    mhc_dt[, netMHCpan := MHC %>% stringr::str_replace_all(stringr::fixed("*"), "")]
    mhc_dt[, netMHC := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "")]
    mhc_dt[, netMHCII := netMHCpan %>% stringr::str_replace(stringr::fixed(":"), "") %>%
                                       stringr::str_replace(fixed("HLA-"), "")]
    mhc_dt[, netMHCIIpan := netMHCII %>% stringr::str_replace("DRB1", "DRB1_")]


  # replace substring with netMHC allele type

    mhc_dt[, netMHCpan := detect_hla(netMHCpan)]
    mhc_dt[, netMHC := detect_hla(netMHC)]
    mhc_dt[, netMHCII := detect_hla(netMHCII)]
    mhc_dt[, netMHCIIpan := detect_hla(netMHCIIpan)]


      dtfn <-
        {
          data.table::rbindlist(
            list(
              mhc_dt[!netMHC %>% is.na &
              pep_length < 15,
                .SD, .SDcols = c("netMHC", "peptide", "pep_length")] %>%
              data.table::copy %>%
              data.table::setkey(netMHC, pep_length) %>%
              write_peps("netMHC"),

              mhc_dt[!netMHCpan %>% is.na &
              pep_length < 15,
                .SD, .SDcols = c("netMHCpan", "peptide", "pep_length")] %>%
              data.table::copy %>%
              data.table::setkey(netMHCpan, pep_length) %>%
              write_peps("netMHCpan"),

              mhc_dt[!netMHCII %>% is.na &
              pep_length == 15,
                .SD, .SDcols = c("netMHCII", "peptide", "pep_length")] %>%
              data.table::copy %>%
              data.table::setkey(netMHCII, pep_length) %>%
              write_peps("netMHCII"),

              mhc_dt[!netMHCIIpan %>% is.na &
              pep_length == 15,
                .SD, .SDcols = c("netMHCIIpan", "peptide", "pep_length")] %>%
              data.table::copy %>%
              data.table::setkey(netMHCIIpan, pep_length) %>%
              write_peps("netMHCIIpan")
                                     ))
        }

  # generate commands

    dtfn[, command :=
      paste(
            type,
            "-p",
            "-l", pep_length,
            "-a", allele,
            "-f", filename)
        ]
    dtfn[type == "netMHCIIpan", command :=
      paste(
            type,
            "-inptype 1",
            "-length", pep_length,
            "-a", allele,
            "-f", filename)
        ]

  # run commands

      # mhcflurry
      system("mhcflurry-predict mhcflurry_input.csv > mhcflurry_output.csv")

      # netMHC
      ag_out_raw <- parallel::mclapply(
         dtfn[, command],
          function(command){
          # run command
           es <- try(system(command, intern = TRUE))

          # if error, return empty dt
            if (es %>% class %>% .[1] == "try-error") {
              return(data.table::data.table(status = 1, command = command))
            }

          # parse results
            # isolate table and header

              dtl <- es %exclude% "^\\#|----|Number|Distance|threshold|version|^$"
              dtn <- dtl %include%
                "[Pp]eptide" %>%
                strsplit("[ ]+") %>%
                unlist %exclude% "^$|^Bind|^Level"

              dt <- dtl[2:length(dtl)] %>%
                stringr::str_replace_all(">", " ") %>%
                stringr::str_replace_all("=", " ") %>%
                stringr::str_replace_all("<", " ") %>%
                stringr::str_replace_all("(SB|WB)", "  ") %>%
                data.table::tstrsplit("[ ]+") %>%
                data.table::as.data.table

          # some output has a null first column
             if (dt$V1 %>% paste(collapse = "") == "") dt$V1 <- NULL

          # apply names to data table
            dt %>% data.table::setnames(dt %>% names, dtn)

          # append command and exit status
            dt$command <- command
            dt$status <- 0

          # set the program type from the command that was run
            ptype <- command %>% stringr::str_extract("net[A-Za-z]+")

          # make netMHC names consistent
            if ("Peptide" %chin% (dt %>% names)) dt %>% data.table::setnames("Peptide", "peptide")
            if ("Affinity(nM)" %chin% (dt %>% names)) dt %>% data.table::setnames("Affinity(nM)", "affinity(nM)")
            if ("Aff(nM)" %chin% (dt %>% names)) dt %>% data.table::setnames("Aff(nM)", "affinity(nM)")
            if ("HLA" %chin% (dt %>% names)) dt %>% data.table::setnames("HLA", "allele")
            if ("Allele" %chin% (dt %>% names)) dt %>% data.table::setnames("Allele", "allele")
            if ("Pos" %chin% (dt %>% names)) dt %>% data.table::setnames("Pos", "pos")
            if ("Icore" %chin% (dt %>% names)) dt %>% data.table::setnames("Icore", "icore")
            if ("iCore" %chin% (dt %>% names)) dt %>% data.table::setnames("iCore", "icore")
            if ("Core" %chin% (dt %>% names)) dt %>% data.table::setnames("Core", "core")

          # fix netMHCpan allele output to match input
            if (command %like% "netMHCpan"){
              dt[, allele := allele %>% stringr::str_replace(fixed("*"), "")]
            }

          # set unique column names based on program
            data.table::setnames(dt, dt %>% names %exclude% "allele|peptide",
                                 paste0((dt %>% names %exclude% "allele|peptide"), "_", ptype))

          # name allele column for merge
          dt %>%
              data.table::setnames("allele", ptype)

            return(dt)
                     })

  saveRDS(ag_out_raw, "ag_out_raw.RDS")

  # merge netMHC by program type
    progl <- lapply(ag_out_raw %>% seq_along, function(dti) {
     ag_out_raw[[dti]]$command[1] %>% stringr::str_extract("net[A-Za-z]+")
     })

    netmprogs <- progl %>% unique %>% unlist

    for (ptype in netmprogs) {

      mhc_dt <- merge(mhc_dt, ag_out_raw[(progl == ptype) %>% which] %>%
            data.table::rbindlist, by = c("peptide", ptype), all.x = TRUE)

    }

  # merge additional rows for mhcflurry
   fdt <- data.table::fread("mhcflurry_output.csv") %>%
      data.table::setnames("allele", "MHC")

   mhc_dt <- merge(mhc_dt, fdt, by = c("peptide", "MHC"), all.x = TRUE)

  # calculatte netMHC consensus score, preferring non-*net tools
    for (col in (mhc_dt %>% names %include% "aff|[Rr]ank|Consensus_scores")) {
      suppressWarnings(set(mhc_dt, j = col, value = mhc_dt[, get(col) %>% as.numeric]))
    }

    mhc_dt[, Consensus_scores := c(c(`affinity(nM)_netMHC`,
                                  `affinity(nM)_netMHCII`,
                                  `affinity(nM)_netMHCpan`,
                                  `affinity(nM)_netMHCIIpan`) %>%
            stats::na.omit, NA) %>% data.table::first, by = 1:nrow(mhc_dt)]

    mhc_dt[, Consensus_scores := c(Consensus_scores,
                                   mhcflurry_prediction) %>%
                                  mean(na.omit = TRUE),
                                  by = 1:nrow(mhc_dt)]

  # remove unpredicted peptides
    mhc_dt <- mhc_dt[!Consensus_scores %>% is.na]

  # calculate DAI
    data.table::setkey(mhc_dt, pep_type, variant_uuid, MHC, pep_length, peptide_index)

    mhc_dt %<>% .[, DAI := .[pep_type == "wt", Consensus_scores] /
                    .[pep_type == "mut", Consensus_scores]]

}

return(mhc_dt)

}


## -------- garnish_predictions_worker
#' Parallelized worker function for garnish_predictions
#'
#' @param basepep_dt Data table. Input data frame from garnish_predictions.
#'
#' @return A data table summary of neoepitope by sample_id.
#' @export garnish_predictions_worker

garnish_predictions_worker <- function(basepep_dt) {


# get unique wt and mut prediction frames (base peptide, mut location)

  basepep_dt <- data.table::rbindlist(list(
                          basepep_dt %>%
                          data.table::copy %>%
                          .[, pep_base := pep_wt] %>%
                          .[, pep_type := "wt"],
                          basepep_dt %>%
                          data.table::copy %>%
                          .[, pep_base := pep_mut] %>%
                          .[, pep_type := "mut"]
                           ))
  basepep_dt[, c("pep_mut", "pep_wt") := NULL]
  basepep_dt %<>% unique

## --- Write IEDB commands and make peptides for each sample

peptide_dt <- parallel::mclapply(1:nrow(basepep_dt), function(x){

  tryCatch({

  ## --- Write peptide fragments

    peptide_dt <- lapply(1:nrow(basepep_dt), function(n) {
      back_trun <- FALSE
        # for every peptide length
            peptide_dt <- lapply((15:8), function(pl) {

            mut_frag_t <- basepep_dt$pep_base[n] %>% strsplit("") %>% unlist
            mut_frag_loc <- basepep_dt$mut_pep_loc[n]

      # if the peptide is long enough
      if ((mut_frag_t %>% length) >= pl) {

          # re-register peptide i the mutant location is not centered due to back truncation
            while (mut_frag_loc > pl){
            mut_frag_t <- mut_frag_t[-1]
            mut_frag_loc <- mut_frag_loc -1
            back_trun <- TRUE
            }

          # re-register peptide if the mutant location is not centered due to front truncation
            while (mut_frag_loc+pl-1 < mut_frag_t %>% length){
            mut_frag_t <- mut_frag_t[-length(mut_frag_t)]
            }

          # reduce the peptide size until it matches the peptide length
            while (length(mut_frag_t) > ((2 * pl) - 1)) {

          # remove extra peptide from front or back depending on whether peptide was back-truncated
          if (!back_trun){
            if (mut_frag_loc != 1) {
                mut_frag_t <- mut_frag_t[-1]
                mut_frag_loc <- mut_frag_loc-1
              }
                if (length(mut_frag_t) == ((2 * pl) - 1)) next
                mut_frag_t <- mut_frag_t[-length(mut_frag_t)]
          }
          if (back_trun){
                  mut_frag_t <- mut_frag_t[-length(mut_frag_t)]
                  if (length(mut_frag_t) == ((2 * pl) - 1)) next
              if (mut_frag_loc != 1) {
                  mut_frag_t <- mut_frag_t[-1]
                  mut_frag_loc <- mut_frag_loc-1
          }
          }
          }

      # create n_mers wrapped in sync to prevent zoo::rollapply stdout
      sink(file = "/dev/null")
      peps <- zoo::rollapply(mut_frag_t, pl, print, partial = FALSE, align = "left") %>%
                            apply(., 1, function(pmr){
                            paste(pmr, sep = "", collapse = "")
                            })
      sink()

  # return a data table of peptides containing
  return(
         data.table::data.table(peptide = peps,
                                peptide_index = 1:length(peps),
                                pep_base = basepep_dt$pep_base[n],
                                pep_type = basepep_dt$pep_type[n],
                                variant_uuid = basepep_dt$variant_uuid[n]
                                ) %>%
         .[, pep_length := peps %>% nchar]
         )

  }}) %>% data.table::rbindlist

  return(peptide_dt %>% unique)

  }) %>% data.table::rbindlist %>% unique

      return(peptide_dt)
  }, error = function(e) {
        print(cat("ERROR :", conditionMessage(e), "\n"))
        return(NULL)
        })

  }) %>% data.table::rbindlist %>% unique
return(peptide_dt)
}



## -------- garnish_summary
#' Summarize epitope prediction.
#'
#' Calculate neoepitope summary statistics over samples.
#'
#' @param dt Data table. Prediction output from garnish_predictions.
#'
#' @export garnish_summary

garnish_summary <- function(dt){

dt <- dt[DAI != Inf & DAI != -Inf & Consensus_scores != Inf & Consensus_scores != -Inf]


  # function to sum top values of a numeric vector
    sum_top_v <- function(x, value = 3){

      x %<>% sort %>% rev
      return(sum(x[1:value]))
    }

dt %>% data.table::setkey(sample_id)


dtn <- parallel::mclapply(dt[, sample_id %>% unique], function(id){

  dt <- dt[sample_id == id]

    return(
        data.table::data.table(
        sample_id = id,
        priority_neos = dt[Consensus_scores < 50 & DAI > 10] %>% nrow,
        classic_neos = dt[Consensus_scores < 50] %>% nrow,
        alt_neos = dt[Consensus_scores < 5000 & DAI > 10] %>% nrow,
        alt_neos_top = dt[Consensus_scores < 5000, DAI %>% sum_top_v],
        classic_neos_top = dt[Consensus_scores < 5000, (1/Consensus_scores) %>% sum_top_v],
        binders = dt[Consensus_scores < 5000] %>% nrow,
        peptides = dt[pep_type == "mut", peptide %>% unique] %>% length,
        predictions = dt[pep_type == "mut"] %>% nrow
        ))

  }) %>% data.table::rbindlist

return(dtn)
}


# global variables

 utils::globalVariables(c(":=", ".", "%<>%", "%>%", "aa_convert", "aa_mutation", "affinity(nM)_netMHC", "affinity(nM)_netMHCII", "affinity(nM)_netMHCIIpan", "affinity(nM)_netMHCpan", "allele", "%chin%", "command", "Consensus_scores", "DAI", "dna_change", "effect_type", "%exclude%", "filename", "FILTER", "frag_begin", "frag_end", "gene_affected", "%include%", "INFO", "%like%", "MHC", "mhcflurry_prediction", "mutant_loc", "mutate_from", "mutate_to", "mut_pep_loc", "netMHC", "netMHCII", "netMHCIIpan", "netMHCpan", "pep_base", "pep_db_aa_from", "pep_db_aa_length", "pep_length", "pep_mut", "peptide", "peptide_index", "pep_type", "pep_uuid", "pep_wt", "protein_change", "protein_coding", "protein_index", "sample_id", ".SD", "se", "set", "transcript_affected", "transcript_affected_v", "type", "uuid", "variant_uuid", "vcf_type"))