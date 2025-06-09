#' This function reads OxBS and BS metadata (csv files) and iDAT files.
#' Performs normalization, CHAMP filtering, and offers as optional analyses sex, smoking, age and cell‐proportion prediction.
#' Estimates 5hmC, and returns only 'filtered_hmC' and 'phenotype_df_bs'.
#' Intermediate objects are removed with 'rm()' and 'gc()' to free memory.
#'
#' @param ox_file Path to OxBS sample‐sheet CSV (requires Sample_Name, Array, Slide, iDAT_PATH, status).
#' @param bs_file Path to BS sample‐sheet CSV (requires same columns).
#' @param annotation_array Array name for minfi annotation (default: "IlluminaHumanMethylationEPICv2").
#' @param annotation_version Genome annotation version (default: "20a1.hg38").
#' @param normalization One of "NOOB", "FUNORM", or "RAW" (default: "NOOB").
#' @param champfilter_arraytype_bs Array type for BS CHAMP filtering (default: "EPICv2").
#' @param champfilter_ProbeCutoff_bs Numeric [0,1] for BS CHAMP filtering (default: 0.01).
#' @param champfilter_arraytype_ox Array type for OxBS CHAMP filtering (default: "EPICv2").
#' @param champfilter_ProbeCutoff_ox Numeric [0,1] for OxBS CHAMP filtering (default: 0.01).
#' @param file_inaccuracies Path to probe inaccuracies CSV or NULL (default: NULL).
#' @param low_variance_threshold_hmc Numeric ≥ 0 for low‐variance 5hmC filtering (default: 0).
#' @param predictSex Logical to predict sex (default: FALSE).
#' @param predictSmoking Logical to predict smoking score (default: FALSE).
#' @param predictAge Logical to predict DNAm age (default: FALSE).
#' @param calculateCellPropPCs Logical to estimate cell proportions-PCs (default: FALSE).
#' @param plotCellProps Logical to plot cell proportions (default: FALSE).
#' @param plotPCA Logical to plot PCA explained variance and SVD (default: FALSE).
#' @param plotSVD Logical to plot SVD (default: FALSE).
#' @param plotHmC Logical to plot hydroxymethylation density (default: FALSE).
#' @param output_dir Directory to write outputs (default: getwd()).
#' @return A named list with elements:
#'   \item{phenotype_df_bs}{Data frame of BS phenotype (with optional columns)}
#'   \item{filtered_hmC}{5hmC long‐format after low‐variance filtering}
#'
#' @export
preprocess_hydroxymethylation_data <- function(
    ox_file,
    bs_file,
    annotation_array           = "IlluminaHumanMethylationEPICv2",
    annotation_version         = "20a1.hg38",
    normalization              = "NOOB",
    champfilter_arraytype_bs   = "EPICv2",
    champfilter_ProbeCutoff_bs = 0.01,
    champfilter_arraytype_ox   = "EPICv2",
    champfilter_ProbeCutoff_ox = 0.01,
    file_inaccuracies          = NULL,
    low_variance_threshold_hmc = 0,
    predictSex                 = FALSE,
    predictSmoking             = FALSE,
    predictAge                 = FALSE,
    calculateCellPropPCs       = FALSE,
    plotCellProps              = FALSE,
    plotPCA                    = FALSE,
    plotSVD                    = FALSE,
    plotHmC                    = FALSE,
    output_dir                 = getwd()
) {
  # Check dependencies
  required_pkgs <- c(
    "minfi", "sesame", "ChAMP", "EpiSmokEr", "wateRmelon",
    "ggplot2", "viridis", "reshape2", "MLML2R", "FlowSorted.Blood.EPIC"
  )
  missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(missing_pkgs) > 0) {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  }
  
  # Validate parameters
  normalization <- toupper(normalization)
  if (!normalization %in% c("NOOB", "FUNORM", "RAW")) {
    stop("`normalization` must be one of 'NOOB', 'FUNORM', or 'RAW'")
  }
  if (!is.numeric(champfilter_ProbeCutoff_bs) ||
      champfilter_ProbeCutoff_bs < 0 || champfilter_ProbeCutoff_bs > 1) {
    stop("`champfilter_ProbeCutoff_bs` must be numeric between 0 and 1")
  }
  if (!is.numeric(champfilter_ProbeCutoff_ox) ||
      champfilter_ProbeCutoff_ox < 0 || champfilter_ProbeCutoff_ox > 1) {
    stop("`champfilter_ProbeCutoff_ox` must be numeric between 0 and 1")
  }
  if (!is.numeric(low_variance_threshold_hmc) || low_variance_threshold_hmc < 0) {
    stop("`low_variance_threshold_hmc` must be numeric ≥ 0")
  }
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  
  # Function: read & validate the input files for the right columns
  .read_and_validate <- function(path, required_cols) {
    if (!file.exists(path)) stop("File not found: ", path)
    df <- read.csv(path, stringsAsFactors = FALSE)
    missing <- setdiff(required_cols, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf("File '%s' missing columns: %s", path, paste(missing, collapse = ", ")))
    }
    df
  }
  
  #  1) Read & validate metadata
  message("[1/14] Reading and validating metadata")
  required_cols <- c("Sample_Name", "Array", "Slide", "iDAT_PATH", "status")
  ox_df <- .read_and_validate(ox_file, required_cols)
  bs_df <- .read_and_validate(bs_file, required_cols)
  
  #  2) Prepare OxBS iDAT files 
  message("[2/14] Reading OxBS iDAT files")
  colnames(ox_df)[colnames(ox_df) == "Sample_Name"] <- "ID"
  ox_df$Sample_Name <- paste0(ox_df$Slide, "_", ox_df$Array)
  ox_df$Basename   <- file.path(ox_df$iDAT_PATH, ox_df$Slide, ox_df$Sample_Name)
  ox_df <- ox_df[order(ox_df$ID), ]
  ox_rg <- minfi::read.metharray.exp(targets = ox_df, force = TRUE, recursive = TRUE)
  ox_rg@annotation <- c(array = annotation_array, annotation = annotation_version)
  ox_rg <- ox_rg[, !is.na(minfi::pData(ox_rg)$Sample_Name)]
  rm(ox_df); gc()
  
  #  3) Prepare BS iDAT files 
  message("[3/14] Reading BS iDAT files")
  colnames(bs_df)[colnames(bs_df) == "Sample_Name"] <- "ID"
  bs_df$Sample_Name <- bs_df$Basename
  bs_df$Basename    <- file.path(bs_df$iDAT_PATH, bs_df$Slide, bs_df$Sample_Name)
  bs_df <- bs_df[order(bs_df$ID), ]
  bs_rg <- minfi::read.metharray.exp(targets = bs_df, force = TRUE, recursive = TRUE)
  bs_rg@annotation <- c(array = annotation_array, annotation = annotation_version)
  bs_rg <- bs_rg[, !is.na(minfi::pData(bs_rg)$Sample_Name)]
  rm(bs_df); gc()
  
  #  4) Optional: Predict sex 
  if (predictSex) {
    message("[4/14] Predicting sex from BS data")
    tryCatch({
      bs_raw    <- minfi::preprocessRaw(bs_rg)
      Rset_raw  <- minfi::ratioConvert(bs_raw, what = "both", keepCN = TRUE)
      GRset_raw <- minfi::mapToGenome(Rset_raw)
      GRset_raw <- minfi::addSnpInfo(GRset_raw)
      predicted_sex <- minfi::getSex(GRset_raw)
      minfi::pData(bs_rg)$predicted_sex <- predicted_sex$predictedSex
      rm(bs_raw, Rset_raw, GRset_raw, predicted_sex); gc()
    }, error = function(e) {
      warning("Sex prediction failed: ", e$message)
    })
  }
  
  #  5) Normalization 
  message("[5/14] Normalizing data using: ", normalization)
  if (normalization == "NOOB") {
    bs_norm <- minfi::preprocessNoob(bs_rg)
    ox_norm <- minfi::preprocessNoob(ox_rg)
  } else if (normalization == "FUNORM") {
    bs_norm <- minfi::preprocessFunnorm(bs_rg)
    ox_norm <- minfi::preprocessFunnorm(ox_rg)
  } else {
    bs_norm <- minfi::preprocessRaw(bs_rg)
    ox_norm <- minfi::preprocessRaw(ox_rg)
  }
  bs_norm <- bs_norm[, order(colnames(bs_norm))]
  ox_norm <- ox_norm[, order(colnames(ox_norm))]
  shared  <- intersect(rownames(bs_norm), rownames(ox_norm))
  bs_norm <- bs_norm[shared, ]
  ox_norm <- ox_norm[shared, ]
  
  #  6a) CHAMP filtering on BS 
  message("[6/14] CHAMP filtering on BS")
  detP_bs <- minfi::detectionP(bs_rg)
  pd_bs   <- as.data.frame(minfi::pData(bs_norm))
  pd_bs$Sample_Name <- rownames(pd_bs)
  detP_bs <- detP_bs[order(rownames(detP_bs)), , drop = FALSE]
  detP_bs <- detP_bs[, order(colnames(detP_bs)), drop = FALSE]
  beta_bs <- minfi::getBeta(bs_norm)
  beta_bs <- beta_bs[order(rownames(beta_bs)), order(colnames(beta_bs))]
  library(ChAMP)
  champRes_bs <- champ.filter(
    beta           = beta_bs,
    pd             = pd_bs,
    arraytype      = champfilter_arraytype_bs,
    ProbeCutoff    = champfilter_ProbeCutoff_bs,
    detP           = detP_bs,
    filterMultiHit = TRUE,
    filterXY       = FALSE
  )
  #  6b) Filter XY for BS
  annot <- minfi::getAnnotation(bs_norm)
  autosomal_probes <- rownames(annot)[!annot$chr %in% c("chrX", "chrY")]
  autosomal_probes <- intersect(autosomal_probes, rownames(champRes_bs$beta))
  champRes_bs_beta <- champRes_bs$beta[autosomal_probes, , drop = FALSE]
  bs_filtered <- bs_norm[rownames(champRes_bs_beta), ]
  rm(detP_bs, bs_rg, pd_bs, beta_bs, champRes_bs,annot,autosomal_probes,champRes_bs_beta); gc()
  
  #  7a) CHAMP filtering on OxBS 
  message("[7/14] CHAMP filtering on OxBS")
  detP_ox <- minfi::detectionP(ox_rg)
  pd_ox   <- as.data.frame(minfi::pData(ox_norm))
  pd_ox$Sample_Name <- rownames(pd_ox)
  detP_ox <- detP_ox[order(rownames(detP_ox)), , drop = FALSE]
  detP_ox <- detP_ox[, order(colnames(detP_ox)), drop = FALSE]
  beta_ox <- minfi::getBeta(ox_norm)
  beta_ox <- beta_ox[order(rownames(beta_ox)), order(colnames(beta_ox))]
  champRes_ox <- ChAMP::champ.filter(
    beta           = beta_ox,
    pd             = pd_ox,
    arraytype      = champfilter_arraytype_ox,
    ProbeCutoff    = champfilter_ProbeCutoff_ox,
    detP           = detP_ox,
    filterMultiHit = TRUE,
    filterXY       = FALSE
  )
  # 7b) Filter XY for OxBS
  annot <- minfi::getAnnotation(ox_norm)
  autosomal_probes <- rownames(annot)[!annot$chr %in% c("chrX", "chrY")]
  autosomal_probes <- intersect(autosomal_probes, rownames(champRes_ox$beta))
  champRes_ox_beta <- champRes_ox$beta[autosomal_probes, , drop = FALSE]
  ox_filtered <- ox_norm[rownames(champRes_ox_beta), ]
  rm(detP_ox, ox_rg, pd_ox, beta_ox, champRes_ox,annot,autosomal_probes,champRes_ox_beta); gc()
  
  #  8) Remove probe inaccuracies 
  if (!is.null(file_inaccuracies) && file.exists(file_inaccuracies)) {
    message("[8/14] Filtering probes via inaccuracies file")
    mapping_inacc <- read.csv(file_inaccuracies, stringsAsFactors = FALSE)
    if ("IlmnID" %in% colnames(mapping_inacc)) {
      bs_filtered <- bs_filtered[!rownames(bs_filtered) %in% mapping_inacc$IlmnID, ]
      ox_filtered <- ox_filtered[!rownames(ox_filtered) %in% mapping_inacc$IlmnID, ]
    } else {
      warning("Inaccuracies file missing 'IlmnID'; skipping that filter.")
    }
    rm(mapping_inacc); gc()
  }
  
  #  9) Build phenotype DataFrame for BS 
  message("[9/14] Building phenotype DataFrame for BS")
  phenotype_df_bs <- as.data.frame(minfi::pData(bs_filtered))
  phenotype_df_bs$Sample_Name <- rownames(phenotype_df_bs)
  phenotype_df_bs$Basename    <- sub("^.*/", "", phenotype_df_bs$Sample_Name)
  
  # 10) Optional: Cell proportions-PCs
  if (calculateCellPropPCs) {
    message("[10/14] Estimating cell proportions and computing PCs")
    library(FlowSorted.Blood.EPIC)
    library(sesame)
    library(minfi)
    
    Beta_c_bs <- sesame::betasCollapseToPfx(minfi::getBeta(bs_filtered))
    Beta_c_ox <- sesame::betasCollapseToPfx(minfi::getBeta(ox_filtered))
    
    idol_probes <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
    comp_table  <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable
    opts_bs     <- intersect(idol_probes, rownames(Beta_c_bs))
    opts_ox     <- intersect(idol_probes, rownames(Beta_c_ox))
    
    cell_bs <- as.data.frame(
      FlowSorted.Blood.EPIC::projectCellType_CP(
        Y           = Beta_c_bs[opts_bs, ],
        coefWBC     = comp_table[opts_bs, ],
        contrastWBC = NULL,
        nonnegative = TRUE,
        lessThanOne = FALSE
      )
    )
    cell_ox <- as.data.frame(
      FlowSorted.Blood.EPIC::projectCellType_CP(
        Y           = Beta_c_ox[opts_ox, ],
        coefWBC     = comp_table[opts_ox, ],
        contrastWBC = NULL,
        nonnegative = TRUE,
        lessThanOne = FALSE
      )
    )
    cell_bs$Basename <- rownames(cell_bs)
    cell_ox$Basename <- rownames(cell_ox)
    
    rm(Beta_c_bs, Beta_c_ox); gc()
    
    pca_bs   <- stats::prcomp(cell_bs[, -ncol(cell_bs)], scale. = TRUE)
    expl_var <- pca_bs$sdev^2 / sum(pca_bs$sdev^2)
    
    if (plotCellProps) {
      long_bs <- reshape2::melt(
        cell_bs,
        id.vars       = "Basename",
        variable.name = "CellType",
        value.name    = "Proportion"
      )
      p_cell <- ggplot2::ggplot(long_bs, ggplot2::aes(x = Basename, y = Proportion, fill = CellType)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        viridis::scale_fill_viridis(discrete = TRUE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::labs(
          title = "Cell Proportions by Sample (BS)",
          x     = "Sample",
          y     = "Proportion"
        )
      plot_title <- p_cell$labels$title
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(plot_title, ".png")),
        plot     = p_cell,
        width    = 8, height = 5
      )
    }
    
    if (plotPCA) {
      plot_title <- "Explained Variance (BS)"
      png(
        filename = file.path(output_dir, paste0(plot_title, ".png")),
        width    = 800, height = 600
      )
      graphics::par(fg = "white")
      graphics::barplot(expl_var[1:6], col = viridis::viridis(6), main = plot_title)
      dev.off()
    }
    
    PC_df <- as.data.frame(pca_bs$x[, 1:6, drop = FALSE])
    PC_df$Sample_Name <- rownames(PC_df)
    phenotype_df_bs <- merge(
      phenotype_df_bs,
      PC_df[, c("Sample_Name", paste0("PC", 1:6))],
      by    = "Sample_Name",
      all.x = TRUE,
      sort  = FALSE
    )
    rm(cell_bs, cell_ox, pca_bs, expl_var, PC_df); gc()
  }
  
  # 11) Optional: Smoking score
  if (predictSmoking) {
    message("[11/14] Predicting smoking scores")
    tryCatch({
      if (!requireNamespace("EpiSmokEr", quietly = TRUE)) {
        stop("Please install EpiSmokEr")
      }
      library(EpiSmokEr)
      data("Illig_data", package = "EpiSmokEr")
      rownames(phenotype_df_bs) <- phenotype_df_bs$Sample_Name
      if (!"Sex" %in% colnames(phenotype_df_bs) &&
          "predicted_sex" %in% colnames(phenotype_df_bs)) {
        phenotype_df_bs$Sex <- phenotype_df_bs$predicted_sex
      }
      Beta_c_bs2 <- sesame::betasCollapseToPfx(minfi::getBeta(bs_filtered))
      smoke_res  <- EpiSmokEr::epismoker(
        Beta_c_bs2,
        samplesheet = phenotype_df_bs,
        method      = "all"
      )
      phenotype_df_bs$smokingScore <- smoke_res$smokingScore
      rm(Beta_c_bs2, smoke_res); gc()
    }, error = function(e) {
      warning("Smoking score prediction failed: ", e$message)
    })
  }
  
  # 12) Optional: DNAm age (Horvath)
  if (predictAge) {
    message("[12/14] Predicting DNAm age (Horvath)")
    tryCatch({
      beta_mat <- minfi::getBeta(bs_filtered)
      rownames(beta_mat) <- sub("_.*$", "", rownames(beta_mat))
      data(coef, package = "wateRmelon")
      library(wateRmelon)
      data("age_coefficients", package = "wateRmelon")
      clock.probes <- names(coef)[-1]
      common       <- intersect(clock.probes, rownames(beta_mat))
      if (length(common) < 1) stop("No Horvath clock probes found")
      beta_sub    <- beta_mat[common, , drop = FALSE]
      DNAm_ages   <- wateRmelon::agep(
        betas   = beta_sub,
        method  = "horvath",
        verbose = TRUE
      )
      phenotype_df_bs$Predicted_Age <- DNAm_ages$horvath.age
      
      rm(beta_mat, DNAm_ages)
      gc()
    }, error = function(e) {
      warning("DNAm age prediction failed: ", e$message)
    })
  }
  
  
  
  #  Write phenotype csv
  phenotype_fname <- file.path(output_dir, "phenotype_table.csv")
  message("Writing phenotype metadata to: ", phenotype_fname)
  write.csv(phenotype_df_bs, phenotype_fname, row.names = FALSE)
  
  #  SVD on BS beta (optional) 
  if (plotSVD) {
    message("Saving SVD plot")
    tryCatch({
      myNorm <- as.data.frame(minfi::getBeta(bs_filtered))
      phenotype_df_bs$Slide <- factor(phenotype_df_bs$Slide)
      ChAMP::champ.SVD(
        beta       = myNorm,
        rgSet      = NULL,
        pd         = phenotype_df_bs,
        resultsDir = paste0(output_dir,"/"),
        PDFplot    = TRUE,
        Rplot      = TRUE
      )
      
      rm(myNorm); gc()
    }, error = function(e) {
      warning("SVD plotting failed: ", e$message)
    })
  }
  
  #  13) 5hmC estimation via MLML
  message("[13/14] Estimating 5hmC via MLML")
  UC_bs <- as.data.frame(minfi::getUnmeth(bs_filtered))
  MC_bs <- as.data.frame(minfi::getMeth(bs_filtered))
  MC_ox <- as.data.frame(minfi::getMeth(ox_filtered))
  UC_ox <- as.data.frame(minfi::getUnmeth(ox_filtered))
  common_probes <- Reduce(
    intersect,
    list(rownames(MC_bs), rownames(UC_bs), rownames(MC_ox), rownames(UC_ox))
  )
  MC_bs <- MC_bs[common_probes, , drop = FALSE]
  UC_bs <- UC_bs[common_probes, , drop = FALSE]
  MC_ox <- MC_ox[common_probes, , drop = FALSE]
  UC_ox <- UC_ox[common_probes, , drop = FALSE]
  MCbs_mat <- as.matrix(MC_bs)
  UCbs_mat <- as.matrix(UC_bs)
  MCox_mat <- as.matrix(MC_ox)
  UCox_mat <- as.matrix(UC_ox)
  library(MLML2R)
  results_mlml <- MLML(
    T.matrix = MCbs_mat,
    U.matrix = UCbs_mat,
    L.matrix = UCox_mat,
    M.matrix = MCox_mat
  )
  rm(UC_bs, MC_bs, MC_ox, UC_ox, MCbs_mat, UCbs_mat, MCox_mat, UCox_mat); gc()
  
  # 5hmC long format
  hmC_mat <- as.data.frame(results_mlml$hmC)
  hmC_mat$CpG_ID <- rownames(results_mlml$hmC)
  long_df <- reshape2::melt(
    hmC_mat,
    id.vars       = "CpG_ID",
    variable.name = "Sample_Name",
    value.name    = "hmC_Value"
  )
  long_df <- merge(long_df, phenotype_df_bs, by = "Sample_Name", all.x = TRUE, sort = FALSE)
  rm(hmC_mat); gc()
  
  if (plotHmC) {
    p_hmC <- ggplot2::ggplot(long_df, ggplot2::aes(x = hmC_Value, colour = Sample_Name)) +
      ggplot2::geom_density() +
      viridis::scale_color_viridis(discrete = TRUE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Hydroxymethylation Density by Sample",
        x     = "hmC Value",
        y     = "Density"
      )
    plot_title <- p_hmC$labels$title
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(plot_title, ".png")),
      plot     = p_hmC,
      width    = 8, height = 5
    )
  }
  
  #  Low‐variance filtering & write csv
  probe_vars  <- apply(results_mlml$hmC, 1, var)
  keep_probe  <- setNames(probe_vars >= low_variance_threshold_hmc, names(probe_vars))
  filtered_hmC <- long_df[keep_probe[long_df$CpG_ID], ]
  rm(results_mlml); gc()
  hmC_fname <- file.path(output_dir, "filtered_hmC.csv")
  message("[14/14] Writing filtered 5hmC to: ", hmC_fname)
  write.csv(filtered_hmC, hmC_fname, row.names = FALSE)
  
  invisible(list(
    phenotype_df_bs = phenotype_df_bs,
    filtered_hmC    = filtered_hmC
  ))
}