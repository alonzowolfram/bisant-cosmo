#' Perform differential expression analysis via negative binomial modeling (no random spatial effect)
#'
#' @param cell_type_i cell type to perform differential expression analysis on
#' @param gene_name name of gene to perform differential expression analysis on
#' @param modeldat object returned by smiDE::xy_kmeans_clusters()
#' @param cts matrix of gene counts
#'
#' @return data frame with differential expression results
process_gene_ic <- function(cell_type_i, gene_name, modeldat, cts) {
  # Capture all output
  log_file_err <- here(output_dir_logs_de, "IC_differential-expression_err.log")
  log_file_out <- here(output_dir_logs_de, "IC_differential-expression_out.log")
  current_info <- paste0("cell type: ", cell_type_i, ", gene: ", gene_name)
  
  tryCatch({
    library(smiDE)
    
    # log_con <- file(log_file_out, open = "a")
    log_message <- paste0("\n", Sys.time(), " - ", current_info)
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    # flush(log_con)  # Force writing to disk immediately
    # close(log_con)  # Close the connection
    
    # Subset `modeldat` for the current `cell_type`
    modeldat_sub <- modeldat %>%
      dplyr::filter(cell_type == cell_type_i) %>%
      dplyr::left_join(otherct, by = "cell_ID")
    modeldat_sub$RankNormOtherct <- RankNorm(modeldat_sub$otherct_exp) #  smide::RankNorm(otherct)
    
    # Iterate over fixed effects
    results <- lapply(fixed_effect_vars, function(fixed_effect_var) {
      y <- cts[gene_name, modeldat_sub$cell_ID]
      fixed_effects <- modeldat_sub %>% dplyr::select(all_of(c(fixed_effect_var, "RankNormOtherct", "nCount_RNA", "k_cluster")))
      fixed_effects$b0 <- 1
      fixed_effects$y <- y
      
      # Define model formula
      mf <- as.formula(paste(
        "y ~ RankNormOtherct +", fixed_effect_var, "+ offset(log(nCount_RNA)) + (1 | k_cluster)"
      ))
      # Define zero-inflated model formula (i.e., what are excess zeros a function of?)
      zi_mf <- as.formula(paste("~ 1")) # If we make the # of 0s a function of fixed-effect variables, we might overparameterize, leading to model convergence issues.
      
      # Fit ZINB model with error handling
      # log_con <- file(log_file_out, open = "a")
      log_message <- paste0("\nFitting model for ", current_info, ", fixed-effect variable: ", fixed_effect_var)
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      # flush(log_con)  # Force writing to disk immediately
      # close(log_con)  # Close the connection
      res <- tryCatch({
        glmmTMB(
          data = fixed_effects,
          formula = mf,
          family = nbinom2 ,
          ziformula = zi_mf
        )
      }, error = function(e) {
        # log_con <- file(log_file_err, open = "a")
        log_message <- paste(Sys.time(), "- IC failed for", current_info, ", fixed-effect variable: ", fixed_effect_var, "\n Error:", e$message, "\n")
        cat(log_message, file = log_file_err, append = TRUE)
        flush.console()
        # flush(log_con)  # Force writing to disk immediately
        # close(log_con)  # Close the connection
        return(NULL)  # Return NULL on failure
      })
      
      if (is.null(res)) return(NULL)  # Skip failed models
      
      # Extract differential expression results
      message("Getting DE results")
      # View model summary
      summary(res)
      # Extract fixed effect estimates
      de_results <- summary(res)$coefficients$cond %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "FixedEffectLevels") %>%
        dplyr::filter(!(FixedEffectLevels %in% c("(Intercept)", "RankNormOtherct")))
      
      # View results
      print(de_results)
      
      # Compute pairwise differential expression
      # Compute estimated marginal means (EMMs)
      emm <- emmeans(res, as.formula(paste("~", fixed_effect_var)))
      
      # Perform pairwise contrasts (Tukey-adjusted)
      pairwise_results <- pairs(emm, adjust = "tukey")
      
      # View results
      print(pairwise_results)
      
      # Store results
      data.frame(
        gene = gene_name,
        cell_type = cell_type_i,
        fixed_effect = fixed_effect_var,
        contrast = summary(pairwise_results)$contrast,
        estimate = summary(pairwise_results)$estimate,
        SE = summary(pairwise_results)$SE,
        z_ratio = summary(pairwise_results)$z.ratio,
        p_val = summary(pairwise_results)$p.value
      )
    })
    
    return(do.call(rbind, results))
    
  }, error = function(e) {
    # log_con <- file(log_file_err, open = "a")
    log_message <- paste(Sys.time(), "- General error in process_gene_ic() for", current_info, ":", e$message, "\n")
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
    # flush(log_con)  # Force writing to disk immediately
    # close(log_con)  # Close the connection
    return(NULL)  # Return NULL on failure
  }, warning = function(w) {
    # log_con <- file(log_file_err, open = "a")
    log_message <- sprintf("\n%s - WARNING in process_gene_ic() for %s, %s,: %s", Sys.time(), cell_type_i, gene_name, conditionMessage(w))
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
    # flush(log_con)  # Force writing to disk immediately
    # close(log_con)  # Close the connection
    
    # Possible causes for "Model convergence problem; non-positive-definite Hessian matrix":
    # 1) Separation Issues / Quasi-Complete Separation
    # If a categorical predictor has levels that are perfectly predictive of the outcome (e.g., all zero or all high values for a certain group), the model struggles to estimate parameters.
    #
    # 2) Too Many Random Effects / Overparameterization
    # If you have too many random effect levels or they have little variation, the Hessian matrix can become singular.
    #
    # 3) Extreme Count Distributions (Sparse Data / Zero Inflation)
    # If a gene has too many zeros or very high dispersion, the model may not converge.
    #
    # 4) Correlation Between Predictors (Multicollinearity)
    # If your fixed effects are highly correlated, this can make it hard for the model to estimate coefficients.
    #
    # 5) Offset Issues 
    # If nCount_RNA has very small values or a skewed distribution, the offset log(nCount_RNA) can cause issues.
    
    # Possible causes for "false convergence":
    # 1) Poor scaling of variables (e.g., log-transformed data with extreme values).
    # 2) Model complexity (too many fixed/random effects).
    # 3) Poor choice of optimizer.
    
    # Possible causes for "singular convergence":
    # 1) A random effect is unnecessary (e.g., (1 | k_cluster) when k_cluster doesn't explain much variation).
    # 2) The data is imbalanced or sparse, making it difficult to estimate variance components.
    # 3) A random effect has very few levels, making it indistinguishable from residual variance.
  })
}

#' Perform differential expression analysis via negative binomial modeling (with GP-INLA random spatial effect)
#'
#' @param cell_type_i cell type to perform differential expression analysis on
#' @param gene_name name of gene to perform differential expression analysis on
#' @param modeldat object returned by smiDE::xy_kmeans_clusters()
#' @param cts matrix of gene counts
#'
#' @return data frame with differential expression results
process_gene_inla <- function(cell_type_i, gene_name, modeldat, cts) {
  inla.setOption(num.threads = 1)  # Limit to 1 thread per parallel worker. This is only for running in RStudio.
  
  # Capture all output
  # log_file_err <- here(output_dir_logs_de, "GP-INLA_differential-expression_err.log")
  # log_file_out <- here(output_dir_logs_de, "GP-INLA_differential-expression_out.log")
  log_file_err <- here(output_dir_logs_de, paste0("GP-INLA_differential-expression_err_", cell_type_i, "_", gene_name, ".log"))
  log_file_out <- here(output_dir_logs_de, paste0("GP-INLA_differential-expression_out_", cell_type_i, "_", gene_name, ".log"))
  # Write each worker's log to a separate file (i.e., serialize logging) to prevent deadlock from multiple workers writing to the same file.
  
  current_info <- paste0("cell type:", cell_type_i, ", gene:", gene_name)
  
  tryCatch({
    library(smiDE)
    
    log_message <- paste0(Sys.time(), " - ", current_info)
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    
    # Subset `modeldat` for the current `cell_type`
    modeldat_sub <- modeldat %>%
      dplyr::filter(cell_type == cell_type_i) %>%
      dplyr::left_join(otherct, by = "cell_ID")
    modeldat_sub$RankNormOtherct <- RankNorm(modeldat_sub$otherct_exp) #  smide::RankNorm(otherct)
    
    # Extract spatial coordinates
    log_message <- paste("\n", Sys.time(), "Extracting spatial coordinates")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    coords <- modeldat_sub[, c("sdimx", "sdimy")] %>% as.data.frame()
    dist_mat <- dist(coords)
    
    # Compute mesh parameters dynamically
    log_message <- paste("\n", Sys.time(), "Computing mesh parameters")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    cutoff_val <- min(as.matrix(dist_mat)[as.matrix(dist_mat) > 0]) / 2  
    bounding_box <- apply(coords, 2, range)
    spatial_extent <- max(bounding_box[2, ] - bounding_box[1, ])
    correlation_range <- spatial_extent * 0.2
    maxedge_val <- c(correlation_range / 5, correlation_range / 2)
    offset_val <- c(correlation_range / 2, correlation_range)
    
    # Create mesh
    log_message <- paste("\n", Sys.time(), "Creating mesh")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    mesh <- inla.mesh.2d(loc = coords, cutoff = cutoff_val, offset = offset_val, max.edge = maxedge_val)
    
    # Projection matrix
    log_message <- paste("\n", Sys.time(), "Creating projection matrix")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords))
    
    # Define SPDE model
    log_message <- paste("\n", Sys.time(), "Defining SPDE model")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    spde <- inla.spde2.pcmatern(
      mesh = mesh,
      alpha = 2, 
      prior.range = c(1, 0.5), 
      prior.sigma = c(0.07, 0.5)
    )
    
    # Iterate over fixed effects
    # results <- lapply(fixed_effect_vars, function(fixed_effect_var) {
    results <- list()
    for (fixed_effect_var in fixed_effect_vars) {
      y <- cts[gene_name, modeldat_sub$cell_ID]
      fixed_effects <- modeldat_sub %>% dplyr::select(all_of(c(fixed_effect_var, "RankNormOtherct", "nCount_RNA")))
      fixed_effects$b0 <- 1
      
      # Create INLA stack
      stk.e <- inla.stack(
        tag = "est",
        data = list(y = y),
        A = list(1, A),
        effects = list(fixed_effects, idx.u = 1:spde$n.spde)
      )
      
      # Define priors
      log_message <- paste("\n", Sys.time(), "Defining priors")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      pcprec <- list(hyper = list(theta = list(prior = 'pc.prec', param = c(1, .1))))
      prior_strength <- ifelse(length(y) < 500, 10, 1)  # Stronger prior for small datasets
      # Small datasets (e.g., <500 cells) need stronger priors to avoid overfitting.
      # Larger datasets allow more data-driven estimates, so weaker priors are better.
      control_fixed <- list(mean = 0, prec = list(default = prior_strength))
      
      # Define model formula
      mf <- as.formula(paste(
        "y ~ 0 + b0 + RankNormOtherct +", fixed_effect_var, "+ offset(log(nCount_RNA)) + f(idx.u, model = spde)"
      ))
      
      # Fit INLA model with error handling
      log_message <- paste("\n", Sys.time(), "Fitting model")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      # system.time({
        res <- tryCatch({
          inla(
            mf,
            family = "zeroinflatednbinomial0",
            control.family = pcprec,
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
            control.fixed = control_fixed,
            data = inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, A = inla.stack.A(stk.e))
          )
        }, error = function(e) {
          log_message <- paste(Sys.time(), "- INLA failed for", cell_type_i, gene_name, fixed_effect_var, "Error:", e$message, "\n")
          cat(log_message, file = log_file_err, append = TRUE)
          flush.console()
          return(NULL)  # Return NULL on failure
        })
      # })
      
      if (is.null(res)) return(NULL)  # Skip failed models
      
      # Extract differential expression results
      log_message <- paste("\n", Sys.time(), "Getting DE results")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      group_effects <- res$summary.fixed %>%
        tibble::rownames_to_column(var = "FixedEffectLevels") %>%
        dplyr::filter(!(FixedEffectLevels %in% c("b0", "RankNormOtherct")))
      rownames(group_effects) <- group_effects$FixedEffectLevels
      
      # Compute pairwise differential expression
      combns <- t(combn(group_effects$FixedEffectLevels, 2))
      de_results <- apply(combns, 1, function(combn) {
        level_A <- combn[1]
        level_B <- combn[2]
        # print(paste(level_A, level_B))
        diff_AB <- group_effects[level_B, "mean"] - group_effects[level_A, "mean"]
        # print(diff_AB)
        diff_AB_CI <- c(
          group_effects[level_B, "0.025quant"] - group_effects[level_A, "0.975quant"],
          group_effects[level_B, "0.975quant"] - group_effects[level_A, "0.025quant"]
        )
        
        # Compute probability of DE
        groupB_marginal <- res$marginals.fixed[[level_B]]
        prob_DE <- 1 - inla.pmarginal(0, groupB_marginal)
        
        # Store results
        data.frame(
          gene = gene_name,
          cell_type = cell_type_i,
          fixed_effect = fixed_effect_var,
          ref_level = level_A,
          comparison_level = level_B,
          direction = paste(level_B, "-", level_A),
          mean_ref = group_effects[level_A, "mean"],
          mean_comp = group_effects[level_B, "mean"],
          diff_mean = diff_AB,
          CI_lower_bound = diff_AB_CI[1],
          CI_upper_bound = diff_AB_CI[2],
          prob_DE = prob_DE
        )
      })
      
      # return(do.call(rbind, de_results))
      results[[fixed_effect_var]] <- do.call(rbind, de_results)
      
      # Clean up
      rm(res, stk.e, A, mesh, spde, fixed_effects, y, group_effects)
      gc()  # Force garbage collection
      
    } # End iteration over fixed effects
    # ) # From the lapply() which has been removed in favor of the for() loop.
    
    return(do.call(rbind, results))
  }, error = function(e) {
    log_message <- paste(Sys.time(), "- General error in process_gene for", cell_type_i, gene_name, ":", e$message, "\n")
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
    return(NULL)  # Return NULL on failure
  }, warning = function(w) {
    log_message <- sprintf("WARNING in process_gene_ic for %s, %s: %s", cell_type_i, gene_name, conditionMessage(w))
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
  })
}

process_gene_inla_test <- function(cell_type_i, gene_name, modeldat, cts) {
  library(R.utils)
  inla.setOption(num.threads = 1)
  
  log_file_err <- here(output_dir_logs_de, paste0("err_", cell_type_i, "_", gene_name, ".log"))
  log_file_out <- here(output_dir_logs_de, paste0("out_", cell_type_i, "_", gene_name, ".log"))
  done_flag <- here(output_dir_logs_de, paste0("done_", cell_type_i, "_", gene_name, ".txt"))
  
  tryCatch({
    cat(Sys.time(), " - Starting gene:", cell_type_i, gene_name, "\n", file = log_file_out, append = TRUE, flush = TRUE)
    
    modeldat_sub <- modeldat %>%
      dplyr::filter(cell_type == cell_type_i) %>%
      dplyr::left_join(otherct, by = "cell_ID")
    modeldat_sub$RankNormOtherct <- RankNorm(modeldat_sub$otherct_exp)
    
    coords <- modeldat_sub[, c("sdimx", "sdimy")] %>% as.data.frame()
    dist_mat <- dist(coords)
    cutoff_val <- min(as.matrix(dist_mat)[as.matrix(dist_mat) > 0]) / 2
    bounding_box <- apply(coords, 2, range)
    spatial_extent <- max(bounding_box[2, ] - bounding_box[1, ])
    correlation_range <- spatial_extent * 0.2
    maxedge_val <- c(correlation_range / 5, correlation_range / 2)
    offset_val <- c(correlation_range / 2, correlation_range)
    
    cat(Sys.time(), " - Creating mesh\n", file = log_file_out, append = TRUE, flush = TRUE)
    mesh <- inla.mesh.2d(loc = coords, cutoff = cutoff_val, offset = offset_val, max.edge = maxedge_val)
    A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords))
    spde <- inla.spde2.pcmatern(
      mesh = mesh, alpha = 2,
      prior.range = c(1, 0.5), prior.sigma = c(0.07, 0.5)
    )
    
    results <- list()
    for (fixed_effect_var in fixed_effect_vars) {
      cat(Sys.time(), " - Fitting model for:", fixed_effect_var, "\n", file = log_file_out, append = TRUE, flush = TRUE)
      
      y <- cts[gene_name, modeldat_sub$cell_ID]
      fixed_effects <- modeldat_sub %>%
        dplyr::select(all_of(c(fixed_effect_var, "RankNormOtherct", "nCount_RNA")))
      fixed_effects$b0 <- 1
      
      stk.e <- inla.stack(
        tag = "est",
        data = list(y = y),
        A = list(1, A),
        effects = list(fixed_effects, idx.u = 1:spde$n.spde)
      )
      
      pcprec <- list(hyper = list(theta = list(prior = 'pc.prec', param = c(1, .1))))
      prior_strength <- ifelse(length(y) < 500, 10, 1)
      control_fixed <- list(mean = 0, prec = list(default = prior_strength))
      mf <- as.formula(paste(
        "y ~ 0 + b0 + RankNormOtherct +", fixed_effect_var, "+ offset(log(nCount_RNA)) + f(idx.u, model = spde)"
      ))
      
      res <- tryCatch({
        withTimeout({
          inla(
            mf,
            family = "zeroinflatednbinomial0",
            control.family = pcprec,
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
            control.fixed = control_fixed,
            data = inla.stack.data(stk.e),
            control.predictor = list(compute = TRUE, A = inla.stack.A(stk.e))
          )
        }, timeout = 300, onTimeout = "silent")  # 5-minute cap
      }, error = function(e) {
        cat(Sys.time(), "- INLA failed:", e$message, "\n", file = log_file_err, append = TRUE)
        return(NULL)
      })
      
      if (is.null(res)) next
      
      group_effects <- res$summary.fixed %>% as.data.frame %>%
        tibble::rownames_to_column(var = "FixedEffectLevels") %>%
        dplyr::filter(!(FixedEffectLevels %in% c("b0", "RankNormOtherct")))
      rownames(group_effects) <- group_effects$FixedEffectLevels
      
      combns <- t(combn(group_effects$FixedEffectLevels, 2))
      de_results <- apply(combns, 1, function(combn) {
        level_A <- combn[1]
        level_B <- combn[2]
        diff_AB <- group_effects[level_B, "mean"] - group_effects[level_A, "mean"]
        diff_AB_CI <- c(
          group_effects[level_B, "0.025quant"] - group_effects[level_A, "0.975quant"],
          group_effects[level_B, "0.975quant"] - group_effects[level_A, "0.025quant"]
        )
        groupB_marginal <- res$marginals.fixed[[level_B]]
        prob_DE <- 1 - inla.pmarginal(0, groupB_marginal)
        
        data.frame(
          gene = gene_name,
          cell_type = cell_type_i,
          fixed_effect = fixed_effect_var,
          ref_level = level_A,
          comparison_level = level_B,
          direction = paste(level_B, "-", level_A),
          mean_ref = group_effects[level_A, "mean"],
          mean_comp = group_effects[level_B, "mean"],
          diff_mean = diff_AB,
          CI_lower_bound = diff_AB_CI[1],
          CI_upper_bound = diff_AB_CI[2],
          prob_DE = prob_DE
        )
      })
      
      results[[fixed_effect_var]] <- do.call(rbind, de_results)
    }
    
    # Cleanup
    rm(res, A, mesh, spde, fixed_effects, y, group_effects, de_results, combns, stk.e)
    gc()
    
    write(Sys.time(), file = done_flag)
    return(do.call(rbind, results))
  }, error = function(e) {
    cat(Sys.time(), "- General error:", e$message, "\n", file = log_file_err, append = TRUE)
    return(NULL)
  })
}

#' Perform differential expression analysis via negative binomial modeling (with GAM random spatial effect)
#'
#' @param cell_type_i cell type to perform differential expression analysis on
#' @param gene_name name of gene to perform differential expression analysis on
#' @param modeldat object returned by smiDE::xy_kmeans_clusters()
#' @param cts matrix of gene counts
#'
#' @return data frame with differential expression results
process_gene_gam <- function(cell_type_i, gene_name, modeldat, cts) {
  # Capture all output
  # log_file_err <- here(output_dir_logs_de, "GAM_differential-expression_err.log")
  # log_file_out <- here(output_dir_logs_de, "GAM_differential-expression_out.log")
  log_file_err <- here(output_dir_logs_de, paste0("GAM_differential-expression_err_", cell_type_i, "_", gene_name, ".log"))
  log_file_out <- here(output_dir_logs_de, paste0("GAM_differential-expression_out_", cell_type_i, "_", gene_name, ".log"))
  done_flag <- here(output_dir_logs_de, paste0("done_", cell_type_i, "_", gene_name, ".txt"))
  # Write each worker's log to a separate file (i.e., serialize logging) to prevent deadlock from multiple workers writing to the same file.
  
  current_info <- paste0("cell type:", cell_type_i, ", gene:", gene_name)
  
  tryCatch({
    library(mgcv)
    library(emmeans)
    
    log_message <- paste0(Sys.time(), " - ", current_info)
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    
    # Subset `modeldat` for the current `cell_type`
    modeldat_sub <- modeldat %>%
      dplyr::filter(cell_type == cell_type_i) %>%
      dplyr::left_join(otherct, by = "cell_ID")
    modeldat_sub$RankNormOtherct <- smiDE::RankNorm(modeldat_sub$otherct_exp) #  smide::RankNorm(otherct)
    
    # Extract spatial coordinates
    log_message <- paste("\n", Sys.time(), "Extracting spatial coordinates")
    cat(log_message, file = log_file_out, append = TRUE)
    flush.console()
    coords <- modeldat_sub[, c("sdimx", "sdimy")] %>% as.data.frame()
    dist_mat <- dist(coords)
    
    # Iterate over fixed effects
    results_list <- list()
    for (fixed_effect_var in fixed_effect_vars) {
      y <- cts[gene_name, modeldat_sub$cell_ID]
      fixed_effects <- modeldat_sub %>% dplyr::select(all_of(c(fixed_effect_var, "RankNormOtherct", "nCount_RNA")))
      fixed_effects$b0 <- 1
      
      # Create the training data frame (`gene_data`).
      gene_data <- cbind(y, modeldat_sub %>% dplyr::select(all_of(c(fixed_effect_var, "sdimx", "sdimy", "RankNormOtherct", "nCount_RNA"))))
      
      # Define model formula
      mf <- as.formula(paste(
        "y ~ RankNormOtherct +", fixed_effect_var, "+ offset(log(nCount_RNA)) + s(sdimx, sdimy)"
      ))
      
      # Fit spatial GAM with error handling
      log_message <- paste("\n", Sys.time(), "Fitting model")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      # system.time({
      res <- tryCatch({
        gam(mf, 
            data = gene_data, 
            family = nb())
      }, error = function(e) {
        log_message <- paste(Sys.time(), "- GAM fitting failed for", cell_type_i, gene_name, fixed_effect_var, "Error:", e$message, "\n")
        cat(log_message, file = log_file_err, append = TRUE)
        flush.console()
        return(NULL)  # Return NULL on failure
      })
      # })
      
      if (is.null(res)) return(NULL)  # Skip failed models
      
      # Extract differential expression results
      summary_res <- summary(res)
      log_message <- paste("\n", Sys.time(), "Getting DE results")
      cat(log_message, file = log_file_out, append = TRUE)
      flush.console()
      group_effects <- summary_res$p.table %>% as.data.frame %>% 
        tibble::rownames_to_column(var = "FixedEffectLevels") %>%
        dplyr::filter(!(FixedEffectLevels %in% c("(Intercept)", "RankNormOtherct")))
      rownames(group_effects) <- group_effects$FixedEffectLevels
      
      # Compute pairwise differential expression
      em <- emmeans(res, specs = as.formula(paste("~", fixed_effect_var)), data = gene_data) # "emmeans() sometimes can't correctly evaluate the reference grid when a GAM has non-parametric terms, offsets, or special data structures [such as smoothed terms]"--hence why we need to include `data = `
      pairs(em)  # gives all pairwise comparisons with estimates, SEs, p-values
      pairwise_results <- as.data.frame(pairs(em)) %>%
        dplyr::select(contrast, estimate, SE, df, p.value)
      # Correct for multiple testing
      pairwise_results$p.adj <- p.adjust(pairwise_results$p.value, method = "BH") 
        
      # Store results
      results_list[[fixed_effect_var]] <- cbind(
        data.frame(
          gene = gene_name,
          cell_type = cell_type_i,
          fixed_effect = fixed_effect_var
        ),
        pairwise_results
      )
    }
      
    results <- do.call(rbind, results_list)
    
    # Clean up
    rm(res, gene_data)
    gc()  # Force garbage collection
    
    write(Sys.time(), file = done_flag)
    return(results)
  }, error = function(e) {
    log_message <- paste(Sys.time(), "- General error in process_gene_gam for", cell_type_i, gene_name, ":", e$message, "\n")
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
    return(NULL)  # Return NULL on failure
  }, warning = function(w) {
    log_message <- sprintf("WARNING in process_gene_gam for %s, %s: %s", cell_type_i, gene_name, conditionMessage(w))
    cat(log_message, file = log_file_err, append = TRUE)
    flush.console()
  })
}