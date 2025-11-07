# corr_tbl_ts.R
suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
  library(lmtest)
  library(nlme)
  library(purrr)
})

# ---- Time-series aware correlation + regression summary per SPECIES × PROJECT ----
corr_tbl_ts <- function(df, x, y, lag_max = 5) {
  df %>%
    group_by(SPECIES, PROJECT) %>%
    group_modify(function(d, key) {
      xv <- dplyr::pull(d, !!enquo(x))
      yv <- dplyr::pull(d, !!enquo(y))
      ok <- is.finite(xv) & is.finite(yv)
      x <- xv[ok]; y <- yv[ok]
      n_ok <- length(x)
      
      out <- tibble(
        n_ok = n_ok,
        lm_slope = NA_real_, lm_se = NA_real_, lm_t = NA_real_, lm_p = NA_real_,
        lm_R2 = NA_real_,
        lm_significant = NA_character_,
        dw_stat = NA_real_, dw_p = NA_real_,
        acf1 = NA_real_, pacf1 = NA_real_,
        acf_sig_lags = NA_character_,
        autocorr_flag = NA_character_,
        gls_slope = NA_real_, gls_se = NA_real_, gls_t = NA_real_, gls_p = NA_real_,
        gls_rho = NA_real_, gls_AIC = NA_real_,
        used_gls = FALSE,
        status = NA_character_
      )
      
      if (n_ok < 2) {
        out$status <- "insufficient pairs for correlation"
        return(out)
      }
      
      # LM
      if (n_ok >= 2) {
        df_lm <- data.frame(y = y, x = x)
        fit_lm <- try(lm(y ~ x, data = df_lm), silent = TRUE)
        if (!inherits(fit_lm, "try-error")) {
          sm <- summary(fit_lm)
          co <- sm$coefficients
          out$lm_slope <- co["x", "Estimate"]
          out$lm_se    <- co["x", "Std. Error"]
          out$lm_t     <- co["x", "t value"]
          out$lm_p     <- co["x", "Pr(>|t|)"]
          out$lm_R2    <- sm$r.squared
          out$lm_significant <- if(out$lm_p <0.05) "significant" else "NS" 
          
          # DW test
          dw <- try(lmtest::dwtest(fit_lm), silent = TRUE)
          if (!inherits(dw, "try-error")) { out$dw_stat <- dw$statistic[1]; out$dw_p <- dw$p.value }
          
          # ACF/PACF
          res <- residuals(fit_lm)
          ac <- try(acf(res, plot = FALSE, na.action = na.pass), silent = TRUE)
          pac <- try(pacf(res, plot = FALSE, na.action = na.pass), silent = TRUE)
          if (!inherits(ac, "try-error") && length(ac$acf) > 1) {
            out$acf1 <- as.numeric(ac$acf[2])
            band <- 1.96 / sqrt(length(res[is.finite(res)]))
            lag_vals <- as.numeric(ac$acf[2:(min(lag_max + 1, length(ac$acf)))])
            sig_lags <- which(abs(lag_vals) > band)
            out$acf_sig_lags <- if (length(sig_lags)) paste(sig_lags, collapse = ",") else "none"
          }
          if (!inherits(pac, "try-error") && length(pac$acf) >= 1) {
            out$pacf1 <- as.numeric(pac$acf[1])
          }
          
          # Autocorr flag
          ac_flag <- FALSE
          if (!is.na(out$dw_p) && out$dw_p < 0.05) ac_flag <- TRUE
          if (!is.null(out$acf_sig_lags) && out$acf_sig_lags != "none") ac_flag <- TRUE
          out$autocorr_flag <- if (ac_flag) "autocorr_detected" else "no_autocorr_detected"
          
          # GLS if autocorr detected
          if (ac_flag && n_ok >= 3) {
            fit_gls <- try(
              nlme::gls(y ~ x, data = df_lm,
                        correlation = nlme::corAR1(form = ~ x),
                        method = "REML"),
              silent = TRUE
            )
            if (!inherits(fit_gls, "try-error")) {
              co_g <- summary(fit_gls)$tTable
              out$gls_slope <- co_g["x","Value"]
              out$gls_se    <- co_g["x","Std.Error"]
              out$gls_t     <- co_g["x","t-value"]
              out$gls_p     <- co_g["x","p-value"]
              rho <- try(as.numeric(coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)), silent = TRUE)
              out$gls_rho <- if (inherits(rho, "try-error")) NA_real_ else rho
              out$gls_AIC <- AIC(fit_gls)
              out$used_gls <- TRUE
            }
          }
        }
      }
      
      out$status <- dplyr::case_when(
        n_ok < 2 ~ "insufficient pairs for correlation",
        isTRUE(out$used_gls) ~ "ok (GLS AR1 used due to autocorrelation)",
        TRUE ~ "ok"
      )
      out
    }) %>% ungroup()
}
