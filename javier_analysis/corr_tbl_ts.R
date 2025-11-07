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
      x  <- xv[ok]; y <- yv[ok]
      n_ok <- length(x)
      
      out <- tibble(
        n_ok = n_ok,
        lm_slope = NA_real_, lm_se = NA_real_, lm_t = NA_real_, lm_p = NA_real_,
        lm_R2 = NA_real_, lm_significant = NA_character_,
        dw_stat = NA_real_, dw_p = NA_real_,
        acf1 = NA_real_, pacf1 = NA_real_,
        acf_sig_lags = NA_character_,
        autocorr_flag = NA_character_,
        gls_slope = NA_real_, gls_se = NA_real_, gls_t = NA_real_, gls_p = NA_real_,
        gls_rho = NA_real_, gls_AIC = NA_real_,
        used_gls = FALSE,
        status = NA_character_,
        note = NA_character_
      )
      
      # Early exits for insufficient data or zero variance (common NA p-value causes)
      if (n_ok < 3) {
        out$status <- "insufficient pairs for correlation"
        out$lm_significant <- "NS"
        out$note <- "n < 3"
        return(out)
      }
      vx <- stats::var(x)
      vy <- stats::var(y)
      if (!is.finite(vx) || vx == 0 || !is.finite(vy) || vy == 0) {
        out$status <- "insufficient pairs for correlation"
        out$lm_significant <- "NS"
        out$note <- "zero variance in x or y"
        return(out)
      }
      
      # LM (safe)
      df_lm <- data.frame(y = y, x = x)
      fit_lm <- try(lm(y ~ x, data = df_lm), silent = TRUE)
      if (!inherits(fit_lm, "try-error")) {
        sm <- summary(fit_lm)
        co <- sm$coefficients
        # Extract safely (could be NA in singular fits)
        out$lm_slope <- co["x", "Estimate"] %||% NA_real_
        out$lm_se    <- co["x", "Std. Error"] %||% NA_real_
        out$lm_t     <- co["x", "t value"] %||% NA_real_
        out$lm_p     <- co["x", "Pr(>|t|)"] %||% NA_real_
        out$lm_R2    <- sm$r.squared %||% NA_real_
        
        # NEVER branch on NA: use is.finite()
        out$lm_significant <- if (is.finite(out$lm_p) && out$lm_p < 0.05) "significant" else "NS"
        
        # DW test (safe)
        dw <- try(lmtest::dwtest(fit_lm), silent = TRUE)
        if (!inherits(dw, "try-error")) {
          out$dw_stat <- as.numeric(dw$statistic[1])
          out$dw_p    <- as.numeric(dw$p.value)
        }
        
        # ACF / PACF of residuals (safe)
        res <- residuals(fit_lm)
        ac  <- try(acf(res, plot = FALSE, na.action = na.pass),  silent = TRUE)
        pac <- try(pacf(res, plot = FALSE, na.action = na.pass), silent = TRUE)
        if (!inherits(ac, "try-error") && length(ac$acf) > 1) {
          out$acf1 <- as.numeric(ac$acf[2])
          band <- 1.96 / sqrt(sum(is.finite(res)))
          lag_vals <- as.numeric(ac$acf[2:(min(lag_max + 1, length(ac$acf)))])
          sig_lags <- which(abs(lag_vals) > band)
          out$acf_sig_lags <- if (length(sig_lags)) paste(sig_lags, collapse = ",") else "none"
        } else {
          out$acf_sig_lags <- "none"
        }
        if (!inherits(pac, "try-error") && length(pac$acf) >= 1) {
          out$pacf1 <- as.numeric(pac$acf[1])
        }
        
        # Autocorr flag (safe checks)
        ac_flag <- (is.finite(out$dw_p) && out$dw_p < 0.05) ||
          (!is.na(out$acf_sig_lags) && out$acf_sig_lags != "none")
        out$autocorr_flag <- if (isTRUE(ac_flag)) "autocorr_detected" else "no_autocorr_detected"
        
        # GLS if autocorr detected & enough data
        if (isTRUE(ac_flag) && n_ok >= 4) {
          fit_gls <- try(
            nlme::gls(y ~ x, data = df_lm,
                      correlation = nlme::corAR1(form = ~ x),
                      method = "REML"),
            silent = TRUE
          )
          if (!inherits(fit_gls, "try-error")) {
            co_g <- summary(fit_gls)$tTable
            out$gls_slope <- co_g["x","Value"]      %||% NA_real_
            out$gls_se    <- co_g["x","Std.Error"]  %||% NA_real_
            out$gls_t     <- co_g["x","t-value"]    %||% NA_real_
            out$gls_p     <- co_g["x","p-value"]    %||% NA_real_
            rho <- try(as.numeric(coef(fit_gls$modelStruct$corStruct, unconstrained = FALSE)), silent = TRUE)
            out$gls_rho <- if (inherits(rho, "try-error")) NA_real_ else rho
            out$gls_AIC <- AIC(fit_gls)
            out$used_gls <- TRUE
          }
        }
      } else {
        out$note <- "lm() failed"
      }
      
      out$status <- dplyr::case_when(
        n_ok < 3 ~ "insufficient pairs for correlation",
        isTRUE(out$used_gls) ~ "ok (GLS AR1 used due to autocorrelation)",
        TRUE ~ "ok"
      )
      if (!is.finite(out$lm_p)) {
        out$lm_significant <- "NS"
        out$note <- coalesce(out$note, "p-value NA (singular/perfect fit?)")
      }
      out
    }) %>% ungroup()
}

# helper for defaulting NULL -> default
`%||%` <- function(a, b) if (is.null(a)) b else a
