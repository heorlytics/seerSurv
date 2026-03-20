#' Plot Blended Survival Curves for LR and DR Populations
#'
#' @description
#' Produces a \code{\link[ggplot2]{ggplot}} comparing up to four survival
#' curves: disease-specific blended LR, disease-specific blended DR, net LR,
#' and net DR.  Designed to give a quick visual QC check after running
#' \code{\link{blend_survival}} and \code{\link{make_background_surv}}.
#'
#' @param S_blend_L A two-column tibble (\code{time}, \code{surv}) — blended
#'   disease-specific LR survival from \code{\link{blend_survival}}.
#' @param S_blend_D A two-column tibble (\code{time}, \code{surv}) — blended
#'   disease-specific DR survival.
#' @param U_L A two-column tibble (\code{time}, \code{surv}) — net LR survival
#'   (disease-specific × background).  Pass \code{NULL} to omit.
#' @param U_D A two-column tibble (\code{time}, \code{surv}) — net DR survival.
#'   Pass \code{NULL} to omit.
#' @param tumour Character scalar used in the plot title (default
#'   \code{"Tumour"}).
#' @param time_unit Character scalar for the x-axis label (default
#'   \code{"Years"}).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' data(lifetable_seer)
#'
#' t <- 0:5
#' s_lr <- c(1, 0.9959, 0.9886, 0.9831, 0.9783, 0.9754)
#' s_dr <- c(1, 0.5491, 0.4396, 0.3945, 0.3630, 0.3410)
#' grid <- seq(0, 39, by = 1)
#'
#' ipd_L <- prep_ipd(t, s_lr)
#' ipd_D <- prep_ipd(t, s_dr)
#' mods_L <- fit_models(ipd_L[-1, ])
#' mods_D <- fit_models(ipd_D[-1, ])
#' wts_L  <- compute_weights(extract_ic(mods_L, "AIC"), 3, "AIC")
#' wts_D  <- compute_weights(extract_ic(mods_D, "AIC"), 3, "AIC")
#' SL     <- blend_survival(mods_L, wts_L, grid)
#' SD     <- blend_survival(mods_D, wts_D, grid)
#'
#' lt     <- lifetable_seer
#' lt$btrate_LR <- lt$Males * 0.58 + lt$Females * 0.42
#' lt$btrate_DR <- lt$Males * 0.69 + lt$Females * 0.31
#' BL <- make_background_surv(lt, 61, "btrate_LR", grid)
#' BD <- make_background_surv(lt, 62, "btrate_DR", grid)
#'
#' UL <- dplyr::mutate(
#'   dplyr::left_join(SL, BL, by = "time"),
#'   surv = surv.x * surv.y
#' )[, c("time", "surv")]
#' UD <- dplyr::mutate(
#'   dplyr::left_join(SD, BD, by = "time"),
#'   surv = surv.x * surv.y
#' )[, c("time", "surv")]
#'
#' p <- plot_survival_curves(SL, SD, UL, UD, tumour = "Melanoma")
#' p
#' }
#'
#' @export
plot_survival_curves <- function(
    S_blend_L,
    S_blend_D,
    U_L       = NULL,
    U_D       = NULL,
    tumour    = "Tumour",
    time_unit = "Years"
) {

  curves <- list(
    "LR (disease-specific)" = S_blend_L,
    "DR (disease-specific)" = S_blend_D
  )
  if (!is.null(U_L)) curves[["LR (net)"]] <- U_L
  if (!is.null(U_D)) curves[["DR (net)"]] <- U_D

  # Combine with a label column
  plot_df <- dplyr::bind_rows(
    lapply(names(curves), function(nm) {
      dplyr::mutate(curves[[nm]], curve = nm)
    })
  )

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = .data$time, y = .data$surv,
                               colour = .data$curve)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      title  = paste0("Survival Curves: ", tumour),
      x      = paste0("Time (", time_unit, ")"),
      y      = "Survival Probability",
      colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
