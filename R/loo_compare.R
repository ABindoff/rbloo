# =====================================================================
# Model comparison on rb_loo objects.
#
# loo::loo_compare() is a generic, so loo_compare(rb1, rb2) dispatches
# here when the first argument is an rb_loo. The method compares the
# embedded RB loo objects (loo_rb) -- the marginal-over-RE predictive --
# after enforcing the comparisons that are actually valid:
#   * every input must be an rb_loo. An rb_loo next to a plain loo would
#     compare a marginal elpd against a conditional one -- different
#     estimands (Merkle, Furr & Rabe-Hesketh 2019) -- so it is refused,
#     not coerced.
#   * PSIS fallbacks carry no RB elpd. An all-fallback comparison
#     demotes, with a loud warning, to the conditional PSIS-LOO
#     comparison (loo_full): same estimand on both sides, just not RB.
#     A mix of real RB results and fallbacks is refused outright.
#   * residual flagged folds (RB k > base_cut) in any model draw a
#     warning naming the per-model counts, because the elpd difference
#     may rest on exactly those folds; rb_loo(fit, reloo = TRUE)
#     replaces them with exact refits first.
# The return value is loo's own "compare.loo" matrix, so printing and
# downstream use match loo_compare() everywhere else.
# =====================================================================

#' @importFrom loo loo_compare
#' @export
loo::loo_compare

#' Compare models on RB-LOO
#'
#' Compares two or more \code{\link{rb_loo}} results on their
#' Rao-Blackwellised elpd (the marginal-over-RE predictive), delegating the
#' arithmetic to \code{loo::loo_compare()} on the embedded \code{loo_rb}
#' objects. Comparisons that would mix estimands are refused: every input
#' must be an rb_loo, and PSIS fallbacks (where RB-LOO was not applied)
#' cannot be compared against real RB results. If every model is a
#' fallback, the comparison demotes to the conditional PSIS-LOO
#' (\code{loo_full}) with a warning. Residual flagged folds in any model
#' also draw a warning, since the elpd difference may rest on exactly
#' those folds; run \code{rb_loo(fit, reloo = TRUE)} first to replace them
#' with exact refits.
#'
#' @param x an \code{rb_loo} object.
#' @param ... further \code{rb_loo} objects (at least one).
#' @return a \code{compare.loo} matrix, as returned by
#'   \code{loo::loo_compare()}: models ordered best-first, with
#'   \code{elpd_diff} and \code{se_diff} in the leading columns.
#' @export
loo_compare.rb_loo <- function(x, ...) {
  objs <- c(list(x), list(...))
  # row labels from the calling expressions, as loo_compare() itself does
  calls <- as.list(match.call(expand.dots = TRUE))[-1]
  nms   <- vapply(calls, function(cl) paste(deparse(cl), collapse = ""),
                  character(1))[seq_along(objs)]

  ok <- vapply(objs, inherits, logical(1), what = "rb_loo")
  if (!all(ok))
    stop("loo_compare: all objects must be rb_loo results (got ",
         paste(unique(vapply(objs[!ok], function(o) class(o)[1L],
                             character(1))), collapse = ", "),
         " alongside rb_loo). RB-LOO estimates the MARGINAL ",
         "(integrated-over-RE) predictive; comparing it against a ",
         "conditional PSIS-LOO mixes two different estimands. Compare like ",
         "with like: rb_loo() every model, or loo::loo_compare() the plain ",
         "loo objects.", call. = FALSE)
  if (length(objs) < 2L)
    stop("loo_compare: needs at least two rb_loo objects.", call. = FALSE)

  fb <- vapply(objs, function(o) !is.null(o$meta$fallback), logical(1))
  if (all(fb)) {
    warning("loo_compare: every model is a PSIS-LOO fallback (RB-LOO was ",
            "not applied to any of them); comparing the plain conditional ",
            "PSIS-LOO (elpd_full) instead. This is a like-for-like ",
            "comparison, but it is not RB-LOO and inherits PSIS-LOO's ",
            "reliability caveats on influential observations.", call. = FALSE)
    loos <- lapply(objs, `[[`, "loo_full")
  } else if (any(fb)) {
    stop("loo_compare: ", paste(nms[fb], collapse = ", "), " is a PSIS-LOO ",
         "fallback (no RB elpd) while the rest carry real RB-LOO results -- ",
         "different estimands, so the comparison would be invalid. Either ",
         "compare the conditional PSIS-LOO for ALL models ",
         "(loo::loo_compare() on the $loo_full objects) or bring every ",
         "model into RB scope.", call. = FALSE)
  } else {
    loos  <- lapply(objs, `[[`, "loo_rb")
    nflag <- vapply(objs, function(o) sum(o$refit_flag, na.rm = TRUE),
                    numeric(1))
    if (any(nflag > 0))
      warning("loo_compare: residual flagged folds remain (",
              paste0(nms[nflag > 0], ": ", nflag[nflag > 0], collapse = ", "),
              "); the elpd difference may rest on exactly these folds. ",
              "rb_loo(fit, reloo = TRUE) replaces them with exact refits.",
              call. = FALSE)
  }
  cmp <- loo::loo_compare(loos)
  # loo labels list input by position ("1".."k" / "model1".."modelk") and
  # ignores list names -- map back to the calling expressions so the table
  # reads like loo_compare(fit1, fit2) would. Newer loo returns a data frame
  # whose `model` column drives print(); older loo used the matrix rownames.
  relabel <- function(lbl) {
    idx <- suppressWarnings(as.integer(sub("^model", "", lbl)))
    if (!anyNA(idx) && setequal(idx, seq_along(loos))) nms[idx] else lbl
  }
  if (is.data.frame(cmp) && "model" %in% names(cmp)) {
    cmp$model <- relabel(cmp$model)
    rownames(cmp) <- cmp$model
  } else {
    rownames(cmp) <- relabel(rownames(cmp))
  }
  cmp
}
