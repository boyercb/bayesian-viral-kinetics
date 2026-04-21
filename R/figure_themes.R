# ──────────────────────────────────────────────────────────────────────────────
# figure_themes.R — Journal-adaptive figure themes and color palettes
#
# Usage:
#   set_journal("pnas")       # switches theme globally
#   p + theme_journal()       # applies current journal theme
#   journal_colors()          # returns named color list
#   save_journal_figure(p, "fig1", width = 7, height = 5)
# ──────────────────────────────────────────────────────────────────────────────

# Module-level state for active journal style
.journal_env <- new.env(parent = emptyenv())
.journal_env$style <- "pnas"

#' Set the active journal style
#'
#' Switches the global theme for all subsequent calls to theme_journal(),
#' journal_colors(), and save_journal_figure().
#'
#' @param style One of "pnas", "plos", "annals"
#' @return Invisible previous style
set_journal <- function(style = c("pnas", "plos", "annals")) {
  style <- match.arg(style)
  old <- .journal_env$style
  .journal_env$style <- style
  invisible(old)
}

#' Get the active journal style
#' @return Character string
get_journal <- function() {
  .journal_env$style
}

#' Journal-specific ggplot2 theme
#'
#' Returns a ggplot2 theme appropriate for the active (or specified) journal.
#' - "pnas": Clean, sans-serif (Helvetica), minimal gridlines, 8pt base
#' - "plos": Modern, sans-serif (Arial), light grid, 10pt base
#' - "annals": Classic academic, serif-adjacent, full axis frames, 10pt base
#'
#' @param style  Override active journal (NULL = use set_journal())
#' @param base_size Override base font size (NULL = journal default)
#' @return ggplot2 theme object
theme_journal <- function(style = NULL, base_size = NULL) {
  if (is.null(style)) style <- get_journal()

  if (style == "pnas") {
    bs <- base_size %||% 8
    ggplot2::theme_minimal(base_size = bs, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          color = "grey92", linewidth = 0.3
        ),
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold", size = bs, hjust = 0
        ),
        strip.background = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = bs),
        axis.text = ggplot2::element_text(size = bs - 1, color = "grey30"),
        legend.title = ggplot2::element_text(size = bs - 1),
        legend.text = ggplot2::element_text(size = bs - 1),
        legend.key.size = grid::unit(0.35, "cm"),
        plot.title = ggplot2::element_text(
          size = bs + 1, face = "bold", hjust = 0
        ),
        plot.subtitle = ggplot2::element_text(
          size = bs, color = "grey40", hjust = 0
        ),
        plot.margin = ggplot2::margin(4, 6, 4, 4, "pt"),
        panel.spacing = grid::unit(0.6, "lines")
      )

  } else if (style == "plos") {
    bs <- base_size %||% 10
    ggplot2::theme_minimal(base_size = bs, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          color = "grey90", linewidth = 0.3
        ),
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold", size = bs, hjust = 0
        ),
        strip.background = ggplot2::element_rect(
          fill = "grey96", color = NA
        ),
        axis.title = ggplot2::element_text(size = bs),
        axis.text = ggplot2::element_text(size = bs - 1, color = "grey25"),
        legend.title = ggplot2::element_text(size = bs - 1, face = "bold"),
        legend.text = ggplot2::element_text(size = bs - 1),
        legend.key.size = grid::unit(0.4, "cm"),
        plot.title = ggplot2::element_text(
          size = bs + 2, face = "bold", hjust = 0
        ),
        plot.subtitle = ggplot2::element_text(
          size = bs, color = "grey40", hjust = 0
        ),
        plot.margin = ggplot2::margin(6, 8, 6, 6, "pt"),
        panel.spacing = grid::unit(0.8, "lines")
      )

  } else if (style == "annals") {
    bs <- base_size %||% 10
    ggplot2::theme_bw(base_size = bs, base_family = "Times") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(
          color = "grey90", linewidth = 0.25
        ),
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "italic", size = bs
        ),
        strip.background = ggplot2::element_rect(
          fill = "grey95", color = "grey70"
        ),
        axis.title = ggplot2::element_text(size = bs),
        axis.text = ggplot2::element_text(size = bs - 1),
        legend.title = ggplot2::element_text(size = bs - 1),
        legend.text = ggplot2::element_text(size = bs - 1),
        plot.title = ggplot2::element_text(
          size = bs + 1, face = "bold", hjust = 0
        ),
        plot.margin = ggplot2::margin(6, 8, 6, 6, "pt"),
        panel.spacing = grid::unit(0.8, "lines")
      )
  }
}

#' Journal-specific color palette
#'
#' Returns a named list of colors with role-based names.
#'
#' @param style Override active journal (NULL = use set_journal())
#' @return Named list: rna, pfu, lfd, sym, sig, nonsig, accent, muted
journal_colors <- function(style = NULL) {
  if (is.null(style)) style <- get_journal()
  if (is.null(style)) style <- "pnas"   # sensible default

  if (style == "pnas") {
    list(
      rna     = "#3C78D8",    # steel blue
      pfu     = "#CC3333",    # muted red
      lfd     = "#6AA84F",    # sage green
      sym     = "#F1C232",    # gold
      sig     = "#3C78D8",    # matches rna
      nonsig  = "#B0B0B0",    # grey
      accent  = "#E69138",    # warm orange
      muted   = "#999999",
      ci_fill = "#3C78D820"   # transparent blue
    )
  } else if (style == "plos") {
    list(
      rna     = "#0072B2",    # accessible blue
      pfu     = "#D55E00",    # vermillion
      lfd     = "#009E73",    # teal
      sym     = "#E69F00",    # amber
      sig     = "#0072B2",
      nonsig  = "#AAAAAA",
      accent  = "#CC79A7",    # pink
      muted   = "#888888",
      ci_fill = "#0072B220"
    )
  } else if (style == "annals") {
    list(
      rna     = "#1F77B4",    # tableau blue
      pfu     = "#D62728",    # tableau red
      lfd     = "#2CA02C",    # tableau green
      sym     = "#FF7F0E",    # tableau orange
      sig     = "#1F77B4",
      nonsig  = "#7F7F7F",
      accent  = "#9467BD",    # purple
      muted   = "#AAAAAA",
      ci_fill = "#1F77B420"
    )
  }
}

#' Journal-specific figure dimensions (in inches)
#'
#' Returns default width/height for single-column, 1.5-column, and
#' full-width figures appropriate for the active journal.
#'
#' @param layout One of "single", "medium", "full"
#' @param style  Override active journal (NULL = use set_journal())
#' @return Named list: width, height
journal_dims <- function(layout = c("single", "medium", "full"),
                         style = NULL) {
  if (is.null(style)) style <- get_journal()
  layout <- match.arg(layout)

  dims <- list(
    pnas = list(
      single = list(width = 3.42, height = 3.0),
      medium = list(width = 5.50, height = 4.5),
      full   = list(width = 7.08, height = 5.5)
    ),
    plos = list(
      single = list(width = 5.20, height = 4.0),
      medium = list(width = 5.20, height = 5.0),
      full   = list(width = 7.50, height = 6.0)
    ),
    annals = list(
      single = list(width = 3.25, height = 3.0),
      medium = list(width = 5.00, height = 4.5),
      full   = list(width = 6.50, height = 5.5)
    )
  )

  dims[[style]][[layout]]
}

#' Save a figure with journal-appropriate defaults
#'
#' Saves the plot as PDF (for LaTeX) with dimensions from journal_dims().
#' Also saves a high-DPI PNG version for review.
#'
#' @param plot     ggplot object
#' @param name     Figure name (without extension)
#' @param layout   One of "single", "medium", "full"
#' @param width    Override width (inches)
#' @param height   Override height (inches)
#' @param out_dir  Output directory
#' @param style    Override active journal (NULL = use set_journal())
#' @return Character vector of saved file paths
save_journal_figure <- function(plot, name,
                                layout = "full",
                                width = NULL, height = NULL,
                                out_dir = "output/figures",
                                style = NULL) {
  if (is.null(style)) style <- get_journal()
  dims <- journal_dims(layout, style)
  w <- width  %||% dims$width
  h <- height %||% dims$height

  dir.create(file.path(out_dir, style), recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(out_dir, style, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, style, paste0(name, ".png"))

  # Use regular pdf device (cairo_pdf may silently fail on some macOS setups)
  ggplot2::ggsave(pdf_path, plot, width = w, height = h, device = "pdf")
  ggplot2::ggsave(png_path, plot, width = w, height = h, dpi = 300)

  c(pdf_path, png_path)
}

#' Generate a figure in all journal styles
#'
#' Calls a figure-generating function with each journal style and saves
#' the output. The function should accept a `style` argument.
#'
#' @param fig_fn   Function that takes (..., style) and returns a ggplot
#' @param name     Figure name (without extension)
#' @param layout   One of "single", "medium", "full"
#' @param ...      Additional arguments passed to fig_fn
#' @return List of saved file paths by journal
save_all_styles <- function(fig_fn, name, layout = "full", ...) {
  styles <- c("pnas", "plos", "annals")
  paths <- list()
  for (s in styles) {
    p <- fig_fn(..., style = s)
    paths[[s]] <- save_journal_figure(p, name, layout = layout, style = s)
  }
  paths
}
