# Module: IMC coordinate parsing.
# Purpose: Convert the IMC `Position` representation to numeric micrometre coordinates.
# Callers: `figures/Figure3/01_cellular_neighborhoods.R` and, through
# `figures/Figure4/01_methods_aligned_orchestration.R`, the shared Figure 4 modules.
# Inputs: Character vectors of `Position` values formatted as `(x, y)`.
# Outputs: Data frames with numeric `x` and `y` columns in input order.
# Ordered use: Source this module before `cellular_neighborhoods.R` or `ptme.R`, then
# call `position_to_xy()` before spatial calculations.

#' Purpose: Parse ordered IMC position strings into numeric micrometre coordinates.
#'
#' @param position Character vector of values formatted as `(x, y)`.
#' @return Data frame with numeric `x` and `y` columns in the order of `position`.
position_to_xy <- function(position) {
    # Accept only the canonical character representation before extracting coordinates.
    if (!is.character(position)) {
        stop("position must be a character vector")
    }

    # Capture exactly two numeric fields while permitting surrounding whitespace.
    number <- "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    pattern <- paste0(
        "^\\s*\\(\\s*(", number, ")\\s*,\\s*(", number,
        ")\\s*\\)\\s*$"
    )
    matches <- regexec(pattern, position, perl = TRUE)
    parts <- regmatches(position, matches)
    valid <- !is.na(position) & lengths(parts) == 3L

    # Preserve the one-input-to-one-coordinate invariant by rejecting any malformed position.
    if (any(!valid)) {
        stop(
            "Malformed Position value(s) at indices: ",
            paste(which(!valid), collapse = ", ")
        )
    }

    # Return numeric micrometre axes in the same row order as the input cells.
    data.frame(
        x = vapply(parts, function(value) as.numeric(value[[2L]]), numeric(1)),
        y = vapply(parts, function(value) as.numeric(value[[3L]]), numeric(1)),
        check.names = FALSE
    )
}
