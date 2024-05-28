#' Define what category to look at
#'
#' @param .array Array to analyze.
#' @param what Either `"voxels"` (default) or `"components"` to look at
#'   individual voxels or connected components.
#' @param of_category Vector of categories to consider. If negative, these are
#'   omitted.
#' @param in_relation_to Vector of categories to consider for cross-category
#'   morphology. If negative, these are omitted.
#' @param kernel_width Width of the kernel for determining connectivity. Can be
#'   a vector (row, column, layer). Not used when looking at voxels.
#' @param kernel_shape Shape of the kernel for determining connectivity. Can be
#'   `"diamond"` (default), `"disc"` or `"box"`. Not used when looking at
#'   voxels.
#' @return Data frame (tibble) with `x`, `y`, `z` and `value` or `component`
#'   columns.
#' @export
look_at <- function(
  .array,
  what = "voxels",
  of_category,
  in_relation_to = NULL,
  kernel_width = 3,
  kernel_shape = "disc"
) {
  .array <- prep_array(.array)
  array_dim <- dim(.array)
  of_category <- prep_category(of_category, .array)
  in_relation_to <- prep_category(in_relation_to, .array)
  array_prop <- mean(.array %in% of_category)
  if (!is.null(in_relation_to)) array_prop <- c(array_prop, mean(.array %in% in_relation_to))
  if (what == "components") {
    if (length(kernel_width) < ndim(.array)) kernel_width <- rep(kernel_width, ndim(.array))
    .array <- mmand::components(
      array(.array %in% of_category, dim = dim(.array)),
      mmand::shapeKernel(kernel_width, type = kernel_shape)
    )
  }
  .df <- .array |>
    reshape2::melt() |>
    tibble::as_tibble() |>
    setNames(c("y", "x", "z", c(voxels = "value", components = "component")[[what]])[c(1:ndim(.array), 4)])
  if (what == "components") return(
    .df[!is.na(.df$component),][, c(2, 1, 3, 4)] |>
      set_attribute("dimensions", array_dim) |>
      set_attribute("proportion", array_prop)
  )
  if (identical(of_category, in_relation_to)) in_relation_to <- NULL
  category <- .df[.df$value %in% of_category,]
  category$value <- 1L
  if (is.null(in_relation_to)) return(
    category[, c(2, 1, 3, 4)] |>
      set_attribute("dimensions", array_dim) |>
      set_attribute("proportion", array_prop)
  )
  in_relation_to <- .df[.df$value %in% in_relation_to,]
  in_relation_to$value <- 2L
  rbind(category, in_relation_to)[, c(2, 1, 3, 4)] |>
    set_attribute("dimensions", array_dim) |>
    set_attribute("proportion", array_prop)
}

#' Define what direction to look in
#'
#' @param .df Data frame with `x`, `y`, `z` and `value` or `component` columns,
#'   typically output from `look_at()`.
#' @param direction Direction to consider. Either 1D (`"x"`, `"y"`, `"z"`), 2D
#'   (`"xy"`, `"xz"`, `"yz"`) or 3D (`"xyz"`).
#' @param every Value indicating how many voxels to consider, e.g. one of two
#'   (`1/2`), one every three (`1/3`), etc. Defaults to every voxel (`1/1`).
#'   This can be a vector (row, column, layer).
#' @return Data frame (tibble) with `x`, `y`, and/or `z`, and `value` or
#'   `component` columns.
#' @export
look_in <- function(.df, direction = "xyz", every = 1) {
  attributes_df <- get_attributes(.df)
  every <- if (length(every) == 1) rep(every, 3) else every
  every <- 1/every
  if (any(every != 1)) {
    if ("x" %in% names(.df)) .df <- .df[.df$x %% every[2] == 0,]
    if ("y" %in% names(.df)) .df <- .df[.df$y %% every[1] == 0,]
    if ("z" %in% names(.df)) .df <- .df[.df$z %% every[3] == 0,]
  }
  if (direction == "xyz") return(
    make3d(.df) |>
      set_attributes(attributes_df) |>
      set_attribute("direction", direction) |>
      set_attribute("every", every)
  )
  if (nchar(direction) == 1) return(
    make1d(.df, direction) |>
      set_attributes(attributes_df) |>
      set_attribute("direction", direction) |>
      set_attribute("every", every)
  )
  if (nchar(direction) == 2) return(
    make2d(.df, direction) |>
      set_attributes(attributes_df) |>
      set_attribute("direction", direction) |>
      set_attribute("every", every)
  )
  std::err("e Wrong direction")
}

#' Define what neighbours to look for
#'
#' @param .df Data frame with `x`, `y`, and/or `z`, and `value` or `component`
#'   columns, typically output from `look_at()` or `look_in()`.
#' @param neighbours Number of nearest neighbours to look for. Can be `Inf` if
#'   `within` is finite, in which case all neighbours within that distance are
#'   looked for.
#' @param within Maximum search distance.
#' @param from_border Look for neighbours of border voxels. Defaults to TRUE. If
#'   FALSE, minus sampling is performed.
#' @return List with the `distance` matrix and the maximum search distance
#'   `within`. When looking at components, additionally a `connected` matrix.
#' @export
look_for <- function(.df, neighbours = 1, within = Inf, from_border = FALSE) {
  attributes_df <- get_attributes(.df)
  if (is.infinite(neighbours)) {
    if (is.infinite(within)) std::err("e Either neighbours or within should be finite.")
    neighbours <- ceiling(1 + 4 / 3 * pi * within^3)
  }
  if ("value" %in% names(.df)) {
    afrom <- .df[.df$value == 1,]
    ato <- .df[.df$value == 2,]
    self <- FALSE
    if (nrow(ato) == 0) {
      ato <- afrom
      self <- TRUE
    }
    dims <- sort(names(.df)[!names(.df) %in% c("value", "border")])
    knn <- nabor::knn(
      ato[, dims],
      if (from_border) afrom[, dims] else afrom[afrom$border > within, dims],
      k = neighbours + if (self) 1 else 0,
      radius = within
    )
    knn$nn.idx[knn$nn.idx == 0] <- NA
    return(
      list(
        distance = if (self) knn$nn.dists[, -1, drop = FALSE] else knn$nn.dists
      ) |>
        set_attributes(attributes_df) |>
        set_attribute("within", within)
    )
  }
  if ("component" %in% names(.df)) {
    dims <- sort(names(.df)[!names(.df) %in% c("component", "border")])
    knn <- nabor::knn(
      .df[, dims],
      if (from_border) .df[, dims] else .df[.df$border > within, dims],
      k = neighbours + 1,
      radius = within
    )
    knn$nn.idx[knn$nn.idx == 0] <- NA
    component <- knn$nn.idx
    component[] <- .df$component[knn$nn.idx]
    return(
      list(
        connected = (component == component[, 1])[, -1, drop = FALSE],
        distance = knn$nn.dists[, -1, drop = FALSE]
      ) |>
        set_attributes(attributes_df) |>
        set_attribute("within", within)
    )
  }
  std::err("e No value or component column found.")
}

#' Describe morphology using summary functions
#'
#' @param .list List with the `distance` matrix and `within`, and potentially
#'   a `connected` matrix, typically output from `look_for()`. Alternatively,
#'   this can be a data frame resulting from a `describe()` including a `name`
#'   argument, to add additional described summaries.
#' @param what Summary function.
#' @param connected Logical. Count connected or disconnected neighbours only.
#'   Defaults to `NULL`, for which all neighbours are considered.
#' @param cumulative Logical. Provide cumulative results (distances <=
#'   threshold) or binned results (distances at center of provided breaks;
#'   default).
#' @param name Name to use for the current described dataset in an extra name
#'   column. If provided, the input list is included as attribute to the output,
#'   and `describe()` can be repeated multiple times in a pipechain.
#' @param at Distances at which to describe the results. Defaults to the integer
#'   sequence from `0L` to the maximum search distance, if available, and
#'   otherwise the maximum distance found.
#' @return Data frame with `distance` and `neighbours` columns. If `name` is
#'   provided, an additional `name` column is included for identifying different
#'   described summaries. In case a `name` column is there, the object has an
#'   extra `morphology.list` attribute, to enable piping into additional
#'   `describe()` calls.
#' @export
describe <- function(
  .list,
  what = mean,
  of = "neighbours",
  at = NULL,
  connected = NULL,
  cumulative = FALSE,
  name = NULL
) {
  # TODO allow returning point cloud, without summary, through function that
  # returns vector rather than summary? then identity() as what? or other funs
  # that compress the data into less samples? eg all percentiles? for boxplot,
  # violin, etc?
  # TODO document `of` argument! to flip the summaries!
  .df <- NULL
  list_attributes <- NULL
  if (.list |> has_attribute("list")) {
    .df <- .list
    .list <- .list |> get_attribute("list")
    .df <- .df |> set_attribute("list") <- NULL
  }
  list_attributes <- get_attributes(.list)
  nn <- .list
  if (!is.null(connected)) {
    if (connected) nn$distance[!nn$connected] <- NA
    if (!connected) nn$distance[nn$connected] <- NA
  }
  if (of == "neighbours") {
    if (is.null(at)) {
      if (is.finite(.list |> get_attribute("within"))) {
        at <- 0:(.list |> get_attribute("within"))
      } else {
        at <- 0:ceiling(max(.list$distance[is.finite(.list$distance)]))
      }
    }
    results <- at * NA
    for (i in seq_along(results)) {
      results[i] <- what(matrixStats::rowSums2(nn$distance <= at[i], na.rm = TRUE))
    }
    if (cumulative) {
      df <- tibble::tibble(distance = at, neighbours = results, name = name)
    } else {
      df <- tibble::tibble(
        distance = (at[1:(length(at) - 1)] + at[2:length(at)]) / 2,
        neighbours = diff(results),
        name = name
      )
    }
  }
  if (of == "distance") {
    if (is.null(at)) {
      at <- 1:ncol(.list$distance)
    }
    results <- at * NA
    if (cumulative) {
      for (i in seq_along(results)) {
        results[i] <- what(nn$distance[, 1:i])
      }
    } else {
      for (i in seq_along(results)) {
        results[i] <- what(nn$distance[, i])
      }
    }
    df <- tibble::tibble(neighbour = at, distance = results, name = name)
  }
  if (!is.null(.df)) df <- rbind(.df, df)
  if (!is.null(name)) df <- df |> set_attribute("list") <- .list
  if (!is.null(list_attributes)) df <- df |> set_attributes(list_attributes)
  df |>
    set_attribute("at", at) |>
    set_attribute("cumulative", cumulative)
}

#' Scale the morphological description
#'
#' @param df Data frame
#' @param ... Things to scale the description by. Can be `"neighbourhood"`, which
#'   automatically chooses between alternative options `"volume"`, `"area"` or `"length"`,
#'   or `"proportion"`, `"inverse proportion"`, another description data frame
#'   with the same distances, ...
#' @return Data frame with scaled `neighbours` column.
#' @export
scale_by <- function(df, ...) {
  dots <- list(...)
  at <- df |> get_attribute("at")
  direction <- df |> get_attribute("direction")
  proportion <- df |> get_attribute("proportion")
  cumulative <- df |> get_attribute("cumulative")
  dimensions <- df |> get_attribute("dimensions")
  every <- df |> get_attribute("every")
  if (!grepl("x", direction)) every[2] <- 1
  if (!grepl("y", direction)) every[1] <- 1
  if (!grepl("z", direction)) every[3] <- 1
  for (i in seq_along(dots)) {
    `%scale%` <- if (grepl("inverse", dots[[i]])) {
      function(a, b) a / b
    } else {
      function(a, b) a * b
    }
    if (grepl("volume", dots[[i]])) {
      if (cumulative) {
        df$neighbours <- df$neighbours %scale% (4 / 3 * pi * df$distance^3)
      } else {
        df$neighbours <- df$neighbours %scale% diff(4 / 3 * pi * at^3)
      }
    }
    if (grepl("area", dots[[i]])) {
      if (cumulative) {
        df$neighbours <- df$neighbours %scale% (pi * df$distance^2)
      } else {
        df$neighbours <- df$neighbours %scale% diff(pi * at^2)
      }
    }
    if (grepl("length", dots[[i]])) {
      if (cumulative) {
        df$neighbours <- df$neighbours %scale% (df$distance * 2)
      } else {
        df$neighbours <- df$neighbours %scale% diff(at * 2)
      }
    }
    if (grepl("neighbourhood", dots[[i]])) {
      if (nchar(direction) == 1) df <- df |> scale_by(gsub("neighbourhood", "length", dots[[i]]))
      if (nchar(direction) == 2) df <- df |> scale_by(gsub("neighbourhood", "area", dots[[i]]))
      if (nchar(direction) == 3) df <- df |> scale_by(gsub("neighbourhood", "volume", dots[[i]]))
    }
    if (grepl("proportion", dots[[i]])) {
      df$neighbours <- df$neighbours %scale% proportion[1]
    }
    if (grepl("every", dots[[i]])) {
      df$neighbours <- df$neighbours %scale% prod(1/every)
    }
    if (grepl("dimensions", dots[[i]])) {
      df$neighbours <- df$neighbours %scale% prod(dimensions)
    }
    if (is.data.frame(dots[[i]])) {
      df$neighbours <- df$neighbours / dots[[i]]$neighbours
    }
  }
  df
}

#' Remove all {morphology} object attributes to allow garbage collection
#'
#' @param .x Object to remove the attributes from.
#' @return Same object, but without the {morphology} attributes.
#' @export
finalise <- function(.x) {
  drop_attributes(.x)
}

#' Visualise the morphological description
#'
#' @param df Data frame resulting from a `describe()` call.
#' @param ... Arguments passed to `ggplot2::labs()` for labelling.
#' @return A ggplot2 plot.
#' @export
visualise <- function(df, ...) {
  ggplot2::ggplot(df) +
    (if ("neighbours" %in% names(df)) {
      ggplot2::aes(distance, neighbours, colour = if (ncol(df) > 2) name else NULL)
    } else {
      ggplot2::aes(neighbour, distance, colour = if (ncol(df) > 2) name else NULL)
    }) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(colour = NULL, ...) +
    ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0)
    )
}
