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
  kernel_shape = "diamond"
) {
  # TODO allow n times dilation/erosion of the indicator arrays, to check effect
  # of size on eg connectivity (ie find the size cut-off), or constrictivity etc.
  # TODO consider delayed_subset argument, to delay the selection of of_category
  # and in_relation_to, to look_for(), so that the main knn query can be cached
  # there? (always full array, so the query only depends on array dimensions,
  # every, direction, etc.). This may drastically speed up repeated evaluation,
  # especially when not considering connectivity or pathlengths.
  .array <- prep_array(.array)
  array_dim <- dim(.array)
  of_category <- prep_category(of_category, .array)
  in_relation_to <- prep_category(in_relation_to, .array)
  array_prop <- mean(.array %in% of_category)
  if (!is.null(in_relation_to))
    array_prop <- c(array_prop, mean(.array %in% in_relation_to))
  if (what %in% c("components", "paths")) {
    if (length(kernel_width) < ndim(.array))
      kernel_width <- rep(kernel_width, ndim(.array))
    .array <- mmand::components(
      array(.array %in% of_category, dim = dim(.array)),
      mmand::shapeKernel(kernel_width, type = kernel_shape)
    )
  }
  .df <- .array |>
    reshape2::melt() |>
    tibble::as_tibble() |>
    setNames(
      c("y", "x", "z",
        c(
          voxels = "value",
          components = "component",
          paths = "component"
        )[[what]]
      )[c(1:ndim(.array), 4)]
    )
  if (what == "paths") {
    yxz_multiply <- 1 + (kernel_width == 1) * 9 # break connection width 1 dims
    n_voxels <- sum(!is.na(.df$component))
    knn <- nabor::knn(
      .df[!is.na(.df$component), c("y", "x", "z")] *
        matrix(yxz_multiply, byrow = TRUE, nrow = n_voxels, ncol = 3),
      k = 7,
      radius = c(diamond = 1, disc = sqrt(2), box = sqrt(3))[kernel_shape]
    )
    el <- knn$nn.idx[, c(1, 2)]
    for (i in 3:7) el <- rbind(el, knn$nn.idx[, c(1, i)])
    el <- el[el[, 2] != 0, ] # remove non-edges
    # from knn ids to original array ids
    el[] <- xyz_to_id(.df[!is.na(.df$component),], array_dim)[c(el)]
    # NOTE {cuRnet} seems to provide a GPU-enabled alternative to {igraph}
    graph <- igraph::graph_from_edgelist(el, directed = FALSE) |>
      igraph::simplify() # removes duplicates
    # FIXME in case kernel shape is not diamond, we need weighted edges!
    # graph[1, 2, attr="weight"]<- 5
    # graph[from=1:3, to=c(2,3,5), attr = "weight"] <- c(1,-1,4)
    return(
      .df[!is.na(.df$component),][, c(2, 1, 3, 4)] |>
        set_attribute("dimensions", array_dim) |>
        set_attribute("proportion", array_prop) |>
        set_attribute("graph", graph)
    )
  }
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
#'   This can be a vector (row, column, layer) to thin the different dimensions
#'   differently. Alternatively it can be a list of two vectors, where the first
#'   one is used for the category of interest (`of_category` in `look_at()`) and
#'   the second one is for the neighbour category (`in_relation_to` in
#'   `look_at()`, or again `of_category` for a self-query).
#' @return Data frame (tibble) with `x`, `y`, and/or `z`, and `value` or
#'   `component` columns.
#' @export
look_in <- function(.df, direction = "xyz", every = 1) {
  attributes_df <- get_attributes(.df)
  if (is.list(every)) {
    if (length(every) != 2)
      std::err("e When {.arg every} is a list, it should have length 2.")
    if (!"value" %in% names(.df)) {
      .df1 <- .df |> transform(value = 1)
      .df2 <- .df |> transform(value = 2)
    } else {
      .df1 <- .df[.df$value == 1, ]
      .df2 <- .df[.df$value == 2, ]
      if (nrow(.df2) == 0) {
        .df2 <- .df1
        .df2$value <- 2
      }
    }
    every[[1]] <- if (length(every[[1]]) == 1) rep(every[[1]], 3) else
      every[[1]]
    every[[1]] <- 1/every[[1]]
    if (any(every[[1]] != 1)) {
      if ("x" %in% names(.df1)) .df1 <- .df1[.df1$x %% every[[1]][2] == 0,]
      if ("y" %in% names(.df1)) .df1 <- .df1[.df1$y %% every[[1]][1] == 0,]
      if ("z" %in% names(.df1)) .df1 <- .df1[.df1$z %% every[[1]][3] == 0,]
    }
    every[[2]] <- if (length(every[[2]]) == 1) rep(every[[2]], 3) else
      every[[2]]
    every[[2]] <- 1/every[[2]]
    if (any(every[[2]] != 1)) {
      if ("x" %in% names(.df2)) .df2 <- .df2[.df2$x %% every[[2]][2] == 0,]
      if ("y" %in% names(.df2)) .df2 <- .df2[.df2$y %% every[[2]][1] == 0,]
      if ("z" %in% names(.df2)) .df2 <- .df2[.df2$z %% every[[2]][3] == 0,]
    }
    .df <- rbind(.df1, .df2) |> set_attributes(attributes_df)
  } else {
    every <- if (length(every) == 1) rep(every, 3) else every
    every <- 1/every
    if (any(every != 1)) {
      if ("x" %in% names(.df)) .df <- .df[.df$x %% every[2] == 0,]
      if ("y" %in% names(.df)) .df <- .df[.df$y %% every[1] == 0,]
      if ("z" %in% names(.df)) .df <- .df[.df$z %% every[3] == 0,]
    }
  }
  .df$id <- xyz_to_id(.df, attributes_df$dimensions)
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
#' @param error Approximate error bound, for approximate nearest neighbour
#'   search.
#' @return List with the `distance` matrix and `xyz` data frame. When looking at
#'   components, additionally a `connected` matrix.
#' @export
look_for <- function(
  .df,
  neighbours = 1,
  within = Inf,
  from_border = FALSE,
  error = 0
) {
  attributes_df <- get_attributes(.df)
  not_dims <- c("id", "value", "border", "component")
  .df <- .df |> finalise()
  if (is.infinite(neighbours)) {
    if (is.infinite(within))
      std::err("e Either neighbours or within should be finite.")
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
    dims <- sort(names(.df)[!names(.df) %in% not_dims])
    # TODO consider caching this query, for repeated evaluation, after the
    # delayed_subset arg is implemented in look_at()
    knn <- nabor::knn(
      ato[, dims],
      if (from_border) afrom[, dims] else afrom[afrom$border > within, dims],
      k = neighbours + if (self) 1 else 0,
      radius = within,
      eps = error
    )
    knn$nn.idx[knn$nn.idx == 0] <- NA
    if ("component" %in% names(.df)) {
      component <- knn$nn.idx
      component[] <- ato$component[knn$nn.idx]
    }
    return(
      list(
        connected = if ("component" %in% names(.df))
          (component == component[, 1])[, -1, drop = FALSE] else NULL,
        distance = if (self) knn$nn.dists[, -1, drop = FALSE] else knn$nn.dists,
        xyz_from = if (from_border) afrom[, c(dims, "id")] else
          afrom[afrom$border > within, c(dims, "id")],
        xyz_to = ato[, c(dims, "id")],
        id = if (self) knn$nn.idx[, -1, drop = FALSE] else knn$nn.idx
      ) |>
        set_attributes(attributes_df) |>
        set_attribute("within", within)
    )
  }
  if ("component" %in% names(.df)) {
    dims <- sort(names(.df)[!names(.df) %in% not_dims])
    knn <- nabor::knn(
      .df[, dims],
      if (from_border) .df[, dims] else .df[.df$border > within, dims],
      k = neighbours + 1,
      radius = within,
      eps = error
    )
    knn$nn.idx[knn$nn.idx == 0] <- NA
    component <- knn$nn.idx
    component[] <- .df$component[knn$nn.idx]
    return(
      list(
        connected = (component == component[, 1])[, -1, drop = FALSE],
        distance = knn$nn.dists[, -1, drop = FALSE],
        xyz_from = if (from_border) .df[, c(dims, "id")] else
          .df[.df$border > within, c(dims, "id")],
        xyz_to = .df[, c(dims, "id")],
        id = knn$nn.idx[, -1, drop = FALSE]
      ) |>
        set_attributes(attributes_df) |>
        set_attribute("within", within)
    )
  }
  # NOTE xyz_from and xyz_to have modified coordinates, but original array id.
  # The id in the list is knn$nn.idx id, so reindexed after applying every.
  std::err("e No value or component column found.")
}

#' Describe morphology using summary functions
#'
#' @param .list List with the `distance` matrix and `within`, and potentially
#'   a `connected` matrix, typically output from `look_for()`. Alternatively,
#'   this can be a data frame resulting from a `describe()` including a `name`
#'   argument, to add additional described summaries.
#' @param what Summary function.
#' @param of Quantity to summarise, and return as `y`. Can be `"neighbours"`
#'   (default), `"distance"` or `"pathlength"`.
#' @param in_function_of Quantity to return as `x`. The default is to use
#'   `"distance"` in case `of` is `"neighbours"` or `"pathlength"`, and
#'   `"neighbours"` otherwise.
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
  in_function_of = switch(
    of,
    neighbours = "distance",
    distance = "neighbours",
    pathlength = "distance"
  ),
  at = NULL,
  connected = NULL,
  cumulative = FALSE,
  name = NULL
) {
  # TODO allow returning point cloud, without summary, through function that
  # returns vector rather than summary? then identity() as what? or other funs
  # that compress the data into less samples? eg all percentiles? for boxplot,
  # violin, etc?
  # TODO allow merging number of neighbours with xyz df, so we can visualise the
  # descriptors on their corresponding voxels!! maybe think of other verb for
  # doing this? This would correspond somehow to the topographical correlation
  # map of @bull2023?
  .df <- NULL
  list_attributes <- NULL
  if (.list |> has_attribute("list")) {
    .df <- .list
    .list <- .list |> get_attribute("list")
    .df <- .df |> set_attribute("list", NULL)
  }
  list_attributes <- get_attributes(.list)
  nn <- .list
  if (!is.null(connected)) {
    if (connected) nn$distance[!nn$connected] <- NA
    if (!connected) nn$distance[nn$connected] <- NA
  }
  if (of == "pathlength" | in_function_of == "pathlength") {
    graph <- .list |> get_attribute("graph")
    dimensions <- .list |> get_attribute("dimensions")
    paths <- igraph::distances(graph, nn$xyz_from$id, nn$xyz_to$id)
    nn$pathlength <- nn$distance
    nn$pathlength[] <- paths[
      data.frame(
        i = rep(1:nrow(nn$xyz_from), ncol(nn$distance)),
        j = c(nn$id)
      ) |> as.matrix()
    ]
    nn$pathlength[is.infinite(nn$pathlength)] <- NA
    # TODO add nn$pathlength to attributes for repeated describe() calls?
  }
  if (of == "neighbours") {
    if (in_function_of == "distance") {
      x <- nn$distance
    } else {
      if (in_function_of == "pathlength") {
        x <- nn$pathlength
      } else {
        std::err("e In function of should be distance or pathlength.")
      }
    }
    if (is.null(at)) {
      if (is.finite(.list |> get_attribute("within"))) {
        at <- 0:(.list |> get_attribute("within"))
      } else {
        at <- 0:ceiling(max(x[is.finite(x)]))
      }
    }
    results <- at * NA
    for (i in seq_along(results)) {
      results[i] <- what(matrixStats::rowSums2(x <= at[i], na.rm = TRUE))
    }
    if (cumulative) {
      df <- tibble::tibble(x = at, y = results, name = name)
    } else {
      df <- tibble::tibble(
        x = (at[1:(length(at) - 1)] + at[2:length(at)]) / 2,
        y = diff(results),
        name = name
      )
    }
  }
  if (in_function_of == "neighbours") {
    if (of == "distance") {
      y <- nn$distance
    } else {
      if (of == "pathlength") {
        y <- nn$pathlength
      } else {
        std::err("e Of should be distance or pathlength.")
      }
    }
    if (is.null(at)) {
      at <- 1:ncol(y)
    }
    results <- at * NA
    if (cumulative) {
      for (i in seq_along(results)) {
        results[i] <- what(y[, 1:i])
      }
    } else {
      for (i in seq_along(results)) {
        results[i] <- what(y[, i])
      }
    }
    df <- tibble::tibble(x = at, y = results, name = name)
  }
  if (of %in% c("distance", "pathlength") &
      in_function_of %in% c("distance", "pathlength")) {
    if (of == "distance") y <- nn$distance
    if (of == "pathlength") y <- nn$pathlength
    if (in_function_of == "distance") x <- nn$distance
    if (in_function_of == "pathlength") x <- nn$pathlength
    if (is.null(at)) {
      if (is.finite(.list |> get_attribute("within"))) {
        at <- 0:(.list |> get_attribute("within"))
      } else {
        at <- 0:ceiling(max(x[is.finite(x)]))
      }
    }
    results <- at * NA
    if (cumulative) {
      for (i in seq_along(results)) {
        results[i] <- what(y[x <= at[i]], na.rm = TRUE)
      }
      df <- tibble::tibble(x = at, y = results, name = name)
    } else {
      for (i in seq_along(results)[-1]) {
        results[i] <- what(y[
          x <= at[i] & x > at[i - 1]
        ], na.rm = TRUE)
      }
      df <- tibble::tibble(
        x = (at[1:(length(at) - 1)] + at[2:length(at)]) / 2,
        y = results[-1],
        name = name
      )
    }
  }
  if (!is.null(.df)) df <- rbind(.df, df)
  if (!is.null(name)) df <- df |> set_attribute("list", .list)
  if (!is.null(list_attributes)) df <- df |> set_attributes(list_attributes)
  df |>
    set_attribute("at", at) |>
    set_attribute("cumulative", cumulative)
}

#' Scale the morphological description
#'
#' @param df Data frame
#' @param ... Things to scale the description by. Can be `"neighbourhood"`,
#'   which automatically chooses between alternative options `"volume"`,
#'   `"area"` or `"length"`, or `"proportion"`, `"inverse proportion"`, another
#'   description data frame with the same distances, ...
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
        df$y <- df$y %scale% (4 / 3 * pi * df$x^3)
      } else {
        df$y <- df$y %scale% diff(4 / 3 * pi * at^3)
      }
    }
    if (grepl("area", dots[[i]])) {
      if (cumulative) {
        df$y <- df$y %scale% (pi * df$x^2)
      } else {
        df$y <- df$y %scale% diff(pi * at^2)
      }
    }
    if (grepl("length", dots[[i]])) {
      if (cumulative) {
        df$y <- df$y %scale% (df$x * 2)
      } else {
        df$y <- df$y %scale% diff(at * 2)
      }
    }
    if (grepl("neighbourhood", dots[[i]])) {
      if (nchar(direction) == 1) df <- df |>
          scale_by(gsub("neighbourhood", "length", dots[[i]]))
      if (nchar(direction) == 2) df <- df |>
          scale_by(gsub("neighbourhood", "area", dots[[i]]))
      if (nchar(direction) == 3) df <- df |>
          scale_by(gsub("neighbourhood", "volume", dots[[i]]))
    }
    if (grepl("proportion", dots[[i]])) {
      df$y <- df$y %scale% proportion[1]
    }
    if (grepl("every", dots[[i]])) {
      df$y <- df$y %scale% prod(1/every)
    }
    if (grepl("dimensions", dots[[i]])) {
      df$y <- df$y %scale% prod(dimensions)
    }
    if (grepl("distance", dots[[i]])) {
      # TODO safer approach for when distance is not x? maybe use "x"?
      df$y <- df$y %scale% df$x
    }
    if (is.data.frame(dots[[i]])) {
      df$y <- df$y / dots[[i]]$y
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
#' @param df Data frame resulting from a `describe()` call, or just a numeric
#'   vector.
#' @param ... Arguments passed to `ggplot2::labs()` for labelling.
#' @return A ggplot2 plot.
#' @export
visualise <- function(df, ...) {
  ggplot2::ggplot(df) +
    ggplot2::aes(x, y, colour = if ("name" %in% names(df)) name else NULL) +
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
