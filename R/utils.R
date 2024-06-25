prep_array <- function(a) {
  if (is.list(a)) return(lapply(a, prep_array))
  if (!is.integer(a)) {
    if (all(is.wholenumber(a))) {
      storage.mode(a) <- "integer"
      std::err("i Array is converted to integers.")
    } else {
      std::err("! Expecting an object of type integer.")
      std::err("e Issue with object type.")
    }
  }
  if (is.array(a)) {
    if (length(dim(a)) == 3) return(a)
    if (length(dim(a)) == 2) return(array(c(a), dim = c(dim(a), 1)))
    if (length(dim(a)) == 1) return(array(c(a), dim = c(1, dim(a), 1)))
  }
  if (is.integer(a) & is.atomic(a)) return(array(a, dim = c(1, length(a), 1)))
  std::err("! Argument {.arg a} is an unsupported object.")
  std::err("e Issue with object type.")
}
prep_category <- function(x, .array) {
  if (all(x >= 0L)) return(x)
  if (all(x < 0L)) {
    .array <- unique(c(.array))
    return(.array[!.array %in% -x])
  }
  std::err("e Negative and positive categories cannot be mixed.")
}
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
make3d <- function(a) {
  a |> add_border(c("x", "y", "z"))
}
make1d <- function(a, direction = "x") {
  a <- a |> add_border(direction)
  value_column <- if ("value" %in% names(a)) "value" else NULL
  component_column <- if ("component" %in% names(a)) "component" else NULL
  if (direction == "x") a$x <- a$x + max(a$x) * (2*a$y - 2) + max(a$x) * max(a$y) * (2*a$z - 2)
  if (direction == "y") a$y <- a$y + max(a$y) * (2*a$x - 2) + max(a$y) * max(a$x) * (2*a$z - 2)
  if (direction == "z") a$z <- a$z + max(a$z) * (2*a$y - 2) + max(a$z) * max(a$y) * (2*a$x - 2)
  return(a[, c(direction, value_column, component_column, "border", "id")])
}
make2d <- function(a, direction = "xy") {
  value_column <- if ("value" %in% names(a)) "value" else NULL
  component_column <- if ("component" %in% names(a)) "component" else NULL
  if (direction %in% c("xy", "yx")) {
    a <- a |> add_border(c("x","y"))
    a$x <- a$x + max(a$x) * (2*a$z - 2)
    return(a[,c("x", "y", value_column, component_column, "border", "id")])
  }
  if (direction %in% c("xz", "zx")) {
    a <- a |> add_border(c("x","z"))
    a$x <- a$x + max(a$x) * (2*a$y - 2)
    return(a[,c("x", "z", value_column, component_column, "border", "id")])
  }
  if (direction %in% c("zy", "yz")) {
    a <- a |> add_border(c("y","z"))
    a$y <- a$y + max(a$y) * (2*a$x - 2)
    return(a[,c("y", "z", value_column, component_column, "border", "id")])
  }
  stop("Wrong direction")
}
add_border <- function(df, dimnames) {
  border <- matrix(nrow = nrow(df), ncol = length(dimnames))
  dims <- get_attribute(df, "dim")
  names(dims) <- c("y", "x", "z")
  for (i in seq_along(dimnames)) {
    dim_i <- df[[dimnames[i]]]
    dim_center <- dims[dimnames[i]] / 2 + 0.5
    border[, i] <- dim_center - 0.5 - abs(dim_i - dim_center)
  }
  df |> transform(border = matrixStats::rowMins(border)) |> tibble::as_tibble()
}
ndim <- function(...) {
  length(dim(...))
}
drop_attributes <- function(x) {
  attr_names <- names(attributes(x))
  attr_names <- attr_names[grepl("^morphology.", attr_names)]
  if (length(attr_names) == 0) return(x)
  for (i in seq_along(attr_names)) {
    attr(x, attr_names[i]) <- NULL
  }
  x
}
set_attributes <- function(x, l) {
  for (i in seq_along(l)) x <- set_attribute(x, names(l)[i], l[[i]])
  x
}
set_attribute <- function(x, name, a) {
  attr(x, paste0("morphology.", name)) <- a
  x
}
get_attributes <- function(x) {
  attr_names <- names(attributes(x))
  attr_names <- attr_names[grepl("^morphology.", attr_names)]
  a <- attributes(x)[attr_names]
  names(a) <- substr(names(a), 12, nchar(names(a)))
  a
}
get_attribute <- function(x, name) {
  attr(x, paste0("morphology.", name))
}
has_attribute <- function(x, name) {
  attr_names <- names(attributes(x))
  paste0("morphology.", name) %in% attr_names
}
xyz_to_id <- function(df, dims) {
  (df$z - 1) * prod(dims[1:2]) + (df$x - 1) * dims[2] + df$y
}
