#' Fully penetrable spheres generator
#'
#' @param dimensions Dimensions of the generated array.
#' @param proportion Proportion of the spheres category in the generated array.
#' @param radius Radius of the spheres.
#' @return Array with fully penetrable spheres
#' @export
gen_fully_penetrable_spheres <- function(
  dimensions = rep(200, 3),
  proportion = 0.5,
  radius = 15
) {
  proportion <- 1 - proportion
  n <- round(-log(proportion)/(4/3*pi*radius^3)*(prod(dimensions - 1 + radius * 2)))
  spheres <- array(0L, dim = dimensions)
  centers <- tibble::tibble(
    x = runif(n, 1 - radius, dimensions[2] + radius),
    y = runif(n, 1 - radius, dimensions[1] + radius),
    z = runif(n, 1 - radius, dimensions[3] + radius)
  )
  df <- reshape2::melt(spheres) |>
    tibble::as_tibble()
  nn <- nabor::knn(
    df[, 1:3],
    centers,
    k = ceiling(4/3*pi*radius^3*1.1),
    radius = radius
  )
  spheres[unique(c(nn$nn.idx))] <- 1L
  spheres
}
