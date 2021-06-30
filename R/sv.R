sv <- function(data, xcoordcol, ycoordcol, residcol, bins = 15, cutoff = NULL, ...) {
  
  # compute spatial distances
  dists <- as.matrix(dist(cbind(data[[xcoordcol]], data[[ycoordcol]]), ...))
  dists <- dists[upper.tri(dists)]
  if (is.null(cutoff)) {
    cutoff <- max(dists) / 2
  }
  dists_index <- dists <= cutoff
  dists <- dists[dists_index]
  
  # compute squared differences
  sqrdiffs <- as.matrix(dist(data[[residcol]]))^2
  sqrdiffs <- sqrdiffs[upper.tri(sqrdiffs)]
  sqrdiffs <- sqrdiffs[dists_index]
  
  # compute semivariogram classes
  dist_classes <- cut(dists, breaks = seq(0, cutoff, length.out = bins + 1))
  
  # compute squared differences within each class
  gamma <- tapply(sqrdiffs, dist_classes, function(x) mean(x) / 2)
  
  # compute pairs within each class
  np <- tapply(sqrdiffs, dist_classes, length)
  
  # compute average distance within each class
  dist <- tapply(dists, dist_classes, mean)
  
  # return output
  sv_out <- data.frame(dist, gamma, np)
  # remove NA
  sv_out <- na.omit(sv_out)
  return(sv_out)
}
