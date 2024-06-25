# {morphology} (development version)


- [ ] check if we can also merge things by eg rbind just before describe, or by
  providing a list of data frames?
  eg useful for isotropic lineal path, for which we need to combine results
  of x, y and z directional lineal path with corresponding 1D kernel for the
  components
- [ ] check if we can add api for uncertainty eval, where original array is
  divided in subarrays, and certain workflow is applied, then we take mean and
  sd over all subarrays? maybe allow input .array to be list of arrays, and
  foresee function to subdivide? or just split, apply, combine, with fn 
  wrapper around workflow?
- [ ] add J function to benchmarks? (1-nearest neighbour distance) /
  (1-empty space), indicates regularity vs clustering
- [ ] replace str() with rui::inspect() after s3 clash has been solved?
- [ ] rename every to keep
- [ ] expose nabor::knn() eps arguments for approximate nn, and speeding up?
- [ ] evaluate whether it would be possible to do incremental updates of
  descriptors based on voxel ids of the subset of modified voxels?
  if straightforward, consider implementation in framework of sa approaches?
  if practically difficult, forget about it
- [ ] add functions to generate model systems? fully penetrable spheres, others?
- [x] Introduction of minus sampling edge correction, see baddeley1993
- [x] Initial set of functions for full morphological description pipechains
