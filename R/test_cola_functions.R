

#' @title  Test cola functions
#' @description Run the main geospatial CoLa tools: Suitability to Resistance, Kernels, Corridors
#' @param run Logical. Run the tests?
#' @param zarr Logical. Run zarr corridors? Default FALSE
#' @return Path with CDPOP results
#' @examples
#' library(cola)
#' test_cola(run = TRUE )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

test_cola <- function(run = FALSE, zarr = FALSE){

  if (run){
    ## Habitat Suitability (HS) to surface resistance (SR)
    require('terra')
    library(terra)

    cat('######\n###### Here some functionalities CoLa demo:\n######\n\n\n')

    ## Intro
    plot(rast(system.file(package = 'cola', 'sampledata/sampleHS.tif')),
         main = 'Sample Habitat Suitability')

    plot(vect(system.file(package = 'cola', 'sampledata/points_sabah_50.shp')),
         add = TRUE)

    cat("plot(rast(system.file(package = 'cola', 'sampledata/sampleHS.tif')),
     main = 'Sample Surface Resistance')

plot(vect(system.file(package = 'cola', 'sampledata/points_sabah_50.shp')),
     add = TRUE)

    ")


    ## SR
    cat("\noutdir <- tempdir()
    \n################## SURFACE RESISTANCE -----------\n
    resistance <- sui2res_py(
      intif = system.file(package = 'cola', 'sampledata/sampleHS.tif'),
      outtif = file.path(outdir, 'resistance.tif'), minval = 0, maxval = 1,
      maxout = 100, shape = 1
    )

    plot(rast(file.path(outdir, 'resistance.tif')), main = 'Created surface resistance')
    ")

    outdir <- tempdir()
    resistance <- sui2res_py(
      intif = system.file(package = 'cola', 'sampledata/sampleHS.tif'),
      outtif = file.path(outdir, 'resistance.tif'), minval = 0, maxval = 1,
      maxout = 100, shape = 1
    )

    plot(rast(file.path(outdir, 'resistance.tif')), main = 'Created surface resistance')


    ## Kernels
    cat("\n\n################## KERNELS -----------\n\n
    kernels <- crk_py(
          inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
          # intif = file.path(outdir, 'resistance.tif'),
          intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
          outtif = file.path(outdir, 'kernels.tif'),
          maxdist = 1000,
          transform = 'linear',
          shape = 'no',
          volume = '1')

        plot(rast(kernels$file), main = 'Created kernels')

        ")

    kernels <- crk_py(
      inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
      # intif = file.path(outdir, 'resistance.tif'),
      intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
      outtif = file.path(outdir, 'kernels.tif'),
      maxdist = 100000,
      transform = 'no',
      shape = 'linear',
      volume = '1')

    plot(rast(kernels$file), main = 'Created kernels')


    cat("\n
    kernels_joblib <- crk_py(
          inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
          # intif = file.path(outdir, 'resistance.tif'),
          intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
          outtif = file.path(outdir, 'kernels_joblib.tif'),
          maxdist = 1000,
          transform = 'linear',
          shape = 'no',
          volume = '1')

        plot(rast(kernels_joblib$file), main = 'Created kernels parallel script')

        ")

    ## kernels
    kernels_joblib <- crk_py(
      inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
      # intif = file.path(outdir, 'resistance.tif'),
      intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
      outtif = file.path(outdir, 'kernels_joblib.tif'),
      maxdist = 100000,
      shape = 'linear',
      transform = 'no',
      volume = '1')

    plot(rast(kernels_joblib$file), main = 'Created kernels parallel script')


    cat("\n\n################## CORRIDORS -----------\n
    corridors <- lcc_py(
          inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
          intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
          outtif = file.path(outdir, 'corridors.tif'),
          maxdist = 100000, smooth = 0, tolerance = 0)

        plot(rast(corridors$file), main = 'Created corridors')

        ")

    corridors <- lcc_py(
      inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
      intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
      outtif = file.path(outdir, 'corridors_short.tif'),
      maxdist = 50000, smooth = 0, tolerance = 0)

    plot(rast(corridors$file), main = 'Created corridors')


    cat("\n\ncorridors_joblib <- lccJoblib_py(
          inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
          intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
          outtif = file.path(outdir, 'corridors_joblib.tif'),
          maxdist = 100000, smooth = 0, tolerance = 0)

    plot(rast(corridors_joblib$file), main = 'Created corridors parallel script')

        ")

    corridors_joblib <- lccJoblib_py(
      inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
      intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
      outtif = file.path(outdir, 'corridors_joblib.tif'),
      maxdist = 100000, smooth = 0, tolerance = 0)

    plot(rast(corridors_joblib$file), main = 'Created corridors parallel script')

    corridors_joblib <- lccJoblib_py(
      inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
      intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
      outtif = file.path(outdir, 'corridors_joblib.tif'),
      maxdist = 100000, smooth = 0, tolerance = 0)

    plot(rast(corridors_joblib$file), main = 'Created corridors parallel script')



    if (zarr){
      cat("\n\ncorridors_zarr <- lccZarr_py(
          inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
          intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
          outtif = file.path(outdir, 'corridors_zarr.tif'),
          maxdist = 100000, smooth = 0, tolerance = 0)

    plot(rast(corridors_zarr$file), main = 'Created corridors zarr parallel script')

        ")

      corridors_zarr <- lccZarr_py(
        inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'),
        intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
        outtif = file.path(outdir, 'corridors_zarr.tif'),
        maxdist = 100000, smooth = 0, tolerance = 0)

      plot(rast(corridors_zarr$file), main = 'Created corridors zarr parallel script')
    }

    ##
    cat("\n\n prioritization <- prio_py(tif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
                              incrk = system.file(package = 'cola', 'sampledata/kernels.tif'),
                              inlcc = system.file(package = 'cola', 'sampledata/corridors.tif'),
                              threshold = 0.7,
                              tolerance = 100000,
                              maskedcsname = file.path(outdir, 'pri_mask.tif'),
                              outshppoint = file.path(outdir, 'pri_points.shp'),
                              outshppol = file.path(outdir, 'pri_polygon.shp'),
                              outshppatch = file.path(outdir, 'pri_patch.shp'),
                              outtifpatch = file.path(outdir, 'pri_patcht.tif'),
                              outtif = file.path(outdir, 'corridors_zarr.tif')
                                )
        ")
    # outdir <- tempdir()
    prioritization <- prio_py(intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
                              incrk = system.file(package = 'cola', 'sampledata/kernels.tif'),
                              inlcc = system.file(package = 'cola', 'sampledata/corridors.tif'),
                              threshold = 0.7,
                              tolerance = 100000,
                              maskedcsname = file.path(outdir, 'pri_mask.tif'),
                              outshppoint = file.path(outdir, 'pri_points.shp'),
                              outshppol = file.path(outdir, 'pri_polygon.shp'),
                              outshppatch = file.path(outdir, 'pri_patch.shp'),
                              outtifpatch = file.path(outdir, 'pri_patcht.tif'),
                              outtif = file.path(outdir, 'corridors_zarr.tif')
                                )


    }
  return( 'Test finished')
}
## test_cola(run = TRUE, zarr = TRUE)
