equalarea_project <- function (path) 
{
  if (grepl(".zip", path, fixed = TRUE)) {
    vsizip <- "/vsizip/"
    out_path <- gsub(".zip", ".tif", path, fixed = TRUE)
  }
  else {
    vsizip <- ""
    out_path <- gsub("_lonlat.tif", ".tif", path, fixed = TRUE)
  }
  stopifnot(out_path != path)
  if (!file.exists(out_path)) {
    r <- raster::raster(paste0(vsizip, path))
    message(paste0("Projecting ", path, " from native WGS84 to Behrmann Equal Areas. This might take a few minutes."))
    if (grepl("/FW_", path, ignore.case = FALSE, fixed = TRUE)) {
      res = 815
    }
    else {
      res = 7000
    }
    out <- raster::projectRaster(r, crs = sdmpredictors::equalareaproj, 
                                 method = "ngb", res = res)
    raster::writeRaster(out, out_path)
  }
  stopifnot(as.character(proj4string(raster(out_path))) == as.character(sdmpredictors::equalareaproj))
  return(out_path)
}

load_layers <- function (layercodes, equalarea = FALSE, rasterstack = TRUE, 
                         datadir = NULL) 
{
  if (is.na(equalarea) || !is.logical(equalarea) || length(equalarea) != 
      1) {
    stop("equalarea should be TRUE or FALSE")
  }
  if (is.data.frame(layercodes)) {
    layercodes <- layercodes$layer_code
  }
  info <- get_layers_info(layercodes)
  counts <- table(info$common$time)
  if (length(unique(info$common$cellsize_equalarea)) > 1) {
    stop("Loading layers with different cellsize is not supported")
  }
  else if (sum(counts) != NROW(layercodes)) {
    layers <- info$common$layer_code
    stop(paste0("Following layer codes where not recognized: ", 
                paste0(layercodes[!(layercodes %in% layers)], collapse = ", ")))
  }
  if (max(counts) != NROW(layercodes)) {
    warning("Layers from different eras (current, future, paleo) are being loaded together")
  }
  if (sdmpredictors:::gdal_is_lower_than_3()) {
    warning("GDAL is lower than version 3. Consider updating GDAL to avoid errors.")
  }
  datadir <- sdmpredictors:::get_datadir(datadir)
  get_layerpath <- function(layercode) {
    layer_url <- subset(info$common, info$common$layer_code == 
                          layercode)$layer_url
    if (grepl(".zip", layer_url, ignore.case = FALSE, fixed = TRUE)) {
      ext <- "_lonlat.zip"
    }
    else {
      ext <- "_lonlat.tif"
    }
    path <- paste0(datadir, "/", layercode, ext)
    if (!file.exists(path)) {
      ok <- -1
      on.exit({
        if (ok != 0 && file.exists(path)) {
          file.remove(path)
        }
      })
      ok <- utils::download.file(layer_url, path, method = "auto", 
                                 quiet = FALSE, mode = "wb")
    }
    ifelse(file.exists(path), path, NA)
  }
  paths <- sapply(layercodes, get_layerpath)
  if (equalarea) {
    customSupressWarning <- function(w) {
      if (any(grepl("point", w))) {
        invokeRestart("muffleWarning")
      }
    }
    paths <- withCallingHandlers(sapply(paths, equalarea_project), 
                                 warning = customSupressWarning)
  }
  if (rasterstack) {
    logical_zip <- grepl(".zip", paths, fixed = TRUE)
    if (all(logical_zip) | !any(logical_zip)) {
      if (all(logical_zip)) {
        vsizip <- "/vsizip/"
      }
      else if (!any(logical_zip)) {
        vsizip <- ""
      }
      st <- raster::stack(paste0(vsizip, paths))
      if ("layer" %in% names(st)) {
        names(st) <- layercodes
      }
      else (names(st) <- sub("_lonlat$", "", names(st)))
      return(st)
    }
    else {
      stop("Rasterstack for zipped and non-zipped files not supported. Try `rasterstack = FALSE`")
    }
  }
  else {
    return(lapply(paths, function(path) {
      if (grepl(".zip", path, fixed = TRUE)) {
        vsizip <- "/vsizip/"
      } else {
        vsizip <- ""
      }
      r <- raster::raster(paste0(vsizip, path))
      if ("layer" %in% names(r)) {
        names(r) <- layercodes
      } else (names(r) <- sub("_lonlat$", "", names(r)))
      r
    }))
  }
}
