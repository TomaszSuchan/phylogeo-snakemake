#!/usr/bin/env Rscript
# Download Natural Earth 10m HR zip (CDN) and extract the global GeoTIFF once per project+basemap kind.
# Params: ne_kind only — not CRS/bbox — so crop/crs changes do not rerun this rule.
# Outputs: hr_bundle.zip, hr_global.tif under results/{project}/maps/naturalearth_hr_global/

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

out_zip <- snakemake@output[["zip"]]
out_tif <- snakemake@output[["tif"]]
ne_kind <- snakemake@params[["ne_kind"]]

message("\n=== FETCH NATURAL EARTH HR GLOBAL (zip + extract) ===\n")
message(sprintf("kind: %s\n", ne_kind))

if (is.null(ne_kind) || identical(ne_kind, "NULL") || !nzchar(as.character(ne_kind)[1L])) {
  stop("fetch_naturalearth_hr_global: ne_kind is missing (map basemap is not ne1/gray_earth?)")
}

ne_bundle <- switch(
  as.character(ne_kind)[1L],
  ne1 = list(
    zip = "NE1_HR_LC_SR_W_DR.zip",
    tif_stem = "NE1_HR_LC_SR_W_DR",
    min_zip_bytes = 310000000L
  ),
  gray_earth = list(
    zip = "GRAY_HR_SR_OB_DR.zip",
    tif_stem = "GRAY_HR_SR_OB_DR",
    min_zip_bytes = 80000000L
  ),
  stop("ne_kind must be ne1 or gray_earth, got: ", ne_kind)
)

ne_cdn_base <- "https://naciscdn.org/naturalearth/10m/raster"
zip_url <- paste0(ne_cdn_base, "/", ne_bundle$zip)

message(sprintf("zip URL: %s\n", zip_url))
message(sprintf("outputs: %s , %s\n", out_zip, out_tif))

dir.create(dirname(out_zip), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)

# If ne1 product or ne1↔gray switched, drop stale hr_bundle.zip (same size threshold is not enough).
bundle_id <- paste(as.character(ne_kind)[1L], ne_bundle$zip, sep = ":")
bundle_marker <- file.path(dirname(out_zip), "bundle_id.txt")
if (file.exists(bundle_marker)) {
  prev <- tryCatch(readLines(bundle_marker, n = 1L, warn = FALSE), error = function(e) "")
  if (!identical(trimws(as.character(prev)), bundle_id)) {
    message("Natural Earth bundle id changed; removing old zip/tif/partials.\n")
    unlink(out_zip, force = TRUE)
    unlink(out_tif, force = TRUE)
    unlink(paste0(out_zip, ".part"), force = TRUE)
    unlink(paste0(out_zip, ".download"), force = TRUE)
  }
}

zip_ok <- file.exists(out_zip) && (file.info(out_zip)$size >= ne_bundle$min_zip_bytes)

#' curl HTTP/1.1 + resume (see cache_naturalearth_basemap history).
download_ne_zip <- function(url, dest_zip, min_bytes) {
  part <- paste0(dest_zip, ".part")
  legacy <- paste0(dest_zip, ".download")
  if (file.exists(legacy)) {
    unlink(legacy, force = TRUE)
  }
  curl_bin <- Sys.which("curl")
  if (nzchar(curl_bin)) {
    has_partial <- isTRUE(file.exists(part)) && (file.info(part)$size > 0L)
    base_args <- c(
      "-L", "-f", "-S", "--http1.1",
      "--connect-timeout", "60",
      "--retry", "12", "--retry-delay", "10", "--retry-all-errors"
    )
    if (has_partial) {
      message("Resuming partial zip download (curl -C -)...\n")
      args <- c(base_args, "-C", "-", "-o", part, url)
    } else {
      if (file.exists(part)) {
        unlink(part, force = TRUE)
      }
      message("Downloading zip with curl (HTTP/1.1, retries)...\n")
      args <- c(base_args, "-o", part, url)
    }
    out <- system2(curl_bin, args, stdout = FALSE, stderr = TRUE)
    st <- attr(out, "status")
    if (!is.null(st) && st != 0L) {
      stop(
        "curl download failed (exit ", st, "): ",
        paste(out, collapse = "\n")
      )
    }
  } else {
    message("curl not on PATH; falling back to download.file.\n")
    options(timeout = max(7200, getOption("timeout")))
    st <- utils::download.file(url, part, mode = "wb", quiet = FALSE, method = "libcurl")
    if (st != 0L) {
      unlink(part, force = TRUE)
      stop("download.file failed: ", url)
    }
  }
  sz <- file.info(part)$size
  if (is.na(sz) || sz < min_bytes) {
    stop(
      "Downloaded zip too small (",
      sz,
      " bytes); expected at least ",
      min_bytes
    )
  }
  if (file.exists(dest_zip)) {
    unlink(dest_zip, force = TRUE)
  }
  if (!file.rename(part, dest_zip)) {
    file.copy(part, dest_zip, overwrite = TRUE)
    unlink(part, force = TRUE)
  }
  invisible(TRUE)
}

if (!zip_ok) {
  part <- paste0(out_zip, ".part")
  legacy_dl <- paste0(out_zip, ".download")
  if (file.exists(legacy_dl) && !file.exists(part)) {
    if (!file.rename(legacy_dl, part)) {
      file.copy(legacy_dl, part, overwrite = TRUE)
      unlink(legacy_dl, force = TRUE)
    }
    message("Renamed legacy .download partial to .part for curl resume.\n")
  }
  if (file.exists(out_zip) && file.info(out_zip)$size > 0L &&
      file.info(out_zip)$size < ne_bundle$min_zip_bytes) {
    if (file.exists(part)) {
      unlink(part, force = TRUE)
    }
    if (!file.rename(out_zip, part)) {
      file.copy(out_zip, part, overwrite = TRUE)
      unlink(out_zip, force = TRUE)
    }
    message("Reusing partial zip for curl resume (-C -).\n")
  }
  message("Fetching HR zip...\n")
  download_ne_zip(zip_url, out_zip, ne_bundle$min_zip_bytes)
} else {
  message("Using existing zip: ", out_zip, "\n")
}

tif_pat <- paste0("^", ne_bundle$tif_stem, "\\.tif$")
need_extract <- !file.exists(out_tif) || file.info(out_tif)$size < 1000L
if (need_extract) {
  message("Extracting global GeoTIFF from zip...\n")
  tmp_unzip <- file.path(dirname(out_tif), "_unzip_tmp")
  if (dir.exists(tmp_unzip)) {
    unlink(tmp_unzip, recursive = TRUE, force = TRUE)
  }
  dir.create(tmp_unzip, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(out_zip, exdir = tmp_unzip, overwrite = TRUE)
  hits <- list.files(
    tmp_unzip,
    pattern = tif_pat,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) == 0L) {
    unlink(tmp_unzip, recursive = TRUE, force = TRUE)
    stop("After unzip, could not find ", ne_bundle$tif_stem, ".tif under ", tmp_unzip)
  }
  src_tif <- hits[[1L]]
  if (file.exists(out_tif)) {
    unlink(out_tif, force = TRUE)
  }
  if (!file.rename(src_tif, out_tif)) {
    ok <- file.copy(src_tif, out_tif, overwrite = TRUE)
    if (!ok) {
      unlink(tmp_unzip, recursive = TRUE, force = TRUE)
      stop("Could not place global tif at ", out_tif)
    }
  }
  unlink(tmp_unzip, recursive = TRUE, force = TRUE)
  message(sprintf("Wrote global tif: %s\n", out_tif))
} else {
  message("Using existing global tif: ", out_tif, "\n")
}

writeLines(bundle_id, bundle_marker, useBytes = FALSE)

sink(type = "message")
sink(type = "output")
close(log_file)
