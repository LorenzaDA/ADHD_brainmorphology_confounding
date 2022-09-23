qdecr_load <- function(path, reload = FALSE) {
  if (!file.exists(path)) stop("Provided path does not exist.")
  if (dir.exists(path)) {
    nc <- nchar(path)
    end <- substring(path, nc, nc)
    if (end == "/") path <- substring(path, 1, nc - 1)
    path <- paste0(path, "/", basename(path), ".rds")
  }
  path <- normalizePath(path)
  x <- readRDS(path)
  x$paths$rds <- path
  if (dirname(path) != x$paths[["final_path"]])
    warning("The project seems to be in a different ",
            "path than where it was originally run! ",
            "Your paths may be outdated. Use `qdecr_update_path` ",
            "to update all the paths.")
  if (reload) x <- reload(x)
  x
}

qdecr_update_path <- function(vw, dir_fshome, dir_subj, dir_project = dirname(vw$paths$rds), mask_path = NULL, qdecr_mask = TRUE, overwrite = FALSE) {
  to_sub <- vw$paths[["dir_out"]]
  path_names <- names(vw$paths)
  path <- dir_project
  vw$paths <- as.list(sub(to_sub, path, vw$paths))
  names(vw$paths) <- path_names
  vw$paths$dir_tmp <- path
  vw$paths$dir_tmp2 <- path
  vw$stack[-1] <- rapply(vw$stack[-1], function(x) sub(to_sub, path, x), how = "replace")
  vw$paths$dir_fshome <- dir_fshome
  vw$paths$dir_subj <- dir_subj
  if(qdecr_mask) mask_path <- system.file("extdata", paste0(vw$input$hemi, ".fsaverage.cortex.mask.mgh"), package = "QDECR")
  if(!is.null(mask_path)) vw$paths$mask_path <- mask_path
  if (overwrite) {
    message("Overwriting .rds file")
    saveRDS(vw, vw$paths$rds)
  } else {
    warning("Note that the paths are updated, but the stored project file is still out of date. Use `overwrite = TRUE` to also update the .rds file.")
  }
  return(vw)
}

quick_compose <- function(temp_mgh_file, hemi, cs, head = FALSE, text = NULL, cs_size = 16, cs_value = c("0", "0"), cs_offset = c("+5+6", "+2+244"), head_size = 150, text_size = 100, alt_text = NULL, font = "") {
  name <- tools::file_path_sans_ext(temp_mgh_file)
  lo <- c("lateral", "medial")
  lo2 <- c("", "")
  lo3 <- c("east", "west")
  
  if (hemi == "rh") {
    lo <- lo[c(2,1)]
    lo2 <- lo2[c(2,1)]
  }
  header <- if (hemi == "lh") "Left" else "Right"
  nl <- paste(name, lo, "tiff", sep = ".")
  csh <- 0
  
  
  for (i in seq_along(nl)) {
    tt <- magick::image_read(nl[i])
    tt <- magick::image_fill(tt, "white")
    w <- magick::image_info(tt)$width
    h <- magick::image_info(tt)$height
    hw <- w/2.2
    hh <- if (i %in% 1:2) h/1.7 else h/3
    fw <- 1.655
    fh <- if (i %in% 1:2) 2.7 else 0.67
    crop_string <- paste0(hw, "x", hh, "+", hw / fw, "+", hh / fh)
    nn <- magick::image_crop(tt, crop_string)
    nn <- magick::image_annotate(nn, lo2[i], size = text_size * 0.75, color = "black", gravity = lo3[i], font = font)
    nn <- magick::image_border(nn, "white", "18x18")
    nn <- magick::image_border(nn, "black", "2x2")
    nn <- magick::image_border(nn, "white", "1x1")
    nn2 <- if (i == 1) nn else c(nn2, nn)
    if (i %in% c(1, 3)) csh <- csh + hh
    if (cs && i == 2) {
      i = 3
      crop_string <- paste0(w / 9, "x", csh, "+", w - w / 9, "+", 0)
      nn <- magick::image_crop(tt, "60x265+1470+10")
      nn <- magick::image_annotate(nn, cs_value[1], size = cs_size, color = "black", location = cs_offset[1], font = font)
      nn <- magick::image_annotate(nn, cs_value[2], size = cs_size, color = "black", location = cs_offset[2], font = font)
      nn <- magick::image_scale(nn, paste0("x", csh))
      nn <-magick::image_border(nn, "white", "0x42")
      nn2 <- c(nn2, nn)
    }
  }
  ia3 <- magick::image_append(c(nn2[1], nn2[2]))
  
  hp_width <- 150
  
  if (!is.null(text)){
    w <- magick::image_info(ia3)$width
    h <- magick::image_info(ia3)$height
    wn <- h / 1.5
    crop_string <- paste0(wn, "x", h, "+", 0, "+", 0)
    hp <- magick::image_crop(ia3, crop_string)
    hp <- magick::image_colorize(hp, 100, "white")
    hp <- magick::image_annotate(hp, text, size = text_size, gravity = "center", color = "black", font = font)
    ia3 <- magick::image_append(c(hp, ia3))
  }
  
  if (head) {
    if (!is.null(text)) wo <- w
    w <- magick::image_info(ia3)$width
    h <- magick::image_info(ia3)$height
    if (!is.null(text)) wt <- w - wo
    hn <- h / 3
    w1 <- if(!is.null(text)) wo else w
    crop_string <- paste0(w1, "x", hn, "+", 0, "+", 0)
    hp <- magick::image_crop(ia3, crop_string)
    hp <- magick::image_colorize(hp, 100, "white")
    hp <- magick::image_annotate(hp, header, size = head_size, color = "black", gravity = "center", font = font)
    if (!is.null(text)) {
      crop_string2 <- paste0(wt, "x", hn, "+", 0, "+", 0)
      hp2 <- magick::image_crop(ia3, crop_string2)
      hp2 <- magick::image_colorize(hp2, 100, "white")
      if (!is.null(alt_text)) hp2 <- magick::image_annotate(hp2, alt_text, size = head_size, color = "black", gravity = "center", font = font)
      hp <- magick::image_append(c(hp2, hp))
    }

    ia3 <- magick::image_append(c(hp, ia3), stack = TRUE)
  }
  
  if (cs) {
    if (head) {
      w <- magick::image_info(nn2[3])$width
      crop_string <- paste0(w, "x", hn, "+", 0, "+", 0)
      hp <- magick::image_crop(ia3, crop_string)
      hp <- magick::image_colorize(hp, 100, "white")
      nn2[3] <- magick::image_append(c(hp, nn2[3]), stack = TRUE)
    }
    ia3 <- magick::image_append(c(ia3, nn2[3]))
  }
  
  ia3
}

quick_snaps <- function(temp_mgh_file, hemi, dir_out, SUBJECTS_DIR, overlay_threshold, inverted = FALSE, overlay_method = "linear", zoom = 1, print_call_only = FALSE) {
  name <- tools::file_path_sans_ext(temp_mgh_file)
  snap_order <- c("lateral", "medial", "superior", "inferior")
  if (hemi == "rh") snap_order[1:2] <- snap_order[2:1]
  snap_names <- paste0(name, ".", snap_order, ".tiff")
  snap_cmd <- c("--viewport 3d",
                QDECR:::qsnap_zoom(zoom),
                "--colorscale", 
                QDECR:::qsnap(snap_names[1]),
                QDECR:::qsnap_a(180),
                QDECR:::qsnap(snap_names[2]),
                if (hemi == "lh") QDECR:::qsnap_a(180),
                QDECR:::qsnap_e(90),
                QDECR:::qsnap(snap_names[3]),
                QDECR:::qsnap_e(180),
                QDECR:::qsnap(snap_names[4]),
                "--quit")
  tfile <-  file.path(dir_out, "tmp_snapshot_qdecr.txt")
  write.table(snap_cmd, tfile, quote = F, row.names = F, col.names = F)
  
  cmdStr <- paste0("freeview --surface ", SUBJECTS_DIR, "/fsaverage/surf/", hemi, ".inflated:overlay=", temp_mgh_file, 
                   if(!is.null(overlay_threshold)){paste0(":overlay_threshold=", paste0(overlay_threshold, collapse = ","))},
                   ":overlay_method=", overlay_method,
                   if (inverted) {":overlay_color=inverse"},
                   " -cmd ", tfile)
  if (print_call_only) return(cmdStr)
  system(cmdStr)
}
