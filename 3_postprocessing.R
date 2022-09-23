library(cowplot)
library(data.table)
library(magick)

main_dir <- "PATH/TO/YOUR/WORKING/DIRECTORY"
out_dir <- file.path(main_dir, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir)

#### LOAD
source("postprocessing_main_genr.R")
source("postprocessing_main_abcd.R")
source("postprocessing_ksads_abcd.R")

genr <- readRDS(file.path(main_dir, "genr.rds"))
abcd <- readRDS(file.path(main_dir, "abcd.rds"))
ksads <- readRDS(file.path(main_dir, "abcd_ksads.rds"))

#### NUMBERS FOR RESULTS SECTION
res_genr <- fread(file.path(main_dir, "regions_extracted_from_genr.txt"))
res_abcd <- fread(file.path(main_dir, "regions_extracted_from_abcd.txt"))
res_ksads <- fread(file.path(main_dir, "regions_extracted_from_ksads.txt"))

outs <- c("area", "volume", "thickness")
res <- lapply(outs, function(x) {
  data.frame(genr = sapply(1:5, function(i) sum(subset(res_genr, model == i & outcome == x)$mm2)),
                     abcd = sapply(1:5, function(i) sum(subset(res_abcd, model == i & outcome == x)$mm2)),
                     ksads = sapply(1:5, function(i) sum(subset(res_ksads, model == i & outcome == x)$mm2))) / 100
})

res_diff <- lapply(res, function(x) {
  y <- sapply(x, function(k) k / c(0, head(k, -1)))
  y[5, ] <- as.matrix(x[5, ] / x[3, ])
  y * 100 - 100
})

#### FIGURE 1

f1 <- magick::image_append(c(fx_abcd$fig_1$row1, fx_genr$fig_1$row2), 
                            stack = TRUE)
                            
magick::image_write(f1, 
                    file.path(out_dir, "f1.jpeg"), 
                    format = "jpeg")

#### FIGURE 3

f3 <- magick::image_append(c(fx_abcd$fig_3$row1, fx_genr$fig_3$row2), 
                            stack = TRUE)
                            
magick::image_write(f3, 
                    file.path(out_dir, "f3.jpeg"), 
                    format = "jpeg")

#### FIGURE 4

f4 <- plot_grid(NULL,
          fx_abcd$fig_4a,
          NULL,
          fx_genr$fig_4b,
          labels = c("", "ABCD", "", "GENR"),
          hjust = 0,
          vjust = -0.5, 
          nrow = 4,
          ncol = 1,
          rel_heights = c(0.1, 1, 0.1, 1))

png(file.path(out_dir, "f4.png"), 
   type = "cairo-png",
   width = 7.5, 
   height = 11, 
   units = "in", 
   res = 600)
   
f4
        
dev.off()

#### FIGURE 6

f6 <- magick::image_append(c(fx_abcd$fig_6$row1, fx_genr$fig_6$row2), 
                            stack = TRUE)
                            
magick::image_write(f6, 
                    file.path(out_dir, "f6.jpeg"), 
                    format = "jpeg")

#### FIGURE 8

f8 <- magick::image_append(c(fx_abcd$fig_8$row1, fx_genr$fig_8$row2), 
                            stack = TRUE)
                            
magick::image_write(f8, 
                    file.path(out_dir, "f8.jpeg"), 
                    format = "jpeg")

##############################################

#### SUPPLEMENTARY FIGURE 2

sf2 <- magick::image_append(c(fx_abcd$sfig_2$row1, fx_genr$sfig_2$row2, fx_abcd$sfig_2$row3, fx_genr$sfig_2$row4), 
                            stack = TRUE)
                            
magick::image_write(sf2, 
                    file.path(out_dir, "sf2.jpeg"), 
                    format = "jpeg")

#### SUPPLEMENTARY FIGURE 3

sf3 <- magick::image_append(c(fx_abcd$sfig_3$row1, fx_genr$sfig_3$row2, fx_abcd$sfig_3$row3, fx_genr$sfig_3$row4), 
                            stack = TRUE)
                            
magick::image_write(sf3, 
                    file.path(out_dir, "sf3.jpeg"), 
                    format = "jpeg")
                    
#### SUPPLEMENTARY FIGURE 4

sf4 <- magick::image_append(c(fx_ksads$sfig_4$row1, fx_ksads$sfig_4$row2, fx_ksads$sfig_4$row3), 
                            stack = TRUE)
                            
magick::image_write(sf4, 
                    file.path(out_dir, "sf4.jpeg"), 
                    format = "jpeg")

#### SUPPLEMENTARY FIGURE 5

sf5 <- magick::image_append(c(fx_ksads$sfig_5$row1, fx_ksads$sfig_5$row2, fx_ksads$sfig_5$row3), 
                            stack = TRUE)
                            
magick::image_write(sf5, 
                    file.path(out_dir, "sf5.jpeg"), 
                    format = "jpeg")

#### SUPPLEMENTARY FIGURE 6

sf6 <- magick::image_append(c(fx_abcd$sfig_6$row1, fx_genr$sfig_6$row2, 
                              fx_abcd$sfig_6$row3, fx_genr$sfig_6$row4,
                              fx_abcd$sfig_6$row5, fx_genr$sfig_6$row6), 
                            stack = TRUE)
                            
magick::image_write(sf6, 
                    file.path(out_dir, "sf6.jpeg"), 
                    format = "jpeg")

#### SUPPLEMENTARY FIGURE 7

sf7 <- plot_grid(fx_abcd$sfig_7a, 
          fx_abcd$sfig_7b, 
          fx_genr$sfig_7c, 
          fx_genr$sfig_7d, 
          labels = LETTERS[1:4],
          nrow = 2,
          ncol = 2)

png(file.path(out_dir, "sf7.png"), 
   type = "cairo-png",
   width = 7, 
   height = 7, 
   units = "in", 
   res = 600)

sf7
          
dev.off()

#### SUPPLEMENTARY FIGURE 8

sf8 <- plot_grid(NULL,
          fx_abcd$sfig_8a,
          NULL,
          fx_genr$sfig_8b,
          labels = c("", "ABCD", "", "GENR"),
          hjust = 0,
          vjust = -0.5, 
          nrow = 4,
          ncol = 1,
          rel_heights = c(0.1, 1, 0.1, 1))

png(file.path(out_dir, "sf8.png"), 
   type = "cairo-png",
   width = 7.5, 
   height = 11, 
   units = "in", 
   res = 600)
   
sf8
        
dev.off()

#### SUPPLEMENTARY FIGURE 9

sf9 <- plot_grid(fx_abcd$sfig_9a, 
          fx_abcd$sfig_9b, 
          fx_genr$sfig_9c, 
          fx_genr$sfig_9d, 
          labels = LETTERS[1:4],
          nrow = 2,
          ncol = 2)

png(file.path(out_dir, "sf9.png"), 
   type = "cairo-png",
   width = 7, 
   height = 7, 
   units = "in", 
   res = 600)

sf9
          
dev.off()

#### SUPPLEMENTARY FIGURE 10

  ii <- list(c(1, 4), c(2, 5))
  lims <- list(c(-0.005, 0), c(-0.01, 0))
  m2 <- 1:3
  models <- paste0("model_", m2)
  outs <- c("Surface area", "Volume")

  # individual
  
  ss_list <- lapply(seq_along(models), function(j) {
  
    mm <- models[j]
    mmc <- paste0("median", j)
    lapply(seq_along(outs), function(i) {
      
      iis <- ii[[i]]
      ss_title <- paste0(list(LETTERS[1:3], LETTERS[4:6])[[i]][j], ": ", outs[i], " (model ", j, ")")
      
      data1 <- abcd$roi_list[[iis[1]]][[3]][, mmc][-35]
      data2 <- genr$roi_list[[iis[1]]][[3]][, mmc][-35]
      
      d_ss <- data.frame(abcd = data1, genr = data2)
      r_ss <- cor(d_ss, method = "spearman")[1, 2]
      r_text <- paste0("r = ", round(r_ss, 2))
      
      ggplot2::ggplot(d_ss, ggplot2::aes(x = abcd, y = genr)) +
        ggplot2::geom_point() +
        cowplot::theme_cowplot() +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
        ggplot2::xlab("Average beta in ABCD") +
        ggplot2::ylab("Average beta in GenR") +
        ggplot2::xlim(lims[[i]]) +
        ggplot2::ylim(lims[[i]]) + 
        ggplot2::annotation_custom(grid::grobTree(grid::textGrob(r_text, x = 0.3,  y = 0.8, gp = grid::gpar(fontsize = 17)))) + 
        ggplot2::ggtitle(ss_title)
      
    })
  })
  
  # compile
  ss_u <- unlist(do.call(function(...) Map(list,...), ss_list), recursive = FALSE)
  sf10 <- do.call(cowplot::plot_grid, c(ss_u, ncol = length(m2), nrow = 2))
  
  png(file.path(out_dir, "sf10.png"), type = "cairo-png", width = 16, height = 10, units = "in", res = 600)
  sf10
  dev.off()

#### SUPPLEMENTARY FIGURE 11

sf11 <- magick::image_append(c(fx_abcd$sfig_11$row1, fx_genr$sfig_11$row2, fx_abcd$sfig_11$row3, fx_genr$sfig_11$row4), 
                            stack = TRUE)
                            
magick::image_write(sf11, 
                    file.path(out_dir, "sf11.jpeg"), 
                    format = "jpeg")

#### SUPPLEMENTARY FIGURE 12

sf12 <- magick::image_append(c(fx_abcd$sfig_12$row1, fx_genr$sfig_12$row2, fx_abcd$sfig_12$row3, fx_genr$sfig_12$row4), 
                            stack = TRUE)
                            
magick::image_write(sf12, 
                    file.path(out_dir, "sf12.jpeg"), 
                    format = "jpeg")
