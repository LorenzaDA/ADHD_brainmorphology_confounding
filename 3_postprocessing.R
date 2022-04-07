library(cowplot)
library(magick)

main_dir <- "PATH/TO/YOUR/WORKING/DIRECTORY"
dir_figures <- file.path(main_dir, "figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures)

#### LOAD
source("postprocessing_main_genr.R")
source("postprocessing_main_abcd.R")
source("postprocessing_ksads_abcd.R")

genr <- readRDS(file.path(main_dir, "genr.rds"))
abcd <- readRDS(file.path(main_dir, "abcd.rds"))
ksads <- readRDS(file.path(main_dir, "abcd_ksads.rds"))

#### FIGURE 1

f1 <- magick::image_append(c(fx_abcd$fig_1$row1, fx_genr$fig_1$row2, fx_abcd$fig_1$row3, fx_genr$fig_1$row4), 
                           stack = TRUE)

magick::image_write(f1, 
                    file.path(dir_figures, "f1.jpeg"), 
                    format = "jpeg")

#### FIGURE 3

f3 <- magick::image_append(c(fx_abcd$fig_3$row1, fx_genr$fig_3$row2, fx_abcd$fig_3$row3, fx_genr$fig_3$row4), 
                           stack = TRUE)

magick::image_write(f3, 
                    file.path(dir_figures, "f3.jpeg"), 
                    format = "jpeg")

#### FIGURE 4

f4_a <- plot_grid(NULL, fx_abcd$fig_4$a, labels = c("A", ""), label_size = 18, label_x = 0, nrow = 2, rel_heights = c(0.1, 1))
f4_b <- plot_grid(NULL, fx_genr$fig_4$b, labels = c("B", ""), label_size = 18, label_x = 0, nrow = 2, rel_heights = c(0.1, 1))
f4_c <- plot_grid(NULL, fx_abcd$fig_4$c, labels = c("C", ""), label_size = 18, label_x = 0, nrow = 2, rel_heights = c(0.1, 1))
f4_d <- plot_grid(NULL, fx_genr$fig_4$d, labels = c("D", ""), label_size = 18, label_x = 0, nrow = 2, rel_heights = c(0.1, 1))

f4 <- plot_grid(f4_a, f4_b, f4_c, f4_d, ncol = 2)

png(file.path(dir_figures, "f4.png"), 
    res = 600,
    width = 4200,
    height = 4200,
    type = "cairo-png",
    family = "sans")

f4

dev.off()


#### FIGURE 6

f6 <- magick::image_append(c(fx_abcd$fig_6$row1, fx_genr$fig_6$row2, fx_abcd$fig_6$row3, fx_genr$fig_6$row4), 
                           stack = TRUE)

magick::image_write(f6, 
                    file.path(dir_figures, "f6.jpeg"), 
                    format = "jpeg")

#### FIGURE SUPPL 2

fs2 <- magick::image_append(c(fx_abcd$fig_s2$row1, fx_genr$fig_s2$row2, fx_abcd$fig_s2$row3, fx_genr$fig_s2$row4),
                            stack = TRUE)

magick::image_write(fs2, 
                    file.path(dir_figures, "sf2.jpeg"), 
                    format = "jpeg")

##############################################

#### KSADS MODEL 1

kf1 <- magick::image_append(c(fx_ksads$fig_1$row1, fx_ksads$fig_1$row3), 
                            stack = TRUE)

magick::image_write(kf1, 
                    file.path(dir_figures, "ksads_m1.jpeg"), 
                    format = "jpeg")

#### KSADS MODEL 1-3

kf3 <- magick::image_append(c(fx_ksads$fig_3$row1, fx_ksads$fig_3$row3), 
                            stack = TRUE)

magick::image_write(kf3, 
                    file.path(dir_figures, "ksads_m1_m3.jpeg"), 
                    format = "jpeg")


