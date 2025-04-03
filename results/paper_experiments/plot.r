#!/usr/bin/env Rscript
library(ggplot2)
library(stringr)
library(dplyr)
library(scales)
library(ggh4x)
library(patchwork)
library(tidyverse)

# change page format (width, height, etc.) in layout.r
source("layout.r")

# for tikz export
library(tikzDevice)
options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
                              "\\usepackage{amsmath}"))
create_save_locations()

d_qaoa <- read.csv("csv-data/qaoa_data.csv", stringsAsFactors=FALSE)
d_maxcut <- d_qaoa %>% filter(problem == "maxcut")
d_vertcover <- d_qaoa %>% filter(problem == "vertcover")
d_max3sat <- d_qaoa %>% filter(problem == "max3sat")

gamma_breaks <- function(x) { if (max(x) < 3) c(-pi/2, 0, pi/2) else c(-pi, 0, pi) }
beta_breaks <- function(x) { if (max(x) < 1.5) c(-pi/4, 0, pi/4) else if (max(x) < 3) c(-pi/2, 0, pi/2) else c(-pi, 0, pi) }
gamma_labels <- function(x) { if (max(x) < 3) c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$") else c("$-\\pi$", "0", "$\\pi$") }
beta_labels <- function(x) { if (max(x) < 1.5) c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$") else if (max(x) < 3) c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$") else c("$-\\pi$", "0", "$\\pi$")}

p_labels <- function(p) {
    paste0("p = ", p)
}

method_labels <- c(
                    `symmetry` = "Restricted",
                    `no_symmetry` = "Unbounded",
                    `sequential_average` = "Sequential",
                    `sequential_standard_deviation` = "Sequential",
                    `sequential`= "Sequential",
                    `optimised_average` = "COBYLA",
                    `optimised` = "COBYLA",
                    `optimised_standard_deviation` = "COBYLA",
                    #`fixed_average` = "Fixed",
                    #`fixed_standard_deviation` = "Fixed",
                    #`fixed` = "Fixed",
                    `linear_ramp_average` = "LR$_{+\\beta}$",
                    `linear_ramp_standard_deviation` = "LR$_{+\\beta}$",
                    `linear_ramp` = "LR$_{+\\beta}$",
                    `linear_ramp_neg_mixer_average` = "LR$_{-\\beta}$",
                    `linear_ramp_neg_mixer_standard_deviation` = "LR$_{-\\beta}$",
                    `linear_ramp_neg_mixer` = "LR$_{-\\beta}$"
)

problem_labels <- c(
  `maxcut` = "MaxCut$_{p=7}$",
  `max3sat` = "Max3SAT",
  `vertcover` = "Vertexcover",
  `higher_depth_maxcut` = "MaxCut$_{p=21}$"
)

param_labels <- c(
  `gamma` = "$\\gamma$",
  `beta` = "$\\beta$"
)

method_colors <- c(
  `sequential_average` =  "#1f78b4",
  `sequential_standard_deviation` =  "#1f78b4",
  `sequential`=  "#1f78b4",
  `optimised_average` =  "#009371",
  `optimised_standard_deviation` =  "#009371",
  `optimised` =  "#009371",
  #`fixed_average` =  "#beaed4",
  #`fixed_standard_deviation` = "#beaed4",
  #`fixed` = "#beaed4",
  `linear_ramp_average` = "#E69F00",
  `linear_ramp_standard_deviation` = "#E69F00",
  `linear_ramp` = "#E69F00",
  `linear_ramp_neg_mixer_average` =  "#1b0a33",
  `linear_ramp_neg_mixer_standard_deviation` =  "#1b0a33",
  `linear_ramp_neg_mixer` =  "#1b0a33"
)

method_shapes <- c(
  `sequential_average` =  21,
  `sequential`=  21,
  `sequential_standard_deviation` =  21,
  `optimised_average` =  22,
  `optimised_standard_deviation` =  22,
  `optimised` =  22,
  #`fixed_average` =  23,
  #`fixed_standard_deviation` = 23,
  #`fixed` = 23,
  `linear_ramp_average` = 24,
  `linear_ramp_standard_deviation` = 24,
  `linear_ramp` = 24,
  `linear_ramp_neg_mixer_average` = 25,
  `linear_ramp_neg_mixer_standard_deviation` = 25,
  `linear_ramp_neg_mixer` = 25
)

######### SYMMETRIES ###################
plot_max_depth <- 3
#### Plot Maxcut Symmetries
d_sym_no_sym_maxcut <- d_maxcut %>% filter(p <= plot_max_depth) %>% filter(method == "no_symmetry" | method == "symmetry")
d_no_sym_maxcut <- d_maxcut %>% filter(p <= plot_max_depth) %>% filter(method == "no_symmetry")

g_maxcut_sym <- ggplot(d_sym_no_sym_maxcut, aes(x = beta, y = gamma, fill = expectation)) +
    theme_paper_legend_right() +
    geom_raster() +
    geom_hline(d_no_sym_maxcut, mapping=aes(yintercept=3*pi/4), colour = COLOURS.LIST[2]) +
    geom_hline(d_no_sym_maxcut, mapping=aes(yintercept=pi/4), colour = COLOURS.LIST[2]) +
    geom_hline(d_no_sym_maxcut, mapping=aes(yintercept=-3*pi/4), colour = COLOURS.LIST[2]) +
    geom_hline(d_no_sym_maxcut, mapping=aes(yintercept=-pi/4), colour = COLOURS.LIST[2]) +
    geom_vline(d_no_sym_maxcut, mapping=aes(xintercept=pi/2), colour = COLOURS.LIST[2]) +
    geom_vline(d_no_sym_maxcut, mapping=aes(xintercept=-pi/2), colour = COLOURS.LIST[2]) +
    facet_grid2(method ~ p, scales="free", independent = "x",
                labeller = labeller(method = as_labeller(method_labels),
                                    p = as_labeller(p_labels))) +
    labs(x="$\\gamma$", y="$\\beta$",title="MaxCut Symmetries") +
    scale_x_continuous(breaks = gamma_breaks, labels = gamma_labels) +
    scale_y_continuous(breaks = beta_breaks, labels = beta_labels) +
    scale_fill_gradientn(
        "$F_p(\\gamma, \\beta)$",
        colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
        breaks = c(-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0),
        guide = guide_colourbar(theme = theme(
            legend.key.width  = unit(0.4, "lines"),
            legend.key.height = unit(10, "lines"))),
    )

#### Plot general Symmetries
d_sym_no_sym_3sat <- d_max3sat %>% filter(p <= plot_max_depth) %>% filter(method == "no_symmetry" | method == "symmetry")
d_no_sym_3sat <- d_max3sat %>% filter(p <= plot_max_depth) %>% filter(method == "no_symmetry")

g_gen_sym <- ggplot(d_sym_no_sym_3sat, aes(x = beta, y = gamma, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  geom_hline(d_no_sym_3sat, mapping=aes(yintercept=pi/2), colour = COLOURS.LIST[2]) +
  geom_hline(d_no_sym_3sat, mapping=aes(yintercept=-pi/2), colour = COLOURS.LIST[2]) +
  facet_grid2(method ~ p, scales="free", independent = "x",
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$",title="General Symmetries") +
  scale_x_continuous(breaks = gamma_breaks, labels = gamma_labels) +
  scale_y_continuous(breaks = beta_breaks, labels = beta_labels) +
  scale_fill_gradientn(
    "$F_p(\\gamma, \\beta)$",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(10, "lines"))),
  )

#combine into one plot

g <- (g_maxcut_sym | g_gen_sym) & theme(legend.position = 'right')

save_name <- "symmetries"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.6*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.6*COLWIDTH)
print(g)
dev.off()

######### MAIN RESULTS ###################

#### Plot main results maxcut
d_maxcut_results <- d_maxcut %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")
d_maxcut_results$method <- factor(d_maxcut_results$method, levels=c("sequential","optimised", "linear_ramp_neg_mixer", "linear_ramp"))

g <- ggplot(d_maxcut_results, aes(x = beta, y = gamma, fill = expectation)) +
    theme_paper_legend_right() +
    geom_raster() +
    facet_grid2(method ~ p,
                labeller = labeller(method = as_labeller(method_labels),
                                    p = as_labeller(p_labels))) +
    labs(x="$\\gamma$", y="$\\beta$") +
    scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
    scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
    scale_fill_gradientn(
        "$F_p(\\gamma, \\beta)$",
        colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
        breaks = c(-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0),
        guide = guide_colourbar(theme = theme(
            legend.key.width  = unit(0.4, "lines"),
            legend.key.height = unit(10, "lines"))),
    )

save_name <- "maxcut_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot main results max3sat
d_max3sat_results <- d_max3sat %>% filter(additional_info == "hard") %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")
d_max3sat_results$method <- factor(d_max3sat_results$method, levels=c("sequential","optimised", "linear_ramp_neg_mixer", "linear_ramp"))

g <- ggplot(d_max3sat_results, aes(x = beta, y = gamma, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$F_p(\\gamma, \\beta)$",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(10, "lines"))),
  )

save_name <- "max3sat_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot main results vertexcover
d_vertcover_results <- d_vertcover %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")
d_vertcover_results$method <- factor(d_vertcover_results$method, levels=c("sequential","optimised", "linear_ramp_neg_mixer", "linear_ramp"))

g <- ggplot(d_vertcover_results, aes(x = beta, y = gamma, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$F_p(\\gamma, \\beta)$",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0,10,50,100,150,200,250,300,350,400),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(10, "lines"))),
  )

save_name <- "vertcover_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

######### AVERAGE AND STANDARD DEVIATION ###################

#### Plot average and standard deviation results for maxcut
## Average
d_maxcut_avg <- d_maxcut %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_maxcut_avg$method <- factor(d_maxcut_avg$method, levels=c("sequential_average","optimised_average","linear_ramp_neg_mixer_average", "linear_ramp_average"))

g <- ggplot(d_maxcut_avg, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
  scale_fill_gradientn(
    "$ \\bar{r} $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )
g
save_name <- "avg_maxcut_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

### Standard Deviation
d_maxcut_std <- d_maxcut %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_maxcut_std$method <- factor(d_maxcut_std$method, levels=c("sequential_standard_deviation","optimised_standard_deviation", "linear_ramp_neg_mixer_standard_deviation", "linear_ramp_standard_deviation"))

g <- ggplot(d_maxcut_std, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
  scale_fill_gradientn(
    "$\\sigma_r$",
    colours = c(COLOURS.LIST[1], COLOURS.LIST[2], "white"),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "std_maxcut_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot average and standard deviation results for max3sat (hard instances)
## Average
d_max3sat_avg <- d_max3sat %>% filter(additional_info == "hard") %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_max3sat_avg$method <- factor(d_max3sat_avg$method, levels=c("sequential_average","optimised_average","linear_ramp_neg_mixer_average","linear_ramp_average"))

g <- ggplot(d_max3sat_avg, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$ \\bar{r} $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "avg_max3sat_results_hard"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

### Standard Deviation
d_max3sat_std <- d_max3sat %>% filter(additional_info == "hard") %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_max3sat_std$method <- factor(d_max3sat_std$method, levels=c("sequential_standard_deviation","optimised_standard_deviation", "linear_ramp_neg_mixer_standard_deviation", "linear_ramp_standard_deviation"))

g <- ggplot(d_max3sat_std, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$\\sigma_r$",
    colours = c(COLOURS.LIST[1], COLOURS.LIST[2], "white"),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "std_max3sat_results_hard"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot average and standard deviation results for max3sat (easy instances)
## Average
d_max3sat_avg <- d_max3sat %>% filter(additional_info == "easy") %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_max3sat_avg$method <- factor(d_max3sat_avg$method, levels=c("sequential_average","optimised_average","linear_ramp_neg_mixer_average","linear_ramp_average"))

g <- ggplot(d_max3sat_avg, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$ \\bar{r} $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "avg_max3sat_results_easy"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

### Standard Deviation
d_max3sat_std <- d_max3sat %>% filter(additional_info == "easy") %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_max3sat_std$method <- factor(d_max3sat_std$method, levels=c("sequential_standard_deviation","optimised_standard_deviation", "linear_ramp_neg_mixer_standard_deviation", "linear_ramp_standard_deviation"))

g <- ggplot(d_max3sat_std, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$\\sigma_r$",
    colours = c(COLOURS.LIST[1], COLOURS.LIST[2], "white"),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "std_max3sat_results_easy"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot average and standard deviation results for vertexcover
## Average
d_vertcover_avg <- d_vertcover %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_vertcover_avg$method <- factor(d_vertcover_avg$method, levels=c("sequential_average","optimised_average","linear_ramp_neg_mixer_average","linear_ramp_average"))

g <- ggplot(d_vertcover_avg, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$ \\bar{r} $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "avg_vertcover_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

### Standard Deviation
d_vertcover_std <- d_vertcover %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_vertcover_std$method <- factor(d_vertcover_std$method, levels=c("sequential_standard_deviation","optimised_standard_deviation", "linear_ramp_neg_mixer_standard_deviation", "linear_ramp_standard_deviation"))

g <- ggplot(d_vertcover_std, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi, 0, pi), labels = c("$-\\pi$", "0", "$\\pi$")) +
  scale_y_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_fill_gradientn(
    "$\\sigma_r$",
    colours = c(COLOURS.LIST[1], COLOURS.LIST[2], "white"),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "std_vertcover_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

############### PLOT PARAMETERS AND APPROX RATIO BY METHOD #############################
d_params <- read.csv("csv-data/param_data.csv", stringsAsFactors=FALSE)
d_maxcut_params <- d_params %>% filter(problem == "maxcut")
d_vertcover_params <- d_params %>% filter(problem == "vertcover")
d_max3sat_params <- d_params %>% filter(problem == "max3sat")

######## Maxcut main ############################################
d_maxcut_param_results <- d_maxcut_params %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")

g_approx <- ggplot(d_maxcut_param_results, aes(x = p, y = approx, color = method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\text{Residual Energy} r $") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))  + geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") + annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)
g_gamma <- ggplot(d_maxcut_param_results, aes(x = p, y = beta, color = method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\gamma$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))
g_beta <- ggplot(d_maxcut_param_results, aes(x = p, y = gamma, color = method, shape=method))+ theme_paper_legend_right() + geom_line() + ylab("$\\beta$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))

g <- (g_beta | g_gamma) / g_approx +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

save_name <- "params_maxcut"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

######## Max3sat main ############################################
d_max3sat_param_results <- d_max3sat_params %>% filter(additional_info=="hard") %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")

g_approx <- ggplot(d_max3sat_param_results, aes(x = p, y = approx, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\text{Residual Energy} r $")  + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) + geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") + annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)
g_gamma <- ggplot(d_max3sat_param_results, aes(x = p, y = beta, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\gamma$")  + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))
g_beta <- ggplot(d_max3sat_param_results, aes(x = p, y = gamma, color=method, shape=method))+ theme_paper_legend_right() + geom_line() + ylab("$\\beta$")  + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))

g <- (g_beta | g_gamma) / g_approx +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

save_name <- "params_max3sat"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

######## Vertcover main ############################################
d_vertcover_param_results <- d_vertcover_params %>% filter(method == "sequential" | method == "linear_ramp" | method == "linear_ramp_neg_mixer" | method == "optimised")

g_approx <- ggplot(d_vertcover_param_results, aes(x = p, y = approx, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\text{Residual Energy} r $")  + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) + geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") + annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)
g_gamma <- ggplot(d_vertcover_param_results, aes(x = p, y = beta, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\gamma$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))
g_beta <- ggplot(d_vertcover_param_results, aes(x = p, y = gamma, color=method, shape=method))+ theme_paper_legend_right() + geom_line() + ylab("$\\beta$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))

g <- (g_beta | g_gamma) / g_approx +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')


save_name <- "params_vertcover"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

######## Parameter quality Average and standard deviation #########################################

############ Maxcut avg and stdev ######################
d_maxcut_param_avg <- d_maxcut_params %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_maxcut_param_avg <- data.frame(lapply(d_maxcut_param_avg, function(x) { gsub("_average", "", x) }))
d_maxcut_param_std <- d_maxcut_params %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_maxcut_param_std <- data.frame(lapply(d_maxcut_param_std, function(x) { gsub("_standard_deviation", "", x) }))
d_maxcut_stats <- d_maxcut_param_avg %>% inner_join(d_maxcut_param_std, by=c("method","p"), suffix=c("_avg","_std"))
d_maxcut_stats <- transform(d_maxcut_stats, p=as.numeric(p), approx_avg = as.numeric(approx_avg), approx_std = as.numeric(approx_std))

g <- ggplot(d_maxcut_stats, aes(x = p, y = approx_avg, color = method, shape=method)) +
  geom_line() +
  scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) +
  geom_errorbar(aes(ymin=(approx_avg-approx_std),ymax=(approx_avg+approx_std),color=method), alpha=0.5, width=0.25) +
  ylab("$\\text{Average Residual Energy } \\bar{r} $") +
  theme_paper_legend_right() +
  scale_fill_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + 
  geom_point() +
  scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) +
  geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") +
  annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)
g
save_name <- "params_maxcut_stats"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()

############ Max3Sat avg and stdev (easy) ######################
d_max3sat_easy <- d_max3sat_params %>% filter(additional_info == "easy")
d_max3sat_param_avg <- d_max3sat_easy %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_max3sat_param_avg <- data.frame(lapply(d_max3sat_param_avg, function(x) { gsub("_average", "", x) }))
d_max3sat_param_std <- d_max3sat_easy %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_max3sat_param_std <- data.frame(lapply(d_max3sat_param_std, function(x) { gsub("_standard_deviation", "", x) }))
d_max3sat_stats <- d_max3sat_param_avg %>% inner_join(d_max3sat_param_std, by=c("method","p"), suffix=c("_avg","_std"))
d_max3sat_stats <- transform(d_max3sat_stats, p=as.numeric(p), approx_avg = as.numeric(approx_avg), approx_std = as.numeric(approx_std))

g <- ggplot(d_max3sat_stats, aes(x = p, y = approx_avg, color = method, shape=method)) +
  geom_line() +
  scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) +
  geom_errorbar(aes(ymin=(approx_avg-approx_std),ymax=(approx_avg+approx_std),color=method), alpha=0.5, width=0.25) +
  ylab("$\\text{Average Residual Energy } \\bar{r} $") +
  theme_paper_legend_right() +
  scale_fill_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + 
  geom_point() +
  scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) +
  geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") +
  annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)

save_name <- "params_max3sat_stats_easy"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()

############ Max3Sat avg and stdev (hard) ######################
d_max3sat_hard <- d_max3sat_params %>% filter(additional_info == "hard")
d_max3sat_param_avg <- d_max3sat_hard %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_max3sat_param_avg <- data.frame(lapply(d_max3sat_param_avg, function(x) { gsub("_average", "", x) }))
d_max3sat_param_std <- d_max3sat_hard %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_max3sat_param_std <- data.frame(lapply(d_max3sat_param_std, function(x) { gsub("_standard_deviation", "", x) }))
d_max3sat_stats <- d_max3sat_param_avg %>% inner_join(d_max3sat_param_std, by=c("method","p"), suffix=c("_avg","_std"))
d_max3sat_stats <- transform(d_max3sat_stats, p=as.numeric(p), approx_avg = as.numeric(approx_avg), approx_std = as.numeric(approx_std))

g <- ggplot(d_max3sat_stats, aes(x = p, y = approx_avg, color = method, shape = method)) +
  geom_line() +
  scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) +
  geom_errorbar(aes(ymin=(approx_avg-approx_std),ymax=(approx_avg+approx_std),color=method), alpha=0.5, width=0.25) +
  ylab("$\\text{Average Residual Energy } \\bar{r} $") +
  theme_paper_legend_right() +
  scale_fill_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + 
  geom_point() +
  scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) +
  geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") +
  annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)

save_name <- "params_max3sat_stats_hard"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
g

############ Vertcover avg and stdev ######################
d_vertcover_param_avg <- d_vertcover_params %>% filter(method == "sequential_average" | method == "linear_ramp_average" | method == "linear_ramp_neg_mixer_average" | method == "optimised_average")
d_vertcover_param_avg <- data.frame(lapply(d_vertcover_param_avg, function(x) { gsub("_average", "", x) }))
d_vertcover_param_std <- d_vertcover_params %>% filter(method == "sequential_standard_deviation" | method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" |  method == "optimised_standard_deviation")
d_vertcover_param_std <- data.frame(lapply(d_vertcover_param_std, function(x) { gsub("_standard_deviation", "", x) }))
d_vertcover_stats <- d_vertcover_param_avg %>% inner_join(d_vertcover_param_std, by=c("method","p"), suffix=c("_avg","_std"))
d_vertcover_stats <- transform(d_vertcover_stats, p=as.numeric(p), approx_avg = as.numeric(approx_avg), approx_std = as.numeric(approx_std))

g <- ggplot(d_vertcover_stats, aes(x = p, y = approx_avg, color = method, shape=method)) +
  geom_line() +
  scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) +
  geom_errorbar(aes(ymin=(approx_avg-approx_std),ymax=(approx_avg+approx_std),color=method), alpha=0.5, width=0.25) +
  ylab("$\\text{Average Residual Energy } \\bar{r} $") +
  theme_paper_legend_right() +
  scale_fill_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + 
  geom_point() +
  scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) +
  geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") +
  annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)

save_name <- "params_vertcover_stats"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()

################# higher depth experiments ###############################################

d_higher_depth <- d_qaoa %>% filter(problem == "higher_depth_maxcut")
d_higher_depth_params <- d_params %>% filter(problem == "higher_depth_maxcut")

# higher depth maxcut

d_21maxcut <- d_higher_depth %>% filter( method == "linear_ramp_neg_mixer" | method == "linear_ramp" | method == "optimised" | method == "sequential")
d_21maxcut$method <- factor(d_21maxcut$method, levels=c("sequential","optimised", "linear_ramp_neg_mixer", "linear_ramp"))

g <- ggplot(d_21maxcut, aes(x = beta, y = gamma, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
  scale_fill_gradientn(
    "$F_p(\\gamma, \\beta)$",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(10, "lines"))),
  )
g
save_name <- "higher_p_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off() 
 
# higher depth maxcut params 
d_21maxcut_params <- d_higher_depth_params %>% filter(method == "linear_ramp_neg_mixer" | method == "linear_ramp" | method == "optimised" | method=="sequential")

g_approx <- ggplot(d_21maxcut_params, aes(x = p, y = approx, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\text{Residual Energy} r $") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) + geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") + annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)
g_gamma <- ggplot(d_21maxcut_params, aes(x = p, y = beta, color=method, shape=method)) + theme_paper_legend_right() + geom_line() + ylab("$\\gamma$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))
g_beta <- ggplot(d_21maxcut_params, aes(x = p, y = gamma, color=method, shape=method))+ theme_paper_legend_right() + geom_line() + ylab("$\\beta$") + scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + geom_point() + scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels))

g <- (g_beta | g_gamma) / g_approx +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

save_name <- "params_higher_p"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

#### Plot average and standard deviation results for maxcut
## Average
d_21maxcut_avg <- d_higher_depth %>% filter(method == "linear_ramp_neg_mixer_average"| method == "linear_ramp_average" | method == "optimised_average" | method=="sequential_average")
d_21maxcut_avg$method <- factor(d_21maxcut_avg$method, levels=c("sequential_average","optimised_average","linear_ramp_neg_mixer_average","linear_ramp_average"))

g <- ggplot(d_21maxcut_avg, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
  scale_fill_gradientn(
    "$ \\bar{r} $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "avg_higher_p_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

### Standard Deviation
d_21maxcut_std <- d_higher_depth %>% filter(method == "linear_ramp_neg_mixer_standard_deviation" |method == "linear_ramp_standard_deviation" | method == "optimised_standard_deviation" | method=="sequential_standard_deviation")
d_21maxcut_std$method <- factor(d_21maxcut_std$method, levels=c("sequential_standard_deviation","optimised_standard_deviation", "linear_ramp_neg_mixer_standard_deviation", "linear_ramp_standard_deviation"))

g <- ggplot(d_21maxcut_std, aes(x = gamma, y = beta, fill = expectation)) +
  theme_paper_legend_right() +
  geom_raster() +
  facet_grid2(method ~ p,
              labeller = labeller(method = as_labeller(method_labels),
                                  p = as_labeller(p_labels))) +
  labs(x="$\\gamma$", y="$\\beta$") +
  scale_x_continuous(breaks = c(-pi/2, 0, pi/2), labels = c("$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$")) +
  scale_y_continuous(breaks = c(-pi/4, 0, pi/4), labels = c("$-\\frac{\\pi}{4}$", "0", "$\\frac{\\pi}{4}$")) +
  scale_fill_gradientn(
    "$\\sigma_r$",
    colours = c(COLOURS.LIST[1], COLOURS.LIST[2], "white"),
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

save_name <- "std_higher_p_results"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.8*COLWIDTH)
print(g)
dev.off()

## higher depthMaxcut avg and stdev params

d_21maxcut_param_avg <- d_higher_depth_params %>% filter(method == "linear_ramp_neg_mixer_average" | method == "linear_ramp_average" | method == "optimised_average"| method=="sequential_average")
d_21maxcut_param_avg <- data.frame(lapply(d_21maxcut_param_avg, function(x) { gsub("_average", "", x) }))
d_21maxcut_param_std <- d_higher_depth_params %>% filter(method == "linear_ramp_neg_mixer_standard_deviation" | method == "linear_ramp_standard_deviation" | method == "optimised_standard_deviation" | method =="sequential_standard_deviation")
d_maxcut_param_std <- data.frame(lapply(d_21maxcut_param_std, function(x) { gsub("_standard_deviation", "", x) }))
d_21maxcut_stats <- d_21maxcut_param_avg %>% inner_join(d_maxcut_param_std, by=c("method","p"), suffix=c("_avg","_std"))
d_21maxcut_stats <- transform(d_21maxcut_stats, p=as.numeric(p), approx_avg = as.numeric(approx_avg), approx_std = as.numeric(approx_std))

g <- ggplot(d_21maxcut_stats, aes(x = p, y = approx_avg, color = method, shape=method)) +
  geom_line() +
  scale_color_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) +
  geom_errorbar(aes(ymin=(approx_avg-approx_std),ymax=(approx_avg+approx_std),color=method), alpha=0.5, width=0.25) +
  ylab("$\\text{Average Residual Energy } \\bar{r} $") +
  theme_paper_legend_right() +
  scale_fill_manual(name = "Method", values = method_colors, labels=as_labeller(method_labels)) + 
  geom_point() +
  scale_shape_manual(name = "Method", values= method_shapes,labels=as_labeller(method_labels)) +
  geom_hline(yintercept=0.0, linetype='dotted', color = "#ed665a") +
  annotate("text", x = 0.5, y = 0.0, size=2, color="#ed665a", label = "optimal", vjust = -1)

save_name <- "params_higher_p_stats"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = 0.5*COLWIDTH)
print(g)
dev.off()

########ALL PARAMS OF OPTIMISED AND sequential #####################
d_allparams <- read.csv("csv-data/all_param_data.csv", stringsAsFactors=FALSE)
names(d_allparams)[5:6] <- c("beta","gamma")
d_allparams$problem <- factor(d_allparams$problem, levels=c("maxcut","higher_depth_maxcut","vertcover","max3sat"))
d_sequential_params <- d_allparams %>% filter(method=="sequential") %>% pivot_longer(cols=gamma:beta,names_to="param",values_to="value")
d_cobyla_params <- d_allparams %>% filter(method=="optimised") %>% pivot_longer(cols=gamma:beta,names_to="param",values_to="value")


g_sequential <- ggplot(d_sequential_params, aes(x = p, y = value, color=num_qubits, group = additional_info)) +
  theme_paper_legend_right() +
  geom_line() +
  facet_grid2(problem ~ param,
              labeller = labeller(problem = as_labeller(problem_labels),
                                  param = as_labeller(param_labels)),
             scales="free",
             independent = "x") +
  labs(title="Sequential",x="depth $p$", y="value") +
  scale_y_continuous(breaks = c(-pi,-pi/2, 0, pi/2, pi), labels = c("$-\\pi$","$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$", "$\\pi$")) +
  scale_color_gradientn(
    "$ n $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(10,15,20,24),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

g_cobyla <- ggplot(d_cobyla_params, aes(x = p, y = value, color=num_qubits, group = additional_info)) +
  theme_paper_legend_right() +
  geom_line() +
  facet_grid2(problem ~ param,
              labeller = labeller(problem = as_labeller(problem_labels),
                                  param = as_labeller(param_labels)),
              scales="free",
              independent = "x") +
  labs(title="COBYLA",x="depth $p$",y="value") +
  scale_y_continuous(breaks = c(-pi,-pi/2, 0, pi/2, pi), labels = c("$-\\pi$","$-\\frac{\\pi}{2}$", "0", "$\\frac{\\pi}{2}$", "$\\pi$")) +
  scale_color_gradientn(
    "$ n $",
    colours = c("#00ffbf", "#1f78b4", "#1b0a33"),
    breaks = c(10,15,20,24),
    guide = guide_colourbar(theme = theme(
      legend.key.width  = unit(0.4, "lines"),
      legend.key.height = unit(3, "lines"))),
  )

g_cobyla$labels$y <- ""

g <- (g_sequential | g_cobyla) +
  plot_layout(guides = "collect",axis_title ="collect") & theme(legend.position = 'right')
g

save_name <- "parameter_arragement"
pdf(str_c(OUTDIR_PDF, save_name, ".pdf"), width = TEXTWIDTH, height = COLWIDTH)
print(g)
dev.off()
tikz(str_c(OUTDIR_TIKZ, save_name, ".tex"), width = TEXTWIDTH, height = COLWIDTH)
print(g)
dev.off()

