suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(rlang) # Necessary for processing command args
})


# Read Command Args ####
if(!exists("argv")) argv = commandArgs(TRUE)
in_file  = argv[1]   %|% "lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.rds"
out_file = argv[2]  %|%  str_replace(in_file, ".rds", ".png")
lower_bound = argv[3] %>% as.integer() %|% 30L
upper_bound = argv[4] %>% as.integer() %|% 200L

density_colors = c("black", #"#480758",  # Both bad
                   # "#30bb50", # Goldilocks
                   "#007BA7", # Good
                   "yellow" # Too high
)

theme_set(theme_classic())

### Make a plot ####
count_data = read_rds(in_file)

tapestry_plot = count_data %>% 
  mutate(density_class = cut(density, 
                             c(0, lower_bound - 1, upper_bound, Inf),
                             labels = c(paste0("<", lower_bound),
                                        paste(lower_bound, upper_bound, sep = " - "),
                                        paste0(upper_bound, "+")),
                             ordered_result = TRUE) )%>% 
  ggplot() + 
  aes(y = -column, x = y, fill = density_class) + 
  geom_raster() + 
  # scale_fill_viridis_d("Read Count Density") + 
  scale_fill_manual("Read Count Density", values = density_colors, 
                    na.translate = FALSE) + 
  # facet_wrap(~column, nrow = 1) + 
  theme(strip.text = element_blank(), strip.background = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank(),
        panel.spacing = unit(1,"mm"))
  # Output the figure
ggsave(out_file, tapestry_plot,
       width = 16, height = 11, dpi = 300)


## Some old plot code:
## 
# dens_range = range(count_data$density)

# squeeze = function(data, min = 1, max = 500) {
#   pmin(data, max) %>% pmax(min)
# }

# figure = str_replace(in_file, "wig", "png")
# figure2 = str_replace(in_file, ".wig", "_extreme.png")
# 
# tst_plot = ggplot(out_data) + 
#   aes(y = -column, x = y, fill = squeeze(density)) + 
#   geom_raster() + 
#   scale_fill_viridis_c("Read Count Density",
#                        trans = "log10",  
#                        # option = 'magma',
#                        breaks = c(1, 10, 100, 500),
#                        labels = c("<1", "10", "100", "500+")) + 
#   # facet_wrap(~column, nrow = 1) + 
#   theme(strip.text = element_blank(), strip.background = element_blank(),
#         axis.text = element_blank(), 
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = "bottom",
#         axis.ticks = element_blank(),
#         panel.spacing = unit(1,"mm"))
# ggsave(file.path("bams/counts", figure), tst_plot,
#        width = 16, height = 11, dpi = 300)

