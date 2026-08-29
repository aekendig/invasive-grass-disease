# load packages
library(RColorBrewer)

# theme
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.margin = margin(-0.3, 0, -0.1, -0.2, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 9,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 9,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 9, hjust = 0.5))

# colors
cb_pal <- palette.colors(n = 9)
# "#000000" "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7" "#999999"
spp_pal <- cb_pal[c(7, 6)] # Mv, Ev
# "#D55E00" "#0072B2" 
outcome_pal <- cb_pal[c(5, 7, 6, 8)] # coex, Mv, Ev, priority
# "#F0E442" "#D55E00" "#0072B2" "#CC79A7"
dis_pal <- cb_pal[c(4, 2)] # disease suppressed, ambient
# "#009E73" "#E69F00"
comp_pal <- c("#687042","#D3B56A") # density, litter (from ppt)
  

# col_pal_grp <- cb_pal[c(4, 6, 3)]
# names(col_pal_grp) <- c("*M. vimineum*", "first-year *E. virginicus*", "adult *E. virginicus*")
# 
# grey_pal <- brewer.pal(n = 9, "Greys")[6:8]
# coral_pal <- c("coral", "coral3", "coral4")
# 
# col_pal4 <- c("black", coral_pal[2], "royalblue", grey_pal[1])
# 
# dis_pal <- c("goldenrod", "#7C9B5B")

# shapes
shape_pal <- c(19, 15, 17, 5)
spp_shape_pal <- shape_pal[c(2, 3)]

# dodge size
dodge_width <- 0.5

# text size
text_size <- 3.5
