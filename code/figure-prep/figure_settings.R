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
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 9,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 9,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))

# colors
col_pal <- palette.colors(n = 6)
col_pal_grp <- col_pal[c(4, 6, 3)]
names(col_pal_grp) <- c("*M. vimineum*", "first-year *E. virginicus*", "adult *E. virginicus*")

grey_pal <- brewer.pal(n = 9, "Greys")[3:9]
grey_pal3 <- shade_pal[c(3, 4, 5)]

coral_pal <- c("coral", "coral3", "coral4")

# dodge size
dodge_width <- 0.5
