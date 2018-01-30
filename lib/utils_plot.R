##---------------------------------------
# Load libraries
##---------------------------------------
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))

#' Define ggplot2 theme
#'
gg_theme <- function(){
    p <- theme(
        plot.title = element_text(size = 20,face = 'bold',
                                  margin = margin(0,0,3,0), hjust = 0.5),
        axis.text = element_text(size = rel(1.05), color = 'black'),
        axis.title = element_text(size = rel(1.45), color = 'black'),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "gainsboro"),
        panel.background = element_blank()
    )
    return(p)
}

#' Plot the predictive distribution
#'
draw_predictive <- function(xs, pred, title = "", ...){
    # Store predictions data
    dt <- data.table(pred) %>% setnames(paste0("C", seq(1:K))) %>%
        melt(variable.name = "Cluster", value.name = "ys") %>% .[, xs := xs]

    p <- ggplot(dt, aes(x = xs, y = ys, color = Cluster)) +
        geom_line(size = 2) +
        # geom_ribbon(aes(ymin=dt$ys_low, ymax=dt$ys_high, fill=Cluster),
        #               alpha=0.23, size = 0.1) +
        scale_x_continuous(limits = c(-1, 1),
                           labels = c("-5kb", "", "TSS", "", "+5kb")) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = title, x = "x", y = "y") + gg_theme()
}
