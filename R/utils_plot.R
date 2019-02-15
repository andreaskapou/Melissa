# Define ggplot2 theme
#
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

# Plot the predictive distribution
#
draw_predictive <- function(xs, pred, K, title = "", ...){
  xs = Cluster = ys = NULL
  # Store predictions data
  dt <- data.table::data.table(pred) %>%
    data.table::setnames(paste0("C", seq(seq_len(K) ))) %>%
    data.table::melt(variable.name = "Cluster", value.name = "ys") %>%
    .[, xs := xs]

  p <- ggplot(dt, aes(x = xs, y = ys, color = Cluster)) +
    geom_line(size = 2) +
    scale_x_continuous(limits = c(-1, 1),
                       labels = c("-5kb", "", "TSS", "", "+5kb")) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = title, x = "x", y = "y") + gg_theme()
  return(p)
}


# Define ggplot2 theme for line plots
line_theme <- function(){
  p <- theme(
    plot.title = element_text(size = 24, face = 'bold',
                              margin = ggplot2::margin(0,0,0,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.15), color = 'black'),
    axis.title = element_text(size = rel(1.55), color = 'black'),
    axis.title.y = element_text(margin = ggplot2::margin(0,15,0,0)),
    axis.title.x = element_text(margin = ggplot2::margin(12,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 18),
    legend.position = "top",
    legend.key.size = unit(1.9, 'lines'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 19),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}

# Define ggplot2 theme for line plots
line_theme_synth <- function(){
  p <- theme(
    plot.title = element_text(size = 28, face = 'bold',
                              margin = margin(0,0,0,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.35), color = 'black'),
    axis.title = element_text(size = rel(1.95), color = 'black'),
    axis.title.y = element_text(margin = margin(0,15,0,0)),
    axis.title.x = element_text(margin = margin(15,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(2.2, 'lines'),
    legend.title = element_text(size = 26, face = 'bold'),
    legend.text = element_text(size = 22),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}

# Define ggplot2 theme for line plots
boxplot_theme <- function(){
  p <- theme(
    plot.title = element_text(size = 24, face = 'bold',
                              margin = ggplot2::margin(0,0,0,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.15), color = 'black'),
    axis.title = element_text(size = rel(1.55), color = 'black'),
    axis.title.y = element_text(margin = ggplot2::margin(0,15,0,0)),
    axis.title.x = element_text(margin = ggplot2::margin(0,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    legend.position = "left",
    legend.key.size = unit(1.9, 'lines'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 19),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}


#' @title Plot predictive methylaation profiles
#'
#' @description This function plots the predictive distribution of the
#'   methylation profiles inferred using the Melissa model. Each colour
#'   corresponds to a different cluster.
#' @param melissa_obj Clustered cell subtypes using Melissa inference functions.
#' @param region Genomic region number.
#' @param title Plot title
#' @param x_axis x axis label
#' @param y_axis x axis label
#' @param x_labels x axis ticks labels
#' @param ... Additional parameters
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
plot_melissa_profiles <- function(melissa_obj, region = 1, title = "Melissa profiles",
                x_axis = "genomic region", y_axis = "met level",
                x_labels = c("Upstream", "", "Centre", "", "Downstream"), ...) {

  W_Sigma <- list()
  for (cl in seq_along(melissa_obj$pi_k)) {
    W_Sigma[[cl]] <- melissa_obj$W_Sigma[[cl]][[region]]
  }
  obj <- list(W = melissa_obj$W[region,,],
              W_Sigma = W_Sigma,
              basis = melissa_obj$basis)
  class(obj) <- c("cluster_profiles_vb_bernoulli", "cluster_profiles_vb")

  return(BPRMeth::plot_cluster_profiles(cluster_obj = obj, title = title,
                                        x_axis = x_axis, y_axis = y_axis,
                                        x_labels = x_labels, ...))
}


# Create AUC errorbar plot
errorbars_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
  x = y = Model = NULL
  sd_y <- 0.001
  pd <- position_dodge(0.1) # move them .1 to the left and right
  p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
    geom_errorbar(aes(ymin = y - sd_y, ymax = y + sd_y),
                  colour = "black", width = .3, position = pd) +
    geom_line(position = pd, size = 1.5) +
    geom_point(position = pd, size = 3, shape = 21, fill = "white") +
    scale_color_manual(values = c("red3", "chocolate2", "dodgerblue4",
                                  "mediumorchid4", "mistyrose4")) +
    labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
  return(p)
}

# Create AUC jitter plot
# TODO: Add this as function
auc_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
  x = y = Model = NULL
  # pd <- position_dodge(0.1) # move them .1 to the left and right
  p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
    geom_jitter(size = 2.1, width = 0.15, height = -0.1, shape = 21) +
    geom_smooth(aes(fill = Model), span = 0.15, method = "loess",
                se = TRUE, size = 1.3, alpha = 0.1) +
    scale_color_manual(values = c("red3", "chocolate2", "dodgerblue4",
                                  "mediumorchid4", "mistyrose4", "darkgreen")) +
    scale_fill_manual(values = c("red3", "chocolate2", "dodgerblue4",
                                 "mediumorchid4", "mistyrose4", "darkgreen")) +
    labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
  return(p)
}

# Create ARI jitter plot
# TODO: Add this as function
ari_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
  x = y = Model = NULL
  # pd <- position_dodge(0.1) # move them .1 to the left and right
  p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
    geom_jitter(size = 2.1, width = 0.15, height = -0.1, shape = 21) +
    geom_smooth(aes(fill = Model), span = 0.15, method = "loess",
                se = FALSE, size = 1.3, alpha = 0.1) +
    scale_color_manual(values = c("red3",  "dodgerblue4", "chocolate2",
                                  "green4", "cornflowerblue")) +
    scale_fill_manual(values = c("red3",  "dodgerblue4", "chocolate2",
                                 "green4", "cornflowerblue")) +
    labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
  return(p)
}


# Create AUC jitter plot
# TODO: Add this as function
eff_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
  x = y = Model = NULL
  # pd <- position_dodge(0.1) # move them .1 to the left and right
  p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
    geom_jitter(size = 2.1, width = 0.2, height = .3, shape = 21) +
    geom_smooth(aes(fill = Model), method = "loess",
                se = FALSE, size = 1.3, alpha = 0.1) +
    scale_color_manual(values = c("red3", "dodgerblue4", "chocolate2",
                                  "green4", "cornflowerblue")) +
    scale_fill_manual(values = c("red3", "dodgerblue4", "chocolate2",
                                 "green4", "cornflowerblue")) +
    labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
  return(p)
}
