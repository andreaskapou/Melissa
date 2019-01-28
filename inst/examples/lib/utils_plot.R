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


ggplot_cluster_profiles <- function(cluster_obj, title = "Clustered profiles",
                                    x_axis = "genomic region", y_axis = "met level",
                                    x_labels = c("-7Kb", "", "Centre", "", "+7Kb"), ...) {
    # Test data
    aes_xs <- seq(from = -1, to = 1, by = 0.01)
    # Number of clusters
    K <- NCOL(cluster_obj$W)
    ys <- matrix(0, ncol = K, nrow = length(aes_xs))
    # For RMD CHECK to pass without NOTEs
    aes_ys = Cluster <- NULL
    ys_low = ys_high <- NULL
    if (methods::is(cluster_obj, "cluster_profiles_mle")) {
        if (methods::is(cluster_obj, "cluster_profiles_mle_gaussian")) {
            for (k in 1:K) {
                ys[, k] <- eval_function(cluster_obj$basis, aes_xs, cluster_obj$W[,k])
            }
        }else{
            for (k in 1:K) {
                ys[, k] <- eval_probit_function(cluster_obj$basis, aes_xs, cluster_obj$W[,k])
            }
        }
    }else if (methods::is(cluster_obj, "cluster_profiles_vb")) {
        tmp <- BPRMeth:::.predictive_cluster_profile(cluster_obj, aes_xs)
        ys <- tmp$W_pred
        if (methods::is(cluster_obj, "cluster_profiles_vb_binomial") ||
            methods::is(cluster_obj, "cluster_profiles_vb_bernoulli")) {
            ys_low <- ys - ys*(1 - ys);
            ys_high <- ys + ys*(1 - ys)
        }else if (methods::is(cluster_obj, "infer_profiles_vb_gaussian")) {
            ys_low <- ys - 2 * tmp$W_sd_pred;
            ys_high <- ys + 2 * tmp$W_sd_pred
        }
    }else{
        stop("No plotting function for this model!")
    }

    dt <- data.table::data.table(aes_xs = numeric(), aes_ys = numeric(),
                                 ys_low = numeric(), ys_high = numeric(), Cluster = numeric())
    if (methods::is(cluster_obj, "cluster_profiles_vb") ||
        methods::is(cluster_obj, "infer_profiles_gibbs") ) {
        for (k in 1:K) {
            dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs, aes_ys = ys[,k], ys_low = ys_low[,k], ys_high = ys_high[,k], Cluster = as.factor(k)))
        }
    }else{
        for (k in 1:K) {
            dt <- rbind(dt, data.table::data.table(aes_xs = aes_xs, aes_ys = ys[,k], ys_low = 0, ys_high = 0, Cluster = as.factor(k)))
        }
    }

    p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Cluster)) +
        geom_line(size = 2)
    if (methods::is(cluster_obj, "cluster_profiles_vb") ||
        methods::is(cluster_obj, "infer_profiles_gibbs") ) {
        p <- p + geom_ribbon(dt, mapping = aes(ymin = ys_low, ymax = ys_high,
                                               fill = Cluster), alpha = 0.2, size = 0.1)
    }
    p <- p + scale_x_continuous(limits = c(-1, 1), labels = x_labels) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        labs(title = title, x = x_axis, y = y_axis) + line_theme()
    return(p)
}



# Create AUC errorbar plot
errorbars_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
    pd <- position_dodge(0.1) # move them .1 to the left and right
    p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
        geom_errorbar(aes(ymin = y - sd_y, ymax = y + sd_y),
                      colour = "black", width = .3, position = pd) +
        geom_line(position = pd, size = 1.5) +
        geom_point(position = pd, size = 3, shape = 21, fill = "white") +
        #scale_color_manual(values = c("red3", "cornflowerblue", "chocolate2", "green4", "dodgerblue4")) +
        scale_color_manual(values = c("red3", "chocolate2", "dodgerblue4", "mediumorchid4", "mistyrose4")) +
        #scale_y_continuous(limits = c(0.78, .96), breaks=pretty_breaks(n=6)) +
        labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
    return(p)
}

# Create AUC jitter plot
auc_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
    # pd <- position_dodge(0.1) # move them .1 to the left and right
    p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
        geom_jitter(size = 2.1, width = 0.15, height = -0.1, shape = 21) +
        # geom_point(data = dt2, aes(x = x, y = y), position = pd, size = 3, shape = 21, fill = "white") +
        geom_smooth(aes(fill = Model), span = 0.15, method = "loess",
                    se = TRUE, size = 1.3, alpha = 0.1) +
        #scale_color_manual(values = c("red3", "cornflowerblue", "chocolate2", "green4", "dodgerblue4")) +
        #scale_fill_manual(values = c("red3", "cornflowerblue", "chocolate2", "green4", "dodgerblue4")) +
        scale_color_manual(values = c("red3", "chocolate2", "dodgerblue4", "mediumorchid4", "mistyrose4", "darkslategray4")) +
        scale_fill_manual(values = c("red3", "chocolate2", "dodgerblue4", "mediumorchid4", "mistyrose4", "darkslategray4")) +
        labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
    return(p)
}

# Create ARI jitter plot
ari_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
    # pd <- position_dodge(0.1) # move them .1 to the left and right
    p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
        geom_jitter(size = 2.1, width = 0.15, height = -0.1, shape = 21) +
        # geom_point(data = dt2, aes(x = x, y = y), position = pd, size = 3, shape = 21, fill = "white") +
        geom_smooth(aes(fill = Model), span = 0.15, method = "loess",
                    se = FALSE, size = 1.3, alpha = 0.1) +
        scale_color_manual(values = c("red3",  "dodgerblue4", "darkslategray4", "chocolate2", "green4", "cornflowerblue")) +
        scale_fill_manual(values = c("red3",  "dodgerblue4", "darkslategray4", "chocolate2", "green4", "cornflowerblue")) +
        labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
    return(p)
}


# Create AUC jitter plot
eff_jitter_plot <- function(dt, title = "", x_lab = "", y_lab = ""){
    # pd <- position_dodge(0.1) # move them .1 to the left and right
    p <- ggplot(dt, aes(x = x, y = y, colour = Model, group = Model)) +
        geom_jitter(size = 2.1, width = 0.2, height = .3, shape = 21) +
        # geom_point(data = dt2, aes(x = x, y = y), position = pd, size = 3, shape = 21, fill = "white") +
        geom_smooth(aes(fill = Model), method = "loess",
                    se = FALSE, size = 1.3, alpha = 0.1) +
        scale_color_manual(values = c("red3", "dodgerblue4", "chocolate2", "green4", "cornflowerblue")) +
        scale_fill_manual(values = c("red3", "dodgerblue4", "chocolate2", "green4", "cornflowerblue")) +
        labs(title = title, x = x_lab, y = y_lab) + line_theme_synth()
    return(p)
}
