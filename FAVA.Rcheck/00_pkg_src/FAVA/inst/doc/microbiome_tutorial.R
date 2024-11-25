## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  fig.width = 8
)

## ----setup--------------------------------------------------------------------
library(FAVA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ----out.width = "1000px", echo = FALSE---------------------------------------
knitr::include_graphics("../man/figures/schematic_data_structure_1.png")

## ----eval=FALSE---------------------------------------------------------------
#  # Example code to read in a data set
#  my_data = read.csv("Path_to/my_data.csv")
#  
#  # If your relative abundances are in a phyloseq object,
#  # make one object combining the sample data (left-hand side)
#  # and the OTU relative abundances (right-hand side)
#  my_data = FAVA::relab_phyloseq(phyloseq_object)
#  
#  # Confirm that your samples each sum to 1
#  # if columns 4 through 10 contain the relevant categories
#  # and columns 1, 2, and 3 contain metadata
#  rowSums(my_data[,c(4:10)])
#  
#  # Example code to convert counts to relative abundances
#  my_data[,c(4:10)] = my_data[,c(4:10)]/rowSums(my_data[,c(4:10)])

## ----eval = FALSE-------------------------------------------------------------
#  # (1)
#  # Here, we assume that you have already generated a phylogenetic tree
#  # and that it is a part of your phyloseq object.
#  tree = phy_tree(phyloseq_object)
#  
#  # (2)
#  distance_matrix = ape::cophenetic.phylo(tree)
#  
#  # (3)
#  # alternative similarity matrices:
#  similarity_matrix = 1/(distance_matrix + 1)
#  similarity_matrix = 1 - distance_matrix/max(distance_matrix)
#  
#  # the similarity matrix we use:
#  similarity_matrix = exp(-distance_matrix)
#  
#  
#  # (4)
#  # Get the names of the species in your relative abundance matrix
#  species_order = colnames(my_data[,c(4:10)])
#  
#  # Confirm that the entries of the similarity matrix
#  # correspond to relative abundance matrix
#  all(species_order == colnames(similarity_matrix))
#  all(species_order == rownames(similarity_matrix))
#  
#  # If they do not, you can re-order the rows and columns of
#  # your similarity matrix to match your data:
#  similarity_matrix_reordered = similarity_matrix[species_order, species_order]
#  
#  # confirm that all diagonal elements are still 1
#  diag(similarity_matrix_reordered)

## ----echo = FALSE, fig.height=1-----------------------------------------------
timeline_data = data.frame(Day = c(1,8,15,22:40, 43, 50, 57, 64)) %>%
  mutate("Day type" = ifelse(Day >28 & Day < 35, "Antiobitic", "Regular"))

timeline <- ggplot() +
  geom_segment(aes(x = 0, xend = 65, y = 0, yend = 0), color = "#666666") +

  geom_segment(aes(x = 1:64, xend = 1:64, y = rep(-1.5, 64), yend = rep(1.5, 64)), color = "#666666") +
  geom_segment(aes(x = seq(from = 1, to = 64, by = 7), xend = seq(from = 1, to = 64, by = 7),
                   y = rep(-2, 10), yend = rep(2, 10)), size = 1, color = "#666666") +
  geom_text(aes(x = seq(from = 1, to = 64, by = 7), y = -4, label = seq(from = 1, to = 64, by = 7))) +

  geom_point(aes(x = Day, y = 0, color = `Day type`), timeline_data, size = 2) +

  geom_text(aes(x = 66, y = -6, label = "Study\nday"), size = 3.5, hjust = 0, vjust = 0, lineheight = 1) +
  geom_text(aes(x = 6, y = 7, label = "Sampling timeline"), size = 5) +

  theme_void() +
  ylim(-10, 10) +
  xlim(0, 70)  +
  scale_color_manual(values = c("#FFB90F", "black")) +
  theme(legend.position = "none") #, axis.title.x = element_text()) + xlab("Study day")
timeline

## ----eval=FALSE---------------------------------------------------------------
#  # open the data set in a new window
#  View(xue_microbiome_sample)
#  
#  # view the structure of the data set
#  str(xue_microbiome_sample)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(xue_microbiome_sample[1:40, 1:20]) %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(xue_species_similarity[1:20, 1:20])  %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----fig.height = 7, echo = FALSE---------------------------------------------
ggplot(xue_species_similarity %>%
         data.frame() %>%
         mutate(name2 = rownames(xue_species_similarity)) %>%
         pivot_longer(cols = 1:524,
                      values_to = "Similarity")) +
  geom_raster(aes(x = name,
                  y = name2,
                  fill = Similarity)) +
  theme_minimal() + scale_fill_viridis_c() +
  theme(axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 2,
                                   angle = -90,
                                   hjust = 0),
        axis.title = element_blank())

## -----------------------------------------------------------------------------
# Make a color palette for all 524 species
set.seed(1)
species_palette = viridis::turbo(524)[sample(1:524)] %>%
  `names<-`(colnames(xue_microbiome_sample)[-c(1:2)])

# Make a ggplot2 stacked bar plot
plot_relabund(xue_microbiome_sample,
              group = "subject",
              time = "timepoint",
              arrange = "vertical",
              K = 524) +
# Specify a custom color scheme
  ggplot2::scale_color_manual(values = species_palette) +
  ggplot2::scale_fill_manual(values = species_palette)

## -----------------------------------------------------------------------------
plot_relabund(xue_microbiome_sample,
              group = "subject",
              arrange = "vertical",
              K = 524) +
  ggplot2::scale_color_manual(values = species_palette) +
  ggplot2::scale_fill_manual(values = species_palette)

## -----------------------------------------------------------------------------
plot_relabund(xue_microbiome_sample,
              arrange = "vertical",
              K = 524) +
  ggplot2::scale_color_manual(values = species_palette) +
  ggplot2::scale_fill_manual(values = species_palette)

## -----------------------------------------------------------------------------
plot_relabund(xue_microbiome_sample,
              arrange = "both",
              K = 524) +
  ggplot2::scale_color_manual(values = species_palette) +
  ggplot2::scale_fill_manual(values = species_palette)

## -----------------------------------------------------------------------------
fava(relab_matrix = xue_microbiome_sample,
     group = "subject",
     K = 524)

## -----------------------------------------------------------------------------
fava(relab_matrix = xue_microbiome_sample,
     K = 524)

## -----------------------------------------------------------------------------
fava(relab_matrix = xue_microbiome_sample,
     group = "subject",
     K = 524, 
     normalized = TRUE)

## -----------------------------------------------------------------------------
antibiotic_data = xue_microbiome_sample %>%
  mutate(Antibiotic = ifelse(timepoint < 29, "Before", 
                             ifelse(timepoint <35, "During", "After")),
         .after = timepoint)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(antibiotic_data[1:20, 1:5])  %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## -----------------------------------------------------------------------------
antibiotic_data = antibiotic_data %>% filter(Antibiotic != "During")

## -----------------------------------------------------------------------------
fava(relab_matrix = antibiotic_data,
     group = c("subject", "Antibiotic"),
     K = 524)

## -----------------------------------------------------------------------------
fava(relab_matrix = xue_microbiome_sample,
     group = "subject",
     K = 524,
     S = xue_species_similarity)

## -----------------------------------------------------------------------------
fava(relab_matrix = xue_microbiome_sample,
     group = "subject",
     K = 524,
     time = "timepoint")

## -----------------------------------------------------------------------------
XMA = filter(xue_microbiome_sample, subject == "XMA")
XMA$timepoint

weights = time_weights(times = XMA$timepoint)
weights
sum(weights)

## -----------------------------------------------------------------------------
ggplot(mapping = aes(x = XMA$timepoint,
                     y = weights)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("Study day") +
  ylab("Weight")

## -----------------------------------------------------------------------------
fava(relab_matrix = XMA,
     K = 524,
     w = weights)

fava(relab_matrix = XMA,
     K = 524,
     time = "timepoint")

## ----fig.width = 5, fig.height = 3--------------------------------------------
fava_out = fava(relab_matrix = xue_microbiome_sample,
     group = "subject",
     K = 524,
     time = "timepoint",
     S = xue_species_similarity)

fava_out

ggplot(fava_out, aes(x = subject, y = FAVA)) + geom_point(size = 4) + theme_bw()

## -----------------------------------------------------------------------------
bootstrap_out = bootstrap_fava(relab_matrix = xue_microbiome_sample,
                               n_replicates = 100,
                               seed = 3,
                               group = "subject",
                               K = 524,
                               # time = "timepoint",
                               S = xue_species_similarity)

## -----------------------------------------------------------------------------
str(bootstrap_out, max.level = 1)

bootstrap_out$P_values


## -----------------------------------------------------------------------------
bootstrap_out$bootstrap_distribution_plot

## -----------------------------------------------------------------------------
window_out = window_fava(relab_matrix = xue_microbiome_sample,
                         window_size = 6, window_step = 1,
                         K = 524,
                         time = "timepoint",
                         S = xue_species_similarity,
                         group = "subject")
head(window_out$window_data)

## -----------------------------------------------------------------------------
window_out$window_plot
window_out$window_plot +
  ggplot2::facet_wrap(~ group)

## -----------------------------------------------------------------------------
str(xue_species_tree)

## -----------------------------------------------------------------------------
ape::plot.phylo(xue_species_tree, cex = 0.2)

## -----------------------------------------------------------------------------
distance_matrix = ape::cophenetic.phylo(xue_species_tree)
str(distance_matrix)

## -----------------------------------------------------------------------------
species_order = colnames(xue_microbiome_sample)[-c(1:2)]
distance_matrix = distance_matrix[species_order, species_order]
str(distance_matrix)

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(distance_matrix[1:40, 1:20]) %>%
    kableExtra::scroll_box(width = "800px", height = "400px")

## ----fig.height = 7, echo = FALSE---------------------------------------------
ggplot(distance_matrix %>%
         data.frame() %>%
         mutate(name2 = rownames(distance_matrix)) %>%
         pivot_longer(cols = 1:524,
                      values_to = "Distance")) +
  geom_raster(aes(x = name,
                  y = name2,
                  fill = Distance)) +
  theme_minimal() + scale_fill_viridis_c() +
  theme(axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 2,
                                   angle = -90,
                                   hjust = 0),
        axis.title = element_blank())

## -----------------------------------------------------------------------------
summary(c(distance_matrix))

## ----echo = FALSE, fig.height = 8---------------------------------------------
pal = c("#157f1f", "#0a2463", "#3891a6")

(data.frame(Distance = seq(from = 0, to = max(distance_matrix), by = 0.05)) %>%
  mutate(Exponential = exp(-Distance), Inverse = 1/(Distance+1), Difference = 1-Distance/max(Distance)) %>%
  pivot_longer(-Distance, names_to = "Transformation", values_to = "Similarity") %>%
  ggplot(aes(x = Distance, y = Similarity, color = Transformation)) + geom_line(size = 3) +
  theme_bw() + scale_color_manual(values = pal)) /
  ggplot() + geom_density(aes(x = c(distance_matrix)), fill = "grey") + theme_bw() +
  xlab("Distance") + ylab("Density") +
  ggtitle("Distribution of pairwise phylogenetic distances between species (distance_matrix)") +
  plot_layout(heights = c(3,2))

## ----echo = FALSE-------------------------------------------------------------
data.frame(Distance = seq(from = 0, to = 100, by = 0.05)) %>%
  mutate(Exponential = exp(-Distance), Inverse = 1/(Distance+1), Difference = 1-Distance/max(Distance)) %>%
  pivot_longer(-Distance, names_to = "Transformation", values_to = "Similarity") %>%
  ggplot(aes(x = Distance, y = Similarity, color = Transformation)) + geom_line(size = 3) +
  theme_bw()+ scale_color_manual(values = pal)

## -----------------------------------------------------------------------------
difference_similarity = 1 - distance_matrix/max(distance_matrix)

exponential_similarity = exp(-distance_matrix)

inverse_similarity = 1/(distance_matrix + 1)

## ----fig.height = 7, echo = FALSE---------------------------------------------
similarity_heatmap <- function(matrix){
  ggplot(matrix %>%
         data.frame() %>%
         mutate(name2 = rownames(matrix)) %>%
         pivot_longer(cols = 1:524,
                      values_to = "Similarity")) +
  geom_raster(aes(x = name,
                  y = name2,
                  fill = Similarity)) +
  theme_minimal() + scale_fill_viridis_c(limits = c(0,1)) +
  theme(axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 2,
                                   angle = -90,
                                   hjust = 0),
        axis.title = element_blank())
}

similarity_heatmap(difference_similarity) + ggtitle("Difference similarity matrix")
similarity_heatmap(exponential_similarity) + ggtitle("Exponential similarity matrix")
similarity_heatmap(inverse_similarity) + ggtitle("Inverse similarity matrix")

## ----echo = FALSE-------------------------------------------------------------
data.frame(Difference = c(difference_similarity),
           Exponential = c(exponential_similarity),
           Inverse = c(inverse_similarity)) %>%
  pivot_longer(cols = 1:3,
               names_to = "Transformation", values_to = "Similarity") %>%
  ggplot(aes(x = Similarity, fill = Transformation)) +
  geom_density(alpha = 0.5) + theme_bw() +
  ylab("Density") +
  # ggtitle("Distributions of pairwise similarities between species for each transformation")+
  scale_fill_manual(values = pal)

## -----------------------------------------------------------------------------
summary(c(difference_similarity))
summary(c(exponential_similarity))
summary(c(inverse_similarity))

## ----echo = FALSE-------------------------------------------------------------
rbind(fava(xue_microbiome_sample, group = "subject", time = "timepoint",
           S = difference_similarity) %>%
        mutate(Transformation = "Difference"),
      fava(xue_microbiome_sample, group = "subject", time = "timepoint",
           S = exponential_similarity) %>%
        mutate(Transformation = "Exponential"),
      fava(xue_microbiome_sample, group = "subject", time = "timepoint",
           S = inverse_similarity) %>%
        mutate(Transformation = "Inverse")) %>%
  ggplot(aes(x = subject, y = FAVA, color = Transformation)) +
  geom_point(size = 6, alpha = 0.5) +
  theme_bw() + geom_line(aes(group = Transformation), size = 1, alpha = 0.7)+ scale_color_manual(values = pal)


## ----echo = FALSE-------------------------------------------------------------
difference_window = window_fava(relab_matrix = xue_microbiome_sample,
                                window_size = 6,
                                group = "subject", time = "timepoint",
                                S = difference_similarity)

exponential_window = window_fava(relab_matrix = xue_microbiome_sample,
                                window_size = 6,
                                group = "subject", time = "timepoint",
                                S = exponential_similarity)

inverse_window = window_fava(relab_matrix = xue_microbiome_sample,
                                window_size = 6,
                                group = "subject", time = "timepoint",
                                S = inverse_similarity)


difference_window$window_plot + ylim(0,0.3) + ggtitle("Difference") +
  theme(legend.position = "none",
        strip.text = element_blank()) +
  facet_grid(group ~ .)  |
  exponential_window$window_plot + ylim(0,0.3)+ ggtitle("Exponential") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank()) +
  facet_grid(group ~ .) |
  inverse_window$window_plot + ylim(0,0.3) + ggtitle("Inverse") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_grid(group ~ .)


## ----echo = FALSE, fig.height=6-----------------------------------------------
fava_window_transformation = rbind(mutate(difference_window$window_data, Transformation = "Difference"),
      mutate(exponential_window$window_data, Transformation = "Exponential"),
      mutate(inverse_window$window_data, Transformation = "Inverse"))

ggplot() +
  geom_polygon(aes(x = x, y= y),
               data.frame(x = rep(c(29, 34, 34, 29), 3),
                          y = rep(c(0, 0, 0.3, 0.3), 3),
                          group = rep(c("XBA", "XDA", "XMA"), each = 4)),
               fill = "grey", alpha = 0.5) +
  geom_point(aes(x = (w6 + w1)/2, y = FAVA, color = Transformation), fava_window_transformation,
             alpha = 0.5, size = 2) +
  geom_line(aes(x = (w6 + w1)/2, y = FAVA, color = Transformation), fava_window_transformation,
            size = 1, alpha = 0.5)  +
  facet_wrap(~ group) +
  theme_bw()+
  scale_color_manual(values = pal) +
  xlab("Study day") +
   ylab("Weighted FAVA") +
  theme(legend.position = c(0.92, 0.8))

## ----echo = FALSE, fig.width = 5, fig.height = 1.8----------------------------
pal = viridis::viridis(n = 5, option = "rocket") %>%
  `names<-`(paste0("Species_", rev(0:4)))

Q = diag(2)%>% 
  `row.names<-`(c("Sample 1: ", 
                  "Sample 2: ")) %>%
  `colnames<-`(c("Species_1", "Species_2"))

patchwork::wrap_plots(plot_relabund(Q) + theme_minimal() +
                        xlab("Sample") +
                        theme(legend.position = "none") +
                        scale_color_manual(values = rep("white", 3))+ 
                        scale_fill_manual(values = pal), 
                     gridExtra::tableGrob(Q),
                       widths = c(1.5,3))

## ----echo = FALSE, fig.width = 7, fig.height = 2------------------------------
Q1 = diag(3) %>% 
  `row.names<-`(c("Sample 1: ", 
                  "Sample 2: ",
                  "Sample 3: ")) %>%
  `colnames<-`(c("Species_1", "Species_2", "Species_3"))

Q2 = matrix(c(1,0,
              0,1, 
              1, 0),
            byrow=TRUE, nrow= 3)%>% 
  `row.names<-`(c("Sample 1: ", 
                  "Sample 2: ",
                  "Sample 3: ")) %>%
  `colnames<-`(c("Species_1", "Species_2"))


patchwork::wrap_plots(plot_relabund(Q1) + theme_minimal() +
                        xlab("Sample") +
                        theme(legend.position = "none") +
                        scale_color_manual(values = rep("white", 3))+ 
                        scale_fill_manual(values = pal) +
                        ggtitle("Matrix A"), 
                      gridExtra::tableGrob(Q1),
                      widths = c(1.7, 3.3))

## ----echo = FALSE, fig.width = 6, fig.height = 2------------------------------

patchwork::wrap_plots(plot_relabund(Q2) + theme_minimal() +
                        xlab("Sample") +
                        theme(legend.position = "none") +
                        scale_color_manual(values = rep("white", 3))+ 
                        scale_fill_manual(values = pal) +
                        ggtitle("Matrix B"), 
                      gridExtra::tableGrob(Q2),
                      widths = c(2,3)) 

## ----echo = FALSE, fig.height = 2, fig.width=11-------------------------------
favamax = function (I, M){
    
    sig1 <- I*M
    J <- ceiling(1/sig1)
    sig1.frac <- sig1 - floor(sig1)
    if (sig1 == I) {
        favaMax <- 0
    }
    else {
        if (sig1 <= 1) {
            favaMax <- ((I - 1) * (1 - sig1 * (J - 1) * (2 - 
                J * sig1)))/(I - (1 - sig1 * (J - 1) * (2 - J * 
                sig1)))
        }
        else {
            favaMax <- (I * (I - 1) - sig1^2 + floor(sig1) - 
                2 * (I - 1) * sig1.frac + (2 * I - 1) * sig1.frac^2)/(I * 
                (I - 1) - sig1^2 - floor(sig1) + 2 * sig1 - sig1.frac^2)
        }
    }
     return(favaMax)
}

x = seq(0.001, 1, by = 0.001)

data.frame(M = x, 
           "I.2" = sapply(x, favamax, I = 2),
           "I.3" = sapply(x, favamax, I = 3),
           "I.6" = sapply(x, favamax, I = 5),
           "I.10" = sapply(x, favamax, I = 10),
           "I.100" = sapply(x, favamax, I = 100)) %>%
  tidyr::pivot_longer(-M, names_to = "I", values_to = "FstMax") %>%
  mutate(I = stringr::str_replace(I, stringr::fixed("."), " = ") %>%
           factor(ordered = TRUE, levels = paste0("I = ", c(2,3,6,10, 100)))) %>%
  ggplot(aes(x = M, y = FstMax)) + geom_line(size = 1) + facet_wrap(~I, nrow = 1) + 
  theme_bw() + ylab(expression(F[ST]^{max}))

## ----echo = FALSE, fig.width = 7, fig.height = 2------------------------------
C = data.frame(Species_1 = c(0.4, 0.75, 0.8),
               Species_2 = c(.55, .05, .05),
               Species_3 = c(0.05, 0.2, 0.15)) %>% 
  `rownames<-`(paste0("Sample ", 1:3))

D = data.frame(Species_1 = c(0.75, 0.85, 0.95),
               Species_2 = c(0.25, 0, 0),
               Species_3 = c(0, 0.15, 0.05)) %>% 
  `rownames<-`(paste0("Sample ", 1:3))

C_plot = plot_relabund(C, arrange = TRUE) + theme_minimal() +
                        xlab("Sample") +
                        theme(legend.position = "none") +
                        scale_color_manual(values = rep("white", 5))+ 
                        scale_fill_manual(values = pal) +
                        ggtitle("Matrix C")

D_plot = plot_relabund(D, arrange = TRUE) + theme_minimal() +
                        xlab("Sample") +
                        theme(legend.position = "none") +
                        scale_color_manual(values = rep("white", 5))+ 
                        scale_fill_manual(values = pal) +
                        ggtitle("Matrix D")


patchwork::wrap_plots(C_plot, 
                      gridExtra::tableGrob(C),
                      widths = c(2,3)) 


patchwork::wrap_plots(D_plot, 
                      gridExtra::tableGrob(D),
                      widths = c(2,3)) 

## -----------------------------------------------------------------------------
fava(C)
fava(D)

## -----------------------------------------------------------------------------
max(colMeans(C))
max(colMeans(D))

## ----fig.height = 4, fig.width = 6, echo = FALSE------------------------------
labels = data.frame(M = c(max(colMeans(C)), max(colMeans(D))),
                    FstMax = c(fava(C), fava(D)),
                    FstMaxMax = c(favamax(I = 3, M = max(colMeans(C))),
                               favamax(I = 3, M = max(colMeans(D)))),
                    Label = c("C", "D"))

data.frame(M = x, 
           "I.3" = sapply(x, favamax, I = 3)) %>%
  tidyr::pivot_longer(-M, names_to = "I", values_to = "FstMax") %>%
  mutate(I = stringr::str_replace(I, stringr::fixed("."), " = ") %>%
           factor(ordered = TRUE, levels = paste0("I = ", c(2,3,6,10, 100)))) %>%
  ggplot(aes(x = M, y = FstMax)) + 
  facet_wrap(~I, nrow = 1) + 
  theme_bw() + ylab(expression(F[ST])) +
  geom_point(data =labels, aes(color = Label), size = 4) + 
  geom_text(data = labels, aes(color = Label, label = Label), 
            size = 8, nudge_x = -.05) +
  geom_segment(data = labels, aes(y = 0, yend = FstMaxMax, 
                                  xend = M, color = Label), size = 2, alpha = 0.4)+
  geom_line(size = 1.5) + 
  theme(legend.position = "none", strip.text = element_text(size = 14)) +
  geom_text(aes(x = .8, y = .9, label = "Upper bound")) 

## -----------------------------------------------------------------------------
fava(C, normalized = TRUE)
fava(D, normalized = TRUE)

