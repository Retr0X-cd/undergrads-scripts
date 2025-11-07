rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
options(rgl.useNULL=TRUE)

#### User Input ####
options(echo = TRUE) 
site_names <- c('Matnog', 'Taluksangay')

links <-

file_path <- paste0('../data_processed/xrays/tps-final/')

wireframe_vis <- paste0 ("../plots/",
                         "truss_wireframes_geomorph.png")

pca_space_file <- paste0('../plots/',
                   'truss_pca_normalized(site)_geomorph.png')

pca_time_file <- paste0('../plots/',
                        'truss_pca_normalized(time)_geomorph.png')

pca_spacextime_file <- paste0 ('../plots/',
                              'truss_pca_normalized(site-time)_geomorph.png')

pca_space_raw_file <- paste0('../plots/',
                             'truss_pca_raw(site)_geomorph.png')

pca_time_raw_file <- paste0('../plots/',
                            'truss_pca_raw(time)_geomorph.png')

pca_spacextime_raw_file <- paste0('../plots/',
                                  'truss_pca_raw(site-time)_geomorph.png')

pca_comparison <- paste0('../plots/',
                         'truss_pca_comparison(site-time)_geomorph.png')

pca_biplot_file <- paste0('../plots/',
                          'truss_pca_biplot(site-time)_geomorph.png')

cva_plot <- paste0('../plots/',
                   'cva_groups.png')

dfa_file <- paste0('../plots/',
                   'dfa.png')

manova_summary_csv <- paste0('../data_processed/xrays/',
                              'truss_manova(site).csv')

truss_anova_csv <- paste0('../data_processed/xrays/',
                           'truss_anova_results(site-time).csv')

diagnostic_plot_file <- paste0('../plots/',
                           'xray_anova_significant_trusses.png')

results_table_csv <- paste0('../data_processed/xrays/',
                            'truss_matrix_normalized(site)_geomorph.csv')

truss_scores_csv <- paste0('../data_processed/xrays/',
                           'truss_matrix_raw(site)_geomorph.csv')

anova_table_csv <- paste0('../data_processed/xrays/',
                          'truss_pca_scores(site)_geomorph.csv')
#### Links ####
truss_links <- matrix(c(
  1,4, 4,5, 5,6, 6,11, 7,11, 7,9, 9,10, 1,10,
  1,2, 2,3, 6,10, 6,9, 11,17, 7,12, 6,12, 7,14,
  9,13, 9,14, 12,14, 13,14, 12,15, 12,13, 14,18,
  13,16, 13,15, 13,18, 14,15, 15,16, 6,15, 16,18, 17,18,
  17,19, 16,19, 18,19, 15,20, 11,20, 17,20
), byrow = TRUE, ncol = 2)

#### Functions ####

get_pca_standard_trusses <- function(truss_matrix, pca, nSD = 3, axis.choice = c(1, 2)) {
  # truss_matrix: specimens × truss distances
  
  # Mean truss vector
  mean_truss <- colMeans(truss_matrix, na.rm = TRUE)
  
  # Prepare list for each axis
  pca_extremes <- vector("list", length(axis.choice))
  names(pca_extremes) <- paste0("Comp", axis.choice)
  
  for (axis in seq_along(axis.choice)) {
    # PCA loading vector for this axis
    loadings <- pca$rotation[, axis.choice[axis]]
    
    # Scale loadings by SD units
    pca_extremes[[axis]] <- map(c(1, -1) * nSD * pca$sdev[axis.choice[axis]], function(scale_val) {
      mean_truss + (loadings * scale_val)
    }) %>%
      set_names(c("max_truss", "min_truss")) %>%
      map(as_tibble_row)
  }
  
  c(list(mean_truss = as_tibble_row(mean_truss)), pca_extremes)
}



make_truss_distance_plot <- function(truss_list,
                                     coord_ref,
                                     links,
                                     comp = "Comp1",
                                     which_extreme = "max_truss",
                                     point_size = 1) {
  # Extract numeric truss vector from list
  truss_values <- as.numeric(truss_list[[comp]][[which_extreme]])
  
  # Build dataframe for plotting
  df_links <- tibble::tibble(
    x = coord_ref[links[, 1], 1],
    y = coord_ref[links[, 1], 2],
    xend = coord_ref[links[, 2], 1],
    yend = coord_ref[links[, 2], 2],
    truss_length = truss_values
  )
  
  ggplot2::ggplot(df_links, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(xend = xend, yend = yend, color = truss_length),
                          linewidth = 1) +
    ggplot2::geom_point(data = as.data.frame(coord_ref),
                        ggplot2::aes(x = V1, y = V2),
                        size = point_size) +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()
}




#### 1. Load Packages ####
library(ggplot2)
library(geomorph) 
library(abind) 
library(tidyverse) 
library(patchwork)
library(dplyr)
library(scales)
library(tibble)
library(stats)
library(RRPP)
library(Morpho)
library(MASS)

#### Read and Combine TPS Data ####
all_tps_raw <- list.files(file_path, pattern = 'TPS$', full.names = TRUE) %>%
  tibble(file = .) %>%
  mutate(site = stringr::str_extract(file, stringr::str_c(site_names, collapse = '|')),
         tps_out = map(file, ~geomorph::readland.tps(.x, specID = 'ID')))
all_tps_raw$n_specimens <- map_int(all_tps_raw$tps_out, ~ifelse(length(dim(.x)) >= 3, dim(.x)[3], 1))

land_all<- abind(all_tps_raw$tps_out, along = 3)
sample_names <- dimnames(land_all)[[3]] %||% paste0("S", seq_len(dim(land_all)[3]))

#### Create Enhanced Metadata ####
metadata <- tibble(
  sample_id = dimnames(land_all)[[3]],
  spatial = case_when(
    str_detect(sample_id, "Mat") ~ "Matnog",
    str_detect(sample_id, "Tlk|Sac") ~ "Taluksangay"
  ),
  temporal = case_when(
    str_detect(sample_id, "CMat") ~ "Contemporary",
    str_detect(sample_id, "AMat") ~ "Historical",
    str_detect(sample_id, "CTlk") ~ "Contemporary",
    str_detect(sample_id, "ASac") ~ "Historical"
  ),
  group = case_when(
    str_detect(sample_id, "CMat") ~ "CMat",
    str_detect(sample_id, "AMat") ~ "AMat",
    str_detect(sample_id, "CTlk") ~ "CTlk",
    str_detect(sample_id, "ASac") ~ "ASac"
  )
) 

#filter samples that are NA 

keep_samples <- metadata %>% 
  filter(!is.na(spatial) & !is.na(temporal) & !is.na(group)) %>% 
  pull(sample_id)

#check whether any specimens will be filtered
cat("Original number of specimens:", dim(land_all)[3], "\n")
cat("Samples with complete metadata:", length(keep_samples), "\n")
cat("Samples being removed:", dim(land_all)[3] - length(keep_samples), "\n")

land_all_filtered <- land_all[,,keep_samples, drop = FALSE]

sample_names <- dimnames(land_all_filtered)[[3]]

metadata <- metadata %>% 
  filter(sample_id %in% keep_samples)

#### 2. Sanity Check ####
sanity_table <- metadata %>%
  count(spatial, temporal, name = "n_specimens") %>%
  pivot_wider(names_from = temporal, values_from = n_specimens, values_fill = 0)

n_landmarks <- dim(land_all)[1]
n_specimens_total <- dim(land_all)[3]

cat("\n===== SANITY CHECK =====\n")
print(sanity_table)
cat("No. of Landmarks: ", n_landmarks, "\n")
cat("No. of truss links (input):", nrow(truss_links), "\n")
cat("Total Specimens: ", n_specimens_total, "\n")

#### 3. Truss Wireframe Visualization ####
plot_truss_wireframe <- function(coords, links) {
  df <- tibble(landmark = 1:nrow(coords), x = coords[,1], y = coords[,2])
  link_df <- tibble(
    x_start = coords[links[,1],1], y_start = coords[links[,1],2],
    x_end   = coords[links[,2],1], y_end   = coords[links[,2],2]
  )
  ggplot(df, aes(x, y)) +
    geom_segment(data=link_df, aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                 color="black", linewidth=0.7) +
    geom_point(color="black", size=3) +
    geom_text(aes(label=landmark), vjust=-1, size=3) +
    coord_fixed(clip="off") +
    expand_limits(x = range(df$x) + c(-2, 2), y = range(df$y) + c(-2, 2)) +
    theme_void(base_size = 14) +
    theme(panel.background = element_rect(fill="white", colour="white"),
          plot.background  = element_rect(fill="white", colour="white"))
}
mean_shape <- geomorph::mshape(land_all_filtered)
ggsave(wireframe_vis,
       plot_truss_wireframe(mean_shape, truss_links),
       height=10, width=10, bg="white")

#### 4. Build Truss Distance Matrix ####

truss_matrix <- map_dfr(seq_len(dim(land_all_filtered)[3]), function(i){
  d <- interlmkdist(land_all_filtered[,,i], truss_links)
  tibble(sample_id = sample_names[i], !!!set_names(as.list(d), paste0("LM", truss_links[,1], "-", truss_links[,2])))
}) %>%
  left_join(metadata, by = "sample_id") %>% 
  relocate(spatial, temporal, group, .after = sample_id)

#mean absolute difference in old custom function vs interlmkdistance: 2.03796e-15 (negligible)

#### Size Adjustment of Measurements ####
# Function to calculate distance between two landmarks
calc_distance <- function(x1, y1, x2, y2) {
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

# Calculate Standard Length (SL) for each specimen (distance between landmarks 1-17)
calculate_SL <- function(coords) {
  x1 <- coords[1, 1]
  y1 <- coords[1, 2]
  x17 <- coords[17, 1]
  y17 <- coords[17, 2]
  calc_distance(x1, y1, x17, y17)
}

# Calculate Standard Lengths for all specimens
Lo <- sapply(seq_len(dim(land_all_filtered)[3]), function(i) {
  calculate_SL(land_all_filtered[,,i])
})

# Calculate grand mean of Standard Lengths (Ls)
Ls <- mean(Lo)

# Function to calculate allometric coefficient (b) for each measurement
calc_b <- function(M, SL) {
  # Avoid log(0) and handle NAs
  valid <- !is.na(M) & !is.na(SL) & M > 0 & SL > 0
  if (sum(valid) < 2) return(NA)
  
  # Calculate logs and center them (as in Excel)
  log_M <- log(M[valid])
  log_SL <- log(SL[valid])
  mean_log_M <- mean(log_M)
  mean_log_SL <- mean(log_SL)
  
  # Center the log values
  log_M_centered <- log_M - mean_log_M
  log_SL_centered <- log_SL - mean_log_SL
  
  # Calculate regression slope manually (as Excel does)
  b <- sum(log_M_centered * log_SL_centered) / sum(log_SL_centered^2)
  return(b)
}

# Calculate size-adjusted measurements
truss_adjusted <- truss_matrix %>%
  mutate(across(starts_with("LM"), function(M) {
    b <- calc_b(M, Lo)
    if (is.na(b)) return(M)  # If can't calculate b, return original
    
    # Apply the adjustment formula with careful handling of precision
    adj <- M * (Ls/Lo)^b
    # Round to match Excel precision (optional, remove if you want full precision)
    round(adj, 10)
  }, .names = "{.col}_adj"))

# Print some values for verification
cat("\nVerification of calculations:\n")
first_adj_col <- grep("_adj$", names(truss_adjusted), value = TRUE)[1]
first_measurement <- as.numeric(truss_adjusted[[first_adj_col]][1])
cat("First specimen (", truss_adjusted$sample_id[1], ") measurements:\n")
cat("First Madj value (", first_adj_col, "):", first_measurement, "\n")
cat("Standard Length (Lo):", Lo[1], "\n")
cat("Mean Standard Length (Ls):", Ls, "\n")

# Print regression info for the first measurement
first_raw_col <- gsub("_adj$", "", first_adj_col)
M <- truss_matrix[[first_raw_col]]
b <- calc_b(M, Lo)
cat("\nRegression info for", first_raw_col, ":\n")
cat("b coefficient:", b, "\n")

#### 5a. Principal Component Analysis (PCA) ####
# Get adjusted columns
adj_cols <- grep("_adj$", names(truss_adjusted), value = TRUE)
pca_df <- as.data.frame(lapply(truss_adjusted[adj_cols], function(x) as.numeric(x)))

# Check for zero variance columns
vars <- sapply(pca_df, var, na.rm = TRUE)
if (any(vars == 0)) {
  zero_cols <- names(vars)[vars == 0]
  warning("Removing zero-variance columns: ", paste(zero_cols, collapse = ", "))
  pca_df <- pca_df[, vars != 0, drop = FALSE]
}

truss_pca <- prcomp(as.matrix(pca_df), center = TRUE, scale. = TRUE)
explained_var <- (truss_pca$sdev^2) / sum(truss_pca$sdev^2)
truss_scores <- as_tibble(truss_pca$x) %>%
  dplyr::bind_cols(dplyr::select(truss_adjusted, sample_id, spatial, temporal, group))

truss_pca_plot_spatial <- ggplot(truss_scores,
                              aes(x = PC1, y = PC2, color = spatial)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  # ggrepel::geom_text_repel(aes(label = sample_id), size = 2, show.legend = TRUE) +
  labs(
    x = paste0("PC1 (", percent(explained_var[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var[2], accuracy = 0.1), ")"),
    colour = "Spatial",
    title = "PCA - Size-adjusted Truss Distances (Spatial)"
  )
print(truss_pca_plot_spatial)
ggsave(pca_space_file, truss_pca_plot_spatial, height = 6, width = 8)

truss_pca_plot_temporal <- ggplot(truss_scores,
                              aes(x = PC1, y = PC2, color = temporal)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  # ggrepel::geom_text_repel(aes(label = sample_id), size = 2, show.legend = TRUE) +
  labs(
    x = paste0("PC1 (", percent(explained_var[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var[2], accuracy = 0.1), ")"),
    colour = "Temporal",
    title = "PCA - Size-adjusted Truss Distances (Temporal)"
  )

#work in progress: PCA extremes visualization
truss_num <- truss_adjusted %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(matches("_adj"))

sample_truss_plots_list <- get_pca_standard_trusses(truss_num, truss_pca, nSD = 5, axis.choice = c(1, 2))

plots_list <- purrr::imap(sample_truss_plots_list[-1], ~
                            make_truss_distance_plot(
                              sample_truss_plots_list,
                              coord_ref = mean_shape,
                              links = truss_links,
                              comp = .y,
                              which_extreme = "max_truss"
                            )
)

design <- "
  4111
  6111
  5111
  7382
"

(truss_pca_plot_temporal + theme(plot.margin = margin())) +
  wrap_plots(plots_list) +
  plot_spacer() + plot_spacer() + plot_spacer() +
  plot_layout(design = design) +
  plot_annotation(title = "PCA of *S. delicatulus* External Morphology")

print(truss_pca_plot_temporal)
ggsave(pca_time_file, truss_pca_plot_temporal, height = 6, width = 8)


truss_pca_plot_spatialandtemporal <- ggplot(truss_scores,
                                     aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  # ggrepel::geom_text_repel(aes(label = sample_id), size = 2, show.legend = TRUE) +
  labs(
    x = paste0("PC1 (", percent(explained_var[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var[2], accuracy = 0.1), ")"),
    colour = "Groups",
    title = "PCA - Size-adjusted Truss Distances (Spatial and Temporal)"
  )
print(truss_pca_plot_spatialandtemporal)
ggsave(pca_spacextime_file, truss_pca_plot_spatialandtemporal, height = 6, width = 8)

#### 5b. PCA on Raw Truss Distances ####
raw_cols <- grep("^LM", names(truss_matrix), value = TRUE)
pca_df_raw <- as.data.frame(lapply(truss_matrix[raw_cols], function(x) as.numeric(x)))
vars_raw <- sapply(pca_df_raw, var, na.rm = TRUE)
if (any(vars_raw == 0)) {
  zero_cols_raw <- names(vars_raw)[vars_raw == 0]
  warning("Removing zero-variance raw columns: ", paste(zero_cols_raw, collapse = ", "))
  pca_df_raw <- pca_df_raw[, vars_raw != 0, drop = FALSE]
}
truss_pca_raw <- prcomp(as.matrix(pca_df_raw), center = TRUE, scale. = TRUE)
explained_var_raw <- (truss_pca_raw$sdev^2) / sum(truss_pca_raw$sdev^2)
truss_scores_raw <- as_tibble(truss_pca_raw$x) %>%
  dplyr::bind_cols(dplyr::select(truss_matrix, sample_id, spatial, temporal, group))

truss_pca_plot_raw_spatial <- ggplot(truss_scores_raw,
                                  aes(x = PC1, y = PC2, color = spatial)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", percent(explained_var_raw[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var_raw[2], accuracy = 0.1), ")"),
    colour = "Spatial",
    title = "PCA - Raw Truss Distances (Spatial)"
  )
print(truss_pca_plot_raw_spatial)
ggsave(pca_space_raw_file, truss_pca_plot_raw_spatial, height = 6, width = 8)

truss_pca_plot_raw_temporal <- ggplot(truss_scores_raw,
                                  aes(x = PC1, y = PC2, color = temporal)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", percent(explained_var_raw[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var_raw[2], accuracy = 0.1), ")"),
    colour = "Temporal",
    title = "PCA - Raw Truss Distances (Temporal)"
  )
print(truss_pca_plot_raw_temporal)
ggsave(pca_time_raw_file, truss_pca_plot_raw_temporal, height = 6, width = 8)

truss_pca_plot_raw_spatialandtemporal <- ggplot(truss_scores_raw,
                                         aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", percent(explained_var_raw[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var_raw[2], accuracy = 0.1), ")"),
    colour = "Groups",
    title = "PCA - Raw Truss Distances (Spatial and Temporal)"
  )
print(truss_pca_plot_raw_spatialandtemporal)
ggsave(pca_spacextime_raw_file, truss_pca_plot_raw_spatialandtemporal, height = 6, width = 8)

comparison_plot <- truss_pca_plot_raw_spatialandtemporal + truss_pca_plot_spatialandtemporal + patchwork::plot_layout(ncol = 2)
print(comparison_plot)
ggsave(pca_comparison, comparison_plot, height = 6, width = 14)

#### 5c. PCA Biplot ####
loadings <- as.data.frame(truss_pca$rotation[, 1:2])
loadings$truss <- rownames(loadings)
loadings$truss <- gsub("_norm$", "", loadings$truss)
arrow_scale <- 3.5 * max(abs(truss_scores$PC1), abs(truss_scores$PC2))
loadings$PC1_arrow <- loadings$PC1 * arrow_scale
loadings$PC2_arrow <- loadings$PC2 * arrow_scale

biplot <- ggplot(truss_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1_arrow, yend = PC2_arrow),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "gray40", linewidth = 1
  ) +
  geom_text(
    data = loadings,
    aes(x = PC1_arrow * 1.08, y = PC2_arrow * 1.08, label = truss),
    size = 3, color = "black", fontface = "italic"
  ) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw(base_size = 15) +
  labs(
    x = paste0("PC1 (", percent(explained_var[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", percent(explained_var[2], accuracy = 0.1), ")"),
    colour = "Groups",
    title = "PCA Biplot - Size-adjusted Truss Distances"
  ) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.title.align = 0.5,
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 13)
  )
print(biplot)
ggsave(pca_biplot_file, biplot, height = 7, width = 9)

#### 6. Canonical Variate Analysis (CVA) ####
cva_data <- as.matrix(pca_df)  # Using the size-adjusted data
cva_groups <- truss_adjusted$group 

#clear incomplete samples
complete_cases <- !is.na(cva_groups)
cva_data <- cva_data[complete_cases, ]
cva_groups <- cva_groups[complete_cases]

cva_result <- CVA(cva_data,cva_groups)
# Debug: Check what the CVA returned
cat("\n--- CVA Debug Info ---\n")
cat("Number of CVs available:", ncol(cva_result$CVs), "\n")
cat("CV names:", colnames(cva_result$CVs), "\n")
cat("Dimensions of CVs:", dim(cva_result$CVs), "\n")

cva_scores <- as_tibble(cva_result$CVs) %>%
  set_names(paste0("CV", 1:ncol(cva_result$CVs))) %>%
  mutate(
    sample_id = truss_adjusted$sample_id[complete_cases],
    group = cva_groups
  ) %>%
  relocate(sample_id, group, .before = CV1)

# variance explained by each CV
cva_var <- cva_result$Var
cva_var_explained <- cva_var / sum(cva_var) * 100

cva_plot_groups <- ggplot(cva_scores, aes(x = CV1, y = CV2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  theme_bw() +
  labs(
    x = paste0("CV1 (", round(cva_var_explained[1], 1), "%)"),
    y = paste0("CV2 (", round(cva_var_explained[2], 1), "%)"),
    color = "Group",
    shape = "Group",
    title = "Canonical Variate Analysis - Groups"
  ) +
  theme(legend.position = "bottom")

print(cva_plot_groups)
ggsave(cva_plot, cva_plot_groups, height = 7, width = 8)


#### 7. Discriminant Function Analysis (DFA) (Size-adjusted Data) ####

pc_keep <- which(cumsum(explained_var) <= 0.95)
if (length(pc_keep) < 2) pc_keep <- 1:2

dfa_data <- truss_scores %>%
  dplyr::select(all_of(paste0("PC", pc_keep)), group, spatial, temporal, sample_id)

# 3. Run DFA (LDA)
dfa_model <- lda(group ~ ., data = dfa_data[, c(paste0("PC", pc_keep), "group")])

# 4. Get DFA scores
dfa_pred <- predict(dfa_model)

dfa_results <- dfa_data %>%
  mutate(DFA1 = dfa_pred$x[,1],
         DFA2 = dfa_pred$x[,2],
         Predicted = dfa_pred$class)

# 5. Plot DFA
dfa_plot <- ggplot(dfa_results, aes(x = DFA1, y = DFA2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  stat_ellipse(type = "norm", level = 0.68, show.legend = FALSE) +
  labs(title = "Discriminant Function Analysis (DFA)",
       x = "DFA Axis 1",
       y = "DFA Axis 2")

dfa_cv <- lda(group ~ ., data = dfa_data[, c(paste0("PC", pc_keep), "group")], CV = TRUE)

confusion <- table(Actual = dfa_data$group, Predicted = dfa_cv$class)
accuracy <- mean(dfa_data$group == dfa_cv$class)

print(confusion)
cat("Classification accuracy:", round(accuracy * 100, 2), "%\n")

ggsave(dfa_file, dfa_plot, height = 7, width = 8)

#### 8. MANOVA on PCs (Size-adjusted Data) ####
pc_keep <- which(cumsum(explained_var) <= 0.95)
if (length(pc_keep) < 2) pc_keep <- 1:2
manova_data <- dplyr::select(truss_scores, group, all_of(paste0("PC", pc_keep))) #MANOVA test of PCA on Size-adjusted data; spatialXtemporal grouping
if (ncol(manova_data) > 1) {
  manova_res <- manova(as.matrix(manova_data[,-1]) ~ group, data=manova_data)
  manova_summary <- summary(manova_res, test="Wilks")
  print(manova_summary)
  write_csv(as.data.frame(manova_summary$stats),
            manova_summary_csv)
}

#post-hoc test
manova_pairwise <- as.matrix(manova_data[,-1])
pairwise_groups <- manova_data$group

fit <- lm.rrpp(manova_pairwise ~ pairwise_groups, iter = 999)

pw <- pairwise(fit, groups = pairwise_groups)
summary(pw)       

#residuals
manova_residuals <- residuals(manova_res)

par(mfrow = c(2, 2))

# 1. Q-Q plot for normality
qq_data <- scale(manova_residuals)  # Standardize residuals
qqnorm(qq_data, main = "Q-Q Plot: MANOVA Residuals")
qqline(qq_data, col = "red")

# 2. Residuals vs Fitted
fitted_values <- fitted(manova_res)
plot(fitted_values, manova_residuals, 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

# 3. Histogram of residuals
hist(manova_residuals, breaks = 20, 
     xlab = "Residuals", main = "Distribution of MANOVA Residuals",
     col = "lightblue", border = "black")

# 4. Scale-Location plot
std_residuals <- sqrt(abs(scale(manova_residuals)))
plot(fitted_values, std_residuals,
     xlab = "Fitted Values", ylab = "√|Standardized Residuals|",
     main = "Scale-Location Plot")

# Reset plot parameters
par(mfrow = c(1, 1))

#### 9. ANOVA on Individual Truss Distances ####
anova_table <- truss_adjusted %>%
  pivot_longer(ends_with("_adj"), names_to="truss", values_to="distance") %>%
  group_by(truss) %>%
  summarise(p_value = summary(aov(distance ~ group, data=cur_data()))[[1]][["Pr(>F)"]][1],
            .groups="drop") %>%
  mutate(p_adj = p.adjust(p_value, method="fdr"))

print(anova_table)
write_csv(anova_table, truss_anova_csv)

#residuals 
significant_trusses <- anova_table %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj) %>%
  head(6)

if(nrow(significant_trusses) > 0) {
  diagnostic_plots <- list()
  
  for(i in 1:min(3, nrow(significant_trusses))) {
    truss_name <- significant_trusses$truss[i]
    
    # Manually extract the data to avoid column name issues
    truss_data <- truss_adjusted %>%
      dplyr::select(sample_id, group, all_of(truss_name)) %>%
      filter(!is.na(group))
    
    # Create a clean data frame for plotting
    plot_df <- data.frame(
      group = truss_data$group,
      distance = truss_data[[truss_name]],  # Use [[ ]] extraction
      truss = truss_name
    )
    
    p <- ggplot(plot_df, aes(x = group, y = distance, fill = group)) +
      geom_boxplot(alpha = 0.7) +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.6) +
      labs(title = paste("Distribution:", truss_name),
           subtitle = paste("p-value:", round(significant_trusses$p_value[i], 4)),
           x = "Group", y = "Size-adjusted Distance") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    diagnostic_plots[[i]] <- p
  }
  
  if(length(diagnostic_plots) > 0) {
    combined_diagnostics <- wrap_plots(diagnostic_plots, ncol = 2)
    ggsave(diagnostic_plot_file, combined_diagnostics, 
           width = 12, height = 8)
    print(combined_diagnostics)
  }
}

# Combine and save diagnostic plots
if(length(diagnostic_plots) > 0) {
  combined_diagnostics <- wrap_plots(diagnostic_plots, ncol = 2)
  ggsave(diagnostic_plot_file, combined_diagnostics, 
         width = 12, height = 8)
  print(combined_diagnostics)
}

#### 10. Export Results ####
# Export measurements with both raw and size-adjusted values
results_table <- truss_matrix %>%
  # Keep metadata and raw measurements
  dplyr::select(sample_id, spatial, temporal, group, starts_with("LM")) %>%
  # Add size-adjusted measurements
  bind_cols(
    truss_adjusted %>%
      dplyr::select(ends_with("_adj")) %>%
      rename_with(~gsub("_adj$", "_Madj", .x))
  ) %>%
  # Rename raw measurements for clarity
  rename_with(~paste0(.x, "_Raw"), .cols = starts_with("LM") & !ends_with("_Madj"))

# Export all results
write_csv(results_table, results_table_csv)
write_csv(truss_scores, truss_scores_csv)
write_csv(anova_table, anova_table_csv)