# Script to create supplementary figures
#
# Author: James S. Santangelo

###################
#### FIGURE S2 ####
###################

## Environmental variable against distance

#' Plot each environmental variable against distance
#' 
#' @param df_all_popMeans Population-mean dataframe for all cities with environmental data
#' @param response Character vector for response variable for which to run robust regression
#' 
#' @return ggplot object
plot_envar_vs_dist <- function(df_all_popMeans, response){
  
  # Pvalues for clines in environmental variable
  clinesSummary <- df_all_popMeans %>% group_split(city) %>% 
    purrr::map_dfr(., rlmStats, response = response)
  
  # Color palette
  pal <- c("#909090", "#FF0000", "#046C9A")
  
  # Y-axis label
  ylab <- case_when(response == 'annualAI_Mean' ~ 'Standardized annual AI',
                    response == 'annualPET_Mean' ~ 'Standardized annual PET',
                    response == 'DEM_Mean' ~ 'Standardized elevation',
                    response == 'GMIS_Mean' ~ 'Standardized % impervious surface',
                    response == 'NDSI_Mean' ~ 'Standardized NDSI',
                    response == 'summerLST_Mean' ~ 'Standardized summer LST',
                    response == 'summerNDVI_Mean' ~ 'Standardized summer NDVI',
                    response == 'winterLST_Mean' ~ 'Standardized winter LST',
                    response == 'winterNDVI_Mean' ~ 'Standardized winter NDVI',)
  
  # Datamae with response variable standardized between 0 and 1
  df <- df_all_popMeans %>%
    left_join(., clinesSummary, by = "city") %>% 
    mutate(significant = ifelse(pvalRLM < 0.05, "Yes", "No"),
           direction = ifelse(betaRLM > 0, "Positive", "Negative"),
           color = case_when(significant == "Yes" & direction == "Positive" ~ "Significantly positive",
                             significant == "Yes" & direction == "Negative" ~ "Significantly negative",
                             TRUE ~ "Not significant")) %>% 
    group_by(city) %>% 
    # mutate(scaled = scale(!!sym(response),
    #                       center = min(!!sym(response)),
    #                       scale = max(!!sym(response)) - min(!!sym(response))))
    mutate(scaled = scale(!!sym(response)))
  
  # Create plot
  plot <- df %>%
    ggplot(., aes(x = std_distance, y = scaled)) +
    geom_line(stat = "smooth", formula = y ~ x,
              method=function(formula,data,weights=weight) rlm(formula,
                                                                 data,
                                                                 weights=weight,
                                                                 maxit=200),
              aes(color = color, alpha = color, group = city),
              size = 0.75, se = FALSE) +
    geom_line(stat = "smooth", formula = y ~ x, method="lm", colour = "black", size = 2.5) +
    xlab("Standardized distance") + ylab(ylab) +
    scale_colour_manual(values = pal, limits = c('Not significant', 'Significantly negative', 'Significantly positive')) +
    scale_alpha_manual(values = c(0.5, 0.5, 0.5), limits = c('Not significant', 'Significantly negative', 'Significantly positive')) +
    scale_y_continuous(breaks = seq(from = -3, to = 3, by = 1)) +
    scale_x_continuous(breaks = seq(from = 0, to = 1.1, by = 0.25)) +
    coord_cartesian(ylim = c(-3.05, 3.05), xlim = c(0, 1), clip = 'off') +
    ng1 + guides(alpha = guide_legend(override.aes = list(alpha = 1)))

  return(plot)
}

a <- plot_envar_vs_dist(df_all_popMeans, 'annualAI_Mean')
b <- plot_envar_vs_dist(df_all_popMeans, 'annualPET_Mean')
c <- plot_envar_vs_dist(df_all_popMeans, 'DEM_Mean')
d <- plot_envar_vs_dist(df_all_popMeans, 'GMIS_Mean')
e <- plot_envar_vs_dist(df_all_popMeans, 'NDSI_Mean')
f <- plot_envar_vs_dist(df_all_popMeans, 'summerLST_Mean')
g <- plot_envar_vs_dist(df_all_popMeans, 'summerNDVI_Mean')
h <- plot_envar_vs_dist(df_all_popMeans, 'winterLST_Mean')
i <- plot_envar_vs_dist(df_all_popMeans, 'winterNDVI_Mean')

### Combine panels for figure 2
figureS2 <- (a | b | c) / (d | e | f) / (h | g | i) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom', 
        legend.direction="horizontal",
        legend.text = element_text(size=20), 
        legend.key = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.spacing.x = unit(0.5, "cm"),
        plot.tag.position = c(0.1, 1.1),
        plot.tag = element_text(size = 20)) 
figureS2

ggsave(filename = "analysis/figures/supplemental/figS2_enVar_vs_distance.pdf", plot = figureS2, 
       device = "pdf", width = 16, height = 15, units = "in", dpi = 600, useDingbats = FALSE) 
  
###################
#### FIGURE S3 ####
###################

### Figure S3A

# Extract PCA object from multi-dispersion analysis
# PCA is performed on the mean multivariate environmental dispersion among urban and rural pops across all cities
enviroVariancePCA <- enviroVariance$pca_stand
PC1_enviroVariance_varEx <- round(enviroVariancePCA$eig[1, 2], 1)  # Percent variance explained by PC1
PC2_enviroVariance_varEx <- round(enviroVariancePCA$eig[2, 2], 1)  # Percent variance explained by PC2

enviroVariancePCA_sites  <- data.frame(enviroVariancePCA$ind$coord) %>%
  dplyr::select(Dim.1, Dim.2) %>% 
  mutate(habitat = enviroVariancePCA$groups,
         city = enviroVariancePCA$city)

# Colors for PCA biplot
pal <- wes_palette('Darjeeling1', 5, type = 'discrete')
urban_col <- pal[4]
rural_col <- pal[2]
cols <- c(urban_col, rural_col)

# Create PCA biplot. grouped by habitat
enviroVariance_PCA_plot <- ggplot(enviroVariancePCA_sites, aes(x = Dim.1, y = Dim.2)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_line(aes(group = city), alpha = 0.7) +
  geom_point(size = 2.75, shape = 21, colour = "black", aes(fill =  habitat)) +
  stat_ellipse(aes(colour = habitat), size = 1.5, level = 0.95) +
  xlab(sprintf("PC1 (%.1f%%)", PC1_enviroVariance_varEx)) + ylab(sprintf("PC1 (%.1f%%)", PC2_enviroVariance_varEx)) +
  scale_x_continuous(breaks = seq(from = -600, to = 1400, by = 400)) +
  scale_y_continuous(breaks = seq(from = -400, to = 800, by = 200)) +
  scale_colour_manual(values = rev(cols)) +
  scale_fill_manual(values = rev(cols)) +
  ng1 + theme(legend.position = "top", 
              legend.direction="horizontal",
              legend.text = element_text(size=15), 
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"))
enviroVariance_PCA_plot

ggsave(filename = "analysis/figures/supplemental/figureS3A_enviroVariancePCA_withLinesAndHulls.pdf", 
       plot = enviroVariance_PCA_plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)


### Figure S3B

## Eigenvectors of how environmental factors are associated with urban/rural habitats

# Extract RDA1 and PC1 species scores
enviroVariancePCA_vars  <- data.frame(enviroVariancePCA$var$coord) %>% 
  dplyr::select(Dim.1, Dim.2) %>% 
  cbind(., enviroVariancePCA$var$contrib %>% 
          as.data.frame() %>% 
          dplyr::select(Dim.1) %>% 
          rename("contrib" = "Dim.1"))

pal <- wes_palette("Darjeeling1", 3, type = "continuous")
enviroVariancePCA_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = enviroVariancePCA_vars, aes(x = 0, xend = Dim.1, y=0, yend = Dim.2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = enviroVariancePCA_vars,
            aes(x = Dim.1, y = Dim.2, label = rownames(enviroVariancePCA_vars),
                hjust = "inward", vjust =  0.5 * (1 - sign(Dim.1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", PC1_enviroVariance_varEx)) + ylab(sprintf("PC2 (%.1f%%)", PC2_enviroVariance_varEx)) +
  scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
  scale_x_continuous(breaks = seq(from = 0, to = 200, by = 50)) +
  scale_y_continuous(breaks = seq(from = -50, to = 100, by = 25)) +
  ng1 + theme(legend.position = "top",
              legend.direction="horizontal",
              # legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text = element_text(size=10)) +
  guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
enviroVariancePCA_variableContrib

ggsave(filename = "analysis/figures/supplemental/figureS3B_enviroVariancePCA_eigenvectorsOnly.pdf", 
       plot = enviroVariancePCA_variableContrib, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

## Figure S3C

# BoxM plot with log determinant
pdf('analysis/figures/supplemental/figureS3C_boxM.pdf', width = 6, height = 6, useDingbats = FALSE)
plot(enviroVariance_boxM)
dev.off()

###################
#### FIGURE S7 ####
###################

# Byplots showing variable loadings on first two PCs for emvironmental mean PCA (A) and 
# environmental slopes PCA (B)

## Figure S7A

# Get percent variance of first two PCs
pca_enviroMeans_eig <- pca_enviroMeans$CA$eig
pca_enviroMeans_percent_var <- pca_enviroMeans_eig * 100 / sum(pca_enviroMeans_eig)
pca_enviroMeans_eig_PC1_varEx <- round(pca_enviroMeans_percent_var[1], 1)  # Percent variance explained by PC1
pca_enviroMeans_eig_PC2_varEx <- round(pca_enviroMeans_percent_var[2], 1)  # Percent variance explained by PC2

# Calculate contribution of each environmental variable to PC2
# https://stackoverflow.com/questions/50177409/how-to-calculate-species-contribution-percentages-for-vegan-rda-cca-objects
pca_enviroMeans_contrib <- round(100*scores(pca_enviroMeans, display = "species", scaling = 0)[,1]^2, 3)

# Extract PC1 and PC2 species scores and bind contributions
pca_enviroMeans_vars  <- scores(pca_enviroMeans, display = 'species', choices = c(1, 2), scaling = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'var') %>% 
  mutate(contrib = pca_enviroMeans_contrib)

# Plot
pal <- wes_palette("Darjeeling1", 3, type = "continuous")
pca_enviroMeans_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = pca_enviroMeans_vars, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = pca_enviroMeans_vars,
            aes(x = PC1, y = PC2, label = var,
                hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", pca_enviroMeans_eig_PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", pca_enviroMeans_eig_PC2_varEx)) +
  scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
  # scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  # scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  ng1 + theme(legend.position = "top",
              legend.direction="horizontal",
              # legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text = element_text(size=10)) +
  guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
pca_enviroMeans_variableContrib

## Figure S7B

# Get percent variance of first two PCs
pca_enviroSlopes_eig <- pca_enviroSlopes$CA$eig
pca_enviroSlopes_percent_var <- pca_enviroSlopes_eig * 100 / sum(pca_enviroSlopes_eig)
pca_enviroSlopes_eig_PC1_varEx <- round(pca_enviroSlopes_percent_var[1], 1)  # Percent variance explained by PC1
pca_enviroSlopes_eig_PC2_varEx <- round(pca_enviroSlopes_percent_var[2], 1)  # Percent variance explained by PC2

# Calculate contribution of each environmental variable to PC2
# https://stackoverflow.com/questions/50177409/how-to-calculate-species-contribution-percentages-for-vegan-rda-cca-objects
pca_enviroSlopes_contrib <- round(100*scores(pca_enviroSlopes, display = "species", scaling = 0)[,1]^2, 3)

# Extract PC1 and PC2 species scores and bind contributions
pca_enviroSlopes_vars  <- scores(pca_enviroSlopes, display = 'species', choices = c(1, 2), scaling = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'var') %>% 
  mutate(contrib = pca_enviroSlopes_contrib)

# Plot
pal <- wes_palette("Darjeeling1", 3, type = "continuous")
pca_enviroSlopes_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = pca_enviroSlopes_vars, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = pca_enviroSlopes_vars,
            aes(x = PC1, y = PC2, label = var,
                hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", pca_enviroSlopes_eig_PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", pca_enviroSlopes_eig_PC2_varEx)) +
  scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
  # scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  # scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  ng1 + theme(legend.position = "top",
              legend.direction="horizontal",
              # legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text = element_text(size=10)) +
  guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
pca_enviroSlopes_variableContrib

# Combine figures
figureS7 <- pca_enviroMeans_variableContrib + pca_enviroSlopes_variableContrib +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag.position = c(0.1, 0.95),
        plot.tag = element_text(size = 20))
figureS7

ggsave(filename = "analysis/figures/supplemental/figS7_enviroPCAs_loadings.pdf", plot = figureS7, 
       device = "pdf", width = 16, height = 7, units = "in", dpi = 600, useDingbats = FALSE)
