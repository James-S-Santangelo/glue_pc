# Script to create main text figures and tables

###############
#### SETUP ####
###############

# Theme used for plotting
ng1 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1),
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black",size=15),
          axis.title=element_text(color="black",size=1),
          axis.title.y=element_text(vjust=2,size=17),
          axis.title.x=element_text(vjust=0.1,size=17),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
          strip.background = element_rect(colour="black"),
          legend.position = "top", legend.direction="vertical",
          legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
          legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))

##################
#### FIGURE 1 ####
##################

# Figure 1 is a map with insets that will be created in QGIS

##################
#### FIGURE 2 ####
##################

### Figure 2A

eig <- enviroPCA$CA$eig
percent_var <- eig * 100 / sum(eig)
PC1_varEx <- round(percent_var[1], 1)  # Percent variance explained by PC1
PC2_varEx <- round(percent_var[2], 1)  # Percent variance explained by PC2

# Get coordinates of urban/rural populations by city. Add habitat and city column for grouping
enviroPCA_sites  <- scores(enviroPCA, display = 'sites', choices = c(1, 2), scaling = 1) %>%
  as.data.frame() %>% 
  dplyr::select(PC1, PC2) %>% 
  mutate(habitat = habitat,
         city = rep(city_names, 2))

# Colors for PCA biplot
pal <- wes_palette('Darjeeling1', 5, type = 'discrete')
urban_col <- pal[4]
rural_col <- pal[2]
cols <- c(urban_col, rural_col)

# Create PCA biplot. grouped by habitat
enviroPCA_plot <- ggplot(enviroPCA_sites, aes(x = PC1, y = PC2)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_line(aes(group = city), alpha = 0.7) +
  geom_point(size = 2.75, shape = 21, colour = "black", aes(fill =  habitat)) +
  stat_ellipse(aes(colour = habitat), size = 1.5, level = 0.95) +
  xlab(sprintf("PC1 (%.1f%%)", PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", PC2_varEx)) +
  scale_colour_manual(values = rev(cols)) +
  scale_fill_manual(values = rev(cols)) +
  scale_x_continuous(breaks = seq(from = -0.6, to = 0.6, by = 0.2), labels = scales::comma) +
  scale_y_continuous(breaks = seq(from = -0.45, to = 0.45, by = 0.15), labels = scales::comma) +
  ng1 + theme(legend.position = "top", 
              legend.direction="horizontal",
              legend.text = element_text(size=15), 
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"))
enviroPCA_plot

ggsave(filename = "analysis/figures/manuscript-panels/figure-2/figure2A_enviroPCA_withLinesAndHulls.pdf",
       plot = enviroPCA_plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)


### Figure 2B

## Eigenvectors of how environmental factors are associated with urban/rural habitats

# Calculate contribution of each environmental variable to PC2
# https://stackoverflow.com/questions/50177409/how-to-calculate-species-contribution-percentages-for-vegan-rda-cca-objects
contrib <- round(100*scores(enviroPCA, display = "sp", scaling = 0)[,2]^2, 3)

# Extract RDA1 and PC1 species scores and bind contributions
enviroPCA_vars  <- scores(enviroPCA, display = 'species', choices = c(1, 2), scaling = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'var') %>% 
  mutate(contrib = contrib)

# Plot
pal <- wes_palette("Darjeeling1", 3, type = "continuous")
enviroPCA_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = enviroPCA_vars, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = enviroPCA_vars,
            aes(x = PC1, y = PC2, label = var,
                hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", PC2_varEx)) +
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
enviroPCA_variableContrib

ggsave(filename = "analysis/figures/manuscript-panels/figure-2/figure2B_enviroPCA_eigenvectorsOnly.pdf", 
       plot = enviroPCA_variableContrib, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, 
       useDingbats = FALSE)

### Figure 2C

# Extract PCA object from multi-dispersion analysis
# PCA is performed on the mean multivariate environmental dispersion among urban and rural pops across all cities
enviroVariancePCA <- enviroVariance$pca_stand
PC1_enviroVariance_varEx <- round(enviroVariancePCA$eig[1, 2], 1)  # Percent variance explained by PC1
PC2_enviroVariance_varEx <- round(enviroVariancePCA$eig[2, 2], 1)  # Percent variance explained by PC2

enviroVariancePCA_sites  <- data.frame(enviroVariancePCA$ind$coord) %>%
  dplyr::select(Dim.1, Dim.2) %>% 
  mutate(habitat = enviroVariancePCA$groups,
         city = enviroVariancePCA$city)

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

ggsave(filename = "analysis/figures/manuscript-panels/figure-2/figure2C_enviroVariancePCA_withLinesAndHulls.pdf", 
       plot = enviroVariance_PCA_plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)


### Figure 2D

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

ggsave(filename = "analysis/figures/manuscript-panels/figure-2/figure2D_enviroVariancePCA_eigenvectorsOnly.pdf", 
       plot = enviroVariancePCA_variableContrib, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

### Combine panels for figure 2
figure2 <- enviroPCA_plot + enviroPCA_variableContrib + enviroVariance_PCA_plot + enviroVariancePCA_variableContrib +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag.position = c(0.05, 0.95),
        plot.tag = element_text(size = 20))
figure2

ggsave(filename = "analysis/figures/manuscript-panels/figure-2/figure2.pdf", plot = figure2, 
       device = "pdf", width = 16, height = 14, units = "in", dpi = 600, useDingbats = FALSE)

##################
#### FIGURE 3 ####
##################

### Figure 3A

## Histogram showing distribution of slopes of clines, with different colors for significantly negative, positive, and no cline.

# Get betas and pvalues from robust regressions
hcnClinesSummary <- df_all_popMeans %>% group_split(city) %>% 
  purrr::map_dfr(., logistic_regression_stats)

pal <- c("#909090", "#FF0000", "#046C9A")
logOddsHistogram <- hcnClinesSummary %>% 
  mutate(Significance = case_when(
    betaLog < 0 & pvalLog < 0.05 ~ "Significantly negative",
    betaLog > 0 & pvalLog < 0.05 ~ "Significantly positive",
    TRUE ~ "Not significant"
  )) %>% 
  ggplot(., aes(x = betaLog, fill = Significance)) +
  geom_histogram(data = . %>% filter(Significance == 'Significantly negative'), 
                 bins = 50,
                 color = 'black') +
  geom_histogram(data = . %>% filter(Significance == 'Significantly positive'), 
                 bins = 50,
                 color = 'black') +
  geom_histogram(data = . %>% filter(Significance == 'Not significant'), 
                 bins = 50,
                 color = 'black',
                 alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = mean(hcnClinesSummary$betaLog, linetype = "dashed", size = 1)) +
  xlab("Standardized slope of cline") + ylab("Count") +
  scale_fill_manual(values = pal) +
  ng1 + theme(legend.position = "top", 
              legend.direction="horizontal",
              legend.text = element_text(size=15), 
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.5, "cm"))
logOddsHistogram

ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3A_logOddsHistogram.pdf", plot = logOddsHistogram,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

### Figure 3B

## Reaction norm plot with HCN vs distance for all cities, colored by whether relationship 
## is positive or negative. Use solid line if the relationship is significant and dashed lines if not. 
## Set alpha = 0.5 on lines to better see overlap and avoid clutter.

HCN_by_city <- df_all_popMeans %>%
  left_join(., hcnClinesSummary, by = "city") %>% 
  mutate(significant = ifelse(pvalLog < 0.05, "Yes", "No"),
         direction = ifelse(betaLog > 0, "Positive", "Negative"),
         color = case_when(significant == "Yes" & direction == "Positive" ~ "Significantly positive",
                           significant == "Yes" & direction == "Negative" ~ "Significantly negative",
                           TRUE ~ "Not significant")) %>%
  ggplot(., aes(x = std_distance, y = freqHCN, weight = total_plants)) +  
  geom_line(stat = "smooth", 
            method="glm", 
            aes(color = color, alpha = color, group = city),
            size = 0.5, 
            show.legend = FALSE,
            method.args = list(family = "binomial")) +
  geom_line(stat = "smooth", 
            method="glm", 
            colour = "black", 
            size = 2.5,
            method.args = list(family = "binomial")) +
  xlab("Standardized distance") + ylab("Frequency of HCN") + 
  scale_colour_manual(values = pal) +
  scale_alpha_discrete(range = c(0.5, 0.8)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
  scale_x_continuous(breaks = seq(from = 0, to = 1.1, by = 0.25)) +
  coord_cartesian(ylim = c(-0.05, 1.05), xlim = c(0, 1), clip = 'off') +
  ng1 
HCN_by_city

ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3B_clineByCity.pdf", plot = HCN_by_city,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

## Figure 3C, D, and E

temuco <- df_all_popMeans %>% filter(city == "Temuco")  # No cline
freehold <- df_all_popMeans %>% filter(city == "Freehold")  # Negative cline
muenster <- df_all_popMeans %>% filter(city == "Muenster")  # Positive cline

#' Function to plot panels C, D, and E for figure 3
#' 
#' @param df Population-mean HCN frequency dataframe
#' 
#' @return ggplot object
fig3_inset_biplot <- function(df){
  
  city <- df %>% pull(city) %>% unique()
  col <- case_when(city == 'Muenster' ~ pal[3],
                   city == 'Freehold' ~ pal[2])
  
  plot <- ggplot(df, aes(x = std_distance, y = freqHCN, weight = total_plants)) +
    geom_point(size = 3.5) +
    geom_smooth(method = 'glm', 
                color = ifelse(city == "Temuco", "black", col),
                fill = ifelse(city == "Temuco", "grey", col),
                method.args = list(family = "binomial"),
                size = 1.5) +
    xlab("Standardized distance") +
    ylab("Frequency of HCN") +
    ng1
  
  return(plot)
}

temuco_plot <- fig3_inset_biplot(temuco)
freehold_plot <- fig3_inset_biplot(freehold)
muenster_plot <- fig3_inset_biplot(muenster)

ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3C_Muenster_HCN_vs_distance.pdf", plot = muenster_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)
ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3D_Freehold_HCN_vs_distance.pdf", plot = freehold_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)
ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3E_Temuco_HCN_vs_distance.pdf", plot = temuco_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

### Combine panels for figure 2
figure3 <- (logOddsHistogram + HCN_by_city) / (muenster_plot | freehold_plot | temuco_plot) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = c(1, 1.2),
        plot.tag.position = c(0.05, 1.05),
        plot.tag = element_text(size = 20))
figure3

ggsave(filename = "analysis/figures/manuscript-panels/figure-3/figure3.pdf", plot = figure3, 
       device = "pdf", width = 16, height = 16, units = "in", dpi = 600, useDingbats = FALSE)

##################
#### FIGURE 4 ####
##################

# Extracts columns from supplementary tables rather than regenerating
df_wNDVImean_NDSI_HCNslope <- final_table %>% 
  dplyr::select(city, betaRLM_freqHCN) %>% 
  left_join(., cityEnviroMeans %>% dplyr::select(city, Mean_winterNDVI, Mean_NDSI),
            by = 'city') %>% 
  filter(!(is.na(Mean_winterNDVI))) %>%  # Remove 2 cities with missing values
  filter(!(city == 'St_Albert'))  # St Albert not in Elastic Net due to outlier
  
### Figure 4A

## HCN Slope vs winter NDVI mean

HCNslope_by_winterNDVI <- ggplot(df_wNDVImean_NDSI_HCNslope, aes(x = Mean_winterNDVI, y = betaRLM_freqHCN)) +
  geom_point(size = 3.5, color = 'black') +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', size = 1.5) +
  xlab('Mean winter NDVI') + ylab('Slope of HCN cline') +
  scale_x_continuous(breaks = seq(from = 0, to = 0.6, by = 0.1)) +
  scale_y_continuous(breaks = seq(from = -0.4, to = 0.8, by = 0.2)) +
  ng1
HCNslope_by_winterNDVI

ggsave(filename = "analysis/figures/manuscript-panels/figure-4/figure4A_HCNslope_wNDVImean.pdf", 
       plot = HCNslope_by_winterNDVI, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)


### Figure 4B

## HCN Slope vs NDSI 

HCNslope_by_NDSI <- ggplot(df_wNDVImean_NDSI_HCNslope, aes(x = Mean_NDSI, y = betaRLM_freqHCN)) +
  geom_point(size = 3.5, color = 'black') +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', size = 1.5) +
  xlab('Mean NDSI') + ylab('Slope of HCN cline') +
  scale_x_continuous(breaks = seq(from = -0.5, to = 0.8, by = 0.2)) +
  scale_y_continuous(breaks = seq(from = -0.4, to = 0.8, by = 0.2)) +
  ng1
HCNslope_by_NDSI

ggsave(filename = "analysis/figures/manuscript-panels/figure-4/figure4B_HCNslope_NDSImean.pdf", 
       plot = HCNslope_by_NDSI, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

# NDVI_by_NDSI <- ggplot(df_wNDVImean_NDSI_HCNslope, aes(x = Mean_NDSI, y = Mean_winterNDVI)) +
#   geom_point(size = 3.5, color = 'black') +
#   geom_smooth(method = 'lm', se = TRUE, color = 'black', size = 1.5) +
#   xlab('Mean NDSI') + ylab('Mean winter NDVI') +
#   scale_x_continuous(breaks = seq(from = -0.5, to = 0.8, by = 0.2)) +
#   scale_y_continuous(breaks = seq(from = 0, to = 0.6, by = 0.1)) +
#   ng1
# NDVI_by_NDSI
# 
# ggsave(filename = "analysis/figures/manuscript-panels/figure-4/figure4B_HCNslope_wNDVImean.pdf", 
#        plot = HCNslope_by_NDSI, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)


### Figure 4C

# Define low, mean, and high AI categories
LSThigh <- round(mean(model_df_withMainEffects$summerLST_Mean) + sd(model_df_withMainEffects$summerLST_Mean), 1)
LSTmean <- round(mean(model_df_withMainEffects$summerLST_Mean), 1)
LSTlow <- round(mean(model_df_withMainEffects$summerLST_Mean) - sd(model_df_withMainEffects$summerLST_Mean), 1)

# Get dataframe with predicted lines from Elastic Net model for effects of winterNDVI_Slope on HCN slopes for each of 
# 3 levels summerLST_Mean
range(model_df_withMainEffects$winterNDVI_Slope) # Used to parameterize x-axis range in list below. 
LST_wNDVI_list <- list(winterNDVI_Slope = seq(from = -2, to = 3, by = 0.1), summerLST_Mean = c(LSTlow, LSTmean, LSThigh))
LST_wNDVI_df <- emmip(predClines_elasticNet_withMainEffects, summerLST_Mean~winterNDVI_Slope, at = LST_wNDVI_list, CIs=TRUE, plotit=FALSE)

# Color palette
LSTlow_col = wes_palette("Zissou1", 5, type = "discrete")[5]
LSTmean_col = wes_palette("Zissou1", 5, type = "discrete")[3]
LSThigh_col = wes_palette("Zissou1", 5, type = "discrete")[1]
cols = c(LSTlow_col, LSThigh_col)

LST_wNDVI_df$fsummerLST_Mean <- factor(LST_wNDVI_df$summerLST_Mean)
LST_wNDVI_plot <- LST_wNDVI_df %>% 
  filter(fsummerLST_Mean != '0') %>% 
  ggplot(., aes(x = winterNDVI_Slope, y = yvar)) + 
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = fsummerLST_Mean), alpha = 0.25) +
  geom_line(aes(color = fsummerLST_Mean), size = 2) +
  ylab("Predicted HCN slope") + xlab("Slope of winter NDVI") +
  scale_color_manual(values = rev(cols), name = "Mean summer LST", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_fill_manual(values = rev(cols), name = "Mean summer LST", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_y_continuous(breaks = seq(from = -0.4, to = 1.0, by = 0.2)) +
  ng1 + theme(legend.position = 'right')
LST_wNDVI_plot

ggsave(filename = "analysis/figures/manuscript-panels/figure-4/figure4C_HCNslope_wNDVIslope_sLST.pdf", 
       plot = LST_wNDVI_plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)



