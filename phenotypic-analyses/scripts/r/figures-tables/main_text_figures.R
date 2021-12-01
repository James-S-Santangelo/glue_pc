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
          axis.ticks=element_line(size = 1, color="black"),
          axis.ticks.length=unit(0.25, 'cm'),
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

# Figure 1 is a map with insets and was created in QGIS

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

ggsave(filename = "analysis/figures/main_text/figure2/figure2A_enviroPCA_withLinesAndHulls.pdf",
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

ggsave(filename = "analysis/figures/main_text/figure2/figure2B_enviroPCA_eigenvectorsOnly.pdf", 
       plot = enviroPCA_variableContrib, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, 
       useDingbats = FALSE)

### Figure 2C

## Histogram showing distribution of slopes of clines, with different colors for significantly negative, positive, and no cline.

pal <- c("#909090", "#FF0000", "#046C9A")
logOddsHistogram <- linearClineTable_mod %>% 
  mutate(Significance = case_when(
    betaLog_Dist < 0 & pvalLog_Dist < 0.05 ~ "Significantly negative",
    betaLog_Dist > 0 & pvalLog_Dist < 0.05 ~ "Significantly positive",
    TRUE ~ "Not significant"
  )) %>% 
  ggplot(., aes(x = betaLog_Dist, fill = Significance)) +
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
  geom_vline(xintercept = mean(linearClineTable_mod$betaLog_Dist), linetype = "dashed", size = 1) +
  xlab("Standardized slope of cline") + ylab("Count") +
  coord_cartesian(ylim = c(0, 21)) +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 5), expand = c(0, 0)) +
  scale_fill_manual(values = pal) +
  ng1 + theme(legend.position = "top", 
              legend.direction="horizontal",
              legend.text = element_text(size=15), 
              legend.key = element_rect(fill = "white"),
              legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.5, "cm"))
logOddsHistogram

ggsave(filename = "analysis/figures/main_text/figure2/figure2C_logOddsHistogram.pdf", plot = logOddsHistogram,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

### Figure 2D

## Binomial regressions with HCN vs distance for all cities, colored by whether relationship 
## is positive or negative. Set alpha = 0.5 on lines to better see overlap and avoid clutter.

HCN_by_city <- df_all_popMeans %>%
  left_join(., linearClineTable_mod, by = "city") %>% 
  mutate(significant = ifelse(pvalLog_Dist < 0.05, "Yes", "No"),
         direction = ifelse(betaLog_Dist > 0, "Positive", "Negative"),
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
  coord_cartesian(ylim = c(-0.01, 1.01), xlim = c(0, 1)) +
  ng1 
HCN_by_city 

ggsave(filename = "analysis/figures/main_text/figure2/figure2D_clineByCity.pdf", plot = HCN_by_city,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

## Figure 2E, F, and G

# These panels are maps for individual cities and were manually generated in QGIS

## Figure 2H, I, and J

temuco <- df_all_popMeans %>% filter(city == "Temuco")  # No cline
freehold <- df_all_popMeans %>% filter(city == "Freehold")  # Negative cline
muenster <- df_all_popMeans %>% filter(city == "Muenster")  # Positive cline

temuco_plot <- fig2_hij(temuco)
freehold_plot <- fig2_hij(freehold)
muenster_plot <- fig2_hij(muenster)

ggsave(filename = "analysis/figures/main_text/figure2/figure2H_Muenster_HCN_vs_distance.pdf", plot = muenster_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)
ggsave(filename = "analysis/figures/main_text/figure2/figure2I_Freehold_HCN_vs_distance.pdf", plot = freehold_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)
ggsave(filename = "analysis/figures/main_text/figure2/figure2J_Temuco_HCN_vs_distance.pdf", plot = temuco_plot,
       device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

##################
#### FIGURE 3 ####
##################

## Figure 3 is a genomics figure showing urban-rural differences in pairwise nucleotide diversity and
##  differentiation. This figure is generated using Jupyter Notebooks in the genomic analyses repo 
## See ../../../genomic-analyses/notebooks/diversity_populationStructure.r.ipynb

##################
#### FIGURE 4 ####
##################

# Convert model matrix to data frame
model_matrix_df <- as.data.frame(model_matrix)
  
## Figure 6A

# PET_Slope by Mean summer NDVI

# Run model fit with OLS for plotting
mod <- lm(betaLog ~ summerNDVI_Mean * annualPET_Slope, data = model_matrix_df)

# Define low, mean, and high categories
high <- round(mean(model_matrix_df$summerNDVI_Mean) + sd(model_matrix_df$summerNDVI_Mean), 1)
mean <- round(mean(model_matrix_df$summerNDVI_Mean), 1)
low <- round(mean(model_matrix_df$summerNDVI_Mean) - sd(model_matrix_df$summerNDVI_Mean), 1)

# Get dataframe with predicted lines from model for effects  
range(model_matrix_df$annualPET_Slope) # Used to parameterize x-axis range in list below. 
val_list <- list(annualPET_Slope = seq(from = -7.5, to = 3, by = 0.1), summerNDVI_Mean = c(low, mean, high))
pred_df <- emmip(mod, summerNDVI_Mean~annualPET_Slope, at = val_list, CIs=TRUE, plotit=FALSE)

# Color palette
low_col = wes_palette("Rushmore1", 5, type = "discrete")[4]
high_col = wes_palette("Darjeeling1", 5, type = "discrete")[3]
cols = c(low_col, high_col)

pred_df$fsummerNDVI_Mean <- factor(pred_df$summerNDVI_Mean)
figure4A <- pred_df %>% 
  filter(fsummerNDVI_Mean != '0') %>% 
  ggplot(., aes(x = annualPET_Slope, y = yvar)) + 
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = fsummerNDVI_Mean), alpha = 0.25) +
  geom_line(aes(color = fsummerNDVI_Mean), size = 2) +
  ylab("Predicted HCN slope") + xlab("Slope of annual PET") +
  scale_color_manual(values = rev(cols), name = "Mean summer NDVI", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_fill_manual(values = rev(cols), name = "Mean summer NDVI", labels = c("low (-1 sd)","high (+1 sd)")) +
  # coord_cartesian(xlim = c(-3, 2)) +
  scale_y_continuous(breaks = seq(from = -2, to = 4, by = 2)) +
  scale_x_continuous(breaks = seq(from = -7.5, to = 3, by = 1.5)) +
  ng1 + theme(legend.position = 'right')
figure4A

ggsave(filename = "analysis/figures/main_text/figure4/figure4A_betaLog_by_sNDVImean_PETslope.pdf", 
       plot = figure4A, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

### Figure 6B

## Log Odds vs winter NDVI mean

logOdds_by_winterNDVI <- ggplot(model_matrix_df, aes(x = winterNDVI_Mean, y = betaLog)) +
  geom_point(size = 3.5, color = 'black') +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', size = 1.5) +
  xlab('Mean winter NDVI') + ylab('Slope of HCN cline') +
  scale_y_continuous(breaks = seq(from = -1, to = 5, by = 2)) +
  # coord_cartesian(xlim = c(-1.7, 2), ylim = c(-1.4, 6)) +
  ng1
logOdds_by_winterNDVI

ggsave(filename = "analysis/figures/main_text/figure4/figure4B_logOdds_by_winterNDVI.pdf", 
       plot = logOdds_by_winterNDVI, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

## Figure 6C

# PET_Slope by Mean GMIS

# Run model fit with OLS for plotting
mod <- lm(betaLog ~ GMIS_Mean * annualPET_Slope, data = model_matrix_df)

# Define low, mean, and high categories
high <- round(mean(model_matrix_df$GMIS_Mean) + sd(model_matrix_df$GMIS_Mean), 1)
mean <- round(mean(model_matrix_df$GMIS_Mean), 1)
low <- round(mean(model_matrix_df$GMIS_Mean) - sd(model_matrix_df$GMIS_Mean), 1)

# Get dataframe with predicted lines from model for effects  
range(model_matrix_df$annualPET_Slope) # Used to parameterize x-axis range in list below. 
val_list <- list(annualPET_Slope = seq(from = -7.5, to = 3, by = 0.1), GMIS_Mean = c(low, mean, high))
pred_df <- emmip(mod, GMIS_Mean~annualPET_Slope, at = val_list, CIs=TRUE, plotit=FALSE)

pred_df$fGMIS_Mean <- factor(pred_df$GMIS_Mean)
figure4C <- pred_df %>% 
  filter(fGMIS_Mean != '0') %>% 
  ggplot(., aes(x = annualPET_Slope, y = yvar)) + 
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = fGMIS_Mean), alpha = 0.25) +
  geom_line(aes(color = fGMIS_Mean), size = 2) +
  ylab("Predicted HCN slope") + xlab("Slope of annual PET") +
  scale_color_manual(values = rev(cols), name = "Mean GMIS", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_fill_manual(values = rev(cols), name = "Mean GMIS", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_y_continuous(breaks = seq(from = -3, to = 3, by = 1.5)) +
  scale_x_continuous(breaks = seq(from = -7.5, to = 3, by = 1.5)) +
  ng1 + theme(legend.position = 'right')
figure4C

ggsave(filename = "analysis/figures/main_text/figure4/figure4C_betaLog_by_GMISmean_PETslope.pdf", 
       plot = plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

## Figure 6D

# Summer NDVI Slope by AI Mean

# Run model fit with OLS for plotting
mod <- lm(betaLog ~ annualAI_Mean * summerNDVI_Slope, data = model_matrix_df)

# Define low, mean, and high categories
high <- round(mean(model_matrix_df$annualAI_Mean) + sd(model_matrix_df$annualAI_Mean), 1)
mean <- round(mean(model_matrix_df$annualAI_Mean), 1)
low <- round(mean(model_matrix_df$annualAI_Mean) - sd(model_matrix_df$annualAI_Mean), 1)

# Get dataframe with predicted lines from model for effects  
range(model_matrix_df$summerNDVI_Slope) # Used to parameterize x-axis range in list below. 
val_list <- list(summerNDVI_Slope = seq(from = -2.5, to = 3, by = 0.1), annualAI_Mean = c(low, mean, high))
pred_df <- emmip(mod, annualAI_Mean~summerNDVI_Slope, at = val_list, CIs=TRUE, plotit=FALSE)

pred_df$fannualAI_Mean <- factor(pred_df$annualAI_Mean)
figure4D <- pred_df %>% 
  filter(fannualAI_Mean != '0') %>% 
  ggplot(., aes(x = summerNDVI_Slope, y = yvar)) + 
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = fannualAI_Mean), alpha = 0.25) +
  geom_line(aes(color = fannualAI_Mean), size = 2) +
  ylab("Predicted HCN slope") + xlab("Slope of summer NDVI") +
  scale_color_manual(values = rev(cols), name = "Mean annual AI", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_fill_manual(values = rev(cols), name = "Mean annual AI", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_y_continuous(breaks = seq(from = -1.2, to = 2, by = 0.4), labels = scales::comma_format(accuracy = 0.1)) +
  ng1 + theme(legend.position = 'right')
figure4D

ggsave(filename = "analysis/figures/main_text/figure4/figure4D_betaLog_by_AImean_sNDVImean.pdf", 
       plot = plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)

