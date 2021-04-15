# Script to create supplementary figures
#
# Author: James S. Santangelo

###################
#### FIGURE SX ####
###################

## BoxM plot with log determinant
pdf('analysis/figures/supplemental/figureSX_boxM.pdf', width = 6, height = 6, useDingbats = FALSE)
plot(enviroVariance_boxM)
dev.off()

###################
#### FIGURE SX ####
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
    mutate(scaled = scale(!!sym(response),
                          center = min(!!sym(response)),
                          scale = max(!!sym(response)) - min(!!sym(response))))
  
  # Create plot
  plot <- df %>%
    ggplot(., aes(x = std_distance, y = scaled)) +
    geom_line(stat = "smooth", method=function(formula,data,weights=weight) rlm(formula,
                                                                 data,
                                                                 weights=weight,
                                                                 maxit=200),
                aes(color = color, alpha = color, group = city),
                size = 0.75, se = FALSE) +
    geom_line(stat = "smooth", method="lm", colour = "black", size = 2.5) +
    xlab("Standardized distance") + ylab(ylab) +
    scale_colour_manual(values = pal, limits = c('Not significant', 'Significantly negative', 'Significantly positive')) +
    scale_alpha_manual(values = c(0.5, 0.5, 0.5), limits = c('Not significant', 'Significantly negative', 'Significantly positive')) +
    scale_y_continuous(breaks = seq(from = -0.4, to = 1.2, by = 0.2)) +
    scale_x_continuous(breaks = seq(from = 0, to = 1.1, by = 0.25)) +
    coord_cartesian(ylim = c(-0.4, 1.2), xlim = c(0, 1), clip = 'off') +
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
figureSX <- (a | b | c) / (d | e | f) / (g | h | i) +
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
figureSX

ggsave(filename = "analysis/figures/supplemental/figSX_enVar_vs_distance.pdf", plot = figureSX, 
       device = "pdf", width = 16, height = 14, units = "in", dpi = 600, useDingbats = FALSE) 
  