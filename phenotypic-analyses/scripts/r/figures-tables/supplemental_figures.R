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
