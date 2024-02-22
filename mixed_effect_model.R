# Import required packages
library("readxl")
library("lme4")
library("tidyverse")
library("xlsx")
library("lmerTest") 
library("sjPlot")
library("dplyr")
library("Matrix")
library("writexl")


# Read database:
npc <- read_excel("input_npc_yield_gap_data.xlsx", sheet = "organic_crop_specific")
#npc <- read_excel("input_npc_yield_gap_data.xlsx", sheet = "previous")


# Remove crops with less than 12 observations:
count_obs <- npc %>% group_by(crop) %>% summarise(total_count=n(),.groups = 'drop')
to_remove <- count_obs %>% filter(total_count < 12)
cat("The following crop variable(s) are being removed due to too little observations: ",paste(to_remove$crop, collapse = ", "), ".","\n")
to_remove <- list(to_remove$crop)
npc_clean <- subset(npc, !(crop %in% to_remove))

# Print the crops remaining, which will be fed into the mixed effect model:
cat("The crops remaining are: ", paste(unique(npc_clean$crop), collapse = ", "), ".", "\n")

# Drop all regions for which no median NPC index was establish due to potentially the extent of the NPC map:
npc_clean <- npc_clean %>% drop_na(median_biogeo_norm)
# All regions have LF-NPC as the total amount of observation remains the same.

#---------------------------- ESTIMATING THE CONTRIBUTION OF NPC TO YIELD-------------------------------------------
# 1.Run normal linear model:
npc_lm = lm(estimate ~ median_biogeo_norm, data = npc_clean)
summary(npc_lm)

#-----------------------------------------------------------------------
# 2.Run mixed effect model (only intercept):
npc_mixed = lmer(estimate ~ median_biogeo_norm + (1 | crop), data = npc_clean)

# Summarize results:
summary(npc_mixed)
#r.squaredGLMM(npc_mixed)
tab_model(npc_mixed)

# Make predictions and add them back to the data frame:
npc_clean$prediction <- predict(npc_mixed)

# Look at and extract random effect (re) of intercept for each crop:
ranef(npc_mixed)$crop %>% head(10)
re <- ranef(npc_mixed)$crop
re <- cbind(crop = rownames(re), re)

# Look at extract intercept and coefficient:
coef(npc_mixed)$crop %>% head(10)
fe <- coef(npc_mixed)$crop
fe <- cbind(crop = rownames(fe), fe)

# Write to disk prediction as well as estimation of fixed and mixed effect:
write_xlsx(list("prediction" = npc_clean,
                "random_effect" = re,
                "intercept_ceofficient" = fe),
                "output_mixed_effect.xlsx")

#-----------------------------------------------------------------------
# 3.Run mixed effect model (full - i.e. intercept and slope are allowed to vary):
npc_mixed_full =  lmer(estimate ~ median_biogeo_norm + (1 + median_biogeo_norm | crop), data = npc_clean)

# Summarize results:
summary(npc_mixed_full)
r.squaredGLMM(npc_mixed_full)
tab_model(npc_mixed_full)

# Make predictions and add them back to the data frame:
predict_with_re = predict(npc_mixed_full)
npc_clean$prediction_full <- predict(npc_mixed_full)

# Extract random effect (for intercept) of each crop:
ranef(npc_mixed_full)$crop %>% head(10)
re_full <- ranef(npc_mixed_full)$crop
re_full <- cbind(crop = rownames(re_full), re_full)

# Look at intercept and coefficient:
coef(npc_mixed_full)$crop %>% head(10)
fe_full <- coef(npc_mixed_full)$crop
fe_full <- cbind(crop = rownames(fe_full), fe_full)

# Write to disk prediction as well as estimation of fixed and mixed effect:
write_xlsx(list("prediction" = npc_clean,
                "random_effect" = re_full,
                "intercept_ceofficient" = fe_full),
                 "output_mixed_effect.xlsx")

#----------Model diagnostics----------------------------------------------

# Print comparison of both mixed effect model:
tab_model(npc_mixed, npc_mixed_full)

# Error calculation:
npc_clean$error<-(npc_clean$estimate-npc_clean$prediction)

# Inspection of error: QQ Plot
qqnorm(npc_clean$error, main="QQ Plot of Errors")
qqline(npc_clean$error)

# Inspection of error: Histogram
hist(npc_clean$error, prob = TRUE, bins = 22, breaks = 15, main = "Distribution of Errors")  
lines(density(npc_clean$error), col = "red") 

#Tukey-Anscombe
plot(npc_mixed, type = c("p", "smooth"), col.line = "black")

#-----------Further analysis------------------------------------------
# Alternative model specification:
#explore_lmer = lmer(estimate ~ I(nobs_organic/nobs)+biogeo +median_biogeo_norm + (1+ median_biogeo_norm | crop), data = npc_clean)
#explore_lmer = lmer(estimate ~ I(nobs_organic/nobs) +median_biogeo_norm + (1+ median_biogeo_norm | crop), data = npc_clean)
explore_lmer =  lmer(estimate ~ biogeo+median_biogeo_norm + (1 + median_biogeo_norm | crop), data = npc_clean)
car::Anova(explore_lmer)
summary(explore_lmer)
tab_model(explore_lmer)

npc_clean$I_nobs_ratio <- npc_clean$nobs_organic / npc_clean$nobs
hist(npc_clean$I_nobs_ratio, main = "Histogram of I(nobs_organic/nobs)", xlab = "I(nobs_organic/nobs)")

plot(npc_clean$median_biogeo_norm ~ npc_clean$I_nobs_ratio, 
     main = "Scatterplot of Residuals vs. I_nobs_ratio",
     xlab = "I_nobs_ratio",
     ylab = "Residuals")

# Boxplot showing yield gap distribution per biogeographic region:
ggplot(npc_clean, aes(x = as.factor(biogeo), y = estimate, fill = as.factor(biogeo))) +
  geom_boxplot() +
  labs(x = "", y = "Yield gap") +
  ggtitle("Yield Gap by Biogeographic Region") +
  scale_fill_discrete(name = "Biogeographic Region") +
  scale_y_continuous(breaks = seq(-60, max(npc_clean$estimate), 10)) +  # Set breaks every 10 steps
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

