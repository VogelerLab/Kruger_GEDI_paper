#------------------------------------------------------------------------------
#  Using EMDI for estimation
#------------------------------------------------------------------------------
library(emdi)
library(arrow)
library(sfarrow)
library(sf)
library(dplyr)
library(survey)
library(stringr)
library(ggplot2)
library(gridExtra)

pop_path <- "J:/projects/ECOFOR/gedi/sae/sae_gedi_pop.parquet"
samp_path <- "J:/projects/ECOFOR/gedi/sae/sae_gedi_samp.gpkg"

pop <- data.frame(read_parquet(pop_path))
samp <- data.frame(st_read(samp_path, int64_as_string = TRUE))

colnames(samp) = gsub("_normal", "", colnames(samp)) # fhd_normal -> fhd

# Select subset of AOIs presented in case studies
subaois <- c("skukuza_se", "plantation_a", "thornybush", "bushbuckridge_a")
pattern <- paste0("^", subaois, collapse = "|")
pop <- pop %>%
  filter(str_starts(domain, pattern))
samp <- samp[samp$aoi %in% subaois,]

# Get summary of sample size across domains
print(summary(data.frame(table(samp$domain))))

# Run unit-level model for all metrics
metrics <- c("cover", "rh98", "fhd")
Xset = "lt.p.s.t"
models <- list()
est_df <- data.frame()

# # Test with small population
# testpop <- subpop %>%
#   group_by(domain) %>%
#   sample_n(100)

for (metric in metrics){
  print(metric)
  pred_col <- paste0("pred_lt.p.s.t_",metric)
  pop_metric <- rename(pop %>% dplyr::select(all_of(c(metric, "domain"))), !!pred_col:=metric)
  formula <- as.formula(paste(metric, "~", pred_col))
  
  if (metric=="rh98"){
    transformation="no"
  } else{
    transformation="box.cox"
  }
  
  model <- ebp(fixed = formula,
               pop_data = pop_metric, pop_domains = "domain",
               smp_data = samp, smp_domains = "domain",
               MSE = TRUE, interval="default",
               transformation=transformation, boot_type = "wild",
               L=200, B=200, seed=123, cpus=20
  )
  
  models[[metric]] <- model
  est <- estimators(model, indicator=c("Mean", "Median"), MSE=TRUE)$ind
  est$metric <- metric
  est_df <- rbind(est_df, est)
}

# Save results
today <- format(Sys.Date(), "%Y%m%d")
model_path <- paste0("J:/projects/ECOFOR/gedi/sae/sae_gedi_models_", today, ".rds")
estimates_path <- paste0("J:/projects/ECOFOR/gedi/sae/sae_gedi_estimates_", today, ".csv")
saveRDS(models, file=model_path)
write.csv(est_df, file=estimates_path)


#-------Make Diagnostic plots--------------------------------------------------

# Load saved models
models <- readRDS(r"(J:\projects\ECOFOR\gedi\sae\sae_gedi_models_20251103.rds)")

theme_8pt <- theme(
  text = element_text(size = 8),
  axis.text = element_text(size = 7),
  plot.title = element_text(size = 9, face = "bold"),
  strip.text = element_text(size = 8)
)

metrics <- c("cover", "rh98", "fhd")

for (metric in metrics){
  qqplots <- qqnorm(models[[metric]], gg_theme=theme_8pt)

  rdf <- data.frame(resids=models[[metric]]$model$residuals[,1],
                    fitted=models[[metric]]$model$fitted[,1])
  resid_fitted <- ggplot(rdf, aes(x = fitted, y = resids)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(se = FALSE, color = "blue") + # Optional: Add a smoothed line
    labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs. Fitted Values Plot") +
    theme_8pt
  
  all_plots <- grid.arrange(resid_fitted, qqplots, ncol=2)
  
  plot_path <- paste0(r"(F:\My Drive\Work\ecofor\manuscript\figs\sae_diagnostics_)", metric, ".svg")
  ggsave(plot_path, all_plots, height=2.5, width=6.5, units="in", dpi=300)
}





# # ----Fill pop domains method ---------------------------------------------------------------
# # Domains in samp need to be in pop even if we don't want to estimate in those domains.
# # copying rows from samp for the missing domains doesn't affect the model fit
# # But it does affect all the estimates (even for the domains of interest)
# # So this should not be done. The full population should be gathered in the years used in the sample.
# domains_notin_pop <- setdiff(unique(samp['domain']), unique(pop['domain']))$domain
# copy_rows <- samp %>%
#               filter(domain %in% domains_notin_pop) %>%
#               select("pred_cover", "pred_rh98", "pred_fhd_normal", "pred_pai", "rain_year", "aoi", "domain") %>%
#               rename(cover=pred_cover, rh98=pred_rh98, fhd=pred_fhd_normal, pai=pred_pai, year=rain_year) %>%
#               mutate(year = type.convert(year, as.is=TRUE))
# pop_filled <- bind_rows(pop, copy_rows)
#
# # Subsetting the copy rows doesn't make a difference in the model stats.
# # It does affect the estimates (even of the domains of interest with full populations)
# copy_rows1 <- copy_rows %>% group_by(domain) %>% slice(1)
# pop_filled1 <- bind_rows(pop, copy_rows1)


