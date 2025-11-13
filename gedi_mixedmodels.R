# Compare alternative random effect structures to the one used for SAE
library(lme4)
library(lmerTest)
library(arrow)
library(sfarrow)
library(sf)
library(dplyr)

path <- "J:/projects/ECOFOR/gedi/sae/sae_gedi_samp.gpkg"
df <- data.frame(st_read(path, int64_as_string = TRUE))

colnames(df) = gsub("_normal", "", colnames(df)) # fhd_normal -> fhd

# Select subset of AOIs used for paper
subaois <- c("skukuza_se", "plantation_a", "thornybush", "bushbuckridge_a")
df <- df[df$aoi %in% subaois,]

# Convert cover to percent
cover_cols <- grep("cover", names(df), value = TRUE)
df[,cover_cols] = df[,cover_cols] * 100

# Get models from emdi to get lambda's used
emdi_models <- readRDS(r"(J:\projects\ECOFOR\gedi\sae\sae_gedi_models_20251031.rds)")


###------------Iterate over metrics---------------------------------
metrics = c("cover", "rh98", "fhd")

reffects <- c("(1|aoi)", "(1|aoi) + (1|rain_year)", "(1|aoi:rain_year)")
for (metric in metrics){
  print(metric)
  models <- list()
  for (reffect in reffects){
    # apply same boxcox transformation used in emdi
    if (metric=="rh98"){
      lambda = 1
    } else {
      lambda <- emdi_models[[metric]]$transform_param$optimal_lambda
    }
    metric_trans <- paste0(metric,"_trans")
    df[metric_trans] <- if (lambda == 0) {
      log(y)
    } else {
      (df[metric]^lambda - 1) / lambda
    }
    
    # Get model fit with random effects
    formula <- as.formula(paste0(metric_trans, " ~ pred_lt.p.s.t_", metric, " + ", reffect))
    mod <- lmer(formula, data=df, REML=F)
    models[[reffect]] <- mod
    print(formula)
    print(VarCorr(mod))
    cat("\n\n")
    
  }
  anova_result <- eval(
    do.call(
      call, 
      c(list("anova"), 
        lapply(names(models), as.symbol),
        test="LRT",
        refit=FALSE), 
      quote = TRUE), 
    models)
  
  # Fix Anova table to get correct DF, LRT(chisq) and p-value for REML likelihood ratio test
  areml <- anova_result
  areml$Df <- max(areml$npar) - areml$npar
  areml$Chisq <- 2*(areml['(1|aoi) + (1|rain_year) + (1|aoi:rain_year)', 'logLik'] - areml$logLik)
  areml$`Pr(>Chisq)` <- pchisq(areml$Chisq, areml$Df, lower.tail = F)
  
  print(areml)
  cat("\n\n\n\n")
}