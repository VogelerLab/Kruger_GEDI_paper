# library(nlme)
library(lme4)
library(lmerTest)
library(arrow)
library(sfarrow)
library(sf)
library(dplyr)

path <- "J:/projects/ECOFOR/gedi/sae/sae_gedi_samp.gpkg"
df <- data.frame(st_read(path, int64_as_string = TRUE))

colnames(df) = gsub("_normal", "", colnames(df)) # fhd_normal -> fhd

# Select subset of AOIs presented in case studies
subaois <- c("skukuza_se", "plantation_a", "thornybush", "bushbuckridge_a")
df <- df[df$aoi %in% subaois,]

# Convert cover to percent
cover_cols <- grep("cover", names(df), value = TRUE)
df[,cover_cols] = df[,cover_cols] * 100

# Get models from emdi to get lambda's used
emdi_models <- readRDS(r"(J:\projects\ECOFOR\gedi\sae\sae_gedi_models_20251031.rds)")


###------------Iterate over metrics---------------------------------
metrics = c("cover", "rh98", "fhd")

reffects <- c("(1|aoi)", "(1|aoi) + (1|rain_year)", "(1|aoi:rain_year)")#, "(1|aoi) + (1|aoi:rain_year)", "(1|aoi) + (1|rain_year) + (1|aoi:rain_year)")
for (metric in metrics){
  print(metric)
  models <- list()
  for (reffect in reffects){
    #TODO: apply same boxcox transformation used in emdi and report on transformed scale
    #      or don't use transformation at all
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
    
    # # Run model comparison with REML which is better for random effects
    # # ranova only drops 1 RE at a time though
    # if (reffect == "(1|aoi) + (1|rain_year) + (1|aoi:rain_year)"){
    #   ranova_result <- ranova(mod)
    # }
    
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
  # This should match output from ranova
  areml <- anova_result
  areml$Df <- max(areml$npar) - areml$npar
  areml$Chisq <- 2*(areml['(1|aoi) + (1|rain_year) + (1|aoi:rain_year)', 'logLik'] - areml$logLik)
  areml$`Pr(>Chisq)` <- pchisq(areml$Chisq, areml$Df, lower.tail = F)
  
  print(areml)
  # print(ranova_result)
  cat("\n\n\n\n")
}

#-----------Test on Cover---------------------------------------------------------
# Model 1
# Not sure how to specify (1|area x time) in Paco's comment
# Using (1|domain) leads to "isSingular"
mod1 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|aoi) + (1|rain_year) + (1|aoi:rain_year), 
             data=df, REML = T) 


# Model 2
mod2 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|aoi) + (1|rain_year), 
             data=df, REML = T)


# Model 3
mod3 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|rain_year), 
             data=df, REML = T)

# Model 4
mod4 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|aoi), 
             data=df, REML = T)

# Model 5
mod5 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|aoi:rain_year), 
             data=df, REML = T)

# # Model 6: To check that aoi:rain_year and domain are the same (it is)
# mod6 <- lmer(cover ~ pred_lt.p.s.t_cover + (1|domain), 
#              data=df, REML = T)

anova(mod1, mod2, mod3, mod4, mod5)

ranova(mod2)
