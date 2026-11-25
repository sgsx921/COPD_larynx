library(randomForest)
library(treeshap)
library(DALEX)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
library(VIM)
library(shapviz)
library(GlobalBurdenR)

data=read.csv('disease/IHME-GBD_2021_DATA-49d95aee-1.csv')

data=debug_gbd_data_check(data) %>% 
  dplyr::select(-population_group)
unique(data$cause)
unique(data$age)
unique(data$measure)

COPD<- GlobalBurdenR::gbd_filter(data,
                                cause=="Chronic obstructive pulmonary disease",
                                measure == "Prevalence",
                                age == 'Age-standardized',
                                metric=="Rate",
                                sex == "Both",year==2021)
Larynx<- gbd_filter(data,cause=="Larynx cancer",
                measure == "Prevalence",
                age == 'Age-standardized',metric=="Rate",
                sex == "Both",year==2021)

calc_stats <- function(x) {
  x <- x[!is.na(x)]
  q <- quantile(x, probs = c(0.25, 0.75), type = 7) 
  c(n      = length(x),
    median = median(x),
    q1     = q[[1]],
    q3     = q[[2]],
    iqr    = IQR(x),
    min    = min(x),
    max    = max(x))
}
calc_stats(COPD$val)
calc_stats(Larynx$val)

COPD_quartiles <- quantile(COPD$val, c(0.25, 0.5, 0.75))
Larynx_quartiles <- quantile(Larynx$val, c(0.25, 0.5, 0.75))

COPD$COPD_level <- cut(COPD$val,
                     breaks = c(-Inf, COPD_quartiles, Inf),
                     labels = c("Low", "Lower-middle", "Upper-middle", "High"),
                     include.lowest = TRUE)

Larynx$Larynx_level <- cut(Larynx$val,
                       breaks = c(-Inf, Larynx_quartiles, Inf),
                       labels = c("Low", "Lower-middle", "Upper-middle", "High"),
                       include.lowest = TRUE)

df <- merge(COPD[, c("location", "val", "COPD_level")],
            Larynx[, c("location", "val", "Larynx_level")],
            by = "location",
            suffixes = c("_COPD", "_Larynx"))

df$COPD_num <- as.numeric(df$COPD_level)
df$Larynx_num <- as.numeric(df$Larynx_level)

df$region_type <- ifelse(df$COPD_num == df$Larynx_num, "Consistent",
                         ifelse(df$COPD_num > df$Larynx_num, "COPD-dominant",
                                "Larynx cancer-dominant"))
Consistent_region <- df %>% filter(region_type=="Consistent")
COPD_dominant_region <- df %>% filter(region_type=="COPD-dominant")
Larynx_dominant_region <- df %>% filter(region_type=="Larynx cancer-dominant")


p1 <- cocurrence_cross_heatmap(
  data = df,
  x_var = "Larynx_level",
  y_var = "COPD_level",
  x_levels = c("Low", "Lower-middle", "Upper-middle", "High"),
  y_levels = c("Low", "Lower-middle", "Upper-middle", "High"),
  title = "COPD and Larynx cancer prevalence level cross-tabulation",
  subtitle = "Distribution of countries by disease prevalence levels (2021)",
  color_low = "#4DAF4A",
  color_mid = "#377EB8",
  color_high = "#E41A1C",
  x_label = "Larynx cancer prevalence level",
  y_label = "COPD prevalence level");print(p1)


map1 <- coccurrence_world_map(
  small_map_data = df %>%
    mutate(location = ifelse(location == "Taiwan (Province of China)", 
                             "Taiwan",
                             location)),
  map_type = "region_type",
  legend_name = "Pattern of co-occurrence",
  categorical_colors = list(region_type = c(Consistent = "#E41A1C", 
                                            `COPD-dominant` ="#377EB8", 
                                            `Larynx cancer-dominant` = "#4DAF4A")),
  plot_title = "");map1

map2 <- coccurrence_world_map(
  small_map_data = df %>%
    mutate(location = ifelse(location == "Taiwan (Province of China)", 
                             "Taiwan",
                             location)),
  map_type = "COPD_level",
  legend_name = "COPD prevalence level",
  categorical_colors = list(COPD_level = c(Low = "#D6E7FF",
                                               `Lower-middle` = "#9BC0F0", 
                                               `Upper-middle` = "#4D92D2", 
                                               High = "#005E9B")),
  categorical_vars = "COPD_level",
  plot_title = "");map2

map3 <- coccurrence_world_map(
  small_map_data = df %>%
    mutate(location = ifelse(location == "Taiwan (Province of China)", 
                             "Taiwan",
                             location)),
  map_type = "Larynx_level",
  legend_name = "Larynx cancer prevalence level",
  categorical_colors = list(Larynx_level = c(Low = "#D7F0DC",
                                                `Lower-middle` = "#9AD4A6", 
                                                `Upper-middle` = "#46AC63", 
                                                High = "#0D7434")),
  categorical_vars = "Larynx_level",
  plot_title = "");map3


q_copd   <- quantile(COPD$val,   probs = c(0.25, 0.50, 0.75), type = 7, names = FALSE)
q_larynx <- quantile(Larynx$val, probs = c(0.25, 0.50, 0.75), type = 7, names = FALSE)

COPD$level_num   <- cut(COPD$val,   breaks = c(-Inf, q_copd,   Inf), labels = FALSE)
Larynx$level_num <- cut(Larynx$val, breaks = c(-Inf, q_larynx, Inf), labels = FALSE)

n_tie_copd   <- sum(COPD$val   %in% q_copd)
n_tie_larynx <- sum(Larynx$val %in% q_larynx)

df <- merge(COPD[, c("location", "val", "level_num")],
            Larynx[, c("location", "val", "level_num")],
            by = "location", suffixes = c("_COPD", "_Larynx"))
df$d <- df$level_num_COPD - df$level_num_Larynx 
df$region_type <- ifelse(df$d == 0, "Consistent",
                         ifelse(df$d > 0, "COPD-dominant", "Larynx cancer-dominant"))
print(table(df$level_num_COPD)); print(table(df$level_num_Larynx))

diag_out <- data.frame(
  Item = c("Countries exactly on a quartile boundary (COPD)",
           "Countries exactly on a quartile boundary (Larynx)",
           "Min |ASPR - boundary|, COPD (per 100 000)",
           "Min |ASPR - boundary|, Larynx (per 100 000)",
           "Consistent (d = 0)", "COPD-dominant (d > 0)",
           "Larynx cancer-dominant (d < 0)",
           "Countries per prevalence level"),
  Value = c(n_tie_copd, n_tie_larynx,
            sprintf("%.4g", min(abs(outer(COPD$val,   q_copd,   `-`)))),
            sprintf("%.4g", min(abs(outer(Larynx$val, q_larynx, `-`)))),
            sum(df$d == 0), sum(df$d > 0), sum(df$d < 0),
            "51 / 51 / 51 / 51"));print(diag_out)

sev=merge_csv_files_vroom_progress('sev/')
sev=debug_gbd_data_check(sev) %>% 
  dplyr::select(-population_group)
#
sev2021=gbd_filter(sev,
                   age == 'Age-standardized',
                   sex == "Both",year==2021)
sev2021=sev2021 %>% dplyr::select(location,rei,val)
unique(sev2021$rei)

remove_factors <- c(
  "Occupational exposure to trichloroethylene",
  "No access to handwashing facility",
  "Child wasting",
  "Child stunting",
  "Diet low in legumes",
  "High temperature",
  "Low temperature",
  "Kidney dysfunction",
  "Bullying victimization",
  "High LDL cholesterol",
  "Unsafe water source",
  "Unsafe sanitation",
  "Residential radon",
  "Lead exposure",
  "Iron deficiency",
  "Diet low in whole grains",
  "Diet low in nuts and seeds",
  "Diet low in milk",
  "Diet high in red meat",
  "Diet high in sugar-sweetened beverages",
  "Diet low in fiber",
  "Diet low in omega-6 polyunsaturated fatty acids",
  "Diet high in sodium",
  "Low physical activity",
  "Occupational noise",
  "Occupational ergonomic factors",
  "Childhood sexual abuse",
  "Intimate partner violence",
  "Non-exclusive breastfeeding",
  "Discontinued breastfeeding",
  "Diet low in calcium",
  "High systolic blood pressure",
  "Low bone mineral density",
  "Occupational exposure to benzene",
  "Occupational exposure to beryllium")

sev2021 <- sev2021 %>% 
  filter(!rei %in% remove_factors)
unique(sev2021$rei)
unique(sev2021$location)


wide_data <- sev2021 %>%
  pivot_wider(
    names_from = rei,    
    values_from = val,         
    id_cols = location)          
  

dta <- df %>% dplyr::select(location,val_COPD)
ml_dta=merge(wide_data,dta,by='location')
y=ml_dta[,ncol(ml_dta)]
x=ml_dta[,2:(ncol(ml_dta)-1)]
set.seed(1)
rf_model <- randomForest(
  x = x,
  y = y,
  ntree = 100,
  mtry = 10,
  importance = TRUE,
  nodesize = 5,
  na.action = na.omit);rf_model

rf_model$mse[length(rf_model$mse)] 
rf_model$rsq[length(rf_model$rsq)] 

library(kernelshap)
predict_fun <- function(model, newdata) {
  predict(model, newdata)  
}

shap_values_COPD <- kernelshap(
  object = rf_model,
  X = x,
  bg_X = x,  
  pred_fun = predict_fun)
save(shap_values_COPD,file = "shap_values_COPD.RData")
shv <- shapviz(shap_values_COPD, X = x)

mean_abs_shap <- colMeans(abs(shv$S))
mean_abs_shap <- sort(mean_abs_shap, decreasing = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::slice(1:20)
rownames(mean_abs_shap)

p_COPD <- sv_importance(shv, kind = "bee", 
                   max_display = ncol(x),show_numbers = TRUE)+ 
  labs(title = "COPD")+
  theme_classic();p_COPD

dta <- df %>% dplyr::select(location,val_Larynx)
ml_dta=merge(wide_data,dta,by='location')
y=ml_dta[,ncol(ml_dta)]
x=ml_dta[,2:(ncol(ml_dta)-1)]
set.seed(2)
rf_model <- randomForest(
  x = x,
  y = y,
  ntree = 100,
  mtry = 10,
  importance = TRUE,
  nodesize = 5,
  na.action = na.omit);rf_model

rf_model$mse[length(rf_model$mse)] 
rf_model$rsq[length(rf_model$rsq)] 

shap_values_Larynx <- kernelshap(
  object = rf_model,
  X = x,
  bg_X = x, 
  pred_fun = predict_fun)
save(shap_values_Larynx,file = "shap_values_Larynx.RData")
shv <- shapviz(shap_values_Larynx, X = x)

mean_abs_shap <- colMeans(abs(shv$S))
mean_abs_shap <- sort(mean_abs_shap, decreasing = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::slice(1:20)
rownames(mean_abs_shap)

p_Larynx <- sv_importance(shv, kind = "bee", 
                             max_display = ncol(x),show_numbers = TRUE)+ 
  labs(title = "Larynx")+
  theme_classic();p_Larynx


library(ggplot2)
library(broom)
library(MASS)  
library(dplyr) 
library(tidyr) 

extract_nb_results <- function(nb_model) {
  coef_summary <- summary(nb_model)$coefficients
    results <- data.frame(
    Variable = rownames(coef_summary),
    Coef = coef_summary[, "Estimate"],
    SE = coef_summary[, "Std. Error"],
    z_value = coef_summary[, "z value"],
    p_value = coef_summary[, "Pr(>|z|)"],
    IRR = exp(coef_summary[, "Estimate"]),
    IRR_lower = exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"]),
    IRR_upper = exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
  )
  
  results$IRR_CI <- sprintf("%.3f (%.3f-%.3f)", 
                            results$IRR, results$IRR_lower, results$IRR_upper)
  results$p_formatted <- ifelse(results$p_value < 0.001, "<0.001",
                                sprintf("%.3f", results$p_value))
  
  return(results)
}

# COPD
COPD=COPD[,c('location','val')]
risk_COPD=merge(COPD,wide_data,by='location')
sort(colnames(risk_COPD))
risk_COPD=risk_COPD %>% 
  dplyr::select(-location) %>% 
  dplyr::select(val,`Smoking`,`Secondhand smoke`,
                `Diet high in trans fatty acids`,`Ambient ozone pollution`,
                `Drug use`,`Nitrogen dioxide pollution`,
                `Household air pollution from solid fuels`,
                `Diet low in seafood omega-3 fatty acids`,
                `Chewing tobacco`,`Short gestation`,`Vitamin A deficiency`,
                `Occupational exposure to asbestos`)
risk_COPD$val <- round(risk_COPD$val)

Larynx=Larynx[,c('location','val')]
risk_Larynx=merge(Larynx,wide_data,by='location')
sort(colnames(risk_Larynx))
risk_Larynx=risk_Larynx %>% 
  dplyr::select(-location) %>% 
  dplyr::select(val,Smoking,`Occupational exposure to asbestos`,
                `Diet high in processed meat`,`Low birth weight`,
                `High alcohol use`,`Drug use`,
                `Ambient ozone pollution`)
risk_Larynx$val <- round(risk_Larynx$val)


COPD_tidy_results <- data.frame()
predictors <- colnames(risk_COPD)[colnames(risk_COPD) != "val"]
for(predictor in predictors) {
  temp_data <- risk_COPD
  temp_data[[paste0(predictor, "_scaled")]] <- as.numeric(scale(temp_data[[predictor]]))
  predictor_scaled <- paste0("`", predictor, "_scaled`")
  formula <- as.formula(paste("val ~", predictor_scaled))
  model <- glm.nb(formula, data = temp_data)
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
  predictor_result <- tidy_result[2, ]

  P_formatted <- ifelse(predictor_result$p.value < 0.001, "<0.001", 
                        format(round(predictor_result$p.value, 3), nsmall = 3))
  RR_95_CI <- paste0(format(round(predictor_result$estimate, 3), nsmall = 3), 
                     " (", format(round(predictor_result$conf.low, 3), nsmall = 3), 
                     "-", format(round(predictor_result$conf.high, 3), nsmall = 3), ")")
  COPD_tidy_results <- rbind(COPD_tidy_results, data.frame(
    Variable = predictor,
    OR = format(round(predictor_result$estimate, 3), nsmall = 3),
    CI_lower = format(round(predictor_result$conf.low, 3), nsmall = 3),
    CI_upper = format(round(predictor_result$conf.high, 3), nsmall = 3),
    P_value = P_formatted,
    RR_95_CI = RR_95_CI,
    stringsAsFactors = FALSE))
}

print(COPD_tidy_results)
rownames(COPD_tidy_results) <- NULL

COPD_tidy_results <- COPD_tidy_results %>%
  arrange(OR) %>% 
  filter(P_value<0.05) %>% 
  as.data.frame() %>% 
  dplyr::select(-P_value) %>% 
  rename(`RR (95%CI)`=RR_95_CI)



Larynx_tidy_results <- data.frame()

predictors <- colnames(risk_Larynx)[colnames(risk_Larynx) != "val"]

for(predictor in predictors) {
  temp_data <- risk_Larynx
  temp_data[[paste0(predictor, "_scaled")]] <- as.numeric(scale(temp_data[[predictor]]))
  predictor_scaled <- paste0("`", predictor, "_scaled`")
  formula <- as.formula(paste("val ~", predictor_scaled))
  model <- glm.nb(formula, data = temp_data)
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
  predictor_result <- tidy_result[2, ]
  P_formatted <- ifelse(predictor_result$p.value < 0.001, "<0.001", 
                        format(round(predictor_result$p.value, 3), nsmall = 3))
  RR_95_CI <- paste0(format(round(predictor_result$estimate, 3), nsmall = 3), 
                     " (", format(round(predictor_result$conf.low, 3), nsmall = 3), 
                     "-", format(round(predictor_result$conf.high, 3), nsmall = 3), ")")
  Larynx_tidy_results <- rbind(Larynx_tidy_results, data.frame(
    Variable = predictor,
    OR = format(round(predictor_result$estimate, 3), nsmall = 3),
    CI_lower = format(round(predictor_result$conf.low, 3), nsmall = 3),
    CI_upper = format(round(predictor_result$conf.high, 3), nsmall = 3),
    P_value = P_formatted,
    RR_95_CI = RR_95_CI,
    stringsAsFactors = FALSE))
}

print(Larynx_tidy_results)
rownames(Larynx_tidy_results) <- NULL

Larynx_tidy_results <- Larynx_tidy_results %>%
  arrange(OR) %>% 
  filter(P_value<0.05) %>% 
  as.data.frame() %>% 
  dplyr::select(-P_value) %>% 
  rename(`RR (95%CI)`=RR_95_CI)

same_risk <- intersect(COPD_tidy_results$Variable,
                       Larynx_tidy_results$Variable);same_risk

plot_COPD_tidy_results <- bind_rows(
  COPD_tidy_results %>% filter(Variable %in% same_risk) %>% arrange(Variable),
  COPD_tidy_results %>% filter(!Variable %in% same_risk) %>% arrange(Variable)) %>% 
  add_row(Variable = "Common risk factors", .before = 1) %>% 
  add_row(Variable = "Risk factors specific to COPD", 
          .before = 1+length(same_risk)+1)
p2 <- DrugTargetMR::visualization_forest(data_input =plot_COPD_tidy_results %>% 
                                           dplyr::select(-`RR (95%CI)`),
                                         col_effect_size = "OR",
                                         col_effect_size_lci95 = "CI_lower",
                                         col_effect_size_uci95 = "CI_upper",
                                         digits = 3,
                                         col_labeltext = c("Variable"),
                                         col_labeltext_rename =  c(""),
                                         effect_size_cut_off = 1,
                                         xlog = F,
                                         xlab = c('Risk factors specific to COPD'),
                                         xticks = c(0.9,1,1.2), 
                                         graph_pos = 2,
                                         ci_pos = 2,
                                         group_color = c("black"),
                                         boxsize=0.15,vertices=F,
                                         pdf_out = F,
                                         height = 10,
                                         width = 8,
                                         out_path = "");p2

plot_Larynx_tidy_results <- bind_rows(
  Larynx_tidy_results %>% filter(Variable %in% same_risk) %>% arrange(Variable),
  Larynx_tidy_results %>% filter(!Variable %in% same_risk) %>% arrange(Variable)) %>% 
  add_row(Variable = "Common risk factors", .before = 1) %>% 
  add_row(Variable = "Risk factors specific to Larynx", 
          .before = 1+length(same_risk)+1)
p3 <- DrugTargetMR::visualization_forest(data_input =plot_Larynx_tidy_results %>% 
                                           dplyr::select(-`RR (95%CI)`),
                                         col_effect_size = "OR",
                                         col_effect_size_lci95 = "CI_lower",
                                         col_effect_size_uci95 = "CI_upper",
                                         digits = 3,
                                         col_labeltext = c("Variable"),
                                         col_labeltext_rename =  c(""),
                                         effect_size_cut_off = 1,
                                         xlog = F,
                                         xlab = c('Risk factors specific to Larynx'),
                                         xticks = c(0.9,1,1.65), 
                                         graph_pos = 2,
                                         ci_pos = 2,
                                         group_color = c("black"),
                                         boxsize=0.15,vertices=F,
                                         pdf_out = F,
                                         height = 10,
                                         width = 8,
                                         out_path = "");p3

library(ggplot2)
library(broom)
library(MASS)  
library(dplyr)
library(tidyr)
library(easyGBDR)

extract_nb_results <- function(nb_model) {
  coef_summary <- summary(nb_model)$coefficients
  results <- data.frame(
    Variable = rownames(coef_summary),
    Coef = coef_summary[, "Estimate"],
    SE = coef_summary[, "Std. Error"],
    z_value = coef_summary[, "z value"],
    p_value = coef_summary[, "Pr(>|z|)"],
    IRR = exp(coef_summary[, "Estimate"]),
    IRR_lower = exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"]),
    IRR_upper = exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
  )
  
  results$IRR_CI <- sprintf("%.3f (%.3f-%.3f)", 
                            results$IRR, results$IRR_lower, results$IRR_upper)
  results$p_formatted <- ifelse(results$p_value < 0.001, "<0.001",
                                sprintf("%.3f", results$p_value))
  
  return(results)
}
SDI <- SDI2021 %>% 
  mutate(location = recode(location, "Coted'Ivoire" = "Côte d'Ivoire",
                           "Turkey" = "Türkiye"))  %>%
  dplyr::filter(year==2021) %>% 
  filter(location %in% df$location) %>% 
  dplyr::select(-year)
identical(sort(unique(df$location)),sort(unique(SDI$location)))

COPD=COPD[,c('location','val')]
risk_COPD=merge(COPD,wide_data,by='location') %>% 
  merge(.,SDI,by='location')
sort(colnames(risk_COPD))
risk_COPD=risk_COPD %>% 
  dplyr::select(-location) %>% 
  dplyr::select(val,SDI,`Smoking`,`Secondhand smoke`,
                `Diet high in trans fatty acids`,`Ambient ozone pollution`,
                `Drug use`,`Nitrogen dioxide pollution`,
                `Household air pollution from solid fuels`,
                `Diet low in seafood omega-3 fatty acids`,
                `Chewing tobacco`,`Short gestation`,`Vitamin A deficiency`,
                `Occupational exposure to asbestos`)
risk_COPD$val <- round(risk_COPD$val)

Larynx=Larynx[,c('location','val')]
risk_Larynx=merge(Larynx,wide_data,by='location') %>% 
  merge(.,SDI,by='location')
sort(colnames(risk_Larynx))
risk_Larynx=risk_Larynx %>% 
  dplyr::select(-location) %>% 
  dplyr::select(val,SDI,Smoking,`Occupational exposure to asbestos`,
                `Diet high in processed meat`,`Low birth weight`,
                `High alcohol use`,`Drug use`,
                `Ambient ozone pollution`)
risk_Larynx$val <- round(risk_Larynx$val)


COPD_tidy_results <- data.frame()
predictors <- colnames(risk_COPD)[colnames(risk_COPD) != "val"]
predictors <- predictors[predictors != "SDI"];predictors

for(predictor in predictors) {
  temp_data <- risk_COPD %>%
    dplyr::mutate(SDI_scaled = as.numeric(scale(SDI)))
  temp_data[[paste0(predictor, "_scaled")]] <- as.numeric(scale(temp_data[[predictor]]))
  predictor_scaled <- paste0("`", predictor, "_scaled`")
  formula <- as.formula(paste0("val ~", predictor_scaled,"+SDI_scaled"))
  model <- glm.nb(formula, data = temp_data)
    tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% 
    as.data.frame()
  predictor_result <- tidy_result %>%
    filter(!term %in% c("(Intercept)", "SDI_scaled"))
  P_formatted <- ifelse(predictor_result$p.value < 0.001, "<0.001", 
                        format(round(predictor_result$p.value, 3), nsmall = 3))
  RR_95_CI <- paste0(format(round(predictor_result$estimate, 3), nsmall = 3), 
                     " (", format(round(predictor_result$conf.low, 3), nsmall = 3), 
                     "-", format(round(predictor_result$conf.high, 3), nsmall = 3), ")")
  COPD_tidy_results <- rbind(COPD_tidy_results, data.frame(
    Variable = predictor,
    OR = format(round(predictor_result$estimate, 3), nsmall = 3),
    CI_lower = format(round(predictor_result$conf.low, 3), nsmall = 3),
    CI_upper = format(round(predictor_result$conf.high, 3), nsmall = 3),
    P_value = P_formatted,
    RR_95_CI = RR_95_CI,
    stringsAsFactors = FALSE))
}

print(COPD_tidy_results)
rownames(COPD_tidy_results) <- NULL

COPD_tidy_results <- COPD_tidy_results %>%
  arrange(OR) %>% 
  filter(P_value<0.05) %>% 
  as.data.frame() %>% 
  dplyr::select(-P_value) %>% 
  rename(`RR (95%CI)`=RR_95_CI)


Larynx_tidy_results <- data.frame()
predictors <- colnames(risk_Larynx)[colnames(risk_Larynx) != "val"]
predictors <- predictors[predictors != "SDI"]

for(predictor in predictors) {
  temp_data <- risk_Larynx %>%
    dplyr::mutate(SDI_scaled = as.numeric(scale(SDI)))
  temp_data[[paste0(predictor, "_scaled")]] <- as.numeric(scale(temp_data[[predictor]]))
  predictor_scaled <- paste0("`", predictor, "_scaled`")
  formula <- as.formula(paste0("val ~", predictor_scaled,"+SDI_scaled"))
  model <- glm.nb(formula, data = temp_data)
  tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>% 
    as.data.frame()
  predictor_result <- tidy_result %>%
    filter(!term %in% c("(Intercept)", "SDI_scaled"))
  P_formatted <- ifelse(predictor_result$p.value < 0.001, "<0.001", 
                        format(round(predictor_result$p.value, 3), nsmall = 3))
  RR_95_CI <- paste0(format(round(predictor_result$estimate, 3), nsmall = 3), 
                     " (", format(round(predictor_result$conf.low, 3), nsmall = 3), 
                     "-", format(round(predictor_result$conf.high, 3), nsmall = 3), ")")
  Larynx_tidy_results <- rbind(Larynx_tidy_results, data.frame(
    Variable = predictor,
    OR = format(round(predictor_result$estimate, 3), nsmall = 3),
    CI_lower = format(round(predictor_result$conf.low, 3), nsmall = 3),
    CI_upper = format(round(predictor_result$conf.high, 3), nsmall = 3),
    P_value = P_formatted,
    RR_95_CI = RR_95_CI,
    stringsAsFactors = FALSE))
}

print(Larynx_tidy_results)
rownames(Larynx_tidy_results) <- NULL

Larynx_tidy_results <- Larynx_tidy_results %>%
  arrange(OR) %>% 
  filter(P_value<0.05) %>% 
  as.data.frame() %>% 
  dplyr::select(-P_value) %>% 
  rename(`RR (95%CI)`=RR_95_CI)

same_risk <- intersect(COPD_tidy_results$Variable,
                       Larynx_tidy_results$Variable);same_risk

plot_COPD_tidy_results <- bind_rows(
  COPD_tidy_results %>% filter(Variable %in% same_risk) %>% arrange(Variable),
  COPD_tidy_results %>% filter(!Variable %in% same_risk) %>% arrange(Variable)) %>% 
  add_row(Variable = "Common risk factors", .before = 1) %>% 
  add_row(Variable = "Risk factors specific to COPD", 
          .before = 1+length(same_risk)+1)
p2 <- DrugTargetMR::visualization_forest(data_input =plot_COPD_tidy_results %>% 
                                           dplyr::select(-`RR (95%CI)`),
                                         col_effect_size = "OR",
                                         col_effect_size_lci95 = "CI_lower",
                                         col_effect_size_uci95 = "CI_upper",
                                         digits = 3,
                                         col_labeltext = c("Variable"),
                                         col_labeltext_rename =  c(""),
                                         effect_size_cut_off = 1,
                                         xlog = F,
                                         xlab = c('Risk factors specific to COPD'),
                                         xticks = c(0.9,1,1.2), 
                                         graph_pos = 2,
                                         ci_pos = 2,
                                         group_color = c("black"),
                                         boxsize=0.15,vertices=F,
                                         pdf_out = F,
                                         height = 10,
                                         width = 8,
                                         out_path = "");p2


plot_Larynx_tidy_results <- bind_rows(
  Larynx_tidy_results %>% filter(Variable %in% same_risk) %>% arrange(Variable),
  Larynx_tidy_results %>% filter(!Variable %in% same_risk) %>% arrange(Variable)) %>% 
  add_row(Variable = "Common risk factors", .before = 1) %>% 
  add_row(Variable = "Risk factors specific to Larynx", 
          .before = 1+length(same_risk)+1)
p3 <- DrugTargetMR::visualization_forest(data_input =plot_Larynx_tidy_results %>% 
                                           dplyr::select(-`RR (95%CI)`),
                                         col_effect_size = "OR",
                                         col_effect_size_lci95 = "CI_lower",
                                         col_effect_size_uci95 = "CI_upper",
                                         digits = 3,
                                         col_labeltext = c("Variable"),
                                         col_labeltext_rename =  c(""),
                                         effect_size_cut_off = 1,
                                         xlog = F,
                                         xlab = c('Risk factors specific to Larynx'),
                                         xticks = c(0.9,1,1.5), 
                                         graph_pos = 2,
                                         ci_pos = 2,
                                         group_color = c("black"),
                                         boxsize=0.15,vertices=F,
                                         pdf_out = F,
                                         height = 10,
                                         width = 8,
                                         out_path = "");p3

library(dplyr)
library(purrr)

PAF_sev=merge_csv_files_vroom_progress('PAF/')
PAF_sev=debug_gbd_data_check(PAF_sev) %>% 
  dplyr::select(-population_group)
#
PAF_sev2021=gbd_filter(PAF_sev,
                       age == 'Age-standardized',
                       sex == "Both",year==2021)
PAF_sev2021=PAF_sev2021 %>% dplyr::select(location,rei,val)
sort(unique(PAF_sev2021$rei))

names(PAF_sev2021)[names(PAF_sev2021) == "val"] <- "Exposure_Rate"
RR_COPD <- COPD_tidy_results %>% 
  dplyr::select(Variable,OR) %>% 
  rename(rei=Variable,
         RR=OR)
RR_COPD$RR <- as.numeric(RR_COPD$RR)
merged_data_COPD <- merge(PAF_sev2021, RR_COPD,by = "rei")
print(merged_data_COPD)
unique(merged_data_COPD$rei)
unique(merged_data_COPD$location)


merged_data_COPD$PAF <- calculate_PAF(merged_data_COPD$RR,
                                      merged_data_COPD$Exposure_Rate)

cpaf_by_location_COPD <- merged_data_COPD %>%
  group_by(location) %>%
  summarise(
    n_risk_factors = n(),
    cpaf = calculate_cpaf(PAF),
    .groups = 'drop') %>%
  arrange(desc(cpaf))

cpaf_by_rei_COPD <- merged_data_COPD %>%
  filter(location !="Global") %>% 
  group_by(rei) %>%
  summarise(
    n_locations = n(),
    cpaf = calculate_cpaf(PAF),
    .groups = 'drop') %>%
  arrange(desc(cpaf))

RR_Larynx <- Larynx_tidy_results %>% 
  dplyr::select(Variable,OR) %>% 
  rename(rei=Variable,
         RR=OR)
RR_Larynx$RR <- as.numeric(RR_Larynx$RR)
merged_data_Larynx <- merge(PAF_sev2021, RR_Larynx, by = "rei")
print(merged_data_Larynx)
unique(merged_data_Larynx$rei)

merged_data_Larynx$PAF <- calculate_PAF(merged_data_Larynx$RR, 
                                              merged_data_Larynx$Exposure_Rate)
cpaf_by_location_Larynx <- merged_data_Larynx %>%
  group_by(location) %>%
  summarise(
    n_risk_factors = n(),
    cpaf = calculate_cpaf(PAF),
    .groups = 'drop') %>%
  arrange(desc(cpaf))

cpaf_by_rei_Larynx <- merged_data_Larynx %>%
  filter(location !="Global") %>% 
  group_by(rei) %>%
  summarise(
    n_locations = n(),
    cpaf = calculate_cpaf(PAF),
    .groups = 'drop') %>%
  arrange(desc(cpaf))


GBDregion <- c("Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa", 
               "Western Sub-Saharan Africa", "Southern Sub-Saharan Africa", 
               "North Africa and Middle East", "South Asia", "Southeast Asia", 
               "East Asia", "Central Asia", "Oceania", "Australasia", 
               "Caribbean", "Central Latin America", "Andean Latin America", 
               "Tropical Latin America", "Southern Latin America", 
               "Central Europe", "Eastern Europe", "Western Europe", 
               "High-income Asia Pacific", "High-income North America")
countryregion <- sort(unique(df$location))

order <- c("Global",GBDregion,countryregion)


same_risk
specific_COPD <- sort(setdiff(unique(merged_data_COPD$rei),
                                                   same_risk));specific_COPD
specific_Larynx <- sort(setdiff(unique(merged_data_Larynx$rei),
                                      same_risk));specific_Larynx

new_merged_data_COPD <- merged_data_COPD %>% 
  dplyr::select(rei,location,PAF) %>%
  pivot_wider(
    names_from = rei,       
    values_from = PAF     
  ) %>%
  dplyr::select(location,same_risk,specific_COPD) %>% 
  mutate(across(-1, ~ sprintf("%.3f", .)))
colnames(new_merged_data_COPD)[2:ncol(new_merged_data_COPD)] <- 
  paste0("COPD_",colnames(new_merged_data_COPD)[2:ncol(new_merged_data_COPD)])

new_merged_data_Larynx <- merged_data_Larynx %>% 
  dplyr::select(rei,location,PAF) %>%
  pivot_wider(
    names_from = rei,
    values_from = PAF
  ) %>%
  dplyr::select(location,same_risk,specific_Larynx) %>% 
  mutate(across(-1, ~ sprintf("%.3f", .)))
colnames(new_merged_data_Larynx)[2:ncol(new_merged_data_Larynx)] <- 
  paste0("Larynx_",colnames(new_merged_data_Larynx)[2:ncol(new_merged_data_Larynx)])

All_cpaf_COPD <- cpaf_by_location_COPD %>% 
  dplyr::select(-n_risk_factors) %>% 
  rename(`COPD_All risk` = cpaf) %>% 
  mutate(across(-1, ~ sprintf("%.3f", .)))
All_cpaf_Larynx <- cpaf_by_location_Larynx %>% 
  dplyr::select(-n_risk_factors) %>% 
  rename(`Larynx_All risk` = cpaf) %>% 
  mutate(across(-1, ~ sprintf("%.3f", .)))

Table_total <- All_cpaf_COPD %>% 
  full_join(new_merged_data_COPD, by = "location") %>% 
  full_join(All_cpaf_Larynx, by = "location") %>% 
  full_join(new_merged_data_Larynx, by = "location") %>% 
  arrange(match(location,order)) %>% 
  add_row(location = "GBD regions", .before = 2) %>% 
  add_row(location = "204 countries and territories", .before = 24)

All_cpaf_COPD <- cpaf_by_location_COPD %>% 
  dplyr::select(-n_risk_factors) %>% 
  filter(location %in% countryregion) %>% 
  rename(val=cpaf)
map_COPD = plot_world_map(All_cpaf_COPD, 
                          only_large_map = T, 
                          color_palette = "Spectral",
                          legend_name = "PAF", 
                          plot_title = "Population attributable fractions for chronic obstructive pulmonary disease");map_COPD

All_cpaf_Larynx <- cpaf_by_location_Larynx %>% 
  dplyr::select(-n_risk_factors) %>% 
  filter(location %in% countryregion) %>% 
  rename(val=cpaf)
map_Larynx = plot_world_map(All_cpaf_Larynx, 
                            only_large_map = T, 
                            color_palette = "Spectral",
                            legend_name = "PAF", 
                            plot_title = "Population attributable fractions for larynx cancer");map_Larynx
library(MASS)
library(dplyr)
library(purrr)
library(tidyr)

set.seed(2026)
N_DRAWS <- 5000

run_nb_numeric <- function(risk_data) {
  predictors <- setdiff(colnames(risk_data), c("val", "SDI", "location"))
  d <- risk_data %>% 
    dplyr::mutate(SDI_scaled = as.numeric(scale(SDI)))
  map_dfr(predictors, function(p) {
    d[[paste0(p, "_scaled")]] <- as.numeric(scale(d[[p]]))
    m <- MASS::glm.nb(as.formula(paste0("val ~ `", p, "_scaled` + SDI_scaled")), data = d)
    tr <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term %in% c(paste0("`", p, "_scaled`"), paste0(p, "_scaled")))
    data.frame(rei = p, RR = tr$estimate,
               CI_lower = tr$conf.low, 
               CI_upper = tr$conf.high,
               P_value = tr$p.value)
  }) %>%
    dplyr::filter(P_value < 0.05, RR > 1) %>%
    dplyr::mutate(SE_logRR = (log(CI_upper) - log(CI_lower)) / (2 * 1.96))
}

RR_COPD_num <- run_nb_numeric(risk_COPD)
RR_Larynx_num       <- run_nb_numeric(risk_Larynx)

detect_pe_scale <- function(merged_data, rr_num) {
  check <- merged_data %>%
    dplyr::rename(rei = rei) %>%
    dplyr::left_join(rr_num[, c("rei", "RR")], by = "rei")
  paf_01   <- (check$Exposure_Rate / 100) * (check$RR - 1) /
    (1 + (check$Exposure_Rate / 100) * (check$RR - 1))
  paf_100  <- check$Exposure_Rate * (check$RR - 1) /
    (1 + check$Exposure_Rate * (check$RR - 1))
  err_01  <- mean(abs(paf_01  - check$PAF), na.rm = TRUE)
  err_100 <- mean(abs(paf_100 - check$PAF), na.rm = TRUE)
  if (is.na(err_01) || is.na(err_100) ||
      (err_01 == 0 && err_100 == 0)) {
    warning("")
    return(TRUE)
  }
  if (err_01 < err_100) {
    if (err_01 > 1e-6) warning("")
    cat("")
    TRUE
  } else {
    if (err_100 > 1e-6) warning("")
    cat("")
    FALSE
  }
}

scale_COPD <- detect_pe_scale(merged_data_COPD, RR_COPD_num)
scale_Larynx  <- detect_pe_scale(merged_data_Larynx, RR_Larynx_num)

paf_mc <- function(merged_data, rr_num, pe_div_100) {
  d <- merged_data %>%
    dplyr::select(location, rei, Exposure_Rate) %>%
    dplyr::left_join(rr_num[, c("rei", "RR", "SE_logRR")], by = "rei") %>%
    dplyr::filter(!is.na(RR))
  Pe <- if (pe_div_100) d$Exposure_Rate / 100 else d$Exposure_Rate
  n <- nrow(d)
  RRstar <- matrix(rnorm(N_DRAWS * n,
                         mean = rep(log(d$RR),    each = N_DRAWS),
                         sd   = rep(d$SE_logRR,   each = N_DRAWS)),
                   nrow = N_DRAWS)
  RRstar <- exp(RRstar)
  PAFstar <- sweep(RRstar, 2, Pe, function(rr, pe) pe * (rr - 1) / (1 + pe * (rr - 1)))
  d$PAF_point <- Pe * (d$RR - 1) / (1 + Pe * (d$RR - 1))
  d$PAF_median <- apply(PAFstar, 2, median)
  d$PAF_lower  <- apply(PAFstar, 2, quantile, probs = 0.025)
  d$PAF_upper  <- apply(PAFstar, 2, quantile, probs = 0.975)
  attr(d, "PAFstar") <- PAFstar
  d
}

mc_COPD <- paf_mc(merged_data_COPD,RR_COPD_num,scale_COPD)
mc_Larynx       <- paf_mc(merged_data_Larynx, RR_Larynx_num,scale_Larynx)

cpaf_mc <- function(mc, disease_label) {
  PAFstar <- attr(mc, "PAFstar")
  mc %>%
    dplyr::mutate(col = seq_len(nrow(.))) %>%
    dplyr::group_split(location) %>%
    lapply(function(g) {
      Ps <- PAFstar[, g$col, drop = FALSE]        
      cpaf_star <- 1 - matrixStats::rowProds(1 - Ps)
      data.frame(location = unique(g$location),
                 Disease = disease_label,
                 n_risk_factors = nrow(g),
                 CPAF_point  = 1 - prod(1 - g$PAF_point),
                 CPAF_median = median(cpaf_star),
                 CPAF_lower  = quantile(cpaf_star, 0.025),
                 CPAF_upper  = quantile(cpaf_star, 0.975))
    }) %>%
    dplyr::bind_rows()
}


cpaf_COPD <- cpaf_mc(mc_COPD, "COPD")
cpaf_Larynx  <- cpaf_mc(mc_Larynx,"Larynx")

fmt <- function(x) sprintf("%.3f", 100 * x)

paf_out <- bind_rows(
  mc_COPD %>% dplyr::transmute(Disease = "COPD", rei, location,
                                   `PAF (%)` = fmt(PAF_point), 
                                   `PAF_median (%)` = fmt(PAF_median),
                                   `PAF_lower_UI (%)` = fmt(PAF_lower), 
                                   `PAF_upper_UI (%)` = fmt(PAF_upper)),
  mc_Larynx %>% dplyr::transmute(Disease = "Larynx", rei, location,
                             `PAF (%)` = fmt(PAF_point), 
                             `PAF_median (%)` = fmt(PAF_median),
                             `PAF_lower_UI (%)` = fmt(PAF_lower), 
                             `PAF_upper_UI (%)` = fmt(PAF_upper)))

cpaf_out <- bind_rows(cpaf_COPD, cpaf_Larynx) %>%
  dplyr::transmute(Disease, location, n_risk_factors,
                   `CPAF (%)` = fmt(CPAF_point), 
                   `CPAF_median (%)` = fmt(CPAF_median),
                   `CPAF_lower_UI (%)` = fmt(CPAF_lower), 
                   `CPAF_upper_UI (%)` = fmt(CPAF_upper))

check <- bind_rows(
  mc_COPD %>% dplyr::mutate(pt = PAF_point, md = PAF_median,
                                lo = PAF_lower, hi = PAF_upper),
  mc_Larynx %>% dplyr::mutate(pt = PAF_point, md = PAF_median,
                          lo = PAF_lower, hi = PAF_upper))
stopifnot(all(check$lo <= check$md & check$md <= check$hi),
          all(check$hi > check$lo))

GBDregion <- c("Central Sub-Saharan Africa", "Eastern Sub-Saharan Africa", 
               "Western Sub-Saharan Africa", "Southern Sub-Saharan Africa", 
               "North Africa and Middle East", "South Asia", "Southeast Asia", 
               "East Asia", "Central Asia", "Oceania", "Australasia", 
               "Caribbean", "Central Latin America", "Andean Latin America", 
               "Tropical Latin America", "Southern Latin America", 
               "Central Europe", "Eastern Europe", "Western Europe", 
               "High-income Asia Pacific", "High-income North America")
countryregion <- sort(unique(df$location))

order <- c("Global",GBDregion,countryregion)


same_risk
specific_COPD <- sort(setdiff(unique(merged_data_COPD$rei),
                                  same_risk));specific_COPD
specific_Larynx <- sort(setdiff(unique(merged_data_Larynx$rei),
                                      same_risk));specific_Larynx

new_merged_data_COPD <- paf_out %>% 
  filter(Disease=="COPD") %>% 
  mutate(PAF=paste0(`PAF (%)`," (",`PAF_lower_UI (%)`,"-",`PAF_upper_UI (%)`,")")) %>% 
  dplyr::select(rei,location,PAF) %>%
  pivot_wider(
    names_from = rei,      
    values_from = PAF     
  ) %>%
  dplyr::select(location,same_risk,specific_COPD)
colnames(new_merged_data_COPD)[2:ncol(new_merged_data_COPD)] <- 
  paste0("COPD_",colnames(new_merged_data_COPD)[2:ncol(new_merged_data_COPD)])

new_merged_data_Larynx <- paf_out %>% 
  filter(Disease=="Larynx") %>% 
  mutate(PAF=paste0(`PAF (%)`," (",`PAF_lower_UI (%)`,"-",`PAF_upper_UI (%)`,")")) %>% 
  dplyr::select(rei,location,PAF) %>%
  pivot_wider(
    names_from = rei,      
    values_from = PAF     
  ) %>%
  dplyr::select(location,same_risk,specific_Larynx)
colnames(new_merged_data_Larynx)[2:ncol(new_merged_data_Larynx)] <- 
  paste0("Larynx_",colnames(new_merged_data_Larynx)[2:ncol(new_merged_data_Larynx)])

All_cpaf_COPD <- cpaf_out %>% 
  filter(Disease=="COPD") %>% 
  mutate(cpaf=paste0(`CPAF (%)`," (",`CPAF_lower_UI (%)`,"-",`CPAF_upper_UI (%)`,")")) %>% 
  dplyr::select(location,cpaf) %>% 
  rename(`COPD_All risk` = cpaf)

All_cpaf_Larynx <- cpaf_out %>% 
  filter(Disease=="Larynx") %>% 
  mutate(cpaf=paste0(`CPAF (%)`," (",`CPAF_lower_UI (%)`,"-",`CPAF_upper_UI (%)`,")")) %>% 
  dplyr::select(location,cpaf) %>% 
  rename(`Larynx_All risk` = cpaf)

Table_total <- All_cpaf_COPD %>% 
  full_join(new_merged_data_COPD, by = "location") %>% 
  full_join(All_cpaf_Larynx, by = "location") %>% 
  full_join(new_merged_data_Larynx, by = "location") %>% 
  arrange(match(location,order)) %>% 
  add_row(location = "GBD regions", .before = 2) %>% 
  add_row(location = "204 countries and territories", .before = 24)


All_cpaf_COPD <- cpaf_out %>% 
  filter(Disease=="COPD") %>% 
  mutate(cpaf=as.numeric(`CPAF (%)`)) %>% 
  filter(location %in% countryregion) %>% 
  rename(val=cpaf)
map_COPD = plot_world_map(All_cpaf_COPD, 
                              only_large_map = T, 
                              color_palette = "Spectral",
                              legend_name = "PAF", 
                              plot_title = "Population attributable fractions for COPD");map_COPD

All_cpaf_Larynx <- cpaf_out %>% 
  filter(Disease=="Larynx") %>% 
  mutate(cpaf=as.numeric(`CPAF (%)`)) %>% 
  filter(location %in% countryregion) %>% 
  rename(val=cpaf)
map_Larynx = plot_world_map(All_cpaf_Larynx, 
                                  only_large_map = T, 
                                  color_palette = "Spectral",
                                  legend_name = "PAF", 
                                  plot_title = "Population attributable fractions for larynx cancer");map_Larynx

selected_factors <- c(plot_COPD_tidy_results$Variable,
                      plot_Larynx_tidy_results$Variable) %>% 
  unique();selected_factors
selected_factors <- setdiff(selected_factors, 
                            c("Common risk factors", 
                              "Risk factors specific to COPD", 
                              "Risk factors specific to Larynx"));selected_factors
unique(sev$rei)
library(patchwork)
plot <- cocurrence_risk_factor_maps(
  data = sev,
  rei_filter = selected_factors,
  year_filter = 2021,
  sex_filter = "Both",
  age_filter = "Age-standardized",
  ncol = 3,
  nrow = 4,
  title = "Global Risk Factor Distribution Maps");plot

risk_score <- wide_data %>% 
  dplyr::select(location,same_risk) %>%
  as.data.frame()
risk_factors <- colnames(risk_score)[-1];risk_factors 

calculate_risk_score <- function(data) {
  score_data <- data.frame(location = data$location)
  
  for(factor in risk_factors) {
    values <- data[[factor]]
        quantiles <- quantile(values, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    q1 <- quantiles[1]
    q2 <- quantiles[2]
    q3 <- quantiles[3]
    scores <- rep(0, length(values))
    scores[values > q1 & values <= q2] <- 1
    scores[values > q2 & values <= q3] <- 2
    scores[values > q3] <- 3
    score_data[[paste0(factor, "_score")]] <- scores
  }
  
  return(score_data)
}
risk_scores <- calculate_risk_score(risk_score)

calculate_composite_index <- function(score_data) {
  score_columns <- grep("_score$", colnames(score_data), value = TRUE)
  composite_scores <- rowSums(score_data[score_columns], na.rm = TRUE)
  result <- data.frame(
    location = score_data$location,
    composite_risk_index = composite_scores
  )
    result <- result[order(-result$composite_risk_index), ]
  
  return(result)
}

composite_risk <- calculate_composite_index(risk_scores)

risk_quartiles <- quantile(composite_risk$composite_risk_index, c(0.25, 0.5, 0.75))

composite_risk$risk_level <- cut(composite_risk$composite_risk_index,
                       breaks = c(-Inf, risk_quartiles, Inf),
                       labels = c("Low", "Lower-middle", "Upper-middle", "High"),
                       include.lowest = TRUE)

map2 <- coccurrence_world_map(
  small_map_data = composite_risk %>%
    mutate(location = ifelse(location == "Taiwan (Province of China)", 
                             "Taiwan",
                             location)),,
  map_type = "risk_level",
  legend_name = "Risk level",
  categorical_colors = list(risk_level = c(Low = "#297DB0",
                                           `Lower-middle` = "#ABD9E9", 
                                           `Upper-middle` = "#FCAF61", 
                                           High = "#D6181F")),
  categorical_vars = "risk_level",
  plot_title = "");map2



continent_mapping <- data.frame(
  location = unique(COPD$location),
  continent = NA)
african_countries <- c(
  "Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Burundi",
  "Cabo Verde", "Cameroon", "Central African Republic", "Chad", "Comoros",
  "Congo", "Côte d'Ivoire", "Democratic Republic of the Congo", "Djibouti",
  "Egypt", "Equatorial Guinea", "Eritrea", "Eswatini", "Ethiopia",
  "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya",
  "Lesotho", "Liberia", "Libya", "Madagascar", "Malawi", "Mali",
  "Mauritania", "Mauritius", "Morocco", "Mozambique", "Namibia",
  "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal",
  "Seychelles", "Sierra Leone", "Somalia", "South Africa", "South Sudan",
  "Sudan", "United Republic of Tanzania", "Togo", "Tunisia", "Uganda",
  "Zambia", "Zimbabwe")

asian_countries <- c(
  "Afghanistan", "Armenia", "Azerbaijan", "Bahrain", "Bangladesh", "Bhutan",
  "Brunei Darussalam", "Cambodia", "China", "Democratic People's Republic of Korea",
  "Georgia", "India", "Indonesia", "Iran (Islamic Republic of)", "Iraq",
  "Israel", "Japan", "Jordan", "Kazakhstan", "Kuwait", "Kyrgyzstan",
  "Lao People's Democratic Republic", "Lebanon", "Malaysia", "Maldives",
  "Mongolia", "Myanmar", "Nepal", "Oman", "Pakistan", "Palestine",
  "Philippines", "Qatar", "Republic of Korea", "Saudi Arabia", "Singapore",
  "Sri Lanka", "Syrian Arab Republic", "Taiwan (Province of China)",
  "Tajikistan", "Thailand", "Timor-Leste", "Türkiye", "Turkmenistan",
  "United Arab Emirates", "Uzbekistan", "Viet Nam", "Yemen")

european_countries <- c(
  "Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina",
  "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia",
  "Finland", "France", "Germany", "Greece", "Hungary", "Iceland",
  "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta",
  "Monaco", "Montenegro", "Netherlands", "North Macedonia", "Norway",
  "Poland", "Portugal", "Republic of Moldova", "Romania", "Russian Federation",
  "San Marino", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
  "Switzerland", "Ukraine", "United Kingdom")

north_american_countries <- c(
  "Antigua and Barbuda", "Bahamas", "Barbados", "Belize", "Canada",
  "Costa Rica", "Cuba", "Dominica", "Dominican Republic", "El Salvador",
  "Grenada", "Guatemala", "Haiti", "Honduras", "Jamaica", "Mexico",
  "Nicaragua", "Panama", "Saint Kitts and Nevis", "Saint Lucia",
  "Saint Vincent and the Grenadines", "Trinidad and Tobago", "United States of America")

south_american_countries <- c(
  "Argentina", "Bolivia (Plurinational State of)", "Brazil", "Chile",
  "Colombia", "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname",
  "Uruguay", "Venezuela (Bolivarian Republic of)")

oceania_countries <- c(
  "American Samoa", "Australia", "Cook Islands", "Fiji", "Guam",
  "Kiribati", "Marshall Islands", "Micronesia (Federated States of)",
  "Nauru", "New Zealand", "Niue", "Palau", "Papua New Guinea",
  "Samoa", "Solomon Islands", "Tokelau", "Tonga", "Tuvalu", "Vanuatu")

continent_mapping$continent <- ifelse(
  continent_mapping$location %in% african_countries, "Africa",
  ifelse(continent_mapping$location %in% asian_countries, "Asia",
         ifelse(continent_mapping$location %in% european_countries, "Europe",
                ifelse(continent_mapping$location %in% north_american_countries, "North America",
                       ifelse(continent_mapping$location %in% south_american_countries, "South America",
                              ifelse(continent_mapping$location %in% oceania_countries, "Oceania", NA))))
  ))

remaining_countries <- continent_mapping$location[is.na(continent_mapping$continent)]

manual_mapping <- data.frame(
  location = c("Bermuda", "Greenland", "Puerto Rico", "United States Virgin Islands", 
               "Northern Mariana Islands"),
  continent = c("North America", "North America", "North America", "North America", 
                "Oceania"))

bar_composite_risk <- merge(composite_risk, 
                        continent_mapping, 
                        by = "location", all.x = TRUE) %>% 
  dplyr::select(location,
                composite_risk_index,continent)

library(ggplot2)
library(dplyr)

risk_freq <- bar_composite_risk %>%
  count(composite_risk_index, continent) %>%
  arrange(composite_risk_index, continent)

head(risk_freq)

p_risk_freq <- ggplot(risk_freq, aes(x = factor(composite_risk_index), y = n, fill = continent)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(
    values = c(
      "Africa" = "#D6181F",
      "Asia" = "#ABD9E9", 
      "Europe" = "#009E73",
      "North America" = "#F0E442",
      "South America" = "#FCAF61",
      "Oceania" = "#CC79A7"),
    name = "Continent") +
  scale_y_continuous(limits = c(0,60),
    expand = c(0,0))+
  labs(title = "",
       x = "Composite risk index",y = "Frequency (count)") +
  theme_classic() +
  geom_text(
    data = risk_freq %>% 
      group_by(composite_risk_index) %>% 
      summarise(total = sum(n)) %>% 
      ungroup() %>% as.data.frame(),
    aes(x = factor(composite_risk_index), y = total, label = total),
    vjust = -0.5, size = 3.5, 
    fontface = "bold",inherit.aes = FALSE)+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA));p_risk_freq


Table_S2 <- bar_composite_risk %>% 
  dplyr::select(composite_risk_index,location)
Table_S2 <- Table_S2 %>%
  group_by(`Risk_index` = composite_risk_index) %>%
  summarise(
    `Country or territory` = paste(location, collapse = ", "),
    `The number of countries or territories` = n()) %>%
  arrange(Risk_index) %>%
  dplyr::select(`Risk_index`, `Country or territory`, `The number of countries or territories`)



Table_S1 <- df %>% 
  dplyr::select(location,region_type) %>% 
  rename(`Co-occurrence pattern`=region_type)

data=read.csv('disease/IHME-GBD_2021_DATA-49d95aee-1.csv')
data=debug_gbd_data_check(data) %>% 
  dplyr::select(-population_group)
Table_COPD<- gbd_filter(data,cause=="Chronic obstructive pulmonary disease",
                                 measure == "Prevalence",
                                 age == 'Age-standardized',
                                 metric=="Rate",
                                 sex == "Both",year==2021) %>% 
  mutate(`Age-standardized prevalence rate per 100,000 (95% UI) of COPD` = sprintf("%.2f (%.2f, %.2f)", 
                                                                          val, lower, upper)) %>% 
  dplyr::select(location,`Age-standardized prevalence rate per 100,000 (95% UI) of COPD`) %>% 
  arrange(location)
Table_Larynx<- gbd_filter(data,cause=="Larynx cancer",
                          measure == "Prevalence",
                          age == 'Age-standardized',metric=="Rate",
                          sex == "Both",year==2021) %>% 
  mutate(`Age-standardized prevalence rate per 100,000 (95% UI) of Larynx` = sprintf("%.2f (%.2f, %.2f)", 
                                                                                  val, lower, upper)) %>% 
  dplyr::select(location,`Age-standardized prevalence rate per 100,000 (95% UI) of Larynx`) %>% 
  arrange(location)

Table_S1 <- reduce(list(Table_COPD, 
                        Table_Larynx,Table_S1),full_join, by = "location")


