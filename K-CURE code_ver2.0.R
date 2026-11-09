
# ============================================================
# Bed-sharing study: consolidated analysis code
# 1. Sankey diagram
# 2. ANOVA + post-hoc
# 3. ACLM
# ============================================================

# ============================================================
# 0. Packages
# ============================================================
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggalluvial)

library(lavaan)
library(mitools)

library(car)
library(agricolae)
library(PMCMRplus)

# ============================================================
# 1. Figure 1: Sankey / alluvial diagram of bed-sharing
# cosleep = bed-sharing
# non-cosleep = non bed-sharing
# Wave2 = Time1
# Wave3 = Time2
# Wave4 = Time3
# Wave5 = Time4
# Wave6 = Time5
# Wave7 = Time6
# ============================================================
coslp <- read.csv("cosleep.csv", header = TRUE)

ggplot(coslp,
       aes(x = wave, stratum = coslp, alluvium = ScrNo, fill = as.factor(coslp), 
           label = ifelse(coslp == 0, "non-cosleep", "cosleep"))) +
  geom_stratum() +
  geom_flow() +
  scale_fill_manual(values = c("0" = "#87CEEB", "1" = "#FFC0CB"),
                    breaks = c(0, 1),
                    labels = c("non-cosleep", "cosleep")) +
  guides(fill = guide_legend(title = "Group")) +
  theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        panel.grid = element_blank(),  
        panel.background = element_blank())

wave_summary <- coslp %>%
  mutate(
    group_label = ifelse(coslp == 0, "non-cosleep", "cosleep")
  ) %>%
  count(wave, group_label, name = "n") %>%
  group_by(wave) %>%
  mutate(
    total = sum(n),
    percent = round(n / total * 100, 1)
  ) %>%
  ungroup()

wave_summary

# ============================================================
# 2. ANOVA
# Bed-sharing trajectory group coding
# Group 1 = Continuous bed-sharing
# Group 2 = Past bed-sharing
# Group 3 = Intermittent bed-sharing
# W2 = Time1; W7 = Time6
# ============================================================
data_anova <- read_excel("df_ANOVA.xlsx") %>%
  mutate(group = as.factor(group))

# Selected all variables beginning with W7
w7_vars <- names(data_anova)[grepl("^W7", names(data_anova))]


# ------------------------------------------------------------
# 2-1. group-difference check (W7)
# ------------------------------------------------------------

anova_results_w7 <- map_dfr(w7_vars, function(var) {
  
  summary_stats <- data_anova %>%
    group_by(group) %>%
    summarise(
      N = sum(!is.na(.data[[var]])),
      Mean = mean(.data[[var]], na.rm = TRUE),
      SD = sd(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(group, N, Mean, SD) %>%
    pivot_wider(
      id_cols = NULL,
      names_from = group,
      values_from = c(N, Mean, SD),
      names_glue = "Group{group}_{.value}",
      names_vary = "slowest"
    )
  
  formula <- as.formula(paste(var, "~ group"))
  aov_result <- aov(formula, data = data_anova)
  anova_summary <- summary(aov_result)[[1]]
  
  result <- tibble(
    Variable = var,
    F_value = anova_summary$`F value`[1],
    p_value = anova_summary$`Pr(>F)`[1]
  )
  
  bind_cols(result, summary_stats)
})

print(anova_results_w7)

write_xlsx(
  anova_results_w7,
  "anova_results_W7.xlsx"
)

# ------------------------------------------------------------
# 2-2. Post-hoc analysis
# ------------------------------------------------------------
posthoc_results <- list()

for (var in w7_vars) {

  formula <- as.formula(paste(var, "~ group"))

  # Omnibus ANOVA
  fit <- aov(formula, data = data_anova)
  omnibus_p <- summary(fit)[[1]]$`Pr(>F)`[1]

  # Skip post-hoc if omnibus ANOVA is not significant
  if (is.na(omnibus_p) || omnibus_p >= .05) {
    next
  }

  # Levene test
  lev <- car::leveneTest(formula, data = data_anova)
  levene_p <- lev$`Pr(>F)`[1]

  if (!is.na(levene_p) && levene_p >= .05) {

    # Equal variances: Duncan
    duncan <- agricolae::duncan.test(
      fit,
      "group",
      alpha = 0.05,
      group = TRUE,
      console = FALSE
    )

    posthoc_results[[var]] <- list(
      variable = var,
      omnibus_p = omnibus_p,
      levene_p = levene_p,
      method = "Duncan",
      result = duncan
    )

    cat("\n==================================================\n")
    cat("Variable:", var, "\n")
    cat("ANOVA p:", omnibus_p, "\n")
    cat("Levene p:", levene_p, "\n")
    cat("Post-hoc: Duncan\n")
    print(duncan$groups)

  } else {

    # Unequal variances: Dunnett T3
    dunnett <- PMCMRplus::dunnettT3Test(
      formula,
      data = data_anova
    )

    posthoc_results[[var]] <- list(
      variable = var,
      omnibus_p = omnibus_p,
      levene_p = levene_p,
      method = "Dunnett T3",
      result = dunnett
    )

    cat("\n==================================================\n")
    cat("Variable:", var, "\n")
    cat("ANOVA p:", omnibus_p, "\n")
    cat("Levene p:", levene_p, "\n")
    cat("Post-hoc: Dunnett T3\n")
    print(dunnett)
  }
}

# ============================================================
# 3. Baseline group-difference check (W2)
# ============================================================
# Selected all variables beginning with W2
w2_vars <- names(data_anova)[grepl("^W2_", names(data_anova))]

data_anova <- data_anova %>%
  mutate(
    across(starts_with("W2_CSHQ"), as.numeric)
  )

anova_results_w2 <- map_dfr(w2_vars, function(var) {

  summary_stats <- data_anova %>%
    group_by(group) %>%
    summarise(
      Mean = mean(.data[[var]], na.rm = TRUE),
      SD = sd(.data[[var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = group,
      values_from = c(Mean, SD),
      names_glue = "Group{group}_{.value}"
    )

  formula <- as.formula(paste(var, "~ group"))
  aov_result <- aov(formula, data = data_anova)
  anova_summary <- summary(aov_result)[[1]]

  tibble(
    Variable = var,
    F_value = anova_summary$`F value`[1],
    p_value = anova_summary$`Pr(>F)`[1]
  ) %>%
    bind_cols(summary_stats)
})

print(anova_results_w2)

# ============================================================
# 4. ACLM
# Analyses were performed using lavaan version 0.6-18.
# ============================================================

data <- read_excel("imp_tdf_rf.xlsx")
imp_list <- split(data, data$.imp)

ordered_vars <- c("coslp3", "coslp4", "coslp5", "coslp6", "coslp7")

run_model <- function(model) {
  
  fits <- lapply(imp_list, function(dat) {
    tryCatch(
      sem(model,
          data = dat,
          ordered = ordered_vars,
          estimator = "WLSMV",
          control = list(iter.max = 10000, rel.tol = 1e-6)),
      error = function(e) NULL
    )
  })
  
  fits <- fits[sapply(fits, function(x) !is.null(x) && lavInspect(x, "converged"))]
  
  fit_avg <- colMeans(
    do.call(rbind, lapply(fits, function(x) {
      fitMeasures(
        x,
        c("cfi", "tli", "rmsea",
          "rmsea.ci.lower", "rmsea.ci.upper", "srmr"),
        fm.args = list(rmsea.ci.level = 0.90)
      )
    })),
    na.rm = TRUE
  )
  
  list(fits = fits, fit = fit_avg)
}


# Table 3. Model fit ------------------------------------------------------------

# Model 1
model1 <- '
  coslp3 ~ coslp2 + W2CSHQ_an
  coslp4 ~ coslp3 + W3CSHQ_an
  coslp5 ~ coslp4 + W4CSHQ_an
  coslp6 ~ coslp5 + W5CSHQ_an
  coslp7 ~ coslp6 + W6CSHQ_an

  W3CSHQ_an ~ W2CSHQ_an + coslp2
  W4CSHQ_an ~ W3CSHQ_an + coslp3
  W5CSHQ_an ~ W4CSHQ_an + coslp4
  W6CSHQ_an ~ W5CSHQ_an + coslp5
  W7CSHQ_an ~ W6CSHQ_an + coslp6
'

# Model 2
model2 <- '
  coslp3 ~ coslp2 + W2CSHQ_an
  coslp4 ~ coslp3 + W3CSHQ_an
  coslp5 ~ coslp4 + W4CSHQ_an
  coslp6 ~ coslp5 + W5CSHQ_an
  coslp7 ~ coslp6 + W6CSHQ_an

  W3CSHQ_an ~ W2CSHQ_an + coslp2
  W4CSHQ_an ~ W3CSHQ_an + coslp3
  W5CSHQ_an ~ W4CSHQ_an + coslp4
  W6CSHQ_an ~ W5CSHQ_an + coslp5
  W7CSHQ_an ~ W6CSHQ_an + coslp6

  coslp3 ~~ coslp5
  coslp3 ~~ coslp6
  coslp3 ~~ coslp7
'

# Model 3
model3 <- '
  coslp3 ~ coslp2 + W2CSHQ_an
  coslp4 ~ coslp3 + W3CSHQ_an
  coslp5 ~ coslp4 + W4CSHQ_an
  coslp6 ~ coslp5 + W5CSHQ_an
  coslp7 ~ coslp6 + W6CSHQ_an

  W3CSHQ_an ~ b * W2CSHQ_an + coslp2
  W4CSHQ_an ~ b * W3CSHQ_an + coslp3
  W5CSHQ_an ~ b * W4CSHQ_an + coslp4
  W6CSHQ_an ~ b * W5CSHQ_an + coslp5
  W7CSHQ_an ~ b * W6CSHQ_an + coslp6

  coslp3 ~~ coslp5
  coslp3 ~~ coslp6
  coslp3 ~~ coslp7
'

# Model 4
model4  <- '
  coslp3 ~ a * coslp2 + W2CSHQ_an
  coslp4 ~ a * coslp3 + W3CSHQ_an
  coslp5 ~ a * coslp4 + W4CSHQ_an
  coslp6 ~ a * coslp5 + W5CSHQ_an
  coslp7 ~ a * coslp6 + W6CSHQ_an
  
  W3CSHQ_an ~ b * W2CSHQ_an + coslp2
  W4CSHQ_an ~ b * W3CSHQ_an + coslp3
  W5CSHQ_an ~ b * W4CSHQ_an + coslp4
  W6CSHQ_an ~ b * W5CSHQ_an + coslp5
  W7CSHQ_an ~ b * W6CSHQ_an + coslp6
  
  coslp3 ~~ coslp5
  coslp3 ~~ coslp6
  coslp3 ~~ coslp7
'

# Model 5
model5 <- '
  coslp3 ~ a * coslp2 + c * W2CSHQ_an
  coslp4 ~ a * coslp3 + c * W3CSHQ_an
  coslp5 ~ a * coslp4 + c * W4CSHQ_an
  coslp6 ~ a * coslp5 + c * W5CSHQ_an
  coslp7 ~ a * coslp6 + c * W6CSHQ_an
  
  W3CSHQ_an ~ b * W2CSHQ_an + coslp2
  W4CSHQ_an ~ b * W3CSHQ_an + coslp3
  W5CSHQ_an ~ b * W4CSHQ_an + coslp4
  W6CSHQ_an ~ b * W5CSHQ_an + coslp5
  W7CSHQ_an ~ b * W6CSHQ_an + coslp6
  
  coslp3 ~~ coslp5
  coslp3 ~~ coslp6
  coslp3 ~~ coslp7
'

# Model 6
model6 <- '
  coslp3 ~ a * coslp2 + c * W2CSHQ_an
  coslp4 ~ a * coslp3 + c * W3CSHQ_an
  coslp5 ~ a * coslp4 + c * W4CSHQ_an
  coslp6 ~ a * coslp5 + c * W5CSHQ_an
  coslp7 ~ a * coslp6 + c * W6CSHQ_an
  
  W3CSHQ_an ~ b * W2CSHQ_an + d * coslp2
  W4CSHQ_an ~ b * W3CSHQ_an + d * coslp3
  W5CSHQ_an ~ b * W4CSHQ_an + d * coslp4
  W6CSHQ_an ~ b * W5CSHQ_an + d * coslp5
  W7CSHQ_an ~ b * W6CSHQ_an + d * coslp6
  
  coslp3 ~~ coslp5
  coslp3 ~~ coslp6
  coslp3 ~~ coslp7
'

fit1 <- run_model(model1)
fit2 <- run_model(model2)
fit3 <- run_model(model3)
fit4 <- run_model(model4)
fit5 <- run_model(model5)
fit6 <- run_model(model6)

Table3 <- rbind(
  Model1 = fit1$fit,
  Model2 = fit2$fit,
  Model3 = fit3$fit,
  Model4 = fit4$fit,
  Model5 = fit5$fit,
  Model6 = fit6$fit
)

round(Table3, 3)


# Table 4 ------------------------------------------------

coef_list <- lapply(fit2$fits, coef)
vcov_list <- lapply(fit2$fits, vcov)

Table4_unstandardized <- summary(
  MIcombine(results = coef_list, variances = vcov_list)
)

Table4_unstandardized

std_list <- lapply(fit2$fits, function(x) {
  standardizedSolution(x) %>%
    select(lhs, op, rhs, est.std, se)
})

Table4_standardized <- bind_rows(std_list) %>%
  group_by(lhs, op, rhs) %>%
  summarize(
    beta = mean(est.std, na.rm = TRUE),
    SE = mean(se, na.rm = TRUE),
    lower_CI = beta - 1.96 * SE,
    upper_CI = beta + 1.96 * SE,
    .groups = "drop"
  )

Table4_standardized %>%
  filter(op == "~")
