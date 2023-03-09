library(ggplot2)
library(dplyr)
library(tidyr)
library(nlme)

myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 legend.key=element_blank())

LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}

urDataPos <- read.csv("Timecourse Urine/Timecourse Urine Positive 20230308.csv")
urDataNeg <- read.csv("Timecourse Urine/Timecourse Urine Negative 20230308.csv")

urMetadata <- read.csv("Timecourse Urine Metadata.csv")

urDataPos <- urDataPos %>%
  select("compound", starts_with(c("A", "W", "Y", "O"))) %>%
  select(-adductName)

urDataNeg <- urDataNeg %>%
  select("compound", starts_with(c("A", "W", "Y", "O"))) %>%
  select(-adductName)

urDataAll <- rbind(urDataPos, urDataNeg)

urDataAll$compound <- make.unique(urDataAll$compound, sep = "_")

urDataWide <- reshape2::melt(urDataAll)
colnames(urDataWide) <- c("compound", "MS.ID", "IC")

urDataWide <- merge(urDataWide,
                    urMetadata,
                    by = "MS.ID",
                    all.x = T, all.y = T)

urDataNorm <- urDataWide %>%
  group_by(MS.ID) %>%
  mutate(med = median(IC)) %>%
  mutate(normIC = IC/med) %>%
  ungroup() %>%
  mutate(min.val = min(abs(normIC[normIC != 0]))/10) %>%
  group_by(MS.ID) %>%
  mutate(normLogIC = log10((normIC + sqrt(normIC^2 + min.val^2))/2))

regression_fxn_age <- function(dat){lm(normLogIC ~ Age, data = dat) %>%
    broom::tidy()}

regression_fxn_ob <- function(dat){lm(normLogIC ~ Genotype, data = dat) %>%
    broom::tidy()}

urDataStats_age <- urDataNorm %>%
  filter(Study == 1) %>%
  nest(data = -compound) %>%
  mutate(df = purrr::map(data, regression_fxn_age)) %>%
  unnest_legacy(df) %>%
  filter(term != "(Intercept)")

urDataStats_ob <- urDataNorm %>%
  filter(Study == 2) %>%
  nest(data = -compound) %>%
  mutate(df = purrr::map(data, regression_fxn_ob)) %>%
  unnest_legacy(df) %>%
  filter(term != "(Intercept)")

DT::datatable(urDataStats_age)
DT::datatable(urDataStats_ob)

urDataNorm$Group <- factor(urDataNorm$Group,
                           levels = c("WT", "OB",
                                      "young", "aged"))

urDataNorm %>%
  filter(compound == "Glucosamine_2") %>%
  ggplot(aes(x = Group, y = normLogIC)) +
  geom_boxplot()

urDataStats_all <- merge(urDataStats_age,
                         urDataStats_ob,
                         by = "compound",
                         suffixes = c(".age", ".ob"))

urDataStats_all$sig <- "Not Sig"
urDataStats_all[urDataStats_all$p.value.age < 0.05,]$sig <- "Sig Age"
urDataStats_all[urDataStats_all$p.value.ob < 0.05,]$sig <- "Sig OB"
urDataStats_all[urDataStats_all$p.value.ob < 0.05 &
                  urDataStats_all$p.value.age < 0.05,]$sig <- "Sig Both"

ggplot(urDataStats_all, aes(x = estimate.age,
                            y = estimate.ob,
                            color = sig)) +
  geom_point() +
  labs(x = "Estimate Age/Young", y = "Estimate OB/WT",
       color = "Significant")

cor(urDataStats_all$estimate.age, urDataStats_all$estimate.ob,
    method = "spearman")
