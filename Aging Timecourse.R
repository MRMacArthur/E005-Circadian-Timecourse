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

ageData_neg <- read.csv("Circadian Timecourse Aging Negative 20230213.csv")
ageData_pos <- read.csv("Circadian Timecourse Aging Positive 20230213.csv")

timecourseMetaData <- read.csv("Circadian Timecourse Sample Order.csv")

ageData_neg <- ageData_neg %>%
  select("compound", starts_with("a"), starts_with("y")) %>%
  select(-c("adductName", "a4_7"))

ageData_pos <- ageData_pos %>%
  select("compound", starts_with("a"), starts_with("y")) %>%
  select(-c("adductName", "a4_7"))

ageDataAll <- rbind(ageData_neg, ageData_pos)

ageDataAllWide <- reshape2::melt(ageDataAll)

colnames(ageDataAllWide) <- c("compound", "MS.ID", "IC")

ageDataAllWide <- merge(ageDataAllWide,
                       timecourseMetaData,
                       by = "MS.ID",
                       all.x = T, all.y = F)

ageDataAllNorm <- ageDataAllWide %>%
  group_by(AnimalID) %>%
  mutate(med = median(IC)) %>%
  mutate(normIC = IC/med) %>%
  ungroup() %>%
  mutate(min.val = min(abs(normIC[normIC != 0]))/10) %>%
  group_by(AnimalID) %>%
  mutate(normLogIC = log10((normIC + sqrt(normIC^2 + min.val^2))/2))

ageDataAllStats <- ageDataAllNorm %>%
  nest(data = -compound) %>%
  rowwise() %>%
  mutate(lme = list(lme(normLogIC ~ Time*Group, 
                        random = ~1|AnimalID,
                        data = data)),
         result = list(broom.mixed::tidy(lme))) %>%
  unnest(result) %>%
  select(-data, -lme) %>%
  filter(term != "(Intercept)",
         term != "sd_(Intercept)",
         term != "sd_Observation")

DT::datatable(ageDataAllStats)  

ageDataAllNorm %>%
  filter(compound == "Adenine") %>%
  ggplot(aes(x = as.numeric(Time), y = normLogIC, color = Group)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25,
                                             dodge.width = 0.5)) +
  stat_smooth(method='lm', formula = y~poly(x,4)) +
  scale_x_continuous(breaks = 0:6, labels = c("0600", "1000", "1400",
                                              "1800", "2200", "0200", "0600")) +
  labs(x = "Time", y = "Median normalized log10 intensity", title = "Glucose") +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=3.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  myTheme
