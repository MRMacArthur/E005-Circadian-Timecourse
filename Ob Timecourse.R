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

obData_neg <- read.csv("Circadian Timecourse Ob Negative 20230213.csv")
obData_pos <- read.csv("Circadian Timecourse Ob Positive 20230213.csv")

timecourseMetaData <- read.csv("Circadian Timecourse Sample Order.csv")

obData_neg <- obData_neg %>%
  select("compound", starts_with("o"), starts_with("w"))

obData_pos <- obData_pos %>%
  select("compound", starts_with("o"), starts_with("w"))

obDataAll <- rbind(obData_neg, obData_pos)

obDataAllWide <- reshape2::melt(obDataAll)

colnames(obDataAllWide) <- c("compound", "MS.ID", "IC")

obDataAllWide <- merge(obDataAllWide,
                       timecourseMetaData,
                       by = "MS.ID",
                       all.x = T, all.y = F)

obDataAllNorm <- obDataAllWide %>%
  group_by(AnimalID) %>%
  mutate(med = median(IC)) %>%
  mutate(normIC = IC/med) %>%
  ungroup() %>%
  mutate(min.val = min(abs(normIC[normIC != 0]))/10) %>%
  group_by(AnimalID) %>%
  mutate(normLogIC = log10((normIC + sqrt(normIC^2 + min.val^2))/2))

obDataAllStats <- obDataAllNorm %>%
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

DT::datatable(obData_neg_stats)  

obData_neg_norm %>%
  filter(compound == "Guanosine") %>%
  ggplot(aes(x = as.numeric(Time), y = normLogIC, color = Group)) +
  geom_point() +
  stat_smooth(method='lm', formula = y~poly(x,4)) +
  scale_x_continuous(breaks = 0:6, labels = c("0600", "1000", "1400",
                                            "1800", "2200", "0200", "0600")) +
  labs(x = "Time", y = "Median normalized log10 intensity", title = "C15:2") +
  annotate("rect", xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  annotate("rect", xmin=3.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.2, fill="gray") +
  myTheme
  


