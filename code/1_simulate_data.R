# simulate data, sodium and systolic blood pressure


# 1. setup ----------------------------------------------------------------

library(here)
library(tidyverse)


# 2. PARAMETERS: usual intake ---------------------------------------------------------

# pull CCHS data
url1 <- "https://health.canada.ca/apps/open-information/sodium-potassium/Sodium-usual-intake-Canada-excl-territories-2015-en.csv"
raw_intake <- read_csv(url1, skip=1) 

# isolate measurements of interest
intake <- raw_intake |> 
  filter(Country=="Canada") |>  # remove data beyond the estimates (includes footnotes etc)
  filter(`Age or Life Stage Group` == "19 years and over") |> # only interested in total adult population, not age subgroups
  transmute(
    Sex,
    Mean,
    sd = (p75-p25) / 1.349 # explain why this number using qnorm(0.75) - qnorm(0.25)
  )

# pull population mean (SD)
mean_sodium <- intake |> filter(Sex=="Both sexes") |> pull(Mean)
sd_sodium <- intake |> filter(Sex=="Both sexes") |> pull(sd)


# 3. PARAMETERS: blood pressure -------------------------------------------------------

# see here for average among Canadian adults: https://www150.statcan.gc.ca/n1/pub/82-625-x/2021001/article/00001-eng.htm
mean_sbp <- 114 # mmHg

# could not easily find population SD, using NHANES estimate (rough guide, does not impact bias in the simulation)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10424908/
sd_sbp <- 16 # mmHg

# 4. PARAMETERS: effect of sodium on systolic blood pressure ----------------------------------------------

# here we need to get data on how sodium impacts systolic blood pressure
# we will use a large Cochrane review, see He et al. BMJ. 2013;346:f1325 (https://www.bmj.com/content/346/bmj.f1325.short)

# abstract in He et. al: "A 100 mmol reduction in 24 hour urinary sodium (6 g/day salt) 
# was associated with a fall in systolic blood pressure of 5.8 mm Hg (2.5 to 9.2, P=0.001) 
# after adjustment for age, ethnic group, and blood pressure status."

# now we have to convert from mmol to mg to align with CCHS data
# 1mmol sodium = 23mg
# 100mmol (which is what He et al. use, see abstract text) = 2300mg
# so, to get the slope in underdstandable units (/1000mg), we do: (5.8/2300)*1000 = 2.52

effect_sodium_sbp <- 2.52 # mmHg per 1000 mg/day


# 5. PARAMETERS: errors in recall and FFQ data, sodium intake -------------------------

# relying primarily on Freedman et al., Am J Epidemiol 2015;181:473–87

# here we have the relative bias when comparing tools to "truth"
# see discussion page 483
relative_bias_recall <- 0.925 # 24HR gets about 92% of actual truth
relative_bias_ffq <- 0.7

# now, to get SD recall, we need to take the attenuation factor, which is how much
# an association shrinks when using a recall vs the "truth", and back-calculate the SD
# you can learn more about this here: https://www.bmj.com/content/340/bmj.c2289.long
# in Freedman (2015, p.478 & Table 5), attenuation factor (AF) for sodium measured w a
# single 24HR is ~0.2
af_recall <- 0.2

# our recall is modelled as: recall = relative_bias_recall * truth + error
# (relative_bias_recall = 0.925, the recall reports ~92.5% of true intake; Freedman p.483)
# because of the scaling, the AF formula is:
# AF = cov(truth, recall) / var(recall)
#    = relative_bias_recall * var(truth) / (relative_bias_recall^2 * var(truth) + var(error))

# var(truth) is sd_sodium^2. now can solve for var(error)

# 0.2 = relative_bias_recall * sd_sodium^2 / (relative_bias_recall^2 * sd_sodium^2 + var(error))
# relative_bias_recall^2 * sd_sodium^2 + var(error) = relative_bias_recall * sd_sodium^2 / 0.2
# var(error) = relative_bias_recall * sd_sodium^2 / 0.2 - relative_bias_recall^2 * sd_sodium^2

# note: if relative_bias_recall = 1 (no proportional bias), this simplifies to
# var(error) = sd_sodium^2 / 0.2 - sd_sodium^2

# using our known values of sd_sodium = 883 and relative_bias_recall = 0.925
var_error <- relative_bias_recall * sd_sodium^2 / af_recall - relative_bias_recall^2 * sd_sodium^2
sd_recall <- sqrt(var_error)   # 1714 mg/day

# finally, we will invent a number for the SD of FFQs. In general, FFQs are considered to have
# lower variance because they ask about foods consumed over a longer period of time, thus
# theoretically averaging out over day-to-day fluctuations. here, for the purpose of the
# activity, we set a low standard deviation for FFQs to make the lesson clear. in reality this
# number is probably not as low.

sd_ffq <- 200


# 6. simulate one study ---------------------------------------------------

# SBP = sodium part + residual. Independent parts add on the variance scale:
# var(SBP) = (effect * sd_sodium/1000)^2 + var(resid)
# sodium part SD = 2.52 * 883/1000 ≈ 2.2; total SBP variance = sd_sbp^2 = 16^2 = 256
# so sd_resid = sqrt(16^2 - 2.2^2) = sqrt(256 - 4.95) ≈ 15.84, keeping total SBP SD at 16
# note: sodium explains only ~2% of SBP variance (4.95/256)
sd_resid <- sqrt(sd_sbp^2 - (effect_sodium_sbp * sd_sodium / 1000)^2)

simulate_study <- function(n){
  tibble(
    truth = rnorm(n, mean_sodium, sd_sodium), # true usual intake
    sbp = mean_sbp + effect_sodium_sbp*(truth - mean_sodium)/1000 + rnorm(n, 0, sd_resid), # outcome
    # here we simulate measured values from recalls and ffqs, and we set a floor of 0 (cannot ahve negative intake)
    recall = pmax(relative_bias_recall * truth + rnorm(n, 0, sd_recall), 0),
    ffq = pmax(relative_bias_ffq * truth + rnorm(n, 0, sd_ffq), 0)  
  )
}

set.seed(12345)
n <- 1000 # 1000 people in the study, for example
study <- simulate_study(n)


# 7. describe exposures ---------------------------------------------------

study |> 
  pivot_longer(
    c(truth, recall, ffq),
    names_to = "measure",
    values_to = "sodium"
  ) |> 
  summarise(
    mean = mean(sodium),
    sd = sd(sodium),
    median = median(sodium),
    p25 = quantile(sodium, 0.25),
    p75 = quantile(sodium, 0.75),
    r_truth = cor(sodium, study$truth),
    .by = measure
  )

study |>
  pivot_longer(
    c(truth, recall, ffq),
    names_to = "measure",
    values_to = "sodium"
  ) |> 
  ggplot(
    aes(sodium, colour = measure, fill = measure)
  )+
  geom_density(alpha = 0.2) +
  geom_vline(
    data = \(d) summarise(d, mean = mean(sodium), .by = measure),
    aes(xintercept = mean, colour = measure),
    linetype = 2, linewidth = 1,
    show.legend = FALSE
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Sodium (mg/day)", 
    title = "Distribution of simulated sodium intake",
    subtitle = "Bias shifts the centre; variance widens the spread"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA)
  )


# 8. estimate association between sodium and SBP, one time --------------------------

# regression
fit_slope <- function(data, measure){
  x <- data[[measure]]/1000 # mmHg per 1000 mg/day
  fit <- lm(data$sbp ~ x) # regression, systolic blood pressure and sodium (for each of the three measures)
  s <- summary(fit)$coefficients[2, ]
  
  tibble(
    measure,
    est = s[[1]],
    se = s[[2]],
    p = s[[4]]
  )
}

# results
map_dfr(c("truth", "recall", "ffq"), \(m) fit_slope(study, m)) |>
  mutate(lo = est-1.96*se, hi = est+1.96*se)

# visualize lm line
study |>
  pivot_longer(c(truth, recall, ffq), names_to = "measure", values_to = "sodium") |>
  ggplot(aes(sodium, sbp)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~measure, scales = "free_x") +
  labs(
    x = "Sodium (mg/day)", 
    y = "SBP (mmHg)",
    title = "Fitted lines, SBP~sodium"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA)
  )


# 9. repeat and get average output ----------------------------------------

n_sims <- 500

sims <- map_dfr(1:n_sims, \(i){
  d <- simulate_study(n)
  map_dfr(c("truth", "recall", "ffq"), \(m) fit_slope(d, m)) |> mutate(sim = i)
})

sims |> 
  group_by(measure) |> 
  summarise(
    mean_est = mean(est),
    bias = mean(est) - effect_sodium_sbp, # systematic error in the slope
    emp_se = sd(est), # spread of the slope across studies
    rmse     = sqrt(mean((est - effect_sodium_sbp)^2)), # bias and variance combined
    coverage = mean(abs(est - effect_sodium_sbp) <= 1.96 * se),  # % of 95% CIs containing the truth
    power    = mean(p < 0.05) # % of studies detecting the effect
  )

sims |>
  ggplot(
    aes(measure, est, fill = measure)
  ) +
  geom_boxplot() +
  geom_hline(yintercept = effect_sodium_sbp, linetype = 2, colour = "red") +
  labs(
    y = "Estimated slope (mmHg per 1000 mg/day)",
    title = "Sampling distribution of the slope by measurement tool"
  ) + 
  coord_flip()+
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA)
  )
