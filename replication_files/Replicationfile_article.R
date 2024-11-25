##############################################
####  REPLICATION FILE FOR MARTINS (2024)  ###
##############################################

# Loading packages

library(tidyverse)
library(haven)
library(survival)
library(survminer)
library(cmprsk)
library(etm)
library(simPH)
library(coxme)
library(texreg)
library(stargazer)
library(broom.mixed)

# Importing dataset

data <- read_dta("data_martins2024.dta")


############### COX PH ANALYSIS ###################


# Creating survival object 

data <- data[order(data$id_figo, data$duration), ]

data$t <- data$duration

data$t0 <- data$t-1

survobj_cox <-Surv(data$t0, data$t, data$death, type = "counting", origin = 0)

survcheck(survobj_cox ~1, id=id_figo, data=data)


##### Figure 1 - Kaplan-Meier Graph #####

# Fitting the Kaplan-Meier model
km_fit <- survfit(survobj_cox ~ 1, data = data)

# Creating the K-M plot
ggsurvplot(
  km_fit,
  data = data,
  censor = TRUE,
  risk.table = TRUE,
  pval = FALSE,
  conf.int = TRUE,
  xlim = c(0, max(data$t)),
  break.time.by = 50,  
  ggtheme = theme_classic()  +
    theme(
      text = element_text(family = "Times", size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.5)
    ),
  palette = "black",
  font.family = "serif",
  linetype = "solid",
  title = "Kaplan-Meier Survival Curve",
  xlab = "Time",
  ylab = "Survival probability",
  risk.table.height = 0.25,
  legend.labs = "Overall"
)


##### Cox PH regression - Table 1 #####

# Cox regression

cox1 <- coxph(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                sum_govtasks + sum_d + nover_live, 
              data = data, id = id_figo, ties = "efron", robust = T, cluster = id_figo)

summary(cox1)

cox1$loglik

# Testing the proportional hazards assumption

cox.zph(cox1)

phtest.id <-cox.zph(cox1,transform="identity", global=T)
phtest.id

phtest.rnk <-cox.zph(cox1,transform="rank", global=T)
phtest.rnk

phtest.km <-cox.zph(cox1,transform="km", global=T)
phtest.km

### Since there is a violation of the PH assumption, some adjustments are due

# Creating a time-dependent covariate

data$sum_govtasks_time <- data$sum_govtasks * log(data$t)

# Fitting the new Cox model with the time-dependent covariate (This one corresponds to Table 1)

cox2 <- coxph(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                sum_govtasks + sum_govtasks_time + sum_d + nover_live, 
              data = data, id = id_figo, ties = "efron", robust = TRUE, cluster = id_figo)

summary(cox2)

num_failures_cox2 <- cox2$nevent 

stargazer(cox2, type = "text", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_cox2)))

##### Figure 2 - Predicted hazard ratios #####


# Creating a sequence of values for the IV
avg_Oij_seq <- seq(0.01, 0.1, by = 0.01)

# Creating a new data frame for predictions
newdata <- data.frame(avg_Oij = avg_Oij_seq)

# Adding mean values for other variables in the model
for(var in c("lnnmem", "sum_issues", "sum_govtasks", "sum_govtasks_time", "sum_d", "nover_live", "niche_width")) {
  newdata[[var]] <- mean(data[[var]], na.rm = TRUE)
}

# Calculating predicted hazard ratios and confidence intervals
pred <- predict(cox2, newdata = newdata, type = "terms", se.fit = TRUE)
hr <- exp(pred$fit[, "avg_Oij"])
lower_ci <- exp(pred$fit[, "avg_Oij"] - 1.96 * pred$se.fit[, "avg_Oij"])
upper_ci <- exp(pred$fit[, "avg_Oij"] + 1.96 * pred$se.fit[, "avg_Oij"])

# Creating a new data frame for plotting the graph
plot_data <- data.frame(avg_Oij = avg_Oij_seq, hr = hr, lower_ci = lower_ci, upper_ci = upper_ci)

# Creating the plot
ggplot(plot_data, aes(x = avg_Oij, y = hr)) +
  geom_line(color = "black") +
  geom_point(color = "black", shape = 16, size = 2) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.002, color = "black") +
  labs(x = "Degree of institutional overlap", y = "Predicted hazard ratio", 
       title = "Adjusted predictions with 95% CIs") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, by = 0.2), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0.01, 0.1, by = 0.01), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, size = 1) +
  geom_vline(xintercept = 0, size = 1)


#### Figure 3 - Survival estimates ####

# Calculating 20th and 80th percentiles of avg_Oij
p20_avg_Oij <- quantile(data$avg_Oij, 0.20, na.rm = TRUE)
p80_avg_Oij <- quantile(data$avg_Oij, 0.80, na.rm = TRUE)

# Creating a new dataset for prediction
newdata <- data.frame(
  avg_Oij = c(p20_avg_Oij, p80_avg_Oij),
  lnnmem = mean(data$lnnmem, na.rm = TRUE),
  sum_issues = mean(data$sum_issues, na.rm = TRUE),
  sum_govtasks = mean(data$sum_govtasks, na.rm = TRUE),
  sum_govtasks_time = mean(data$sum_govtasks_time, na.rm = TRUE),
  sum_d = mean(data$sum_d, na.rm = TRUE),
  nover_live = mean(data$nover_live, na.rm = TRUE)
)

# Generating new survival curves
fit <- survfit(cox2, newdata = newdata)

# Plotting the new survival curves at different percentiles of the IV
ggsurvplot(
  fit,
  data = newdata,
  conf.int = FALSE,
  legend.labs = c("Institutional overlap = Low", "Institutional overlap = High"),
  palette = c("blue", "red"),
  title = "Survival estimates at different levels of institutional overlap",
  xlab = "Analysis time",
  ylab = "Survival",
  ylim = c(0.50, 1.00),
  break.time.by = 20,
  risk.table = FALSE,
  ggtheme = theme_classic() +
    theme(
      text = element_text(family = "Times", size = 12),
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.5)
    ),
  font.family = "serif"
)



########## COMPETING RISKS ANALYSIS ############



# Removing the IOs without information about the type of death

remove_IOs <- c("AATPO", "AFRAND", "AfrOPDA", "AFTE", "AP", "APPA", "ATO", "BALTBAT", 
                "CAAD", "EIPA", "EURATOM", "IAFC", "IALong", "IBPMP", "IPedI", 
                "IPhyL", "ITRO", "LGIDA", "MCPTTC", "MWN", "PCB", "RCAELA", "SAMI", 
                "SugU", "VALDIVIA")

data_cr <- data[!data$ioname %in% remove_IOs, ]


# Creating the variable containing different types of death

data_cr$status <- as.factor(data_cr$termination_type_num)
class(data_cr$status)
table(data_cr$status)

# Creating the survival object

survobj_cr <- Surv(data_cr$t0, data_cr$t, data_cr$status, type = "counting", origin = 0)

ck <- survcheck(survobj_cr ~ 1, data = data_cr, id = id_figo)

ck$transitions


#### Cumulative Incidence Function plot - Figure 4 ####

# Fitting the Aalen-Johansen estimator with a competing risks model
curves_cr <- survfit(Surv(t0, t, status, type = "counting", origin = 0) ~ 1, 
                     id = id_figo, data = data_cr)

# Extracting time and state probabilities
fig_time <- rep(curves_cr$time, 3)  
fig_ps1 <- curves_cr$pstate[,2]  # Contract probability
fig_ps2 <- curves_cr$pstate[,3]  # Transformation probability
fig_ps3 <- curves_cr$pstate[,4]  # Inertia probability
f_allp <- c(fig_ps1, fig_ps2, fig_ps3)

# Creating termination labels
termination_type <- c(rep("Contract", length(fig_ps1)), 
                      rep("Transformation", length(fig_ps2)), 
                      rep("Inertia", length(fig_ps3)))

# Plotting the graph

fig_data <- data.frame(time = fig_time, probability = f_allp, Termination = termination_type)

ggplot(data = fig_data, aes(x = time, y = probability, color = Termination)) +
  geom_step(direction = "hv", size = 1) +
  labs(x = "Time", y = "Probability in State", title = "Cumulative Incidence Function for Competing Risks") +
  theme_classic() +
  theme(
    text = element_text(family = "Times", size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  ) +
  scale_color_manual(values = c("blue", "red", "green")) +  
  xlim(0, 200) +
  ylim(0, 0.3)



##### Fine-Gray subdistribution hazards regression #####


# Creating Fine-Gray datasets for each competing risk

fg_contract <- finegray(Surv(t, status) ~ avg_Oij + lnnmem + sum_issues + 
                          sum_govtasks + sum_d + nover_live + id_figo,
                        data = data_cr, id = id_figo, etype = "1", na.action = na.omit, count = "count")

fg_transformation <- finegray(Surv(t, status) ~ avg_Oij + lnnmem + sum_issues + 
                                sum_govtasks + sum_d + nover_live + id_figo,
                              data = data_cr, id = id_figo, etype = "2", na.action = na.omit, count = "count")

fg_inertia <- finegray(Surv(t, status) ~ avg_Oij + lnnmem + sum_issues + 
                         sum_govtasks + sum_d + nover_live + id_figo,
                       data = data_cr, id = id_figo, etype = "3", na.action = na.omit, count = "count")


# Adjusting start-stop times for each dataset

correct_times <- function(fg_data) {
  fg_data$fgstart <- ifelse(fg_data$count == 0, fg_data$fgstop - 1, fg_data$fgstart)
  return(fg_data)
}

cor_fg_contract <- correct_times(fg_contract)
cor_fg_transformation <- correct_times(fg_transformation)
cor_fg_inertia <- correct_times(fg_inertia)


# Weighted Cox proportional hazards model for each competing risk

contract_fg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                           sum_govtasks + sum_d + nover_live, 
                         data = cor_fg_contract, ties = "efron", weights = fgwt, 
                         id = id_figo, robust = TRUE, x = TRUE)

transformation_fg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                                 sum_govtasks + sum_d + nover_live, 
                               data = cor_fg_transformation, ties = "efron", weights = fgwt, 
                               id = id_figo, robust = TRUE, x = TRUE)

inertia_fg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                          sum_govtasks + sum_d + nover_live, 
                        data = cor_fg_inertia, ties = "efron", weights = fgwt, 
                        id = id_figo, robust = TRUE, x = TRUE)

summary(contract_fg)
summary(transformation_fg)
summary(inertia_fg)


# Testing the PH assumption

cox.zph(contract_fg)
cox.zph(transformation_fg)
cox.zph(inertia_fg)

# Correcting offending variables

cor_fg_transformation$fgstop_time <- log(cor_fg_transformation$fgstop)
cor_fg_transformation$sum_issues_time <- cor_fg_transformation$fgstop * cor_fg_transformation$sum_issues
cor_fg_transformation$sum_govtasks_time <- cor_fg_transformation$fgstop * cor_fg_transformation$sum_govtasks

transformation_fg_time <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + sum_issues_time +
                                      sum_govtasks + sum_govtasks_time + sum_d + nover_live, 
                                    data = cor_fg_transformation, ties = "efron", weights = fgwt, 
                                    id = id_figo, robust = TRUE, x = TRUE)

summary(transformation_fg_time)

num_failures_cfg <- contract_fg$nevent
num_failures_tfg <- transformation_fg_time$nevent
num_failures_ifg <- inertia_fg$nevent


#### Predicted CIF of Contract-based terminations by level of institutional overlap - Figure 5 ####

# Calculating 20th and 80th percentiles of avg_Oij

p20_avg_Oij_cr <- quantile(data_cr$avg_Oij, 0.20, na.rm = TRUE)
p80_avg_Oij_cr <- quantile(data_cr$avg_Oij, 0.80, na.rm = TRUE)

# Create new dataset for predictions
newdata_cr <- data.frame(
  avg_Oij = c(p20_avg_Oij_cr, p80_avg_Oij_cr),
  lnnmem = mean(data_cr$lnnmem, na.rm = TRUE),
  sum_issues = mean(data_cr$sum_issues, na.rm = TRUE),
  sum_govtasks = mean(data_cr$sum_govtasks, na.rm = TRUE),
  sum_d = mean(data_cr$sum_d, na.rm = TRUE),
  nover_live = mean(data_cr$nover_live, na.rm = TRUE)
)


fg_model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + 
                    sum_issues + sum_govtasks + sum_d + nover_live,
                  data = cor_fg_contract, weights = fgwt)

# Get the baseline hazard
bh <- basehaz(fg_model)

# Calculate linear predictors for new data
lp <- predict(fg_model, newdata = newdata_cr, type = "lp")

# Calculate CIF for each covariate pattern
cif_low <- 1 - exp(-bh$hazard * exp(lp[1]))
cif_high <- 1 - exp(-bh$hazard * exp(lp[2]))

# Create data frame for plotting
plot_data <- data.frame(
  time = c(bh$time, bh$time),
  cif = c(cif_low, cif_high),
  overlap = factor(rep(c("Institutional overlap = low", 
                         "Institutional overlap = high"), 
                       each = length(bh$time)))
)

# Create the plot using ggplot2
ggplot(plot_data, aes(x = time, y = cif, color = overlap)) +
  geom_line(size = 1) +
  labs(x = "Time",
       y = "Probability in State",
       title = "Predicted Cumulative Incidence of Contract-based \n terminations by level of institutional overlap",
       color = "Institutional Overlap") +
  scale_color_manual(values = c("blue", "red")) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 0.25), 
                     breaks = seq(0, 0.25, by = 0.05)) +
  scale_x_continuous(breaks = seq(0, max(bh$time), by = 50))


##### Cause-specific hazards ######


data_cr$status <- as.factor(data_cr$termination_type_num)

# Creating the survival object


survobj_csh <- Surv(data_cr$t0, data_cr$t, data_cr$status, type = "counting", origin = 0)
ck_csh <- survcheck(survobj_csh ~ 1, data = data_cr, id = id_figo)
print(ck_csh$transitions)

# Cox models for each competing risk with cause-specific hazards

contract_csh <- coxph(Surv(t0, t, status == "1") ~ avg_Oij + lnnmem + sum_issues + 
                            sum_govtasks + sum_d + nover_live,
                          data = data_cr, ties = "efron", id = id_figo, robust = TRUE)

transformation_csh <- coxph(Surv(t0, t, status == "2") ~ avg_Oij + lnnmem + sum_issues + 
                                  sum_govtasks + sum_d + nover_live,
                                data = data_cr, ties = "efron", id = id_figo, robust = TRUE)

inertia_csh <- coxph(Surv(t0, t, status == "3") ~ avg_Oij + lnnmem + sum_issues + 
                           sum_govtasks + sum_d + nover_live,
                         data = data_cr, ties = "efron", id = id_figo, robust = TRUE)

summary(contract_csh)
summary(transformation_csh)
summary(inertia_csh)

# Testing the PH assumption

cox.zph(contract_csh)
cox.zph(transformation_csh)
cox.zph(inertia_csh)

# Correcting offending variables

data_cr$log_time <- log(data_cr$t) 
data_cr$sum_issues_time <- data_cr$sum_issues * data_cr$log_time


transformation_csh_time <- coxph(Surv(t0, t, status == "2") ~ avg_Oij + lnnmem + sum_issues + sum_issues_time +
                              sum_govtasks + sum_d + nover_live,
                            data = data_cr, ties = "efron", id = id_figo, robust = TRUE)

summary(transformation_csh_time)

num_failures_ccs <- contract_csh$nevent
num_failures_tcs <- transformation_csh_time$nevent
num_failures_ics <- inertia_csh$nevent

##
stargazer(contract_csh, contract_fg, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_ccs, num_failures_cfg)))

stargazer(transformation_csh_time, transformation_fg_time, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_tcs, num_failures_tfg)))

stargazer(inertia_csh, inertia_fg, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_ics, num_failures_ifg)))

##########################################
########     ONLINE APPENDIX      ########
##########################################



# Creating summary statistics (Table A1)

summary_table <- data.frame(
  Variable = c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d"),
  Mean = sapply(data[, c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d")], mean, na.rm = TRUE),
  SD = sapply(data[, c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d")], sd, na.rm = TRUE),
  Median = sapply(data[, c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d")], median, na.rm = TRUE),
  Min = sapply(data[, c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d")], min, na.rm = TRUE),
  Max = sapply(data[, c("death", "duration", "avg_Oij", "nmem", "nover_live", "sum_issues", "sum_govtasks", "sum_d")], max, na.rm = TRUE)
)

stargazer(summary_table, type = "text", summary = FALSE, title = "Summary Statistics")

# Creating datasets for Tables A2 to A4

cow <- read_dta("C:/Users/joao_/Downloads/igo_year_format_3.dta")

contract <- data_cr |> filter(termination_type_num == 1)
transformation <- data_cr |> filter(termination_type_num == 2)
inertia <- data_cr |> filter(termination_type_num == 3)

cow_reduced <- cow %>%
  group_by(ionum) %>%
  slice(1) %>%
  ungroup() %>%
  select(ionum, longorgname)

contract <- contract %>%
  left_join(cow_reduced, by = c("id_figo" = "ionum"))

contract <- contract |> select(ioname, longorgname, strt_yr,end_yr, Fate_death_reason)

contract <- contract %>%
  mutate(Fate_death_reason = case_when(
    Fate_death_reason == "DISS" ~ "Dissolution",
    Fate_death_reason == "EXP" ~ "Expiry",
    TRUE ~ Fate_death_reason  
  ))

transformation <- transformation %>%
  left_join(cow_reduced, by = c("id_figo" = "ionum"))

transformation <- transformation |> select(ioname, longorgname, strt_yr,end_yr, Fate_death_reason)

transformation <- transformation %>%
  mutate(Fate_death_reason = case_when(
    Fate_death_reason == "SUCC" ~ "Succession",
    Fate_death_reason == "ABS" ~ "Absorption",
    Fate_death_reason == "MER" ~ "Merger",
    Fate_death_reason == "EXP/ABS" ~ "Absorption",
    TRUE ~ Fate_death_reason  
  ))

inertia <- inertia %>%
  left_join(cow_reduced, by = c("id_figo" = "ionum"))

inertia <- inertia |> select(ioname, longorgname, strt_yr,end_yr, Fate_death_reason)


inertia <- inertia  %>%
  mutate(Fate_death_reason = case_when(
    Fate_death_reason == "DES" ~ "Desuetude",
    TRUE ~ Fate_death_reason  
  ))


# Plotting the Schoenfeld residuals graphs for variables in the uncorrected model 1 (Figure A1)

par(mfrow = c(3, 2)) 

plot(phtest.km[1], ylab = "Residuals for Institutional overlap")
plot(phtest.km[2], ylab = "Residuals for Number of members (logged)")
plot(phtest.km[3], ylab = "Residuals for Policy scope")
plot(phtest.km[4], ylab = "Residuals for Governance mandate")
plot(phtest.km[5], ylab = "Residuals for Institutionalization")
plot(phtest.km[6], ylab = "Residuals for Number of IOs alive")

par(mfrow = c(1, 1))

# Results from Table 1 using the Breslow method for ties (Table A5)

cox2_breslow <- coxph(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                sum_govtasks + sum_govtasks_time + sum_d + nover_live, 
              data = data, id = id_figo, ties = "breslow", robust = T, cluster = id_figo)

summary(cox2_breslow)

cox.zph(cox2_breslow)

stargazer(cox2_breslow, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_cox2)))

# Results from Table 2 using the Breslow method for ties (Table A6)


# Weighted Cox proportional hazards model for each competing risk with Breslow ties

contract_fg_breslow <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                               sum_govtasks + sum_d + nover_live, 
                             data = cor_fg_contract, ties = "breslow", weights = fgwt, 
                             id = id_figo, robust = TRUE, x = TRUE)

transformation_fg_breslow <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                                     sum_govtasks + sum_d + nover_live, 
                                   data = cor_fg_transformation, ties = "breslow", weights = fgwt, 
                                   id = id_figo, robust = TRUE, x = TRUE)

inertia_fg_breslow <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + 
                              sum_govtasks + sum_d + nover_live, 
                            data = cor_fg_inertia, ties = "breslow", weights = fgwt, 
                            id = id_figo, robust = TRUE, x = TRUE)

summary(contract_fg_breslow)
summary(transformation_fg_breslow)
summary(inertia_fg_breslow)


# Testing the PH assumption for Breslow models

cox.zph(contract_fg_breslow)
cox.zph(transformation_fg_breslow)
cox.zph(inertia_fg_breslow)


# Correcting offending variables in transformation model

cor_fg_transformation$fgstop_time <- log(cor_fg_transformation$fgstop)
cor_fg_transformation$sum_issues_time <- cor_fg_transformation$fgstop * cor_fg_transformation$sum_issues
cor_fg_transformation$sum_govtasks_time <- cor_fg_transformation$fgstop * cor_fg_transformation$sum_govtasks

transformation_fg_breslow_time <- coxph(Surv(fgstart, fgstop, fgstatus) ~ avg_Oij + lnnmem + sum_issues + sum_issues_time +
                                          sum_govtasks + sum_govtasks_time + sum_d + nover_live, 
                                        data = cor_fg_transformation, ties = "breslow", weights = fgwt, 
                                        id = id_figo, robust = TRUE, x = TRUE)

summary(transformation_fg_breslow_time)


# Cause-specific hazards with Breslow ties

contract_csh_breslow <- coxph(Surv(t0, t, status == "1") ~ avg_Oij + lnnmem + sum_issues + 
                                sum_govtasks + sum_d + nover_live,
                              data = data_cr, ties = "breslow", id = id_figo, robust = TRUE)

transformation_csh_breslow <- coxph(Surv(t0, t, status == "2") ~ avg_Oij + lnnmem + sum_issues + 
                                      sum_govtasks + sum_d + nover_live,
                                    data = data_cr, ties = "breslow", id = id_figo, robust = TRUE)

inertia_csh_breslow <- coxph(Surv(t0, t, status == "3") ~ avg_Oij + lnnmem + sum_issues + 
                               sum_govtasks + sum_d + nover_live,
                             data = data_cr, ties = "breslow", id = id_figo, robust = TRUE)

summary(contract_csh_breslow)
summary(transformation_csh_breslow)
summary(inertia_csh_breslow)


# Testing the PH assumption for Breslow models in competing risks

cox.zph(contract_csh_breslow)
cox.zph(transformation_csh_breslow)
cox.zph(inertia_csh_breslow)


# Correcting offending variables in transformation model for Breslow

data_cr$log_time <- log(data_cr$t) 
data_cr$sum_issues_time <- data_cr$sum_issues * data_cr$log_time

transformation_csh_breslow_time <- coxph(Surv(t0, t, status == "2") ~ avg_Oij + lnnmem + sum_issues + sum_issues_time +
                                           sum_govtasks + sum_d + nover_live,
                                         data = data_cr, ties = "breslow", id = id_figo, robust = TRUE)

summary(transformation_csh_breslow_time)

num_failures_ccs <- contract_csh_breslow$nevent
num_failures_tcs <- transformation_csh_breslow_time$nevent
num_failures_ics <- inertia_csh_breslow$nevent


# Generating tables with stargazer

stargazer(contract_csh_breslow, contract_fg_breslow, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_ccs, num_failures_cfg)))

stargazer(transformation_csh_breslow_time, transformation_fg_breslow_time, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_tcs, num_failures_tfg)))

stargazer(inertia_csh_breslow, inertia_fg_breslow, type = "latex", omit.stat = c("rsq", "max.rsq"),
          add.lines = list(c("Number of failures", num_failures_ics, num_failures_ifg)))



# Conditional linear coefficient for governance tasks variable (Figure A2)

beta_govtasks <- cox2$coefficients[4]  
beta_interaction <- cox2$coefficients[5] 

sim.t <- seq(1, max(data$t), 1)  
sim.lnt <- log(sim.t)

cond_coef <- beta_govtasks + beta_interaction * sim.lnt

plot_data <- data.frame(
  Time = sim.t,
  Coefficient = cond_coef
)

ggplot(data = plot_data, aes(x = Time, y = Coefficient)) +
  geom_line(color = "purple4", size = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Conditional Linear Coefficient",
    title = "Time-Varying Effect of Government Tasks"
  ) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )


# Simulated first differences for governance tasks (Figure A3)

sim_govtasks <- coxsimtvc(
  obj = cox2,                      
  b = "sum_govtasks",             
  btvc = "sum_govtasks_time",    
  qi = "First Difference",        
  Xj = 1,                         
  tfun = "log",                   
  from = 0,                       
  to = max(data$t),             
  by = 30                         
)

simGG(
  sim_govtasks, 
  type = "ribbons",              
  lsize = 1,                     
  legend = FALSE,                
  alpha = 0.3,                   
  xlab = "Time",                 
  ylab = "First Difference",     
  title = "Simulated First Differences in Government Tasks Effect"
) +
  theme_bw() +                   
  theme(text = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )
    

# Model with a shared frailty on function (Table A7)

cox_function_frailty <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                                sum_govtasks + sum_d + nover_live + (1 | function_c),
                              data = data)

summary(cox_function_frailty)

cox.zph(cox_function_frailty) 

cox_function_frailty_time <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                                sum_govtasks + sum_govtasks_time + sum_d + nover_live + (1 | function_c),
                              data = data)

summary(cox_function_frailty_time)


# Model with a shared frailty on scope (Table A8)

cox_scope_frailty <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                             sum_govtasks + sum_d + nover_live + (1 | scope_c),
                           data = data)

summary(cox_scope_frailty)

cox.zph(cox_scope_frailty) 

cox_function_scope_time <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                                     sum_govtasks + sum_govtasks_time + sum_d + nover_live + (1 | scope_c),
                                   data = data)

summary(cox_function_scope_time)


# Model with a shared frailty on region (Table A9)

cox_region_frailty <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                              sum_govtasks + sum_d + nover_live + (1 | region_c),
                            data = data)

summary(cox_region_frailty)

cox.zph(cox_region_frailty) 

cox_function_region_time <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                                   sum_govtasks + sum_govtasks_time + sum_d + nover_live + (1 | region_c),
                                 data = data)

summary(cox_function_region_time)

# Model with a shared frailty on membership format (Table A10)

cox_memb_form_frailty <- coxme(survobj_cox ~ avg_Oij + lnnmem + sum_issues + 
                                 sum_govtasks + sum_d + nover_live + (1 | memb_form),
                               data = data)

summary(cox_memb_form_frailty)

cox.zph(cox_memb_form_frailty) 
