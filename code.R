
library(tidyverse)
library(survival)
library(flexsurv)
library(msm)
library(survMisc)
library(muhaz)
library(kableExtra)
library(pander)
library(survminer)
library(data.table)
library(mstate)
library(cmprsk)

bmt_df <- read_csv(file = "~/BIOST537/BIOST537-Project/bmt.csv")

source("fitparametric.R")

# create age at diagnosis variable in years
bmt_df$agediagnosis <- bmt_df$age - ((bmt_df$waittime)/365)

# create age at event (death, relapse or censoring) in years
bmt_df$ageevent <- bmt_df$agediagnosis + ((bmt_df$tdfs)/365)


# Creating Survival Object w/ delayed entry
s_bmt_de <- with(bmt_df, Surv(agediagnosis, ageevent, deltadfs==1))

# Creating Suvival Object w/o delayed entry
s_bmt <- with(bmt_df, Surv(tdfs, deltadfs == 1))

sfit_bmt <- survfit(s_bmt ~ 1, data = bmt_df, conf.type = "log-log")

sfit_weibull <- flexsurvreg(s_bmt ~ 1, dist = "weibull", data = bmt_df)

fitparametric(s_bmt, dist = "weibull")

sfit_gg <- flexsurvreg(s_bmt ~ 1, dist = "gengamma", data = bmt_df)

# mean follow-up time
mean(bmt_df$tdfs)

# proportion of censored of observations
1 - mean(bmt_df$deltadfs)

# Kaplan Meier Plot
km_dir1 <- ggsurvplot(sfit_bmt,
                      conf.int = TRUE,
                      surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate")

weibull_dir1_plot <- ggflexsurvplot(sfit_weibull, conf.int = F,
                                    surv.median.line = "hv") + 
  labs(title = "Weibull Estimate w/ Kaplan-Meier Survival Estimate")

ggamma_dir1_plot <- ggflexsurvplot(sfit_gg, conf.int = F,
                                   surv.median.line = "hv") + 
  labs(title = "GGamma Estimate w/ Kaplan-Meier Survival Estimate")

# Median Survival time, C.I. for median surivival time, and other summary stats
summary(sfit_bmt)$table

# Parametric fits
fit_weibull <- fitparametric(s_bmt, dist = "weibull")
fit_ggamma <- fitparametric(s_bmt, dist = "gengamma")

# Log likelihood p-value to compare weibull and gengamma
1 - pchisq(q = 2*(-650.19 - (-657.77 )), df = 1)

# p-value from likelihood ratio test comparing the weibull to ggamma
1 - pchisq(2 * (fit_ggamma$loglik - fit_weibull$loglik), df=1)

bmt_df_summary <- bmt_df %>%
  summarise(mean_age = round(mean(age), 3),
            sd_age = round(sd(age), 3),
            count_males = sum(male),
            prop_males = round(sum(male) / length(male), 3),
            count_females = length(male) - sum(male),
            prop_females = round((length(male) - sum(male)) / length(male), 3),
            count_cmv = sum(cmv),
            prop_cmv = round(sum(cmv) / length(cmv), 3),
            count_mtx = sum(mtx),
            prop_mtx = round(sum(mtx) / length(mtx), 3),
            count_hospital = sum(hospital),
            mean_donor_age = round(mean(donorage), 3),
            sd_donor_age = round(sd(donorage), 3),
            count_donor_males = sum(donormale),
            prop_donor_males = round(sum(donormale) / length(donormale), 3),
            count_donor_cmv = sum(donorcmv),
            prop_donor_cmv = round(sum(donorcmv) / length(donorcmv), 3))

# Transposed plot for easier viewing
bmt_df_summary_t <- transpose(bmt_df_summary)

kable(bmt_df_summary)

# Survival Object based on disease subgrouping
sfit_bmt_byDisgroup <- survfit(s_bmt ~ disgroup,
                               data = bmt_df,
                               conf.type = "log-log")

print(sfit_bmt_byDisgroup)


pander(survdiff(s_bmt ~ disgroup, data = bmt_df))

# Kaplan Meier curve
disgroup_dir2_plot <- ggsurvplot(sfit_bmt_byDisgroup,
                                 pval = TRUE,
                                 pval.coord = c(2100, 1),
                                 conf.int = F,
                                 surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by Disease Group")

#table 1: columns = disease groups, rows= baseline characteristics
bmt_df_dis <- bmt_df %>%
  group_by(disgroup) %>%
  summarise(mean_age = round(mean(age), 3),
            sd_age = round(sd(age), 3),
            count_males = sum(male),
            prop_males = round(sum(male) / length(male), 3),
            count_females = length(male) - sum(male),
            prop_females = round((length(male) - sum(male)) / length(male), 3),
            count_cmv = sum(cmv),
            prop_cmv = round(sum(cmv) / length(cmv), 3),
            count_mtx = sum(mtx),
            prop_mtx = round(sum(mtx) / length(mtx), 3),
            count_hospital = sum(hospital),
            mean_donor_age = round(mean(donorage), 3),
            sd_donor_age = round(sd(donorage), 3),
            count_donor_males = sum(donormale),
            prop_donor_males = round(sum(donormale) / length(donormale), 3),
            count_donor_cmv = sum(donorcmv),
            prop_donor_cmv = round(sum(donorcmv) / length(donorcmv), 3))

# Transposed plot for easier viewing
bmt_df_dis_t <- transpose(bmt_df_dis)
colnames(bmt_df_dis_t) <- c("Disease Group 1",
                            "Disease Group 2",
                            "Disease Group 3")
rownames(bmt_df_dis_t) <- colnames(bmt_df_dis)
bmt_df_dis_t <- bmt_df_dis_t[c(-1),]
rownames(bmt_df_dis_t) <- str_replace(rownames(bmt_df_dis_t), "_", " ")
rownames(bmt_df_dis_t) <- str_replace(rownames(bmt_df_dis_t), "_", "")

kableExtra::kable(bmt_df_dis_t)

# For test statistics
test_stats_disgroup <- comp(ten(sfit_bmt_byDisgroup))


# Survival Object based on "When we was FAB" classification
sfit_bmt_byFAB <- survfit(s_bmt ~ fab,
                          data = bmt_df,
                          conf.type = "log-log")

print(sfit_bmt_byFAB)

pander(survdiff(s_bmt ~ fab, data = bmt_df))

FAB_dir2_plot <- ggsurvplot(sfit_bmt_byFAB,
                            pval = TRUE,
                            pval.coord = c(2100, 1),
                            conf.int = F,
                            surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by FAB Group")

#table 2: columns = FAB classifications, rows = baseline characteristics
bmt_df_fab <- bmt_df %>%
  group_by(fab) %>%
  summarise(mean_age = round(mean(age), 3),
            sd_age = round(sd(age), 3),
            count_males = sum(male),
            prop_males = round(sum(male) / length(male), 3),
            count_females = length(male) - sum(male),
            prop_females = round((length(male) - sum(male)) / length(male), 3),
            count_cmv = sum(cmv),
            prop_cmv = round(sum(cmv) / length(cmv), 3),
            count_mtx = sum(mtx),
            prop_mtx = round(sum(mtx) / length(mtx), 3),
            count_hospital = sum(hospital),
            mean_donor_age = round(mean(donorage), 3),
            sd_donor_age = round(sd(donorage), 3),
            count_donor_males = sum(donormale),
            prop_donor_males = round(sum(donormale) / length(donormale), 3),
            count_donor_cmv = sum(donorcmv),
            prop_donor_cmv = round(sum(donorcmv) / length(donorcmv), 3))

# Transposed plot for easier viewing
bmt_df_fab_t <- transpose(bmt_df_fab)
colnames(bmt_df_fab_t) <- c("FAB Classification 0",
                            "FAB Classification 1")
rownames(bmt_df_fab_t) <- colnames(bmt_df_fab)
bmt_df_fab_t <- bmt_df_fab_t[c(-1),]
rownames(bmt_df_fab_t) <- str_replace(rownames(bmt_df_fab_t), "_", " ")
rownames(bmt_df_fab_t) <- str_replace(rownames(bmt_df_fab_t), "_", "")

kable(bmt_df_fab_t)

test_stats_fabgroup <- comp(ten(sfit_bmt_byFAB))

pchisq(q = (2.6559)^2, df=1, lower.tail=FALSE)

chisq.test(bmt_df$fab, bmt_df$age)
chisq.test(bmt_df$fab, bmt_df$male)
chisq.test(bmt_df$fab, bmt_df$cmv)
chisq.test(bmt_df$fab, bmt_df$hospital)
chisq.test(bmt_df$fab, bmt_df$mtx)
chisq.test(bmt_df$fab, bmt_df$donorage)
chisq.test(bmt_df$fab, bmt_df$donormale)
chisq.test(bmt_df$fab, bmt_df$donorcmv)


# Survival Object based on sex subgrouping
sfit_bmt_byMale <- survfit(s_bmt ~ male,
                           data = bmt_df,
                           conf.type = "log-log")

# For log-rank test p-value
pander(survdiff(s_bmt ~ male, data = bmt_df))

bymale_dir3_plot <- ggsurvplot(sfit_bmt_byMale,
                               pval = TRUE,
                               pval.coord = c(2100, 1),
                               conf.int = F,
                               surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by Sex")


# Survival Object based on CMV subgrouping
sfit_bmt_byCMV <- survfit(s_bmt ~ cmv,
                          data = bmt_df,
                          conf.type = "log-log")

# For log-rank test p-value
pander(survdiff(s_bmt ~ cmv, data = bmt_df))

byCMV_dir3_plot <- ggsurvplot(sfit_bmt_byCMV,
                              pval = TRUE,
                              pval.coord = c(2100, 1),
                              conf.int = F,
                              surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by CMV")


# Survival Object based on sex subgrouping of donor
sfit_bmt_byDonerMale <- survfit(s_bmt ~ donormale,
                                data = bmt_df,
                                conf.type = "log-log")

# For log-rank test p-value
pander(survdiff(s_bmt ~ donormale, data = bmt_df))

byDonerMale_dir3_plot <- ggsurvplot(sfit_bmt_byDonerMale,
                                    pval = TRUE,
                                    pval.coord = c(2100, 1),
                                    conf.int = F,
                                    surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by Donor Sex")


# Survival Object based on CMV subgrouping of donor
sfit_bmt_byDonerCMV <- survfit(s_bmt ~ donorcmv,
                               data = bmt_df,
                               conf.type = "log-log")

# For log-rank test p-value
pander(survdiff(s_bmt ~ donorcmv, data = bmt_df))

byDonerCMV_dir3_plot <- ggsurvplot(sfit_bmt_byDonerCMV,
                                   pval = TRUE,
                                   pval.coord = c(2100, 1),
                                   conf.int = F,
                                   surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by Donor CMV")


# Survival Object based on hospital subgrouping
sfit_bmt_byHospital <- survfit(s_bmt ~ hospital,
                               data = bmt_df,
                               conf.type = "log-log")

sfit_bmt_byHospital

# For log-rank test p-value
pander(survdiff(s_bmt ~ hospital, data = bmt_df))

byHosptial_dir3_plot <- ggsurvplot(sfit_bmt_byHospital,
                                   pval = TRUE,
                                   pval.coord = c(2100, 1),
                                   conf.int = F,
                                   surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by Hospital")

byHosptial_dir3_plot

bmt_df %>% group_by(hospital) %>% summarise(mtx = sum(mtx)) %>% kable()

# Survival Object based on mtx subgrouping
sfit_bmt_byMTX <- survfit(s_bmt ~ mtx,
                          data = bmt_df,
                          conf.type = "log-log")

# For log-rank test p-value
pander(survdiff(s_bmt ~ mtx, data = bmt_df))

byMTX_dir3_plot <- ggsurvplot(sfit_bmt_byMTX,
                              pval = TRUE,
                              pval.coord = c(2100, 1),
                              conf.int = F,
                              surv.median.line = "hv") + 
  labs(title = "Kaplan-Meier survival estimate, by MTX")

#create long dataset with time-varying GVHD variable
bmt.tvc <- tmerge(data1=bmt_df, data2=bmt_df, id=id,
                  dfs=event(tdfs, deltadfs), gvhd.tv=tdc(ta))

#create survival object for disease-free survival with time-varying GVHD variable
surv.bmt.tv = with(bmt.tvc,
                   Surv(tstart, tstop, dfs))

summary(coxph(surv.bmt.tv ~ gvhd.tv + age + cmv + donorcmv + strata(hospital),
              data=bmt.tvc))

#Is it associated with a decreased risk of relapse?
#create long dataset for relapse with time-varying GVHD variable
bmtrelapse.tvc <- tmerge(data1 = bmt_df, data2 = bmt_df, id=id,
                         dfs=event(tdfs, deltar), gvhd.tv=tdc(ta))

#create survival object for relapse with time varying GVHD variable
surv.relapse.tv = with(bmtrelapse.tvc, Surv(tstart, tstop, dfs))

#competing risk
attach(bmtrelapse.tvc)

cov.matrix <- model.matrix(~gvhd.tv + age + cmv + donorcmv + strata(hospital))

head(cov.matrix)

#relapse risk: relapse
result.relapse.crr <- crr((bmtrelapse.tvc$tstop-bmtrelapse.tvc$tstart),
                          deltar, cov1 = cov.matrix[, -1], failcode=1)

summary(result.relapse.crr)

#competing risk: death
result.death.crr <- crr((bmt.tvc$tstop - bmt.tvc$tstart), deltadfs,
                        cov1=cov.matrix[,-1], failcode=1)

summary(result.death.crr)


#survival object for those who develop aGVHD
s.dfsgvhd <- with(bmt_df[bmt_df$deltaa==1, ], Surv(tdfs, deltadfs==1))

#by patient sex
survfit.bymale <- survfit(s.dfsgvhd ~ male,
                          data = bmt_df[bmt_df$deltaa==1, ],
                          conf.type = "log-log")

survdiff(formula=s.dfsgvhd~male, data=bmt_df[bmt_df$deltaa==1, ])

bymale_dir5_plot <- ggsurvplot(survfit.bymale,
                               conf.int = F,
                               pval = TRUE,
                               pval.coord = c(1700, 1),
                               xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier GVHD survival est,",
       subtitle = "by sex in patients with aGVHD")

#by donor sex
survfit.bydonormale <- survfit(s.dfsgvhd ~ donormale,
                               data = bmt_df[bmt_df$deltaa==1, ],
                               conf.type = "log-log")

survdiff(formula=s.dfsgvhd~donormale, data=bmt_df[bmt_df$deltaa==1, ])

bydonormale_dir5_plot <- ggsurvplot(survfit.bydonormale,
                                    conf.int = F,
                                    pval = TRUE,
                                    pval.coord = c(1700, 1),
                                    xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier GVHD survival est,",
       subtitle = "by donor sex in patients with aGVHD")

#by cmv
survfit.bycmv <- survfit(s.dfsgvhd ~ cmv,
                         data = bmt_df[bmt_df$deltaa==1, ],
                         conf.type = "log-log")

survdiff(formula=s.dfsgvhd~cmv, data=bmt_df[bmt_df$deltaa==1, ])

bycmv_dir5_plot <- ggsurvplot(survfit.bycmv,
                              conf.int = F,
                              pval = TRUE,
                              pval.coord = c(1700, 1),
                              xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier disease-free survival est,",
       subtitle = "by patient CMV status in patients who develop aGVHD")

#by donorcmv
survfit.bydonorcmv <- survfit(s.dfsgvhd ~ donorcmv,
                              data = bmt_df[bmt_df$deltaa==1, ],
                              conf.type = "log-log")

survdiff(formula=s.dfsgvhd~donorcmv, data=bmt_df[bmt_df$deltaa==1, ])

bydonorcmv_dir5_plot <- ggsurvplot(survfit.bydonorcmv,
                                   conf.int = F,
                                   pval = TRUE,
                                   pval.coord = c(1700, 1),
                                   xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier GVHD survival est,",
       subtitle = "by donor CMV status in patients with aGVHD")

#by disease group
survfit.bydg <- survfit(s.dfsgvhd ~ disgroup,
                        data = bmt_df[bmt_df$deltaa==1, ],
                        conf.type = "log-log")

survdiff(formula=s.dfsgvhd~disgroup, data=bmt_df[bmt_df$deltaa==1, ])

bydis_dir5_plot <- ggsurvplot(survfit.bydg,
                              conf.int = F,
                              pval = TRUE,
                              pval.coord = c(1700, 1),
                              xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier GVHD survival est,",
       subtitle = "by hospital in patients who develop aGVHD")

#by FAB classification
survfit.byfab <- survfit(s.dfsgvhd ~ fab,
                         data = bmt_df[bmt_df$deltaa==1, ],
                         conf.type = "log-log")

survdiff(formula=s.dfsgvhd~fab, data=bmt_df[bmt_df$deltaa==1, ])

byfab_dir5_plot <- ggsurvplot(survfit.bydg,
                              conf.int = F,
                              pval = TRUE,
                              pval.coord = c(1700, 1),
                              xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier GVHD survival est,",
       subtitle = "by FAB classification in patients with aGVHD")

#by mtx
survfit.bymtx <- survfit(s.dfsgvhd ~ mtx,
                         data = bmt_df[bmt_df$deltaa==1, ],
                         conf.type = "log-log")

survdiff(formula=s.dfsgvhd~mtx, data=bmt_df[bmt_df$deltaa==1, ])

bymtx_dir5_plot <- ggsurvplot(survfit.bymtx,
                              conf.int = F,
                              pval = TRUE,
                              pval.coord = c(1700, 1),
                              xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier disease-free survival est,",
       subtitle = "by MTX status in patients who develop aGVHD")

#by hospital
survfit.byhospital <- survfit(s.dfsgvhd ~ hospital,
                              data = bmt_df[bmt_df$deltaa==1, ],
                              conf.type = "log-log")

survdiff(formula=s.dfsgvhd~hospital, data=bmt_df[bmt_df$deltaa==1, ])

byhospital_dir5_plot <- ggsurvplot(survfit.byhospital,
                                   conf.int = F,
                                   pval = TRUE,
                                   pval.coord = c(1700, 1),
                                   xlab = "Time (in days)") + 
  labs(title = "Kaplan-Meier disease-free survival est,",
       subtitle = "by hospital in patients who develop aGVHD")

byhospital_dir5_plot

library(corrplot)
library(RColorBrewer)
M <-cor(bmt_df)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

