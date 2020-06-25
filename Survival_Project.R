library(survival)
library(RcmdrPlugin.survival)

data <- read.table('http://socserv.mcmaster.ca/jfox/Books/Companion/data/Rossi.txt', header = T)


# The employment variables are time-dependent
data.time <- unfold(data, 'week', 'arrest', 11:62, 'employed', lag = 1)
attach(data.time)

fit <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + wexp + mar + paro + prio + factor(educ) + employed, data = data.time)
summary(fit)

plot(survfit(fit), ylim = c(0.7, 1), xlab = 'Weeks', ylab = 'Proportion not Rearrested')


# Find the "best" model

fit <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + wexp + mar + paro + prio + factor(educ) + employed, method = 'breslow', data = data)
summary(fit)

# Remove wexp from model

fit.1 <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + mar + paro + prio + factor(educ) + employed, method = 'breslow', data = data.time)
summary(fit.1)

# Remove paro from model

fit.2 <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + mar + prio + factor(educ) + employed, method = 'breslow', data = data.time)
summary(fit.2)

# Remove educ from model

fit.3 <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + mar + prio + employed, method = 'breslow', data = data.time)
summary(fit.3)

# Remove mar from model

fit.4 <- coxph(Surv(start, stop, arrest.time) ~ fin + age + race + prio, method = 'breslow', data = data.time)
summary(fit.4)

# Remove race from model

fit.5 <- coxph(Surv(start, stop, arrest.time) ~ fin + age + prio + employed, method = 'breslow', data = data.time)
summary(fit.5)

# Fin is slightly significant, while age ang prio are highly significant


# For dianostics, we create a new subset of the data

data.new <- subset(data.time, stop == 52 | arrest.time == 1)
attach(data.new)

# Check Linearity

fit1 <- coxph(Surv(week, arrest) ~ fin + prio + employed , method = 'breslow', data = data.new)

plot(age, resid(fit1), xlab = 'Age', ylab = 'Residual')
lines(lowess(age, resid(fit1), f = .5))
abline(lm(resid(fit1) ~ age), lty = 3)

# Check Age^2

plot(age^2, resid(fit1), xlab = 'Age', ylab = 'Residual')
lines(lowess(age^2, resid(fit1), f = .5))
abline(lm(resid(fit1) ~ age^2), lty = 3)

# Check sqrt(Age)

plot(sqrt(age), resid(fit1), xlab = 'Age', ylab = 'Residual')
lines(lowess(sqrt(age), resid(fit1), f = .5))
abline(lm(resid(fit1) ~ sqrt(age)), lty = 3)

# Age can stay as a linear term in the model

fit2 <- coxph(Surv(week, arrest) ~ fin + age + employed, method = 'breslow', data = data.new)

plot(prio, resid(fit2), xlab = 'Number of Prior Convictions', ylab = 'Residual')
lines(lowess(prio, resid(fit2), f = .5))
abline(lm(resid(fit2) ~ prio), lty = 3)

# Check Transformations

plot(prio^2, resid(fit2), xlab = 'Number of Prior Convictions', ylab = 'Residual')
lines(lowess(prio^2, resid(fit2), f = .5))
abline(lm(resid(fit2) ~ prio^2), lty = 3)

plot(log(prio+1), resid(fit2), xlab = 'Number of Prior Convictions', ylab = 'Residual')
lines(lowess(log(prio+1), resid(fit2), f = .5))
abline(lm(resid(fit2) ~ log(prio+1)), lty = 3)

plot(sqrt(prio), resid(fit2), xlab = 'Number of Prior Convictions', ylab = 'Residual')
lines(lowess(sqrt(prio), resid(fit2), f = .5))
abline(lm(resid(fit2) ~ sqrt(prio)), lty = 3)


# Prio can stay as a linear term in the model


# PH Assumptions

temp <- cox.zph(fit.5)
print(temp)

# Age does not follow proportional hazards

#Check interaction between time and age

fit.6 <- coxph(Surv(start, stop, arrest.time) ~fin + age + age:stop + prio + employed, data=data.time)
summary(fit.6)

temp2 <- cox.zph(fit.6)
print(temp2)

# Stratify Age

data.time$age.cat <- recode(data.time$age, " lo:19=1; 20:25=2; 26:30=3; 31:hi=4 ")
xtabs(~ age.cat, data=data.time)


fit.7 <- coxph(Surv(start, stop, arrest.time) ~fin + strata(age.cat) + prio + employed, data=data.time)
summary(fit.7)

temp3 <- cox.zph(fit.7)
print(temp3)

# Check AIC for best model

AIC(fit.6)
AIC(fit.7)

# Use stratified model

fit.final <- fit.7

# Outliers 

data.new$age.cat <- recode(data.new$age, " lo:19=1; 20:25=2; 26:30=3; 31:hi=4 ")
xtabs(~ age.cat, data=data.new)

fit.final.sub <- coxph(Surv(week, arrest) ~fin + strata(age.cat) + prio + employed, data=data.new)
summary(fit.final.sub)

fin.yes <- ifelse(data.new$fin == 'yes', 1, 0)
employed.yes <- ifelse(data.new$employed == 'yes', 1, 0)

rs <- -.29085 * fin.yes + .08686 * data.new$prio - .75241 * employed.yes
rs <- as.vector(rs)

plot(rs, resid(fit.final.sub, 'dev'), xlab = 'Risk Score', ylab = 'Deviance Residual')
abline(0,0)
identify(rs, resid(fit.final.sub, 'dev'))

# Influential

sresid <- resid(fit.final.sub, 'score')
dim(sresid)

par(mfrow = c(2, 2))
plot(1:431, sresid[, 1], xlab = 'Financial Aid', ylab = 'Residual')
identify(1:431, sresid[, 1], 1:431)

plot(1:431, sresid[, 2], xlab = 'Prior Arrests', ylab = 'Residual')
identify(1:431, sresid[, 2], 1:431)

plot(1:431, sresid[, 3], xlab = 'Employment Status', ylab = 'Residual')
identify(1:431, sresid[, 3], 1:431)


# Lets look at the estimated survival curves of fin

data.fin <- with(data.new, data.frame(fin = c(0, 1), age = rep(mean(age), 2), prio = rep(mean(prio), 2)))
plot(survfit(td.model.5, newdata = data.fin), conf.int = T, lty = c(1, 2), ylim = c(0.6, 1), xlab = 'Weeks', ylab = 'Proportion not Rearrested')
legend('bottomleft', legend = c('fin = no', 'fin = yes'), lty = c(1, 2), inset = .02)











