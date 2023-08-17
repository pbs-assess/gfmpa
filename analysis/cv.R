M <- 300
n <- 100

surv_dat <- readRDS("data-generated/dat_to_fit.rds")
dat <- filter(surv_dat, species_common_name == "pacific cod") |>
  filter(survey_abbrev == "SYN WCHG") |> mutate(density = density_kgpm2 * 10000)

dat |> group_by(year) |> summarise(var = var(density)) |> summarise(var = mean(var))

s2 <- 0.3

# eq. 4 Schnute and Haigh CJFAS
Vbar <- ((M - n) / M) / (s2 / n)
Vbar



library(dplyr)
N <- 20
CV <- 0.2
second_mean = 0.9
cv <- function(x) sd(x) / mean(x)
# # cv(x1)

test_diff <- function(CV = 0.5, N = 100, decline = 0.1) {
  SD <- sqrt(log(CV^2 + 1))
  x1 <- rlnorm(N, meanlog = log(1) - (SD^2)/2 , sdlog = SD)
  x2 <- rlnorm(N, meanlog = log(1 - decline) - (SD^2)/2, sdlog = SD)
  d <- data.frame(x = c(rep(1, N), rep(2, N)), y = c(x1, x2))
  m <- lm(log(y) ~ as.factor(x), data = d)
  est <- coef(m)[[2]]
  se <- summary(m)$coef[2,2]
  lwr <- est - qnorm(0.975) * se
  upr <- est + qnorm(0.975) * se
  data.frame(lwr = lwr, est = est, upr = upr)
}


# 0.25
# 0.32
REPS <- 1000
set.seed(1)
out1 <- purrr::map_dfr(seq_len(REPS), function(i) test_diff(CV = 0.25, N = 100, decline = 0.1))
out1$cv <- 0.25
out2 <- purrr::map_dfr(seq_len(REPS), function(i) test_diff(CV = 0.32, N = 100, decline = 0.1))
out2$cv <- 0.32
bind_rows(out1, out2) |> mutate(sig = upr < 0) |>
  group_by(cv) |> summarise(m = mean(sig))


set.seed(1)
out1 <- purrr::map_dfr(seq_len(REPS), function(i) test_diff(CV = 0.25 + 0.2, N = 100, decline = 0.2))
out1$cv <- 0.25
out2 <- purrr::map_dfr(seq_len(REPS), function(i) test_diff(CV = 0.32 + 0.2, N = 100, decline = 0.2))
out2$cv <- 0.32
bind_rows(out1, out2) |> mutate(sig = upr < 0) |>
  group_by(cv) |> summarise(m = mean(sig))


a <- 100
b <- -0.08
CV <- 0.4

test_slope <- function(CV = 0.4, N = 100, change_per_decade = -0.7, a = 100, proc_sd = 0.2) {
  yrs <- seq(0, 9)
  SD <- sqrt(log(CV^2 + 1))
  decades <- yrs / 10
  eta <- exp(log(a) + change_per_decade * decades)
  eta_true <- eta
  # add 0.2 process error:
  eta <- rlnorm(length(eta), log(eta) - (SD^2) / 2, proc_sd)
  # plot(yrs, eta)
  # eta
  dat <- purrr::map_dfr(yrs, function(y)
    data.frame(year = y + 1, obs = rlnorm(N, log(eta[y + 1]) - (SD^2) / 2, SD)))
  dat$decade <- dat$year / 10
  m <- lm(log(obs) ~ decade, data = dat)
  est <- coef(m)[[2]]
  se <- summary(m)$coef[2,2]
  lwr <- est - qnorm(0.975) * se
  upr <- est + qnorm(0.975) * se
  data.frame(lwr = lwr, est = est, upr = upr, fraction_remain = eta_true[length(eta_true)] / eta_true[1], true_change = change_per_decade, proc_sd = proc_sd, cv = CV)
}

# decline = 1 - exp(-coef)
# decline - 1 = - exp(-coef)
# 1 - decline = exp(-coef)
# -log(1 - decline) = coef
decl <- log(1 - 0.3)
decl
exp(decl)

REPS <- 500
set.seed(1)
out1 <- purrr::map_dfr(seq_len(REPS), function(i) test_slope(CV = 0.2, N = 100,
  change_per_decade = decl, proc_sd = 0.2))
set.seed(1)
out2 <- purrr::map_dfr(seq_len(REPS), function(i) test_slope(CV = 0.4, N = 100,
  change_per_decade = decl, proc_sd = 0.2))

bind_rows(out1, out2) |>
  mutate(detect_decline = upr < 0, sign_pos = lwr > 0) |>
  group_by(cv) |>
  summarise(
    power = mean(detect_decline),
    m_error = mean(est / true_change),
    s_error = mean(sign_pos),
    fraction_remain = mean(fraction_remain),
    proc_sd = mean(proc_sd)
  )

