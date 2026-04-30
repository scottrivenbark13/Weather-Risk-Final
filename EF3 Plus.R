

library(tidyverse)
library(sf)
library(maps)
library(extRemes)
library(ggplot2)
library(viridis)
library(patchwork)
library(lubridate)
library(MASS)         


# Downloading Archive Data

spc_url    <- "https://www.spc.noaa.gov/wcm/data/1950-2023_actual_tornadoes.csv"
local_file <- "spc_tornadoes.csv"

if (!file.exists(local_file)) {
  message("Downloading SPC tornado database...")
  download.file(spc_url, destfile = local_file, mode = "wb")
} else {
  message("Using cached SPC tornado database.")
}

raw <- read_csv(local_file, show_col_types = FALSE)
message("Columns: ", paste(names(raw), collapse = ", "))


# Data orgnaization, filtering for the fields we need for the analysis

tornadoes <- raw %>%
  rename(
    year   = yr,
    month  = mo,
    day    = dy,
    state  = st,
    rating = mag,
    slat   = slat,
    slon   = slon,
    elat   = elat,
    elon   = elon,
    length = len,
    width  = wid
  ) %>%
  filter(
    !is.na(slat), !is.na(slon),
    slat > 20, slat < 55,
    slon > -130, slon < -65,
    rating >= 0
  ) %>%
  mutate(date = make_date(year, month, day))

message(sprintf("Total tornadoes in database: %d", nrow(tornadoes)))


# Split times periods by EF and F ratings

ef_cutover <- as.Date("2007-02-01")

tornadoes <- tornadoes %>%
  mutate(
    scale     = if_else(date < ef_cutover, "F", "EF"),
    ef_rating = rating
  )

message("\n--- Tornado counts by rating scale era ---")
tornadoes %>%
  group_by(scale) %>%
  summarise(n = n(), min_rating = min(ef_rating), max_rating = max(ef_rating)) %>%
  print()


# Finding Stong EF3+ and F3+ torndaoes

violent <- tornadoes %>%
  filter(ef_rating >= 3)

message(sprintf("\nEF3+ tornadoes in dataset: %d", nrow(violent)))


# Define Peridos of study, Early and Recent Period

period1_range <- c(1954, 1988)
period2_range <- c(1989, 2023)

violent <- violent %>%
  mutate(
    period = case_when(
      year >= period1_range[1] & year <= period1_range[2] ~ "1954-1988",
      year >= period2_range[1] & year <= period2_range[2] ~ "1989-2023",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period))

message("\n--- EF3+ count by period ---")
violent %>% count(period) %>% print()


# Spatial grid setup

# use 2-degree cells (~220 km) for the Poisson rate ratio map.
# coarser 2-degree cells have enough counts to make the Poisson test stable.

grid_res <- 2.0

violent <- violent %>%
  mutate(
    lat_bin = floor(slat / grid_res) * grid_res + grid_res / 2,
    lon_bin = floor(slon / grid_res) * grid_res + grid_res / 2
  )

grid_counts <- violent %>%
  group_by(period, lat_bin, lon_bin) %>%
  summarise(n_tornadoes = n(), .groups = "drop")

grid_wide <- grid_counts %>%
  pivot_wider(names_from = period, values_from = n_tornadoes, values_fill = 0) %>%
  rename(n_early = `1954-1988`, n_recent = `1989-2023`) %>%
  mutate(diff = n_recent - n_early)


# US map setup for multiple maps

us_states <- map_data("state") %>%
  filter(!(region %in% c("alaska", "hawaii")))

map_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    axis.text       = element_blank(),
    axis.title      = element_blank(),
    axis.ticks      = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", size = 13)
  )

xlim <- c(-105, -75)
ylim <- c(25, 45)

# KDE Smoothing: run 2D kernel density estimation on a subset of tornado points
# and return a tidy data frame of smoothed density values clipped to the map bounds.
kde_for_period <- function(df, period_label, n_grid = 150) {
  pts <- df %>% filter(period == period_label)
  if (nrow(pts) < 5) return(tibble())
  k <- kde2d(pts$slon, pts$slat,
             n = n_grid,
             lims = c(xlim[1], xlim[2], ylim[1], ylim[2]))
  expand.grid(lon = k$x, lat = k$y) %>%
    mutate(density = as.vector(k$z),
           period  = period_label)
}

kde_early  <- kde_for_period(violent, "1954-1988")
kde_recent <- kde_for_period(violent, "1989-2023")

# Normalized each period's density to [0, 1] so the color scales are comparable
kde_early  <- kde_early  %>% mutate(density_norm = density / max(density))
kde_recent <- kde_recent %>% mutate(density_norm = density / max(density))

# dfference surface: recent density minus early density (both on [0,1] scale)
# grids are identical so they can be can subtract directly.
# use dplyr::select() explicitly
kde_diff <- kde_early %>%
  dplyr::select(lon, lat) %>%
  mutate(
    d_early  = kde_early$density_norm,
    d_recent = kde_recent$density_norm,
    diff     = d_recent - d_early
  )


# EF3+ Early Period

p_early <- ggplot() +
  geom_raster(data = kde_early,
              aes(x = lon, y = lat, fill = density_norm),
              interpolate = TRUE) +
  scale_fill_gradientn(
    colors   = c("white", "#fff7bc", "#fc8d59", "#d73027"),
    name     = "Relative\nDensity",
    limits   = c(0, 1)
  ) +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "grey30", linewidth = 0.35) +
  coord_fixed(ratio = 1.3, xlim = xlim, ylim = ylim) +
  labs(title = "EF3+ Tornado Density - Early Period (1954-1988)") +
  map_theme

print(p_early)


#EF3+ Recent Period

p_recent <- ggplot() +
  geom_raster(data = kde_recent,
              aes(x = lon, y = lat, fill = density_norm),
              interpolate = TRUE) +
  scale_fill_gradientn(
    colors   = c("white", "#fff7bc", "#fc8d59", "#d73027"),
    name     = "Relative\nDensity",
    limits   = c(0, 1)
  ) +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "grey30", linewidth = 0.35) +
  coord_fixed(ratio = 1.3, xlim = xlim, ylim = ylim) +
  labs(title = "EF3+ Tornado Density - Recent Period (1989-2023)") +
  map_theme

print(p_recent)

#smoothed difference map

p_diff <- ggplot() +
  geom_raster(data = kde_diff,
              aes(x = lon, y = lat, fill = diff),
              interpolate = TRUE) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#d73027",
    midpoint = 0,
    name     = "Recent - Early\n(normalized)"
  ) +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "grey30", linewidth = 0.35) +
  coord_fixed(ratio = 1.3, xlim = xlim, ylim = ylim) +
  labs(title    = "Change in EF3+ Tornado Density (Recent - Early)",
       subtitle = "Blue = fewer recently  |  Red = more recently") +
  map_theme

print(p_diff)

# poisson rate ratio figure

n_years_early  <- diff(period1_range) + 1
n_years_recent <- diff(period2_range) + 1

poisson_summary <- grid_wide %>%
  mutate(
    lambda_early  = n_early  / n_years_early,
    lambda_recent = n_recent / n_years_recent,
    rate_ratio    = (n_recent + 0.5) / (n_early + 0.5),
    pvalue = mapply(function(x, y) {
      tryCatch(
        poisson.test(c(x, y), c(n_years_early, n_years_recent))$p.value,
        error = function(e) NA_real_
      )
    }, n_early, n_recent),
    significant = pvalue < 0.05
  )

sig_cells <- sum(poisson_summary$significant, na.rm = TRUE)
message(sprintf("\nGrid cells with significant change in EF3+ rate: %d", sig_cells))

active_poisson <- poisson_summary %>%
  filter(n_early + n_recent > 0)

# computed national rate ratio across all active cells combined
total_early_all  <- sum(active_poisson$n_early,  na.rm = TRUE)
total_recent_all <- sum(active_poisson$n_recent, na.rm = TRUE)
national_rr      <- (total_recent_all / n_years_recent) /
  (total_early_all  / n_years_early)

message(sprintf("National rate ratio (recent/early): %.3f  (log2 = %.3f)",
                national_rr, log2(national_rr)))

active_poisson <- active_poisson %>%
  mutate(log2_std_rr = log2(rate_ratio / national_rr))

# stippled points for statistically significant cells
sig_points <- active_poisson %>% filter(significant)

p_rate_ratio <- ggplot() +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey92", color = "white", linewidth = 0.3) +
  geom_tile(data = active_poisson,
            aes(x = lon_bin, y = lat_bin, fill = log2_std_rr),
            alpha = 0.9) +
  scale_fill_gradient2(low      = "#2166ac",
                       mid      = "white",
                       high     = "#d73027",
                       midpoint = 0,
                       name     = "log2(Standardized\nRate Ratio)",
                       na.value = "transparent") +
  # stippled dots mark cells with statistically significant absolute rate change
  geom_point(data = sig_points,
             aes(x = lon_bin, y = lat_bin),
             shape = 3, size = 1.5, color = "black", alpha = 0.7) +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "grey30", linewidth = 0.35) +
  coord_fixed(ratio = 1.3, xlim = xlim, ylim = ylim) +
  labs(title    = "Poisson Rate Ratio - EF3+ Tornadoes (Relative to National Trend)",
       subtitle = paste0("Red = relative gain  |  Blue = relative loss  |  ",
                         "+ marks p < 0.05 absolute change")) +
  map_theme

print(p_rate_ratio)

# Poisson Rate and Return Period Model

target_states <- c("TX", "OK", "KS", "MO", "AR", "LA", "MS", "AL",
                   "TN", "KY", "GA", "IL", "IN", "NE", "IA")

n_years_total <- 2023 - 1954 + 1   # 70 years of record

ef_thresholds <- c(3, 4, 5)

poisson_rp <- map_dfr(ef_thresholds, function(ef_thr) {
  n_events <- tornadoes %>%
    filter(state %in% target_states,
           year  >= 1954,
           ef_rating >= ef_thr) %>%
    nrow()
  
  lambda_hat <- n_events / n_years_total
  
  # Exact Poisson CI on the rate (not the count)
  pci <- poisson.test(n_events, T = n_years_total, conf.level = 0.95)
  lambda_low  <- pci$conf.int[1]
  lambda_high <- pci$conf.int[2]
  
  # For EF3+ (lambda >> 1), report as events/year instead
  tibble(
    ef_threshold      = ef_thr,
    n_events          = n_events,
    lambda_per_yr     = round(lambda_hat,  3),
    lambda_ci_low     = round(lambda_low,  3),
    lambda_ci_high    = round(lambda_high, 3),
    # Return period: meaningful only when < 1 event/yr; NA otherwise
    return_period_yr  = if_else(lambda_hat < 1,
                                round(1 / lambda_hat, 1), NA_real_),
    rp_ci_low         = if_else(lambda_hat < 1,
                                round(1 / lambda_high, 1), NA_real_),
    rp_ci_high        = if_else(lambda_hat < 1,
                                round(1 / lambda_low,  1), NA_real_),
    # Exceedance probability: correct Poisson formula for ALL lambda values
    exceed_prob_pct   = round((1 - exp(-lambda_hat)) * 100, 1),
    ep_ci_low_pct     = round((1 - exp(-lambda_high)) * 100, 1),
    ep_ci_high_pct    = round((1 - exp(-lambda_low))  * 100, 1)
  )
})

message("\n--- Poisson Rate and Return Period by EF Threshold (1954-2023) ---")
message("lambda = avg events/year | return_period = avg years between events (EF4+/EF5+ only)")
message("exceed_prob = P(at least 1 event in a given year) using 1 - exp(-lambda)")
print(poisson_rp)

# EAD proxy: expected annual count weighted by EF severity
lambda_vec <- poisson_rp$lambda_per_yr
ead_approx <- sum(diff(c(0, rev(lambda_vec))) * rev(ef_thresholds))
message(sprintf("\nEAD Proxy (Poisson intensity-weighted): %.3f EF-tornado units/year", ead_approx))
message("(avg annual sum of EF ratings for violent tornadoes in the study region)")


# =============================================================================
# SECTION 13: PLOT 5 -- Two-panel: rate chart (EF3+) and return period (EF4+/5+)
# =============================================================================
# EF3+ has lambda >> 1 so we plot events/year with CI.
# EF4+ and EF5+ have lambda < 1 so return period (years between events) is shown.

# Panel A: Annual rate for all three thresholds
p_lambda <- ggplot(poisson_rp,
                   aes(x = factor(ef_threshold,
                                  labels = c("EF3", "EF4", "EF5")))) +
  geom_col(aes(y = lambda_per_yr),
           fill = "#d73027", alpha = 0.75, width = 0.5) +
  geom_errorbar(aes(ymin = lambda_ci_low, ymax = lambda_ci_high),
                width = 0.15, linewidth = 0.9, color = "grey20") +
  geom_text(aes(y = lambda_per_yr,
                label = sprintf("%.2f/yr\n(%.1f%%)", lambda_per_yr, exceed_prob_pct)),
            vjust = -0.4, size = 3.3, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Annual Poisson Rate by EF Threshold",
    subtitle = "Southern US + Great Plains, 1954-2023  |  Error bars = 95% exact Poisson CI",
    x        = "Minimum EF Rating",
    y        = "Average Events per Year (lambda)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# Panel B: Return period for EF4+ and EF5+ only (where T = 1/lambda is meaningful)
rp_data <- poisson_rp %>% filter(!is.na(return_period_yr))

p_rp <- ggplot(rp_data,
               aes(x = factor(ef_threshold,
                              labels = paste0("EF5")))) +
  geom_col(aes(y = return_period_yr),
           fill = "#fc8d59", alpha = 0.75, width = 0.4) +
  geom_errorbar(aes(ymin = rp_ci_low, ymax = rp_ci_high),
                width = 0.12, linewidth = 0.9, color = "grey20") +
  geom_text(aes(y = return_period_yr,
                label = sprintf("%.1f yr\n(%.1f%%/yr)", return_period_yr, exceed_prob_pct)),
            vjust = -0.4, size = 3.3, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Return Period for Rare EF Thresholds",
    subtitle = "T = 1/lambda  (valid only when lambda < 1 event/yr)",
    x        = "Minimum EF Rating",
    y        = "Return Period (years between events)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_lambda)
print(p_rp)

message("\n--- Annual Exceedance Probabilities (with 95% CI) ---")
poisson_rp %>%
  dplyr::select(ef_threshold, n_events, lambda_per_yr,
                lambda_ci_low, lambda_ci_high,
                return_period_yr, rp_ci_low, rp_ci_high,
                exceed_prob_pct, ep_ci_low_pct, ep_ci_high_pct) %>%
  print()



# annual EF3+ Count Time Series

annual_ef3plus <- tornadoes %>%
  filter(state %in% target_states, ef_rating >= 3, year >= 1954) %>%
  count(year, name = "n_ef3plus")

p_trend <- ggplot(annual_ef3plus, aes(x = year, y = n_ef3plus)) +
  geom_col(fill = "#fc8d59", alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.4,
              color = "#d73027", linewidth = 1.2, se = TRUE) +
  geom_vline(xintercept = 2007, linetype = "dashed", color = "grey40") +
  annotate("text", x = 2008, y = max(annual_ef3plus$n_ef3plus) * 0.95,
           label = "EF-scale\nadopted", hjust = 0, size = 3.5, color = "grey30") +
  labs(
    title    = "Annual EF3+ Tornado Count - Southern US & Great Plains",
    subtitle = "Trend curve with 95% confidence band",
    x        = "Year",
    y        = "Number of EF3+ Tornadoes"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_trend)


# Summary Statistics Table (printed to Console)


summary_stats <- violent %>%
  filter(state %in% target_states) %>%
  group_by(period) %>%
  summarise(
    total_ef3plus      = n(),
    annual_avg         = round(n() / 35, 1),
    pct_ef3            = round(mean(ef_rating == 3) * 100, 1),
    pct_ef4            = round(mean(ef_rating == 4) * 100, 1),
    pct_ef5            = round(mean(ef_rating == 5) * 100, 1),
    mean_path_len_mi   = round(mean(length, na.rm = TRUE), 2),
    median_path_len_mi = round(median(length, na.rm = TRUE), 2),
    .groups = "drop"
  )

message("\n--- Summary Statistics: EF3+ Tornadoes in Target States ---")
print(summary_stats)