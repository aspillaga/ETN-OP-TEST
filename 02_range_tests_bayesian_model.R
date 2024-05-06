################################################################################
##
## OP-TESTS
## Model detection range using Bayesian inference
##
## Author: 
## Eneko Aspillaga (aspillaga@imedea.uib-csic.es)
## Instituto Mediterraneo de Estudios Avanzados (IMEDEA, CSIC)
##
## Last update: November 3, 2023
##
################################################################################

# Load libraries
library(data.table)
library(R2jags)       # Bayesian analysis
library(mcmcplots)    # Plotting MCMC

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load detection data in 2h bins
detect_2h <- fread("./data/range_tests/detections_2h_intervals.csv")

# Logistic function. Same function adjusted in the model.
logFun <- function(L, inflexion, b, x) {
  a = -inflexion * b
  # https://en.wikipedia.org/wiki/Logistic_function
  return(L / (1 + exp(-(a + b*x))))
}


# 1. PREPARE DATA FOR THE ANALYSIS =============================================

# Create new variable with receiver and transmitter manufacturer combination
detect_2h[, manufacturer := paste(receiver_manufacturer, tag_manufacturer,
                                  sep = "_")]

# Add date as a separate variable
detect_2h[, date := as.Date(date_time)]

# Delete non-compatible manufacturer combinations (combinations without
# any detections)
comp_tab <- tapply(detect_2h$n, detect_2h$manufacturer, sum)
rm_man <- names(comp_tab[comp_tab == 0])
detect_2h <- detect_2h[!manufacturer %in% rm_man, ]

# Delete the first and last dates of each test (they have a lower number of
# observations)
rm_indx <- unlist(lapply(unique(detect_2h$test_id), function(t) {
  indx_t <- which(detect_2h$test_id == t)
  return(indx_t[detect_2h$date[indx_t] %in% range(detect_2h$date[indx_t])])
}))
detect_2h <- detect_2h[-rm_indx, ]

# Delete Lotek receiver data at Yser from 2021-05-22 on (receiver run out of 
# battery)
detect_2h <- detect_2h[!(date >= as.Date("2021-05-22") & 
                           detect_2h$test_id == "Yser-OPs"), ]


# # 2. FIT BAYESIAN MODEL ======================================================

# Skip this section if models have been already fitted and saved in
# "./data/range_tests/model_fits/". Models are loaded in step 3.

# # 2.1. Pack data for JAGS ----------------------------------------------------
# 
# # Extract unique test locations
# tests <- unique(detect_2h$test_location)
# 
# # Generate a list with the packed data for each location
# data_list <- lapply(tests, function(t) {
# 
#   # Subset data corresponding to each location
#   data_sub <- detect_2h[test_location == t]
# 
#   # Inititalize object to store all the data
#   data <- list()
# 
#   # Create method and trial factors.
#   # "method": method to be compared, that is, each manufacturer combination
#   # "trial": each unique date and method combination
#   data_sub$ID_method <- factor(data_sub$manufacturer)
#   data_sub$ID_trial <- factor(paste(data_sub$date, data_sub$ID_method))
# 
#   # Observation level (2h intervals)
#   data$n_obs <- nrow(data_sub)
#   data$detections_obs <- data_sub$n
#   data$emissions_obs <- data_sub$expected_n
#   data$ID_method_obs <- as.numeric(data_sub$ID_method)
#   data$ID_trial_obs <- as.numeric(data_sub$ID_trial)
#   data$x <- data_sub$distance
# 
#   # Method level (manufacturer combination)
#   data$labels_method <- levels(data_sub$ID_method)
#   data$n_method <- length(data$labels_method)
# 
#   # Trial level (days)
#   data$label_trial <- levels(data_sub$ID_trial)
#   data$n_trial <- length(data$label_trial)
#   data$ID_method_trial <- match(sub("(.*) (.*)", "\\2", data$label_trial),
#                                 data$labels_method)
#   data$label_day_trial <- sub("(.*) (.*)", "\\1", data$label_trial)
# 
#   return(data)
# })
# 
# names(data_list) <- tests
# 
# 
# # 2.2. Model definition ------------------------------------------------------
# 
# sink("./model.txt")
# cat("model {
# 
#   # Observation level
#   for (i in 1:n_obs) {
# 
#     # Binomial distribution
#     detections_obs[i] ~ dbin(p[i], emissions_obs[i])
# 
#     # Logistic model for probability
#     L_obs[i] <- L_t[ID_trial_obs[i]]
#     inflexion_obs[i] <- inflexion_t[ID_trial_obs[i]]
#     b_obs[i] <- b_t[ID_trial_obs[i]]
#     a[i] <- -inflexion_obs[i] * b_obs[i]
# 
#     p[i] <- L_obs[i] / (1 + exp(-(a[i] + b_obs[i] * x[i])))
#   }
# 
#   # Trial level, random factor: Separate days
#   for (t in 1:n_trial) {
#     L_t[t] ~ dnorm(L_m[ID_method_trial[t]], tol_L_m[ID_method_trial[t]])
#     inflexion_t[t] ~ dnorm(inflexion_m[ID_method_trial[t]],
#                            tol_inflexion_m[ID_method_trial[t]])
#     b_t[t] ~ dnorm(b_m[ID_method_trial[t]], tol_b_m[ID_method_trial[t]])
# 
#   }
# 
#   # Method level, fixed factor: combination of receiver and tag manufacturers
#   for (m in 1:n_method) {
#     L_m[m] ~ dunif(0, 1)
#     sd_L_m[m] ~ dgamma(1, 1)
#     tol_L_m[m] <- 1/sd_L_m[m]^2
# 
#     inflexion_m[m] ~ dunif(0, 400)
#     sd_inflexion_m[m] ~ dgamma(1, 0.01)
#     tol_inflexion_m[m] <- 1/sd_inflexion_m[m]^2
# 
#     b_m[m] ~ dnorm(0, 100)T(-1, 0) # Must be negative
#     sd_b_m[m] ~ dgamma(1, 1)
#     tol_b_m[m] <- 1/sd_b_m[m]^2
#   }
# 
# }", fill = TRUE)
# 
# sink()
# 
# 
# # 2.3. Loop to fit a separate model to each test -----------------------------
# 
# # Create directory for model outputs
# dir.create("./data/range_tests/model_fits/", showWarnings = FALSE)
# 
# # Loop to adjust the model to each test period
# set.seed(42)
# 
# for (i in seq_along(data_list)) {
# 
#   test <- names(data_list)[i]
#   data <- data_list[[i]]
# 
#   cat("Adjusting test", test, "\n")
# 
#   # Initial values
#   inits <- function() list()
# 
#   # Parameters to be monitored
#   params <- c("L_m", "inflexion_m", "b_m",
#               "sd_L_m", "sd_inflexion_m", "sd_b_m",
#               "L_t", "inflexion_t", "b_t")
# 
#   # For initial checks
#   # out <- jags(data, inits, params, "model.txt", n.chains = 3, n.thin = 1,
#   #             n.iter = 1000, n.burnin = 500)
#   # out <- update(out, n.iter = 10000, n.thin = 100)
# 
#   # Definitive run
#   time0 <- Sys.time()
#   out <- jags(data, inits, params, "model.txt", n.chains = 3, n.thin = 100,
#               n.iter = 30000, n.burnin = 20000)
#   Sys.time() - time0
#   cat("\n\n")
# 
#   # Export the model
#   saveRDS(list(data, out), 
#           file = paste0("./data/range_tests/model_fits/", test, ".rds"),
#           compress = "xz")
# }


# 3. LOAD FITTED MODELS ========================================================

# List model files
model_files <- list.files("./data/range_tests/model_fits/", pattern = ".rds", 
                          full.names = TRUE)

# Load models
out_list <- lapply(model_files, readRDS)
names(out_list) <- sub(".*//(.*).rds", replacement = "\\1", model_files)


# 4. PLOT MODEL CONVERGENCE CHECKS =============================================

# Generate directory to save the plots
dir.create("./plots/range_tests/model_convergence/", showWarnings = FALSE, 
           recursive = TRUE)

# Function to transform tolerance to sd
tol2sd <- function(tol) 1 / sqrt(tol)

# Function to create positive normal priors for "b_m"
rnormTrunc <- function(n, mean, sd) {
  result <- c()
  while(length(result) < n) {
    result <- c(result, rnorm(n - length(result), mean, sd))
    result <- result[result > 0]
  }
  return(result)
}

# Define prior parameters (sames as in the model)
n <- 10000
priors <- list(L_m = runif(n, 0, 1),
               sd_L_m = rgamma(n, 1, 1),
               inflexion_m = runif(n, 0, 400),
               sd_inflexion_m = rgamma(n, 1, 0.01),
               b_m = rnormTrunc(n, 0, 100),
               sd_b_m = rgamma(n, 1, 1))

# Function to plot prior vs posterior distributions
plotPost <- function(priors, out, var) {
  prior_d <- density(priors[[var]])
  post <- out$BUGSoutput$sims.list[[var]]
  post_d <- lapply(1:ncol(post), function(i) density(post[, i]))
  max_p <- max(unlist(lapply(post_d, function(x) max(x$y))))
  plot(prior_d, ylim = range(0, max(prior_d$y), max_p), main = var,
       xlab = "", lwd = 2)
  for (i in seq_along(post_d)) {
    lines(post_d[[i]], col = 2)
  }
}

# Plot model convergence outputs
for (i in seq_along(out_list)) {
  t <- names(out_list)[i]
  cat(t, "\n")
  
  pdf(paste0("./plots/range_tests/model_convergence/", t, ".pdf"), width = 8, 
      height = 8)
  out <- out_list[[i]][[2]]
  
  # Plot chain convergence graphs
  traplot(out, parms = c("L_m"))
  traplot(out, parms = c("inflexion_m"))
  traplot(out, parms = c("b_m"))
  traplot(out, parms = c("sd_L_m"))
  traplot(out, parms = c("sd_inflexion_m"))
  traplot(out, parms = c("sd_b_m"))

  # Plot prior vs posterior distributions
  par(mfrow = c(3, 2), mar = c(3.1, 4.1, 2.1, 1.1), oma = c(0, 0, 3, 0))
  plotPost(priors, out, "L_m")
  plotPost(priors, out, "sd_L_m")
  plotPost(priors, out, "inflexion_m")
  plotPost(priors, out, "sd_inflexion_m")
  plotPost(priors, out, "b_m")
  plotPost(priors, out, "sd_b_m")
  title(main = "Priors vs posterior values", outer = TRUE)
  dev.off()
}


# 5. PLOT ACOUSTIC RANGE RESULTS ===============================================

dir.create("./plots/range_tests/acoustic_range_models/", showWarnings = FALSE)

# Lists of receiver and transmitter manufacturers
rec_man <- c("Thelma", "Lotek", "Sonotronics",
             "Innovasea-MAP114", "Innovasea-MAP115")
tag_man <- c("Thelma-OPi", "Thelma-OPs", "Thelma-R64K", 
             "Lotek-OPi", "Lotek-OPs", 
             "Sonotronics-OPi", "Sonotronics-OPs",
             "Innovasea-A69-1602", "Innovasea-A69-9007")

# Maximum distance to represent in the plots
xmax <- 500

# Loop for each test
for (i in seq_along(out_list)) {
  cat(names(out_list)[i], "\n")
  
  # Extract data for the especific test
  data <- out_list[[i]][[1]]
  out <- out_list[[i]][[2]]
  methods <- data$labels_method
  
  # Extract receiver and transmitter manufacturers used in the test
  receivers <- rec_man[rec_man %in% gsub("(.*)_(.*)", "\\1", methods)]
  tags <- tag_man[tag_man %in%  gsub("(.*)_(.*)", "\\2", methods)]
  
  # Initialize plot device
  jpeg(paste0("./plots/range_tests/acoustic_range_models/", names(out_list)[i], 
              ".jpg"), units = "in", res = 200, pointsize = 14,
       width = 1.5 * length(receivers) + 2.8, height = 1.3 * length(tags) + 1.8)
  par(mfrow = c(length(tags), length(receivers)), mai = rep(0.05, 4), 
      omi = c(0.8, 1.8, 1, 1))
  
  for (t in tags) {
    for (r in receivers) {
      cat("\t", r, t, "\n")
      
      # Extract data corresponding to a receiver and manufacturer combination
      m <- which(methods == paste0(r, "_", t))
      indx <- which(data$ID_method_obs == m)
      data$detect_prob <- data$detections_obs / data$emissions_obs
      
      if (length(m) == 0) {
        # If there is no data, generate an empty plot
        plot(1, type = "n", axes = FALSE, ann = FALSE, ylim = c(0, 1), 
             xlim = c(0, xmax))
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
             col = "#ffe9c1", border = FALSE)
        points(xmax/2, 0.5, pch = 4, cex = 5,  col = "#eec477", lwd = 3)
        box()
        
      } else {
        # If there is data, plot the data points and the estimated logistic
        # model from posterior values
        
        # Extract posterior values
        L <- out$BUGSoutput$sims.list$L_m[, m]
        b <- out$BUGSoutput$sims.list$b_m[, m]
        inflexion <- out$BUGSoutput$sims.list$inflexion_m[, m]
        
        # Apply logistic model to each iteration
        dist <- -50:(xmax + 100)
        val <- sapply(seq_along(L), function(j) {
          curve <- logFun(L[j], inflexion[j], b[j], dist)
        })
        
        # Extract the 95% limits of all the iterations
        v_min95 <- apply(val, 1, quantile, 0.025)
        v_max95 <- apply(val, 1, quantile, 0.975)
        
        # Apply logistic model to the average estimated parameters
        val_median <- logFun(out$BUGSoutput$median$L_m[m], 
                             out$BUGSoutput$median$inflexion_m[m], 
                             out$BUGSoutput$median$b_m[m], dist)
        
        # Start with an empty plot        
        plot(1, xlim = c(0, xmax), ylim = c(0, 1), type = "n", ann = FALSE, 
             axes = FALSE)
        
        # Add observed detection probability points with some jitter
        points(jitter(data$detect_prob[indx], 1) ~ jitter(data$x[indx], 2), 
               pch = 21, cex = 1, bg = adjustcolor("#b8dff5", 0.4), 
               col = adjustcolor("#6aa3db", 0.4), lwd = 0.4)
        
        # Add 95% posterior distribution
        polygon(c(dist, rev(dist)), y = c(v_min95, rev(v_max95)), 
                col = adjustcolor("#dcf022", 0.7), border = "#80a901", 
                lwd = 0.6)
        
        # Add median value
        lines(dist, val_median, col = "#4d6601", lwd = 2.2)
        
        box()
        
      }
      
      # Add the name of the tag manufacturer in the right
      if (r == receivers[1]) {
        label <- sub("-A69","\nA69",  sub("-O","\nO",  sub("-R","\nR", t)))
        text(-50, 0.5, label, xpd = NA, pos = 2, font = 1, cex = 1.6)
      }
      
      # Plot Y axis at the left (detection probability)
      if (r == receivers[length(receivers)]) {
        axis(4, at = seq(0, 1, 0.5), lwd = 0, lwd.ticks = 1, las = 1, 
             cex.axis = 1.2)
      }
      
      # Add X axis at the bottom (distance)
      if (t == tags[length(tags)]) {
        axis(1, at = pretty(c(0, xmax), 3),lwd = 0, lwd.ticks = 1, 
             cex.axis = 1.2)
      }
    }
  }
  
  # Add receiver manufacturer at the top
  mtext(outer = TRUE, side = 3,  sub("-MAP","\nMAP", receivers), line = 0, 
        font = 1, at = ((1:length(receivers)*2) - 1) / (length(receivers)*2), 
        cex = 1)
  
  # Add axis titles
  title(xlab = "Distance (m)", outer = TRUE, line = 2.8, cex.lab = 1.4)
  mtext(text = "Detection probability", side = 4, outer = TRUE, cex.lab = 1.4,
        line = 3.8)
  mtext(outer = TRUE, side = 3,  "Receivers", line = 4, font = 2, cex = 1.2)
  mtext(outer = TRUE, side = 2,  "Transmitter protocols", line = 9, font = 2, 
        cex = 1.2)
  
  dev.off()
  
}


# 6. PLOT MODEL PARAMETERS =====================================================

# 6.1. Extract posteriors ------------------------------------------------------

post_list <- lapply(seq_along(out_list), function(i) {
  
  test <- names(out_list)[[i]]
  cat(test, "\n")
  
  data <- out_list[[i]][[1]]
  out <- out_list[[i]][[2]]

  methods <- data$labels_method
  
  post_tmp <- lapply(seq_along(methods), function(j) {
    
    rec_man <- gsub("(.*)_(.*)", methods[j], replacement = "\\1")
    tag_man <- gsub("(.*)_(.*)", methods[j], replacement = "\\2")
    
    # Extract posteriors
    L <- out$BUGSoutput$sims.list$L_m[, j]
    b <- out$BUGSoutput$sims.list$b_m[, j]
    inflexion <- out$BUGSoutput$sims.list$inflexion_m[, j]
    
    # Calculate maximum detection probability (P(x) at distance = 0 m)
    P_0 <-  L / (1 + exp(b * inflexion))
    
    # Calculate acoustic range (distance at which probability is 50% of P_0)
    x_r <- inflexion - (log(1+ 2*exp(b * inflexion)))/b
    
    # Remove outliers in poor convergence chains (x_r > 1000)
    x_r[x_r > 1000] <- NA
    
    post_j <- data.frame(test, method = methods[j], rec_man, tag_man,
                         L = L, b = b, inflexion = inflexion,
                         P_0 = P_0, x_r = x_r)
    return(post_j)
    
  })
  
  return(rbindlist(post_tmp))
  
})

post <- rbindlist(post_list)


# 6.2. Plot posterior values of each parameter ---------------------------------

# Test labels
test_lab <- c(CP = "Open sea (BE)", Yser = "River (BE)", 
              Formosa = "Coastal lagoon (PT)", 
              Palma = "Coastal habitat (SP)")

# Transmitter manufacturer and protocol labels
tag_lab <-  c("Thelma-OPi", "Thelma-OPs", "Thelma-R64K", 
              "Lotek-OPi", "Lotek-OPs", 
              "Sonotronics-OPi", "Sonotronics-OPs",
              "Innovasea-A69-1602", "Innovasea-A69-9007")
tag_lab_plot <- gsub("^(.*?)-", "\\1\n", tag_lab, perl = TRUE)
names(tag_lab_plot) <- tag_lab

# Receiver manufactuer labels and colors
rec_lab <- c("Thelma", "Lotek", "Sonotronics", "Innovasea-MAP114", 
             "Innovasea-MAP115")
rec_col <-  c("#77D0F7", "#9DED7C", "#F6E472", "#F7A46B", "#B69BD3")
rec_col_l <-  c("#4A6DB4", "#024F18", "#B2A31B", "#CC4B12", "#5D2E88")


# Parameters
param_list <- c("L", "b", "inflexion", "P_0", "x_r")
param_lab <- c("L (detection prob.)", "b (slope)", "Inflection point (m)", 
               expression("P"[0] ~ "(Detect. prob. at 0 m)"),
               expression("x"[r]~"(50% range, m)"))
param_yline <- c(2.6, 3.4, 3.2, 2.6, 3.2)

# Create directory for plots
dir <- "./plots/range_tests/model_posteriors/"
dir.create(dir, showWarnings = FALSE)

# Loop for each variable
for (p in param_list) {
  
  # Subset posterior values and calculate median and 95% range
  post[, var := data.frame(post)[, p]]
  post_sub <- post[, .(min = quantile(var, 0.025, na.rm = TRUE),
                       median = quantile(var, 0.5, na.rm = TRUE),
                       max = quantile(var, 0.95, na.rm = TRUE)),
                   by = list(test, method, rec_man, tag_man)]
  
  # Order data frame
  post_sub[, rec_man := factor(rec_man, levels = rec_lab)]
  post_sub[, tag_man := factor(tag_man, levels = tag_lab)]
  post_sub <- post_sub[order(tag_man, rec_man), ]
  
  jpeg(paste0(dir, p, ".jpeg"), width = 2100, height = 2500, res = 300, 
       pointsize = 12)
  
  # Define plot layout
  layout(matrix(c(rep(1:(length(test_lab) - 1), each = 2), 
                  length(test_lab) + 0:1), nrow = length(test_lab), ncol = 2, 
                byrow = TRUE), widths = c(6, 4))
  
  # Outer margins
  par(mar = c(3.8, 5.1, 2.6, 0.2), oma = c(1.1, 0, 0, 1.1))
  
  # Loof for each test
  for (t in seq_along(test_lab)) {
    
    # Subset data
    post_t <- post_sub[test == names(test_lab)[t]]
    post_t$method <- factor(post_t$method, levels = unique(post_t$method))
    
    # Set positions for data points
    w <- 1.5 # Separation between different protocols
    post_t[, at := rep(1, nrow(post_t))]
    indx <- c(FALSE, post_t$tag_man[-1] != post_t$tag_man[-nrow(post_t)])
    post_t$at[indx] <- w
    post_t[, at := cumsum(at)]
    
    # Limit for the Y axis
    ylim <- c(range(c(post_t$min, post_t$max)))
    if (p == "L" | p == "P_0") ylim <- c(0, 1)
    if (p == "inflexion") ylim <- c(0, 400)
    if (p == "x_r") ylim <- c(0, 850)
    
    # Empty plot
    plot(median ~ at, data = post_t, ylim = ylim, type = "n", xaxs = "i",
         xlim = c(min(post_t$at) - w/2, max(post_t$at) + w/2), axes = FALSE,
         ann = FALSE)
    
    # Add posterior 95% intervals
    arrows(x0 = post_t$at, x1 = post_t$at, y0 = post_t$min, y1 = post_t$max, 
           code = 3, angle = 90, length = 0.07, lwd = 1.2, 
           col = rec_col_l[post_t$rec_man])
    
    # Plot median values
    points(median ~ at, pch = 21, data = post_t, cex = 1.4,
           bg = rec_col[post_t$rec_man], col = rec_col_l[post_t$rec_man])
    
    # Y axis
    axis(2, at = pretty(ylim, 4), lwd = 0, lwd.ticks = 1, las = 1)
    title(ylab = param_lab[param_list == p], 
          line = param_yline[param_list == p])
    
    # Separation lines for tag manufacturers
    abline(v = post_t$at[indx] - w/2, col = "gray80")
    
    # X axis - Protocols
    lab_at <- tapply(post_t$at, post_t$tag_man, mean)
    prot <- gsub("^(.*?)-(.*)", "\\2", names(lab_at), perl = TRUE)
    for (i in 1:2) {
      axis(1, at = lab_at[1:2 == i], labels = prot[1:2 == i], line = -1.1, 
           lwd = 0, cex.axis = 0.8)
    }
    axis(1, at =  post_t$at[indx] - w/2, lwd = 0, lwd.ticks = 1, tcl = -1,
         labels = FALSE, col.ticks = "gray80")
    
    # X axis - Tag manufacturers
    manuf <- gsub("^(.*?)-(.*)", "\\1", names(lab_at), perl = TRUE)
    lab_at_m <- tapply(lab_at, manuf, mean, na.rm = T)
    par(xpd = TRUE)
    rect(xleft = par("usr")[1], xright = par("usr")[2], 
         ytop = par("usr")[3] - 1*diff(par("usr")[3:4])/10,
         ybottom = par("usr")[3] - 2.2* diff(par("usr")[3:4])/10, 
         col = "gray90")
    par(xpd = FALSE)
    axis(1, at = lab_at_m, labels = names(lab_at_m), line = 0, lwd = 0, 
         cex.axis = 0.9, font = 2)
    line_indx <- c(FALSE, substr(post_t$tag_man[-1], 1, 1) != 
                     substr(post_t$tag_man[-nrow(post_t)], 1, 1))
    axis(1, labels = FALSE, lwd = 0, lwd.ticks = 1, tcl = -2,
         at = c(par("usr")[c(1, 2)], post_t$at[line_indx] - w/2))
    box()
  
    title(xlab = c("Transmitter manufacturers / Protocols"), line = 2.5,
          font.lab = 2)
    
    # Panel title
    title(main = test_lab[t], font.main = 2, adj = 0, line = 0.5)
    mtext(side = 2, text = LETTERS[t], font = 2, padj = -5.2, line = 3, las = 1, 
          cex = 1.2)
      
  }
  
  # Add legend
  par(mar = c(2.8, 0.1, 2.6, 0.2))
  plot(1, type = "n", ann = FALSE, axes = FALSE)
  rec_sub <- levels(post_t$rec_man)
  legend("bottomleft", legend = rec_sub, cex = 1.1, pch = 21, pt.cex = 1.4,
         pt.bg = rec_col, col = rec_col_l, inset = c(0.2, 0.01), 
         box.col = "gray50", title = " Receiver manufacturers ", title.font = 2, 
         title.adj = 0.25)
  
  dev.off()
}


# 6.3. Median values and ranges for the manuscript -----------------------------

# Median values for each test and manufacturer combination
post_m <- post[, .(x_r_median = median(x_r, na.rm = TRUE),
                   P_0_median = median(P_0, na.rm = TRUE),
                   L_median = median(L, na.rm = TRUE),
                   inflexion_median = median(inflexion, na.rm = TRUE),
                   b_median = median(b, na.rm = TRUE)),
               by = list(test, method, rec_man, tag_man)]

# Medians and 95% range by test
x_r_summ <- post_m[, .(median = quantile(x_r_median, 0.5),
                       min = quantile(x_r_median, 0.025),
                       max = quantile(x_r_median, 0.975)),
                   by = list(test)]

P_0_summ <- post_m[, .(median = quantile(P_0_median, 0.5),
                       min = quantile(P_0_median, 0.025),
                       max = quantile(P_0_median, 0.975)),
                   by = list(test)]
P_0_summ
