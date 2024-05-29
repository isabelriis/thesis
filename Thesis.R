# Additions to the extremogram package written by Cribben. 

# Functions right and left tail extremograms

Extremogram1_Extended <- function(x, quant, maxlag, type, plotting = 0, cutoff = 0.3, start = 0) {
  level = quantile(x, prob = quant)
  n = length(x)
  rhohat = rep(0, maxlag)
  
  if (type == 1) { # pos-pos
    rhohat[1] = 1
    for (i in 1:(maxlag - 1)) {
      rhohat[i + 1] = length((1:(n - i))[x[1:(n - i)] > level & x[(i + 1):n] > level])
      rhohat[i + 1] = rhohat[i + 1]/length((1:(n - i))[x[1:(n - i)] > level])
    }
  } else if (type == 2) { # neg-neg
    rhohat[1] = 1
    for (i in 1:(maxlag - 1)) {
      rhohat[i + 1] = length((1:(n - i))[x[1:(n - i)] < level & x[(i + 1):n] < level])
      rhohat[i + 1] = rhohat[i + 1]/length((1:(n - i))[x[1:(n - i)] < level])
    }
  } else if (type == 3) { # pos-neg
    rhohat[1] = 0
    for (i in 1:(maxlag - 1)) {
      rhohat[i + 1] = length((1:(n - i))[x[1:(n - i)] > level & x[(i + 1):n] < level])
      rhohat[i + 1] = rhohat[i + 1]/length((1:(n - i))[x[1:(n - i)] > level])
    }
  } else if (type == 4) { # neg-pos
    rhohat[1] = 0
    for (i in 1:(maxlag - 1)) {
      rhohat[i + 1] = length((1:(n - i))[x[1:(n - i)] < level & x[(i + 1):n] > level])
      rhohat[i + 1] = rhohat[i + 1]/length((1:(n - i))[x[1:(n - i)] < level])
    }
  }
  
  if (plotting == 1) {
    plot((start:(maxlag - 1)), rhohat[(start + 1):maxlag], type = "n", xlab = "lag", ylab = "extremogram", ylim = c(0, cutoff))
    lines((start:(maxlag - 1)), rhohat[(start + 1):maxlag], col = 1, lwd = 1, type = "h")
    abline((0:(maxlag - 1)), 0, col = 1, lwd = 1)
  }
  
  return(rhohat)
}

bootconf1_extended <- function (x, R, l, maxlag, quant, type, par=0, start = 1, cutoff = 0.3, alpha = 0.05){
  if (par == 1) {
    n = parallel::detectCores()
    boot = boot::tsboot(x, Extremogram1_Extended, R, l = l, sim = "geom", 
                        endcorr = TRUE, maxlag = maxlag, quant = quant, type = type, 
                        parallel = "snow", ncpus = n)
  } else {
    boot = boot::tsboot(x, Extremogram1_Extended, R, l = l, sim = "geom", 
                        endcorr = TRUE, maxlag = maxlag, quant = quant, type = type)
  }
  
  tmp = boot[[2]]
  mat = tmp[, (2:maxlag)]
  k = dim(mat)[2]
  pocket = matrix(0, ncol = 3, nrow = k)
  for (i in 1:k) {
    pocket[i, 1] = quantile(mat[, i], prob = (alpha/2))
    pocket[i, 2] = mean(mat[, i])
    pocket[i, 3] = quantile(mat[, i], prob = (1 - alpha/2))
  }
  
  results_df <- data.frame(
    lag = start:(k - 1 + start),
    lower_confidence = pocket[start:(k - 1 + start), 1],
    mean = pocket[start:(k - 1 + start), 2],
    upper_confidence = pocket[start:(k - 1 + start), 3]
  )
  
  return(results_df)
}

permfn1_extended <- function (x, p, m, type, exttype, maxlag, start = 1, alpha = 0.05) {
  results <- list() 
  
  if (type == 1) {
    
    cc = matrix(0, ncol = m, nrow = maxlag)
    for (i in 1:m) {
      pBACC = sample(x)
      cc[, i] = Extremogram1_Extended(pBACC, p, maxlag = maxlag, type = exttype)
    }
    results$type1 = cc
  }
  
  if (type == 2) {
    cc = matrix(0, ncol = m, nrow = maxlag)
    for (i in 1:m) {
      pBACC = sample(x)
      cc[, i] = Extremogram1_Extended(pBACC, p, maxlag = maxlag, type = exttype)
    }
    k = dim(cc)[1]
    pocket = matrix(0, ncol = 3, nrow = k)
    for (i in 1:k) {
      pocket[i, 1] = quantile(cc[i, ], prob = (alpha/2))
      pocket[i, 2] = mean(cc[i, ]) 
      pocket[i, 3] = quantile(cc[i, ], prob = (1 - alpha/2))
    }
    results$type2 = pocket 
  }
  
  if (type == 3) {
    
    cc = matrix(0, ncol = m, nrow = 2) 
    for (i in 1:m) {
      pBACC = sample(x)
      
      cc[, i] = Extremogram1_Extended(pBACC, p, maxlag = 2, type = exttype)
    }
    
    dde_value = quantile(cc[2, ], prob = alpha/2)
    gge_value = quantile(cc[2, ], prob = 1 - alpha/2)
    
    
    results$type3 = data.frame(
      lag = 1:(maxlag - 1),
      lower = rep(dde_value, maxlag - 1),
      upper = rep(gge_value, maxlag - 1)
    )
  }
  
  return(results)
}

plot_extremogram1 <- function( x, quant, maxlag, type, typeper, R = 10, l = 30, exttype = 1, m=99, cutoff, title, color) {
  
  ext <- Extremogram1_Extended(x = x, quant=quant, maxlag = maxlag, type = type, plotting = 0)
  per <-permfn1_extended(x=x, p = quant, m = m, type = typeper, exttype = exttype, maxlag, 1, 0.05)
  boot <- bootconf1_extended(x = x, R = R, l = l, maxlag = maxlag, quant = quant, type = type, par = 0)
  
  data_df <- data.frame(
    extremogram = ext[-1],
    lag = boot$lag,
    per_lower = per$type3$lower,
    per_upper = per$type3$upper,
    boot_upper = boot$upper_confidence,
    boot_lower = boot$lower_confidence, 
    boot_mean = boot$mean
  )
  
  plot <- ggplot(data_df, aes(x = lag)) +
    geom_line(aes(y = boot_mean), color = color, size=0.75) +
    geom_line(aes(y = extremogram), color = "black") +
    geom_segment(aes(xend = lag, y = 0, yend = extremogram), color = "grey60") +
    geom_line(aes(y = per_lower), lty=2) +
    geom_line(aes(y = per_upper), lty=2) +
    geom_line(aes(y = boot_lower), lty=2, alpha=0.8, color = color, size=0.75) +
    geom_line(aes(y = boot_upper), lty=2, alpha=0.8, color = color, size=0.75) +
    ylim(c(0,cutoff))+
    labs(title = title, x = "Lag", y = "Extremogram") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18, hjust=0.75))
  
  return(plot)
  
}


# Functions for Cross-extremograms 

extremogram2_extended <- function (a, quant1, quant2, maxlag, type, ploting = 0, cutoff = 1, start = 0) {
  x = a[, 1]
  y = a[, 2]
  level1 = quantile(a[, 1], prob = quant1)
  level2 = quantile(a[, 2], prob = quant2)
  n = length(a[, 1])
  rhohat = rep(0, maxlag)
  if (type == 1) {
    for (i in 1:maxlag) {
      rhohat[i] = length((1:(n - i))[x[1:(n - i + 1)] > 
                                       level1 & y[i:n] > level2])
      rhohat[i] = rhohat[i]/length((1:(n - i))[x[1:(n - 
                                                      i + 1)] > level1])
    }
  }
  else if (type == 2) {
    for (i in 1:maxlag) {
      rhohat[i] = length((1:(n - i))[x[1:(n - i + 1)] < 
                                       level1 & y[i:n] < level2])
      rhohat[i] = rhohat[i]/length((1:(n - i))[x[1:(n - 
                                                      i + 1)] < level1])
    }
  }
  else if (type == 3) {
    for (i in 1:maxlag) {
      rhohat[i] = length((1:(n - i))[x[1:(n - i + 1)] < 
                                       level1 & y[i:n] > level2])
      rhohat[i] = rhohat[i]/length((1:(n - i))[x[1:(n - 
                                                      i + 1)] < level1])
    }
  }
  else if (type == 4) {
    for (i in 1:maxlag) {
      rhohat[i] = length((1:(n - i))[x[1:(n - i + 1)] > 
                                       level1 & y[i:n] < level2])
      rhohat[i] = rhohat[i]/length((1:(n - i))[x[1:(n - 
                                                      i + 1)] > level1])
    }
  }
  if (ploting == 1) {
    plot((start:(maxlag - 1)), rhohat[(start + 1):maxlag], 
         type = "n", xlab = "lag", ylab = "extremogram", ylim = c(0, 
                                                                  cutoff))
    lines((start:(maxlag - 1)), rhohat[(start + 1):maxlag], 
          col = 1, lwd = 1, type = "h")
    abline((0:(maxlag - 1)), 0, col = 1, lwd = 1)
  }
  return(rhohat)
}

bootconf2_extended <- function (x, R, l, maxlag, quant1, quant2, type, par, start = 1, cutoff = 1, alpha = 0.05) {
  if (par == 1) {
    n = parallel::detectCores()
    boot = boot::tsboot(x, extremogram2_extended, R, l = l, sim = "geom", 
                        endcorr = TRUE, maxlag = maxlag, quant1 = quant1, 
                        quant2 = quant2, type = type, parallel = "snow", 
                        ncpus = n)
    tmp = boot[[2]]
    mat = tmp[, (2:maxlag)]
    k = dim(mat)[2]
    pocket = matrix(0, ncol = 3, nrow = k)
    for (i in 1:k) {
      pocket[i, 1] = quantile(mat[, i], prob = (alpha/2))
      pocket[i, 2] = mean(mat[, i])
      pocket[i, 3] = quantile(mat[, i], prob = (1 - alpha/2))
    }
  }
  else {
    boot = boot::tsboot(x, extremogram2_extended, R, l = l, sim = "geom", 
                        endcorr = TRUE, maxlag = maxlag, quant1 = quant1, 
                        quant2 = quant2, type = type)
    tmp = boot[[2]]
    mat = tmp[, (2:maxlag)]
    k = dim(mat)[2]
    pocket = matrix(0, ncol = 3, nrow = k)
    for (i in 1:k) {
      pocket[i, 1] = quantile(mat[, i], prob = (alpha/2))
      pocket[i, 2] = mean(mat[, i])
      pocket[i, 3] = quantile(mat[, i], prob = (1 - alpha/2))
    }
  }
  results_df <- data.frame(
    lag = start:(k - 1 + start),
    lower_confidence = pocket[start:(k - 1 + start), 1],
    mean = pocket[start:(k - 1 + start), 2],
    upper_confidence = pocket[start:(k - 1 + start), 3]
  )
  
  return(results_df)
}

permatrix <- function (x) {
  x = as.matrix(x)
  nrow1 = dim(x)[1]
  ncol1 = dim(x)[2]
  mat = rep(0, ncol1 * nrow1)
  sequence = seq(1, nrow1)
  junk = sample(sequence)
  for (i in 1:nrow1) {
    mat[(i * ncol1 - (ncol1 - 1)):(i * ncol1)] = x[junk[i], 
    ]
  }
  mat = matrix(mat, ncol = ncol1, byrow = T)
  return(mat)
}

permfn2_extended <- function (x, p1, p2, m, type, exttype, maxlag, start = 1, alpha = 0.05) {
  results <- list() 
  if (type == 1) {
    for (i in 1:m) {
      pBACC = permatrix(x)
      cc = extremogram2_extended(pBACC, p1, p2, maxlag, exttype)
      lines((start:(maxlag - 1)), cc[(start + 1):maxlag], 
            col = 1, lwd = 1)
    }
  }
  if (type == 2) {
    cc = matrix(0, ncol = m, nrow = maxlag)
    for (i in 1:m) {
      pBACC = permatrix(x)
      cc[, i] = extremogram2_extended(pBACC, p1, p2, maxlag, exttype)
    }
    k = dim(cc)[1]
    pocket = matrix(0, ncol = 3, nrow = k)
    for (i in 1:k) {
      pocket[i, 1] = quantile(cc[i, ], prob = (alpha/2))
      pocket[i, 3] = quantile(cc[i, ], prob = (1 - alpha/2))
    }
    lines(start:(k - 1 + start), pocket[start:(k - 1 + start), 
                                        2], col = 1, lwd = 2)
    lines(start:(k - 1 + start), pocket[start:(k - 1 + start), 
                                        3], col = 1, lwd = 2)
  }
  if (type == 3) {
    cc = matrix(0, ncol = m, nrow = (start + 1))
    for (i in 1:m) {
      pBACC = permatrix(x)
      cc[, i] = extremogram2_extended(pBACC, p1, p2, (start + 1), exttype)
    }
    dde = as.numeric()
    gge = as.numeric()
    for (i in 1:maxlag) {
      dde[i] = quantile(cc[(start + 1), ], prob = (alpha/2))
      gge[i] = quantile(cc[(start + 1), ], prob = (1 - 
                                                     alpha/2))
    }
    
    dde_value = quantile(cc[2, ], prob = alpha/2)
    gge_value = quantile(cc[2, ], prob = 1 - alpha/2)
    
    
    results$type3 = data.frame(
      lag = 1:(maxlag - 1), 
      lower = rep(dde_value, maxlag - 1),
      upper = rep(gge_value, maxlag - 1)
    )
  }
}

plot_extremogram2 <- function(x, quant1, quant2, maxlag, type, typeper, R = 1000, l = 30, exttype = 1, m=99, cutoff, title, color) {
  
  ext <- extremogram2_extended(a = x, quant1 = quant1, quant2 = quant2, maxlag = maxlag, type = type)
  per <-permfn2_extended(x=x, p1 = quant1, p2=quant2, m = m, type = typeper, exttype = exttype, maxlag, 1, 0.05)
  boot <- bootconf2_extended(x = x, R = R, l = l, maxlag = maxlag, quant1 = quant1, quant2 = quant2, type = type, par = 0)
  
  data_df <- data.frame(
    extremogram = ext[-1],
    lag = boot$lag,
    per_lower = per$lower,
    per_upper = per$upper,
    boot_upper = boot$upper_confidence,
    boot_lower = boot$lower_confidence, 
    boot_mean = boot$mean
  )
  
  plot <- ggplot(data_df, aes(x = lag)) +
    geom_line(aes(y = boot_mean), color = color, size=0.75) +
    geom_line(aes(y = extremogram), color = "black") +
    geom_segment(aes(xend = lag, y = 0, yend = extremogram), color = "grey60") +
    geom_line(aes(y = per_lower), lty=2) +
    geom_line(aes(y = per_upper), lty=2) +
    geom_line(aes(y = boot_lower), lty=2, alpha=0.8, color = color, size=0.75) +
    geom_line(aes(y = boot_upper), lty=2, alpha=0.8, color = color, size=0.75) +
    ylim(c(0,cutoff))+
    labs(title = title, x = "Lag", y = "Extremogram") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size=18, hjust = 0.75))
  
  return(plot)
  
}




# Return times extremograms 
Extremogramr_extended <- function (x, type, maxlag, uplevel = 1, lowlevel = 0, histogram = 0, cutoff = 1) {
  n = length(x)
  mat = mat = rep(0, n)
  if (type == 1) {
    uplevel = quantile(x, prob = uplevel)
    for (i in 1:n) {
      mat[i] = ifelse(x[i] > uplevel, 1, 0)
    }
  }
  else if (type == 2) {
    lowlevel = quantile(x, prob = lowlevel)
    for (i in 1:n) {
      mat[i] = ifelse(x[i] < lowlevel, 1, 0)
    }
  }
  else if (type == 3) {
    uplevel = quantile(x, prob = uplevel)
    lowlevel = quantile(x, prob = lowlevel)
    for (i in 1:n) {
      mat[i] = ifelse(x[i] > uplevel || x[i] < lowlevel, 
                      1, 0)
    }
  }
  sequence = seq(1, n)
  gar = cbind(sequence, mat)
  junk = matrix(0, ncol = 2, nrow = n)
  for (i in 1:n) {
    if (mat[i] == 1) {
      junk[i, ] = gar[i, ]
    }
  }
  ind <- rowSums(junk == 0) != ncol(junk)
  junk = junk[ind, ]
  n = dim(junk)[1]
  return_time = rep(0, n - 1)
  for (i in 1:n - 1) {
    return_time[i] = (junk[i + 1, 1] - junk[i, 1])
  }
  if (histogram == 1) {
    MASS::truehist(return_time, nbins = max(return_time), 
                   xlim = c(0, (maxlag - 1)), col = 0, prob = TRUE, 
                   ylim = c(0, cutoff))
  }
  aa = as.matrix(table(return_time))
  for (j in 1:dim(aa)[1]) {
    aa[j, 1] = aa[j, 1]/length(return_time)
  }
  aa = as.double(aa)[1:maxlag]
  aa[is.na(aa)] <- 0
  return(list(aa, return_time, mean(return_time)))
}

bootconfr_extended <- function (x, R, l, maxlag, uplevel = 1, lowlevel = 0, type, par, start = 1, cutoff = 1, alpha = 0.05) {
  library(boot)
  
  if (par == 1) {
    n = parallel::detectCores()
    boot <- boot::tsboot(x, permbootr, R, l = l, sim = "geom", 
                         endcorr = TRUE, maxlag = maxlag, uplevel = uplevel, 
                         lowlevel = lowlevel, type = type, parallel = "snow", 
                         ncpus = n)
  } else {
    boot <- boot::tsboot(x, permbootr, R, l = l, sim = "geom", 
                         endcorr = TRUE, maxlag = maxlag, uplevel = uplevel, 
                         lowlevel = lowlevel, type = type)
  }
  
  tmp <- boot$t0
  mat <- boot$t
  k <- min(ncol(mat), maxlag)
  pocket <- data.frame(
    Lag = start:(k + start - 1),
    Lower = apply(mat[, 1:k, drop = FALSE], 2, quantile, prob = alpha/2),
    Mean = colMeans(mat[, 1:k, drop = FALSE]),
    Upper = apply(mat[, 1:k, drop = FALSE], 2, quantile, prob = 1 - alpha/2)
  )
  
  return(pocket)
}

permb_extended <- function (x, m, type, exttype, maxlag, uplevel = 1, lowlevel = 0, start = 1, alpha = 0.05) {
  results_list <- list()  
  
  if (type == 1) {
    results <- matrix(0, nrow = maxlag, ncol = m) 
    for (i in 1:m) {
      pBACC = sample(x)
      results[, i] <- permbootr(pBACC, exttype, uplevel, lowlevel, maxlag)
    }
    results_list$bootstrap_results = results
  }
  
  if (type == 2) {
    cc = matrix(0, ncol = m, nrow = maxlag)
    for (i in 1:m) {
      pBACC = sample(x)
      cc[, i] <- permbootr(pBACC, exttype, uplevel, lowlevel, maxlag)
    }
    pocket = matrix(0, ncol = 3, nrow = maxlag)
    for (i in 1:maxlag) {
      pocket[i, 1] <- quantile(cc[i, ], prob = alpha/2)
      pocket[i, 2] <- mean(cc[i, ])
      pocket[i, 3] <- quantile(cc[i, ], prob = 1 - alpha/2)
    }
    results_list$pocket = pocket
  }
  
  if (type == 3) {
    cc = matrix(0, ncol = m, nrow = maxlag)
    for (i in 1:m) {
      pBACC = sample(x)
      cc[, i] <- permbootr(pBACC, exttype, uplevel, lowlevel, maxlag)
    }
    dde = quantile(cc[start, ], prob = alpha/2)
    gge = quantile(cc[start, ], prob = 1 - alpha/2)
    results_list$lower_bounds = dde
    results_list$upper_bounds = gge
  }
  
  return(results_list)
}




# Spectral measure functions and plot
Spectral_measure <- function(data, uplevel, lowlevel, no_bins, lag = 0, type, bootstrap=TRUE)
{
  #first column is lagged 
  data <- data.frame(cbind(x=data[1:(dim(data)[1]- lag), 1], 
                           y=data[(1+lag):dim(data)[1], 2]))
  
  scatter_dataextreme<- data %>%
    mutate(
      extreme_x = x > quantile(x, uplevel) | x < quantile(x, lowlevel),
      extreme_y = y > quantile(y, uplevel) | y < quantile(y, lowlevel),
      any_extreme = extreme_x | extreme_y,
      both_extreme = extreme_x & extreme_y,
      justone_extreme = extreme_x & extreme_y==FALSE | extreme_y & extreme_x==FALSE,
      pos_extreme = (extreme_x & x >= 0) | (extreme_y & y >= 0),
      neg_extreme = (extreme_x & x < 0) | (extreme_y & y < 0)
    )
  
  if (type == 1){
    data <- scatter_dataextreme %>%
      filter(any_extreme)
    
  } else if (type == 2){
    data <- scatter_dataextreme %>%
      filter(both_extreme)
  }
  
  else if (type == 3){
    data <- scatter_dataextreme %>%
      filter(justone_extreme)
  }
  
  angle_radians <- atan2(data[,2], data[,1]) #atan2 takes care of the different quadrants 
  
  data <- data %>%
    mutate(
      angle_circle = angle_radians   
    )
  
  bins <- seq(-pi, pi, by = pi/no_bins)
  Count <- cut(data$angle_circle, bins)
  df <- data.frame(table(Count),
                   midpoints  <- (bins[-length(bins)] + bins[-1]) / 2)
  
  rownames(df) <- NULL
  names(df) <- c('direction', 'magnitude', 'midpoints')
  df$magnitude <- df$magnitude/length(data$angle_circle)
  
  
  return( if (bootstrap == TRUE) {df$magnitude} else {df})
}


bootconf_spectral <- function (data, R, l, lag, uplevel = 0.95, lowlevel = 0.05, no_bins, type, 
                               start = 1, cutoff = 1, alpha = 0.05) 
{
  library(boot)
  
  boot = boot::tsboot(tseries = data, statistic = Spectral_measure, R = R, l = l, sim = "geom", 
                      endcorr = TRUE, uplevel = uplevel, lowlevel = lowlevel, no_bins = no_bins, 
                      lag = lag, type = type)
  
  tmp = boot[[2]]
  k = dim(tmp)[2]
  pocket = matrix(0, ncol = 3, nrow = k)
  for (i in 1:k) {
    pocket[i, 1] = quantile(tmp[, i], prob = (alpha/2))
    pocket[i, 2] = mean(tmp[, i])
    pocket[i, 3] = quantile(tmp[, i], prob = (1 - alpha/2))
  }
  colnames(pocket) <- c("lower", "mean", "upper")
  return(pocket)
}

plot_spectral_measure <- function(data, lag, type, title, R, l, color=color, no_bins=50)
{
  df <- Spectral_measure(data, uplevel =0.95, no_bins=no_bins, lowlevel = 0.05, type = type, 
                         lag=lag, bootstrap = FALSE)
  conf <- bootconf_spectral(data, R = R, l = l, uplevel = 0.95, lowlevel=0.05, no_bins=no_bins, 
                            type = type, lag=lag)
  df = data.frame(df, conf)
  
  plot <- ggplot(df) +
    geom_bar(aes(x = midpoints, y = magnitude), width = 0.07, stat="identity",
             colour = "black", fill=color) +
    coord_polar(theta = "x", start =pi/2, direction=-1) + 
    labs(y = "Magnitude", x = "Angle", title = title) +
    scale_x_continuous(
      breaks = c(-pi/2, 0, pi/2, pi),
      labels = c(expression(-pi/2), "0", expression(pi/2), expression(pi)))+
    geom_col(aes(x = midpoints, y = magnitude), position = position_dodge(width=0.2), 
             colour = "black", fill=color)+
    geom_errorbar(aes(x = midpoints, y = magnitude, ymin = lower, ymax = upper), width = 0.05)+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1, 1, 1, 1), "mm"))
  
  
  return(plot)
  
}







