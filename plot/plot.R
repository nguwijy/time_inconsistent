rm(list = ls())

library(utils)
library(dplyr)

# 1 for quadratic
# 2 for quartic
which_prob <- 1

# set working directory as the file directory
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
setwd(dirname(script.name))

# for the exact solution...
x0 <- 0.0
t0 <- 0.0
Upper_bound_coeff <- 2
Upper_bound_func <- Upper_bound_coeff^10
Tt = 0.1
a <- 1.0
b <- 1.0
sigma <- 1.0
gamma <- 1.0

beta = function(s) {
  return(gamma*exp(a*(Tt-s)))
}

integrand_1 = function(s) {
  if (a==0) {
    ans <- exp(b*b*gamma*(Tt-s))*gamma*exp(2*a*(Tt-s))
  } else {
    ans <- exp(b*b*gamma/a*(exp(a*(Tt-s)) - 1))*gamma*exp(2*a*(Tt-s))
  }

  return(ans)
}

alpha = function(s) {
  integrated_1 <- rep(0, length(s))
  for (iter in 1:length(s)) {
    integrated_1[iter] <- integrate(integrand_1, s[iter], Tt)$value
  }
  if (a==0) {
    w <- exp(-b*b*gamma*(Tt-s))*(1+ b*b*integrated_1)
  } else {
    w <- exp(-b*b*gamma/a*(exp(a*(Tt-s)) - 1)) * (1+ b*b*integrated_1)
  }
  return(gamma*exp(2*a*(Tt-s))/w)
}

integrand_2 = function(s) {
  return(a+b*b*(beta(s)-alpha(s)))
}

gamma_fun = function(t1, t2) {
  if (length(t1) != 1) {
    integrated_2 <- rep(0, length(t1))
    for (iter in 1:length(t1)) {
      integrated_2[iter] <- integrate(integrand_2, t1[iter], t2)$value
    }
  } else if (length(t2) != 1) {
    integrated_2 <- rep(0, length(t2))
    for (iter in 1:length(t2)) {
      integrated_2[iter] <- integrate(integrand_2, t1, t2[iter])$value
    }
  } else {
    integrated_2 <- integrate(integrand_2, t1, t2)$value
  }
  return(exp(integrated_2))
}

gamma_fun_squared = function(t1, t2) {
  return(gamma_fun(t1, t2)*gamma_fun(t1, t2))
}

integrand_3 = function(u, s) {
  return(gamma_fun(u,s)*gamma_fun(u,s))
}

integrand_4 = function(s, t, x) {
  integrated_3 <- rep(0, length(s))
  for (iter in 1:length(s)) {
    integrated_3[iter] <- integrate(integrand_3, t, s[iter], s = s[iter])$value
  }
  return((beta(s) - alpha(s))*(beta(s) - alpha(s))*
           (gamma_fun(t,s)*x*gamma_fun(t,s)*x + sigma*sigma*integrated_3))
}

exact_sol_fun = function(t, x) {
  integrated_4 <- integrate(integrand_4, t, Tt, t = t, x = x)$value
  integrated_3 <- integrate(integrand_3, t, Tt, s = Tt)$value
  return(integrated_4/2 + sigma*sigma*integrated_3*gamma/2)
}

mean_u = function(t, x) {
  return(b*(beta(t) - alpha(t))*(x*gamma_fun(0,t)))
}

sig_u = function(t) {
  var_x <- integrate(gamma_fun_squared, lower = 0, upper = t, t2 = t)$value*sigma*sigma
  return(abs(b*(beta(t) - alpha(t)))*sqrt(var_x))
}

min_u <- -0.3
max_u <- 0.3
group_width <- 0.001
uu_grid <- seq(min_u, max_u, by = group_width)
uu_grid <- ceiling((uu_grid-group_width/2)/group_width)*group_width

tt_N <- 400
del_t <- (Tt - t0)/tt_N
tt_grid <- seq(t0, Tt, by = del_t)

cut_off_dense <- 500
cut_off_sd <- 4*1e-4

exact_gen = function() {
  merge_grid <- expand.grid(control = uu_grid, tt_now = tt_grid)
  merge_grid$probability <- 0
  for (titer in 1:length(tt_grid)) {
    print(titer)
    sd_now <- sig_u(tt_grid[titer])
    for (uiter in 1:length(uu_grid)) {
      index <- (titer - 1)*length(uu_grid) + uiter
      if (sd_now > cut_off_sd) {
        # these two ways are the same by substitution, both it seems that working with standard normal distribution has less error
        # merge_grid$probability[iter] <- pnorm((merge_grid$control[iter]+group_width/2), sd = sd_now) - pnorm((merge_grid$control[iter]-group_width/2), sd = sd_now)
        merge_grid$probability[index] <- pnorm((merge_grid$control[index]+group_width/2)/sd_now) - pnorm((merge_grid$control[index]-group_width/2)/sd_now)
        merge_grid$density[index] <- dnorm(merge_grid$control[index], sd = sd_now)
        merge_grid$cdf[index] <- pnorm(merge_grid$control[index], sd = sd_now)
      } else {
        merge_grid$probability[index] <- ifelse(merge_grid$control[index]==0, 1, 0)
        merge_grid$density[index] <- ifelse(merge_grid$control[index]==0, cut_off_dense, 0)
        merge_grid$cdf[index] <- ifelse(merge_grid$control[index]>=0, 1, 0)
      }
    }
  }
  return(merge_grid)
}


# for the N-person differential game
Nperson_terminal_alpha_gen = function(N) {
  ans <- rep(0, N)
  deldel <- Tt/N
  ans[N] <- gamma
  for(ii in seq(N-1, 1, by = -1)) {
    integrated <- integrate(Nperson_drift, lower = ii*deldel, upper = (ii+1)*deldel, t2 = (ii+1)*deldel, N = N, terminal_vect = ans)$value
    ans[ii] <- ans[ii+1]*(Nperson_gamma_fun(N, ii*deldel, (ii+1)*deldel, ans) + b*integrated)*exp(a*deldel)
  }
  return(ans)
}

Nperson_square = function(N, s, terminal_vect) {
  ans <- rep(0, length(s))
  deldel <- Tt/N
  for (iter in 1:length(s)) {
    round_N <- floor(s[iter]/deldel+1e-10)
    if (round_N > 0) {
      for (iiter in 1:round_N) {
        ans[iter] <- ans[iter]*(Nperson_gamma_fun(N, (iiter-1)*deldel, iiter*deldel, terminal_vect) +
                        b*integrate(Nperson_drift, (iiter-1)*deldel, iiter*deldel, N = N, t2 = iiter*deldel, terminal_vect = terminal_vect)$value)^2 +
                        integrate(Nperson_gamma_fun_squared, (iiter-1)*deldel, iiter*deldel, N = N, t2 = iiter*deldel, terminal_vect = terminal_vect)$value*
                        sigma*sigma
      }
    }
    ans[iter] <- ans[iter]*(Nperson_gamma_fun(N, round_N*deldel, s[iter], terminal_vect) +
                                b*integrate(Nperson_drift, round_N*deldel, s[iter], N = N, t2 = s[iter], terminal_vect = terminal_vect)$value)^2 +
                                integrate(Nperson_gamma_fun_squared, round_N*deldel, s[iter], N = N, t2 = s[iter], terminal_vect = terminal_vect)$value*
                                sigma*sigma
  }
  return(ans)
}

Nperson_cross = function(N, s, terminal_vect) {
  ans <- rep(0, length(s))
  deldel <- Tt/N
  for (iter in 1:length(s)) {
    round_N <- floor(s[iter]/deldel+1e-10)
    if (round_N > 0) {
      for (iiter in 1:round_N) {
        ans[iter] <- ans[iter]*(Nperson_gamma_fun(N, (iiter-1)*deldel, iiter*deldel, terminal_vect) +
                                  b*integrate(Nperson_drift, (iiter-1)*deldel, iiter*deldel, N = N, t2 = iiter*deldel, terminal_vect = terminal_vect)$value)^2 +
                                  integrate(Nperson_gamma_fun_squared, (iiter-1)*deldel, iiter*deldel, N = N, t2 = iiter*deldel, terminal_vect = terminal_vect)$value*
                                  sigma*sigma
      }
    }
    ans[iter] <- ans[iter]*(Nperson_gamma_fun(N, round_N*deldel, s[iter], terminal_vect) +
                              b*integrate(Nperson_drift, round_N*deldel, s[iter], N = N, t2 = s[iter], terminal_vect = terminal_vect)$value)
  }
  return(ans)
}

Nperson_u = function(N, s, terminal_vect) {
  ans <- rep(0, length(s))
  deldel <- Tt/N
  for (iter in 1:length(s)) {
    round_N <- floor(s[iter]/deldel+1e-10)
    ans[iter] <- b*b*Nperson_beta(N, s[iter], terminal_vect)*Nperson_beta(N, s[iter], terminal_vect)*Nperson_square(N, round_N*deldel, terminal_vect) +
                  b*b*Nperson_alpha(N, s[iter], terminal_vect)*Nperson_alpha(N, s[iter], terminal_vect)*Nperson_square(N, s[iter], terminal_vect) -
                  2*b*b*Nperson_alpha(N, s[iter], terminal_vect)*Nperson_beta(N, s[iter], terminal_vect)*Nperson_cross(N, s[iter], terminal_vect)
  }
  return(ans)
}

Nperson_value_func = function(N, terminal_vect) {
  ans <- integrate(Nperson_u, 0, Tt, N = N, terminal_vect = terminal_vect)$value/2 + gamma/2*Nperson_square(N, Tt, terminal_vect)
  return(ans)
}

Nperson_alpha = function(N, s, terminal_vect) {
  ans <- rep(0, length(s))
  for (iter in 1:length(s)) {
    if (s[iter] == Tt) {
      ans[iter] <- gamma
    } else {
      deldel <- Tt/N
      round_t <- floor(s[iter]/deldel)+1
      tt_ter <- round_t*deldel
      temp_ter <- terminal_vect[round_t]
      C <- tt_ter + log(temp_ter/(2*a-b*b*temp_ter))/(2*a)
      ans[iter] <- (2*a*exp(2*a*C)/(b*b*exp(2*a*C)+exp(2*a*s[iter])))
    }
  }
  return(ans)
}

Nperson_beta = function(N, s, terminal_vect) {
  ans <- rep(0, length(s))
  for (iter in 1:length(s)) {
    if (s[iter] == Tt) {
      ans[iter] <- gamma
    } else {
      deldel <- Tt/N
      round_t <- floor(s[iter]/deldel)+1
      tt_ter <- round_t*deldel
      temp_ter <- gamma*exp(a*(Tt - tt_ter))
      integrated <- integrate(Nperson_alpha, s[iter], tt_ter, N = N, terminal_vect = terminal_vect)$value
      ans[iter] <- temp_ter*exp(a*(tt_ter - s[iter]) - b*b*integrated)
    }
  }
  return(ans)
}

Nperson_gamma_fun = function(N, t1, t2, terminal_vect) {
  if (length(t1) != 1) {
    integrated_2 <- rep(0, length(t1))
    for (iter in 1:length(t1)) {
      integrated_2[iter] <- integrate(Nperson_alpha, t1[iter], t2, N = N, terminal_vect = terminal_vect)$value
    }
  } else if (length(t2) != 1) {
    integrated_2 <- rep(0, length(t2))
    for (iter in 1:length(t2)) {
      integrated_2[iter] <- integrate(Nperson_alpha, t1, t2[iter], N = N, terminal_vect = terminal_vect)$value
    }
  } else {
    integrated_2 <- integrate(Nperson_alpha, t1, t2, N = N, terminal_vect = terminal_vect)$value
  }
  return(exp(a*(t2-t1) - b*b*integrated_2))
}

Nperson_gamma_fun_squared = function(N, t1, t2, terminal_vect) {
  return(Nperson_gamma_fun(N, t1, t2, terminal_vect)*Nperson_gamma_fun(N, t1, t2, terminal_vect))
}

Nperson_drift = function(N, t1, t2, terminal_vect) {
  return(Nperson_gamma_fun(N, t1, t2, terminal_vect)*Nperson_beta(N, t1, terminal_vect))
}

Nperson_mean_u = function(t, x) {
  return(0)
}

Nperson_sig_u = function(N, t, terminal_vect) {
  deldel <- Tt/N
  round_t <- floor(t/deldel)
  temp_var <- 0
  for (i in 0:round_t) {
    cur_round_t <- i * deldel
    if (i > 0) {
      drift_integrate <- integrate(Nperson_drift, lower = cur_round_t-deldel, upper = cur_round_t, t2 = cur_round_t, N = N, terminal_vect = terminal_vect)$value
      var_integrate <- integrate(Nperson_gamma_fun_squared, lower = cur_round_t-deldel, upper = cur_round_t, t2 = cur_round_t, N = N, terminal_vect = terminal_vect)$value
      temp_var <- temp_var * (Nperson_gamma_fun(N, cur_round_t-deldel, cur_round_t, terminal_vect) + b*drift_integrate)^2 +
                  var_integrate*sigma*sigma
    }
  }
  cur_round_t <- round_t * deldel
  drift_integrate <- integrate(Nperson_drift, lower = cur_round_t, upper = t, t2 = t, N = N, terminal_vect = terminal_vect)$value
  var_integrate <- integrate(Nperson_gamma_fun_squared, lower = cur_round_t, upper = t, t2 = t, N = N, terminal_vect = terminal_vect)$value
  total_var <- (b*Nperson_beta(N, t, terminal_vect) - b*Nperson_alpha(N, t, terminal_vect)*
                (Nperson_gamma_fun(N, cur_round_t, t, terminal_vect = terminal_vect) + b*drift_integrate))^2 * temp_var +
                (b*Nperson_alpha(N, t, terminal_vect))^2*var_integrate*sigma*sigma
  return(sqrt(total_var))
}

N_person_gen = function(N_interest) {
  N_interest <- 20
  Nperson_merge_grid <- expand.grid(control = uu_grid, tt_now = tt_grid)
  Nperson_merge_grid$probability <- 0
  terminal_alpha <- Nperson_terminal_alpha_gen(N_interest)
  for (titer in 1:length(tt_grid)) {
    print(titer)
    sd_now <- Nperson_sig_u(N_interest, tt_grid[titer], terminal_alpha)
    for (uiter in 1:length(uu_grid)) {
      index <- (titer - 1)*length(uu_grid) + uiter
      if (sd_now > cut_off_sd) {
        # these two ways are the same by substitution, both it seems that working with standard normal distribution has less error
        # merge_grid$probability[iter] <- pnorm((merge_grid$control[iter]+group_width/2), sd = sd_now) - pnorm((merge_grid$control[iter]-group_width/2), sd = sd_now)
        Nperson_merge_grid$probability[index] <- pnorm((Nperson_merge_grid$control[index]+group_width/2)/sd_now) - pnorm((Nperson_merge_grid$control[index]-group_width/2)/sd_now)
        Nperson_merge_grid$density[index] <- dnorm(Nperson_merge_grid$control[index], sd = sd_now)
        Nperson_merge_grid$cdf[index] <- pnorm(Nperson_merge_grid$control[index], sd = sd_now)
      } else {
        Nperson_merge_grid$probability[index] <- ifelse(Nperson_merge_grid$control[index]==0, 1, 0)
        Nperson_merge_grid$density[index] <- ifelse(Nperson_merge_grid$control[index]==0, cut_off_dense, 0)
        Nperson_merge_grid$cdf[index] <- ifelse(Nperson_merge_grid$control[index]>=0, 1, 0)
      }
    }
  }
  return(Nperson_merge_grid)
}


# for quadratic_cdf.txt and quartic_cdf.txt
# read files
main_file <- ifelse(which_prob==1,  "../final_log/T_0.1quadratic_optimal_U.txt", "../final_log/T_0.1quartic_optimal_U.txt")
datafr <- read.table( main_file, sep=",", header=T )

# select only 0, 0.02, 0.04 etc.
datafr <- datafr[(datafr$M == 20) & (((datafr$tt_now*5*10)%%1) == 0), ]

# data manipulations
# since when doing mesh in python, some errors might be introduced
datafr$tt_now <- signif(datafr$tt_now, digits = 8)
datafr$control_round <- ceiling((datafr$control-group_width/2)/group_width)*group_width
if ( which_prob == 2 ) {
    # ONLY for quartic, for better graph
    datafr[datafr$tt_now==0.02,]$control_round <- ceiling((datafr[datafr$tt_now==0.02,]$control-0.005/2)/0.005)*0.005
    datafr[datafr$tt_now==0.04,]$control_round <- ceiling((datafr[datafr$tt_now==0.04,]$control-0.002/2)/0.002)*0.002
    datafr[datafr$tt_now==0.06,]$control_round <- ceiling((datafr[datafr$tt_now==0.06,]$control-0.001/2)/0.001)*0.001
    datafr[datafr$tt_now==0.08,]$control_round <- ceiling((datafr[datafr$tt_now==0.08,]$control-0.0001/2)/0.0001)*0.0001
}

# we use the old way of grouping control
# to avoid points too close to each other
# after grouping, we calculate the cdf
datafr <- aggregate(probability ~ control_round + tt_now + N + M, datafr, sum)
datafr$cum_prob <- ave(datafr$probability, datafr$tt_now, FUN=cumsum)

# only do this for quadratic since it has analytic solution
if ( which_prob == 1 ) {
    # calculate Nperson cdf
    terminal_alpha <- Nperson_terminal_alpha_gen(20)
    datafr$Nperson_cdf <- 100
    for (tt in unique(datafr$tt_now)) {
      print(tt)
      sd_now <- Nperson_sig_u(20, tt, terminal_alpha)
      if (sd_now > cut_off_sd) {
        datafr[datafr$tt_now==tt, ]$Nperson_cdf <- pnorm(datafr[datafr$tt_now==tt, ]$control_round, sd = sd_now)
      } else {
        datafr[datafr$tt_now==tt, ]$Nperson_cdf <- ifelse(datafr[datafr$tt_now==tt, ]$control_round>=0, 1, 0)
      }
    }

    # calculate cdf and the diff compared to N-person
    datafr$cdf_diff <- abs(datafr$cum_prob - datafr$Nperson_cdf)

    # calculate exact sol
    datafr$exact <- 100
    for (tt in unique(datafr$tt_now)) {
      print(tt)
      sd_now <- sig_u(tt)
      if (sd_now > cut_off_sd) {
        datafr[datafr$tt_now==tt, ]$exact <- pnorm(datafr[datafr$tt_now==tt, ]$control_round, sd = sd_now)
      } else {
        datafr[datafr$tt_now==tt, ]$exact <- ifelse(datafr[datafr$tt_now==tt, ]$control_round>=0, 1, 0)
      }
    }
}

if (which_prob == 1) {
    file.remove("quadratic_cdf.txt")
    write.table(datafr, file = "quadratic_cdf.txt", sep = " ",
                row.names = F, col.names = F, append = T)
}else if (which_prob == 2) {
    file.remove("quartic_cdf.txt")
    write.table(datafr, file = "quartic_cdf.txt", sep = " ",
                row.names = F, col.names = F, append = T)
}

# for quadratic_value_func.txt and quartic_value_func.txt
# read files
main_file <- ifelse(which_prob==1, "../final_log/T_0.1_quadratic.txt", "../final_log/T_0.1_quartic.txt")
datafr <- read.table( main_file, sep=",", header=T )

# data manipulation
# since when doing mesh in python, some errors might be introduced
datafr$tt_now <- signif(datafr$tt_now, digits = 8)
final <- aggregate(runtime ~ N + M, datafr, sum)
final <- merge(datafr[datafr$tt_now==0, -6], final, by = c('N', 'M'))
final <- final[order(final$N, final$M), ]

# only do this for quadratic since quartic has no analytic solution
if ( which_prob == 1 ) {
    exact_sol <- exact_sol_fun(0, x0)
    final$exact_sol <- exact_sol
    final$error <- abs(final$numerical_sol-final$exact_sol)/final$exact_sol
}
if (which_prob == 1) {
    file.remove("quadratic_value_func.txt")
    write.table(final, file = "quadratic_value_func.txt", sep = " ",
                row.names = F, col.names = F, append = T)
}else if (which_prob == 2) {
    file.remove("quartic_value_func.txt")
    write.table(final, file = "quartic_value_func.txt", sep = " ",
                row.names = F, col.names = F, append = T)
}


# for quadratic_exact_sol.txt and quadratic_20person.txt
if (which_prob == 1) {
    file.remove("quadratic_exact_sol.txt")
    exact_grid <- exact_gen()
    grid_ordered <- sort(unique(exact_grid$tt_now))
    for (i in 1:length(grid_ordered)) {
      tt <- grid_ordered[i]
      curdata <- cbind(exact_grid[exact_grid$tt_now == tt, ]$control, exact_grid[exact_grid$tt_now == tt, ]$tt_now, exact_grid[exact_grid$tt_now == tt, ]$density)
      write.table(curdata, file = "quadratic_exact_sol.txt", sep = " ",
                  row.names = F, col.names = F, append = T)
      if (i != length(grid_ordered))  cat("",file="quadratic_exact_sol.txt",sep="\n", append = T)
    }

    N_interested <- 20
    N_name <- paste("quadratic_", N_interested, "person.txt", sep = "")
    file.remove(N_name)
    Nperson_grid <- N_person_gen(N_interested)
    grid_ordered <- sort(unique(Nperson_grid$tt_now))
    for (i in 1:length(grid_ordered)) {
      tt <- grid_ordered[i]
      curdata <- cbind(Nperson_grid[Nperson_grid$tt_now == tt, ]$control, Nperson_grid[Nperson_grid$tt_now == tt, ]$tt_now, Nperson_grid[Nperson_grid$tt_now == tt, ]$density)
      write.table(curdata, file = N_name, sep = " ",
                  row.names = F, col.names = F, append = T)
      if (i != length(grid_ordered))  cat("",file=N_name,sep="\n", append = T)
    }
}
