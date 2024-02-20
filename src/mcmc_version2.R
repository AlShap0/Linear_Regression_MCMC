---
title: An R Markdown document converted from "mcmc_version2.ipynb"
output: html_document
---

```{r}
library(ggplot2)
library(purrr)
library(gridExtra)
library(latex2exp)
options(repr.plot.width=12, repr.plot.height=8)
```

```{r}
# Define all models 1-4
linear_model = function(params, x, z, noise){
    output = params[1] + params[2]*x 
    if (noise == TRUE){
        output = output + rnorm(length(x), mean=0, sd=params[3])
    }
    return(output)}

quad_model = function(params, x, z, noise){
    output = params[1] + params[2]*x + params[3]*x^2
    if (noise == TRUE){
        output = output + rnorm(length(x), mean=0, sd=params[4])
    }
    return(output)}

mix_model = function(params, x, z, noise){
    output = params[1] + params[2]*x + params[3]*z
    if (noise == TRUE){
        output = output + rnorm(length(x), mean=0, sd=params[4])
    }
    return(output)}

mix_interaction_model = function(params, x, z, noise){
    output = params[1] + params[2]*x + params[3]*z + params[4]*x*z
    if (noise == TRUE){
        output = output + rnorm(length(x), mean=0, sd=params[5])
    }
    return(output)}
```

```{r}
log_likelihood = function(x, y, z, params, model){
    y_pred = model(params, x, z, noise = FALSE)
    like = c()
    for (i in 1:length(y)){
        like = c(like, dnorm(y[i], mean = y_pred[i], sd = tail(params, n=1), log = TRUE))
    } 
    return(sum(like))
}

log_prior = function(params, unif_min, unif_max, std_max){
    len = length(params)
    p = c()
    for (i in 1:len-1){
       p = c(p, dunif(params[i], unif_min, unif_max, log = TRUE))
    }
    p = c(p, dunif(params[len], 0, std_max, log = TRUE))
    return(sum(p))
}

log_prior_2 = function(params, guess_params, std_prior){
    len = length(params)
    p = c()
    for (i in 1:len){
       p = c(p, dnorm(params[i], mean = guess_params[i], sd = std_prior[i], log = TRUE))
    }
    return(sum(p))
}

log_post = function(x, y, z, params, model, unif_min, unif_max, std_max){
    result = log_likelihood(x, y, z, params, model) + log_prior(params, unif_min, unif_max, std_max)
    return(result)
}

log_post_2 = function(x, y, z, params, model, guess_params, std_prior){
    result = log_likelihood(x, y, z, params, model) + log_prior_2(params, guess_params, std_prior)
    return(result)
}
```

```{r}
proposal_dist = function(params_now, std_prop){
    len = length(params_now)
    params_next = numeric(len)
    for (i in 1:len){
        params_next[i] = rnorm(1, params_now[i], std_prop[i])}
    return(params_next)
}


find_next_step = function(x, y, z, params_now, model, unif_min, unif_max, std_max, std_prop){
    # Draw next step from proposal dist
    params_next = proposal_dist(params_now, std_prop)
    
    r = log(runif(1))
    ratio = log_post(x, y, z, params_next, model, unif_min, unif_max, std_max) - 
            log_post(x, y, z, params_now, model, unif_min, unif_max, std_max)

    if (ratio > r){
    # Accept new stepe
        return(params_next)}
    else {
    # Keep current step
        return(params_now)}
}

find_next_step_2 = function(x, y, z, params_now, model, guess_params, std_prior, std_prop){
    # Draw next step from proposal dist
    params_next = proposal_dist(params_now, std_prop)
    
    r = log(runif(1))
    ratio = log_post_2(x, y, z, params_next, model, guess_params, std_prior) - 
            log_post_2(x, y, z, params_now, model, guess_params, std_prior)

    if (ratio > r){
    # Accept new stepe
        return(params_next)}
    else {
    # Keep current step
        return(params_now)}
}


MCMC = function(x, y, z, initial_params, model, unif_min, unif_max, std_max, std_prop,
                num_walkers, num_of_steps){

    walkers = list() # initialize vector to hold all chains
    
    for(i in 1:num_walkers){
        chain = list(initial_params + rnorm(length(initial_params), 0, 0.1*max(initial_params)))
        
        for (j in 2:num_of_steps){
          par = tail(chain, n=1)[[1]]
          chain[[j]] = find_next_step(x, y, z, par, model, unif_min, unif_max, std_max, std_prop)
        }
        
        walkers[[i]] = chain

    }
    #chain = np.delete(chain, [0,1])
    return(walkers)
}

MCMC_2 = function(x, y, z, initial_params, model, guess_params, std_prior, std_prop,
                num_walkers, num_of_steps){

    walkers = list() # initialize vector to hold all chains
    
    for(i in 1:num_walkers){
        chain = list(initial_params + rnorm(length(initial_params), 0, 0.1*max(initial_params)))
        
        for (j in 2:num_of_steps){
          par = tail(chain, n=1)[[1]]
          chain[[j]] = find_next_step_2(x, y, z, par, model, guess_params, std_prior, std_prop)
        }
        
        walkers[[i]] = chain

    }
    #chain = np.delete(chain, [0,1])
    return(walkers)
}

lm_generator = function(model, true_params, x_min, x_max, n_data, n_rep, if_z){
    df = data.frame()
    for (i in 1:n_rep){
        set.seed(as.integer(runif(1, 1, 10^6)))
        x = runif(n_data, x_min, x_max)
        if (if_z == TRUE){
            z = sample(0:1, n_data, replace = TRUE)
        }
        else {
            z = numeric(n_data)
        }
        y = model(true_params, x, z, noise = TRUE)

        lr = lm(y~x+I(x^2))
        df = rbind(df, c(rbind(lr$coefficients), summary(lr)$sigma))
    }
    return(df)
}

mcmc_generator = function(initial_params, true_params, model, unif_min, unif_max, std_max, std_prop,
                       n_walkers, n_iter, n_rep, n_data, x_min, x_max, burn_in, thin){
    df_main = data.frame()
    for (i in 1:n_rep){
        set.seed(as.integer(runif(1, 1, 10^6)))
        x = runif(n_data, x_min, x_max)
        if (if_z == TRUE){
            z = sample(0:1, n_data, replace = TRUE)
        }
        else {
            z = numeric(n_data)
        }
        y = model(true_params, x, z, noise = TRUE)

        walkers = MCMC(x, y, z, initial_params, model, unif_min, unif_max, std_max, std_prop,
                       n_walkers, n_iter)

        df = data.frame(nrow = length(true_params))
        for (j in 1:n_walkers){
            walk = walkers[[j]] %>% transpose() %>% map(unlist) 
            dff = data.frame()
            for (z in 1:length(true_params)){
                walk_p = walk[[z]]
                walk_p = walk_p[burn_in:n_iter]
                len_walk = length(walk_p)
                ind = c()
                for (a in 1:len_walk){
                    if (a %% thin == 0){
                        ind = c(ind, a)
                    }
                }
                walk_p = walk_p[ind]
                dff = rbind(dff, walk_p)
            }
            df = cbind(df, dff)
        }
        df_main = rbind(df_main, rowMeans(df)) 
    }
    return(df_main)
}

mcmc_generator_2 = function(initial_params, model, true_params, std_prior, std_prop,
                n_walkers, n_iter, n_rep, n_data, x_min, x_max, burn_in, thin){
    df_main = data.frame()
    for (i in 1:n_rep){
        set.seed(as.integer(runif(1, 1, 10^6)))
        x = runif(n_data, x_min, x_max)
        if (if_z == TRUE){
            z = sample(0:1, n_data, replace = TRUE)
        }
        else {
            z = numeric(n_data)
        }
        y = model(true_params, x, z, noise = TRUE)

        walkers = MCMC_2(x, y, z, initial_params, model, true_params, std_prior, std_prop,
                n_walkers, n_iter)

        df = data.frame(nrow = length(true_params))
        for (j in 1:n_walkers){
            walk = walkers[[j]] %>% transpose() %>% map(unlist) 
            dff = data.frame()
            for (z in 1:length(true_params)){
                walk_p = walk[[z]]
                walk_p = walk_p[burn_in:n_iter]
                len_walk = length(walk_p)
                ind = c()
                for (a in 1:len_walk){
                    if (a %% thin == 0){
                        ind = c(ind, a)
                    }
                }
                walk_p = walk_p[ind]
                dff = rbind(dff, walk_p)
            }
            df = cbind(df, dff)
        }
        df_main = rbind(df_main, rowMeans(df)) 
    }
    return(df_main)
}
```

```{r}
# hyperparameters for generating the datasets
n_data = 50
true_params = c(2, 5, 10, 6)
model = quad_model
x_min = -5
x_max = 5

# hyperparameters for MCMC
initial_params = c(3, -5, 7, 4)
n_walkers = 5
n_iter = 5000
n_rep = 100
unif_min = -50
unif_max = 50
std_max = 10
if_z = FALSE
std_prior = 0.1*true_params
std_prop = c(0.05, 0.01, 0.01, 0.03)*true_params

burn_in = 2000
thin = 10
```

```{r}
set.seed(60603)
x = runif(n_data, min = x_min, max = x_max)
y = model(true_params, x, z, noise = TRUE)

ggplot(data.frame(x=x,y=y), aes(x=x,y=y)) + geom_point() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
                      plot.title = element_text(size=26), legend.text = element_text(size=18), 
                      legend.title = element_text(size=20)) + labs(x = 'x',y = 'y')
```

```{r}
#walkers = MCMC(x, y, z, initial_params, model, unif_min, unif_max, std_max, std_prop,
#                   n_walkers, n_iter)
walkers = MCMC_2(x, y, z, initial_params, model, true_params, std_prior, std_prop,
                n_walkers, n_iter)
```

```{r}
walk = walkers[[1]] %>% transpose() %>% map(unlist) 
p1 = walk[[1]]
p2 = walk[[2]]
p3 = walk[[3]]
p4 = walk[[4]]

y_pred = c()
lines = c()

mod = lm(y~x+I(x^2))
lm_params = as.numeric(mod$coefficients)
y_lm_pred = model(lm_params, x, z, FALSE)

for (i in 1:n_iter){
    y_pred = c(y_pred, model(c(p1[i], p2[i], p3[i], p4[i]), x, z, noise=FALSE))
    lines = c(lines, paste('line',i))
}
df = data.frame(x = rep(x, n_iter), y = y_pred, group = rep(lines, each = length(x)))

ggplot() + geom_point(data.frame(x=x,y=y), mapping=aes(x=x,y=y,size=10), color='blue') +
                geom_line(data = df, mapping = aes(x=x, y=y, color = 'blue', group=group),alpha=0.02) +
                geom_line(data = data.frame(x=x,y=y_lm_pred), aes(x=x,y=y), color='black') +
                theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
                      plot.title = element_text(size=26), legend.text = element_text(size=18), 
                      legend.title = element_text(size=20)) + labs(x = 'x',y = 'y')
```

```{r}
line_colors = c('darkgreen', 'red', 'royalblue', 'black', 'orange', 'pink', 'tan1', 'violet',
                'purple', 'yellow')
plot1 = ggplot() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
        plot.title = element_text(size=26), legend.text = element_text(size=18), 
        legend.title = element_text(size=20)) + labs(x = 'n_iteration',y = TeX(r'{$a$}'))
for (i in 1:n_walkers){
    walk = walkers[[i]] %>% transpose() %>% map(unlist) 
    walk_p1 = walk[[1]]
    plot1 = plot1 + geom_line(data = data.frame(x=1:n_iter, y=walk_p1), 
                            aes(x = x, y = y), color = line_colors[i])
}

plot2 = ggplot() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
        plot.title = element_text(size=26), legend.text = element_text(size=18), 
        legend.title = element_text(size=20)) + labs(x = 'n_iteration',y = TeX(r'{$b$}'))
for (i in 1:n_walkers){
    walk = walkers[[i]] %>% transpose() %>% map(unlist) 
    walk_p2 = walk[[2]]
    plot2 = plot2 + geom_line(data = data.frame(x=1:n_iter, y=walk_p2), 
                            aes(x = x, y = y), color = line_colors[i])
}

plot3 = ggplot() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
        plot.title = element_text(size=26), legend.text = element_text(size=18), 
        legend.title = element_text(size=20)) + labs(x = 'n_iteration',y = TeX(r'{$c$}'))
for (i in 1:n_walkers){
    walk = walkers[[i]] %>% transpose() %>% map(unlist) 
    walk_p3 = walk[[3]]
    plot3 = plot3 + geom_line(data = data.frame(x=1:n_iter, y=walk_p3), 
                            aes(x = x, y = y), color = line_colors[i])
}

plot4 = ggplot() + theme(axis.text = element_text(size=18), axis.title = element_text(size=18), 
        plot.title = element_text(size=26), legend.text = element_text(size=18), 
        legend.title = element_text(size=20)) + labs(x = 'n_iteration',y = TeX(r'{$\sigma$}'))
for (i in 1:n_walkers){
    walk = walkers[[i]] %>% transpose() %>% map(unlist) 
    walk_p4 = walk[[4]]
    plot4 = plot4 + geom_line(data = data.frame(x=1:n_iter, y=walk_p4), 
                            aes(x = x, y = y), color = line_colors[i])
}

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
```

```{r}
df_mcmc = mcmc_generator_2(initial_params, model, true_params, std_prior, std_prop, 
                      n_walkers, n_iter, n_rep, n_data, x_min, x_max, burn_in, thin)

esti_params = colMeans(df_mcmc)
x = runif(n_data, x_min, x_max)
z = numeric(n_data)
y = model(true_params, x, z, noise = TRUE)
y_pred = model(esti_params, x, z, noise = FALSE)

rmse_mcmc = sqrt(mean((y - y_pred)^2))
#perf = sum(abs(esti_params - true_params)) + sum(apply(df_mcmc, 2, sd)) 


for (i in 1:ncol(df_mcmc)){
    cat('parameter', i, ': \n', 'mean:', mean(df_mcmc[,i]), ' \ ', 'std:', sd(df_mcmc[,i]), '\n \n')
}
#cat('performance measure (smaller the better):', perf, '\n')
cat('RMSE:', rmse_mcmc)
```

```{r}
df_lm = lm_generator(model, true_params, x_min, x_max, n_data, n_rep, FALSE)

esti_params = colMeans(df_lm)
y_pred = model(esti_params, x, z, noise = FALSE)

rmse_lm = sqrt(mean((y - y_pred)^2))
#perf = sum(abs(colMeans(df_lm) - true_params)) + sum(apply(df_lm, 2, sd)) 
for (i in 1:ncol(df_lm)){
    cat('parameter', i, ': \n', 'mean:', mean(df_lm[,i]), ' \ ', 'std:', sd(df_lm[,i]), '\n \n')
}
#cat('performance measure (smaller the better):', perf, '\n')
cat('RMSE:', rmse_lm)
```

```{r}
cat('RMSE ratio:', rmse_mcmc/rmse_lm)
```

```{r}
bins = 20
plot = vector("list", length = length(true_params))
label = c(r'{$a$}', r'{$b$}', r'{$c$}', r'{$\sigma$}')
for (i in 1:length(true_params)){
  plot[[i]] = ggplot() + geom_histogram(data.frame(x=df_mcmc[,i]), 
                                                          mapping=aes(x=x, fill='MCMC'), bins=bins, alpha=0.3) + 
    geom_histogram(data.frame(x=df_lm[,i]), mapping=aes(x=x, fill='OLS'), bins=bins, alpha=0.3) + 
    geom_vline(xintercept = true_params[i], color = 'green') + 
    # std lines for MCMC
    geom_vline(xintercept = mean(df_mcmc[,i]) + sd(df_mcmc[,i]), color = 'blue', linetype="longdash") + 
    geom_vline(xintercept = mean(df_mcmc[,i]) - sd(df_mcmc[,i]), color = 'blue', linetype="longdash") + 
    # std lines for OLS
    geom_vline(xintercept = mean(df_lm[,i]) + sd(df_lm[,i]), color = 'red', linetype="longdash") + 
    geom_vline(xintercept = mean(df_lm[,i]) - sd(df_lm[,i]), color = 'red', linetype="longdash") + 
    scale_fill_manual(name="Legend", values = c("MCMC" = "blue", "OLS" = "red")) +
  theme(legend.position = "top", axis.text = element_text(size=10), axis.title = element_text(size=10), 
          plot.title = element_text(size=20), legend.text = element_text(size=10), 
          legend.title = element_text(size=10), plot.margin=unit(c(.05,.1,.05,.1), "cm")) + labs(x = TeX(label[i]),y = 'Count')
}
plot[[3]] = plot[[3]] + theme(legend.position="none")
plot[[4]] = plot[[4]] + theme(legend.position="none")

grid.arrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2, ncol=2)
```

```{r}
set_n_data = c(10, 100, 300)
set_n_rep = c(5, 20, 50)

df_ratio = data.frame()
df_p1 = data.frame()
df_p2 = data.frame()
df_p3 = data.frame()
df_p4 = data.frame()
for (n_data in set_n_data){
    ratio = c()
    rel_err_ratio_p1 = c()
    rel_err_ratio_p2 = c()
    rel_err_ratio_p3 = c()
    rel_err_ratio_p4 = c()
    for (n_rep in set_n_rep){
        df_mcmc = mcmc_generator_2(initial_params, model, true_params, std_prior, std_prop, 
                      n_walkers, n_iter, n_rep, n_data, x_min, x_max, burn_in, thin)
        esti_params_mcmc = colMeans(df_mcmc)
        sd_mcmc = apply(df_mcmc, 2, sd)
        rel_err_mcmc = sd_mcmc/esti_params_mcmc
        
        df_lm = lm_generator(model, true_params, x_min, x_max, n_data, n_rep, FALSE)
        esti_params_lm = colMeans(df_lm)
        sd_lm = apply(df_lm, 2, sd)
        rel_err_lm = sd_lm/esti_params_lm
        
        rel_err_ratio = rel_err_mcmc/rel_err_lm
        rel_err_ratio_p1 = c(rel_err_ratio_p1, rel_err_ratio[1])
        rel_err_ratio_p2 = c(rel_err_ratio_p2, rel_err_ratio[2])
        rel_err_ratio_p3 = c(rel_err_ratio_p3, rel_err_ratio[3])
        rel_err_ratio_p4 = c(rel_err_ratio_p4, rel_err_ratio[4])
        
        x = runif(n_data, x_min, x_max)
        z = numeric(n_data)
        y = model(true_params, x, z, noise = TRUE)
        
        y_pred_mcmc = model(esti_params_mcmc, x, z, noise = FALSE)
        rmse_mcmc = sqrt(mean((y - y_pred_mcmc)^2))
        
        y_pred_lm = model(esti_params_lm, x, z, noise = FALSE)
        rmse_lm = sqrt(mean((y - y_pred_lm)^2))
        
        ratio = c(ratio, rmse_mcmc/rmse_lm)
    }
    df_ratio = rbind(df_ratio, ratio)
    df_p1 = rbind(df_p1, rel_err_ratio_p1)
    df_p2 = rbind(df_p2, rel_err_ratio_p2)
    df_p3 = rbind(df_p3, rel_err_ratio_p3)
    df_p4 = rbind(df_p4, rel_err_ratio_p4)
}
```

```{r}
df_ratio
```

```{r}
df_p1
```

```{r}
df_p2
```

```{r}
df_p3
```

```{r}
df_p4
```

```{r}
trans_df = function(df, set_n_data, set_n_rep){
    dff = data.frame()
    for (i in 1:length(set_n_data)){
        for (j in 1:length(set_n_rep)){
            dff = rbind(dff, c(set_n_data[i], set_n_rep[j], df[i,j]))
        }
    }
    colnames(dff) = c("n_data", "n_rep", "target")
    return(dff)
}
```

```{r}
dff_ratio = trans_df(df_ratio, set_n_data, set_n_rep)
dff_ratio

plot = ggplot() + 
    stat_summary_2d(data = dff_ratio, aes(x = n_data, y = n_rep, z = target), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'RMSE ratio') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
dff = trans_df(df_p1, set_n_data, set_n_rep)
dff

plot = ggplot() + 
    stat_summary_2d(data = dff, aes(x = n_data, y = n_rep, z = target), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'rel. error ratio of a') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
dff = trans_df(df_p2, set_n_data, set_n_rep)
dff

plot = ggplot() + 
    stat_summary_2d(data = dff, aes(x = n_data, y = n_rep, z = target), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'rel. error ratio of b') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
dff = trans_df(df_p3, set_n_data, set_n_rep)
dff

plot = ggplot() + 
    stat_summary_2d(data = dff, aes(x = n_data, y = n_rep, z = target), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'rel. error ratio of c') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
dff = trans_df(df_p4, set_n_data, set_n_rep)
dff

plot = ggplot() + 
    stat_summary_2d(data = dff, aes(x = n_data, y = n_rep, z = target), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'rel. error ratio of standard deviation') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
x = c(10,10,10,100,100,100,300,300,300)
y = rep(c(5, 20, 50), 3)

ratio = c(0.9129818,1.0134971,1.1132446,0.9936071,1.0018118,1.0031808,0.9964045,
          1.0002807,0.9997543)
p1 = c(0.014459426,0.005523706,0.006222232,0.034223981,0.041654588,0.049353021,0.099955432,
        0.111116587,0.124064510)
p2 = c(0.1568092,0.3734230,0.2823070,0.3470520,0.7806239,1.0716771,1.5104333,
       0.7362043,0.7275589)
p3 = c(0.4339801,0.5995456,0.4320102,0.3103714,0.4913456,0.6817522,0.5188061,
       0.5747249,0.7498710)
p4 = c(0.07620723,0.10932993,0.13258436,0.59916250,0.59957432,0.69163076,
       1.04474336,1.27549092,0.87032248)

df = data.frame(x=x,y=y,ratio=ratio,p1=p1,p2=p2,p3=p3,p4=p4)
df
```

```{r}
df[,5]
```

```{r}
plot = ggplot() + 
    stat_summary_2d(data = df, aes(x = x, y = y, z = ratio), bins=9) + 
    scale_fill_gradient(low = 'blue', high = 'red') + 
    labs(x = 'n_data', y = 'n_rep', fill = 'RMSE ratio') + 
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
          plot.title = element_text(size=26), legend.text = element_text(size=18), 
          legend.title = element_text(size=20))
plot
```

```{r}
plot = vector("list", length = length(true_params))

plot1 = ggplot() + stat_summary_2d(data = df, aes(x = x, y = y, z = p1), bins=9) + 
                    scale_fill_gradient(low = 'blue', high = 'red') + 
                    labs(x = 'n_data', y = 'n_rep', fill = TeX(r'{$a$}')) + 
                    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
                          plot.title = element_text(size=26), legend.text = element_text(size=18), 
                          legend.title = element_text(size=20))

plot2 = ggplot() + stat_summary_2d(data = df, aes(x = x, y = y, z = p2), bins=9) + 
                    scale_fill_gradient(low = 'blue', high = 'red') + 
                    labs(x = 'n_data', y = 'n_rep', fill = TeX(r'{$b$}')) + 
                    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
                          plot.title = element_text(size=26), legend.text = element_text(size=18), 
                          legend.title = element_text(size=20))

plot3 = ggplot() + stat_summary_2d(data = df, aes(x = x, y = y, z = p3), bins=9) + 
                    scale_fill_gradient(low = 'blue', high = 'red') + 
                    labs(x = 'n_data', y = 'n_rep', fill = TeX(r'{$c$}')) + 
                    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
                          plot.title = element_text(size=26), legend.text = element_text(size=18), 
                          legend.title = element_text(size=20))

plot4 = ggplot() + stat_summary_2d(data = df, aes(x = x, y = y, z = p4), bins=9) + 
                    scale_fill_gradient(low = 'blue', high = 'red') + 
                    labs(x = 'n_data', y = 'n_rep', fill = TeX(r'{$\sigma$}')) + 
                    theme(axis.text = element_text(size=20), axis.title = element_text(size=20), 
                          plot.title = element_text(size=26), legend.text = element_text(size=18), 
                          legend.title = element_text(size=20))
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol=2)
```

