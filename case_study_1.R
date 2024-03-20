# This file has the implementation of our proposed method for case study 1

# file configuration 
root_path = getwd()
print(root_path)

# ------------------------- data collection and processing ----------------------------------
library(fields)
raw_data = read.csv("data.csv",header=T)
str(raw_data)

x_original = raw_data$CTOD
y_original = raw_data$Stress
data = data.frame(x_original, y_original)
colnames(data) = c("x","y")
data = na.omit(data)
str(data)

# adjust the data so that it converges to zero as CTOD increases
data$y = data$y - 0.05
plot(data$x, data$y, ylim=c(0,1), col='blue')
abline(h=0)

right_data = data

rem_id = which(right_data$x < 0)
right_data = right_data[-rem_id,]

y_max_minus_min = max(right_data$y) - min(right_data$y)
y_min = min(right_data$y)
right_data$y = (right_data$y - y_min)/y_max_minus_min

data = right_data # ensure that the two data sets are exactly the same

seed_current = 1

# ------------------------- train/test split ----------------------------------
# here we perform disproportionate stratified sampling (DSS) to obtain the testing data

df_test = data.frame(matrix(ncol = 2, nrow = 0))
colnames(df_test) = c("x","y")

edges = seq(min(data$x), max(data$x), length.out=201) # 200 strata

list_rows_sampled = c()
xs_mid_strata = c()
variances_strata = c()

for (i in 1:(length(edges)-1)) {
  
  index_this_strata = which(data$x > edges[i] & data$x < edges[i+1])
  data_strata = data[index_this_strata,]
  
  # this is to plot the variances across different strata
  xs_mid_strata = append(xs_mid_strata, (edges[i] + edges[i+1])/2)
  if (nrow(data_strata) > 1) {
    variances_strata = append(variances_strata, var(data_strata$y))
  } else {
    variances_strata = append(variances_strata, 0)
  }
  
  # this is to create the testing data set
  if (length(index_this_strata) > 0) {
    
    set.seed(seed_current * 20)
    rows_sampled = sample(1:nrow(data_strata), size=1, replace=F) # index relative to data within strata
    list_rows_sampled = append(list_rows_sampled, index_this_strata[rows_sampled]) # index relative to all the data
    
    df_test[nrow(df_test) + 1,] = c(data_strata[rows_sampled,1], data_strata[rows_sampled,2])
  }
}

plot(xs_mid_strata, variances_strata, type='l', col='red', xlim=c(0,20))

test_data = df_test
train_data = data[-list_rows_sampled,]
str(train_data)
str(test_data)

plot(train_data$x, train_data$y, pch=20, col='grey')
points(test_data$x, test_data$y, pch=20, col='blue')

# ------------------------- remove the replications ----------------------------------
# if a CTOD value has replications, we take the average of the replications

x_unique = train_data$x[duplicated(train_data$x) == FALSE]

indexes_this_x = replacements = list(length = length(x_unique))
indexes_replications = 0

for (i in 1:length(x_unique)) {
  
  indexes_this_x[[i]] = which(train_data$x == x_unique[i])
  replacements[[i]] = rep(mean(train_data$y[indexes_this_x[[i]]]), length(indexes_this_x[[i]]))
  
  if (length(indexes_this_x[[i]]) >= 2) { # if there are replications
    indexes_replications = c(indexes_replications, indexes_this_x[[i]][2:length(indexes_this_x[[i]])])
  }
}
indexes_replications = indexes_replications[-1]

train_data_new = train_data
train_data_new$y = as.vector(unlist(replacements))
train_data_new = train_data_new[-indexes_replications,]

# ------------------------- physics mean function ----------------------------------
mf = function(x_in, sigma_0, delta_0, delta_c, alpha) {
  
  # This define the physics mean function that takes in CTOD values and returns traction values
  # inputs:
  #   - x_in: all the CTOD values
  #   - sigma_0, delta_0, delta_c, alpha: function parameters
  # output:
  #   - y_out: all the corresponding traction values
  
  y_out = vector(length = length(x_in))
  id_small = which(x_in <= delta_0)
  
  if(length(id_small) != 0){
    
    y_out[id_small] = x_in[id_small] / delta_0
    
    numero = 1 - exp(-alpha * (x_in[-id_small] - delta_0)/(delta_c - delta_0))
    denom = 1 - exp(-alpha)
    y_out[-id_small] = pmax(1-(numero/denom), 0)
    
  } else {
    
    numero = 1 - exp(-alpha * (x_in - delta_0)/(delta_c - delta_0))
    denom = 1 - exp(-alpha)
    y_out = pmax(1-(numero/denom), 0)
  }
  
  y_out = sigma_0 * y_out
  return(y_out)
}

# ------------------------- optimization for parameters ----------------------------------

#install.packages("Rfast")
library(Rfast)
matrix_dist = rdist(train_data_new$x) # the absolute difference between a pair of CTOD values

mle = function(pars) {
  
  # This is a parametric objective function to be minimized later
  # Given both mean and GP parameters, 
  # this function calculates the likelihood of the observed data,
  # more specifically, it calculates the value that would be equivalent to 
  # negative log-likelihood in terms of optimization
  # inputs:
  #   - pars[1]: signal variance
  #   - pars[2]: length scale
  #   - pars[3]: noise variance
  #   - pars[5]-pars[8]: mean function parameters
  # outputs:
  #   - NLL: a value equivalent (in optimization) to negative log-likelihood 
  
  # the covariance part
  pars=exp(pars)
  a = pars[1]
  b=pars[2]
  c=pars[3]
  
  S = a * exp(-0.5*(matrix_dist/b)^2) # kernel matrix for training data
  diag(S) = diag(S)  + c + 10^-10
  
  # the mean part
  index_base = 5
  
  # pars[index_base+2] is delta_c - delta_0
  pars[index_base+2] = pars[index_base+1] + pars[index_base+2]
  mf_fit = mf(train_data_new$x,
              pars[index_base],
              pars[index_base+1],
              pars[index_base+2],
              pars[index_base+3])
  
  r = train_data_new$y - mf_fit
  
  # calculate the value equivalent to negative log-likelihood
  mu = 0
  U = cholesky(S, parallel=T) # S = U^TU
  
  term1 = 2 * sum(log(abs(diag(U))))
  
  # solve the linear system of equation A * x = b
  # A is triangular so the computation can be efficient 
  temp = forwardsolve(t(U), r-mu) # for lower triangular matrix
  alpha = backsolve(U,temp) # for upper triangular matrix
  term2 = t(r-mu) %*% alpha
  
  NLL = (term1 + term2) / 2
  
  return(NLL)
}

start_pt = c(-3.7,0.43,-0.05,1,log(c(1,1,30,10)))
optimizing = nlm(mle,start_pt,stepmax=3,steptol=10^-10,iterlim=300)
optimizing

# the covariance part
a = exp(optimizing$estimate[1])
b = exp(optimizing$estimate[2])
c = exp(optimizing$estimate[3])

K = a * exp(-0.5*(matrix_dist/b)^2)
diag(K)= diag(K) + c + 10^-10
K[1:5,1:5]
invcov = solve(K)

# the mean part
index_base = 5

# print out the optimized parameters for mean function
exp(optimizing$estimate[index_base])
exp(optimizing$estimate[index_base+1])
exp(optimizing$estimate[index_base+2]) + exp(optimizing$estimate[index_base+1])
exp(optimizing$estimate[index_base+3])

mf_fit = mf(train_data_new$x,
            exp(optimizing$estimate[index_base]),
            exp(optimizing$estimate[index_base+1]),
            exp(optimizing$estimate[index_base+2]) + exp(optimizing$estimate[index_base+1]),
            exp(optimizing$estimate[index_base+3]))
r = train_data_new$y - mf_fit

# ------------------------- obtaining posterior using Bayes rule ----------------------------------
# now that we have the optimized parameters, the prior GP has been determined
# we can now simply apply Bayes rule to obtain the posterior

krigging = function (x_in, y_in, test_mat, Kinv) {
  
  # This function produces the predictive mean for any given CTOD value
  # inputs:
  #   - x_in: the CTOD values used in training as a column matrix
  #   - y_in: the traction values used in training as a column matrix
  #   - test_mat: the testing data as a matrix, but we will only use its CTOD column (the x column)
  #   - Kinv: inverse of kernel matrix for training data
  # outputs:
  #   - out: a matrix whose first column is the predictive means for the testing CTOD values
  
  out = matrix(nrow=dim(test_mat)[1], ncol=2)
  out_mean = out_var = vector(length=dim(test_mat)[1])
  
  for (i in 1:dim(test_mat)[1]) { 
    
    dist_vec = abs(x_in[,1] - test_mat[i,1])
    k_test = a * exp(-0.5*(dist_vec/b)^2) 

    term1 = Kinv %*% as.matrix(k_test)
    out_mean[i] = t(term1) %*% y_in
    
    out[i,] = c(out_mean[i],out_var[i])
  }
  
  return(out)
}

gp_out = krigging(x_in=as.matrix(train_data_new$x),
                  y_in=as.matrix(r),
                  test_mat=as.matrix(test_data),
                  Kinv=invcov)

gp_train_out = krigging(x_in=as.matrix(train_data_new$x),
                        y_in=as.matrix(r),
                        test_mat=as.matrix(train_data_new),
                        Kinv=invcov)

gp_preds = gp_out[,1]
gp_fits = gp_train_out[,1]

preds = mf(test_data$x,
           exp(optimizing$estimate[index_base]),
           exp(optimizing$estimate[index_base+1]),
           exp(optimizing$estimate[index_base+2]) + exp(optimizing$estimate[index_base+1]),
           exp(optimizing$estimate[index_base+3])) + gp_preds

fits = mf(train_data_new$x,
          exp(optimizing$estimate[index_base]),
          exp(optimizing$estimate[index_base+1]),
          exp(optimizing$estimate[index_base+2]) + exp(optimizing$estimate[index_base+1]),
          exp(optimizing$estimate[index_base+3])) + gp_fits

# convert it back to the original scale
preds = preds * y_max_minus_min + y_min
fits = fits * y_max_minus_min + y_min

# ------------------------- plotting and error calculation ----------------------------------

# plot the original training data as our background
x_plot = train_data$x
y_plot = train_data$y * y_max_minus_min + y_min
plot(x_plot, y_plot, panel.first=grid(), col='grey')

x_plot_new = train_data_new$x
y_plot_new = train_data_new$y * y_max_minus_min + y_min

# plot the mean function fitted on the actual training data
mf_fit_plot = mf(x_plot_new,
                exp(optimizing$estimate[index_base]),
                exp(optimizing$estimate[index_base+1]),
                exp(optimizing$estimate[index_base+2]) + exp(optimizing$estimate[index_base+1]),
                exp(optimizing$estimate[index_base+3]))
mf_fit_plot = mf_fit_plot * y_max_minus_min + y_min
lines(x_plot_new,mf_fit_plot,col="blue")

# handle negative predictions
preds_real = preds
fits_real = fits

# plot the fit on the actual training data
lines(x_plot_new,fits_real,col="red")

# error calculation
test_y_used = test_data$y
train_y_used = train_data$y # not used
test_y_used = test_data$y * y_max_minus_min + y_min
train_y_used = train_data$y * y_max_minus_min + y_min # not used

RMSD = sqrt(mean((test_y_used - preds_real)^2))
NRMSD = 100*RMSD/(max(test_y_used)-min(test_y_used))

print(paste('current RMSD:', RMSD))
print(paste('current NRMSD:', NRMSD))

# plot the prediction on the testing data
plot(test_data$x, test_y_used, panel.first=grid(), col='grey')
lines(test_data$x, preds_real, col="red")

