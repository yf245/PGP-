# This file has the implementation of our proposed method for case study 2

ROOT_PATH = getwd()
print(ROOT_PATH) # view the current directory

library(Rfast)
library(readxl)

# ---------------------------------- importing data ------------------------------------
rawdata = read_excel("nucleation_data.xlsx")
data = as.data.frame(rawdata)
str(data) # just to view
colnames(data) = c('x1', 'x2', strsplit(toString(3:12),", ")[[1]])
str(data) # just to view

collect_all_x1 = c()
collect_all_x2 = c()
collect_all_y = c()
for (i in 3:12) {
  collect_all_x1 = append(collect_all_x1, data$x1)
  collect_all_x2 = append(collect_all_x2, data$x2)
  collect_all_y = append(collect_all_y, data[[i]])
}

actual_data = data.frame(collect_all_x1, collect_all_x2, collect_all_y)
str(actual_data) # just to view
actual_data = na.omit(actual_data)
str(actual_data) # just to view

# data update
data = actual_data
colnames(data) = c('x1','x2','y')
str(data) # just to view

#plot the original data
vector_considered = data$x1
unique(vector_considered) # only 10 so 10 colors are needed
colors = c('hotpink','green','blue','orange','red','cyan','purple','pink','yellow','brown')

vector_for_color = c()
for (i in vector_considered) {
  for (j in 1:length(unique(vector_considered))) {
    if (i == unique(vector_considered)[j]) {
      vector_for_color = append(vector_for_color,colors[j])
    }
  }
}
plot(data$x2,data$y,col=vector_for_color,xlim=c(),ylim=c(),pch=20)

# ---------------------------------- train-test split ----------------------------------

all_10_temperatures = sort(unique(data$x1))
testdata = data[data$x1 == all_10_temperatures[6],]
traindata = data[data$x1 %in% all_10_temperatures[c(5,7)],]

str(traindata)
str(testdata)

testdata_NRMSD = testdata # this is just in case that testdata gets modified later
traindata_NRMSD = traindata # this is just in case that traindata gets modified later

# initialize newtraindata, which will be modified later
newtraindata = traindata
  
# ----------------------------- first obtain the fitted mean function ----------------------
  
mf = function(x_in,t_0,E_0,sigma_c_part,mean_m,mean_n){
  
  # Physics-informed mean function
  # Input:
  #  - the entire data (albeit we are only interested in x)
  #  - the parameters (keep in mind that the fourth argument here represents max(x_in$x2) - sigma_c)
  # Output: 
  #  - the corresponding output of the function

  sigma_c = sigma_c_part + max(x_in$x2)
  
  zhat_term1 = 1 - (x_in$x2 / sigma_c)^mean_m
  zhat_term2 = E_0 * zhat_term1^mean_n
  kB = 8.617 * 10^(-5)
  y_out = t_0 * exp(zhat_term2 / (kB * x_in$x1))

  return(y_out)
}

start.pt = c(0.1,1.37,0.750049,22026.47,3.720076e-10) # choose an appropriate starting point based on domain knowledge

ls_fit = function(pars,x_in,y_in){
  
  # Use least squared error (LSE) to fit the parameters of the mean function
  # Input:
  #  - pars: the current parameters
  #  - x_in: the training data, but we are only interested in the x
  #  - y_in: the y of the training data
  # Output: 
  #  - the resulting mean squared error (MSE)
  
  pars = exp(pars) # all parameters must be positive (domain knowledge)

  pars[4] = (1 - 0) * pars[4] / (pars[4] + 1) + 0 
  pars[5] = (2 - 1) * pars[5] / (pars[5] + 1) + 1
  
  y_out = mf(x_in,pars[1],pars[2],pars[3],pars[4],pars[5])
  obj = mean((y_in - y_out)^2)
  return(obj)
}


ls_fit(pars=log(start.pt),newtraindata,newtraindata$y) # just checking

# now perform the fitting using LSE
fopt = nlm(ls_fit,p=log(start.pt),x_in=newtraindata,y_in=newtraindata$y,stepmax=1,steptol = 10^-10,iterlim=400)

save.start.pt = start.pt
fopt # check convergence, code 1 and code 2 are more acceptable
opt_pars = exp(fopt$estimate)
opt_pars_3 = opt_pars[3]
opt_pars[4] = (1 - 0) * opt_pars[4] / (opt_pars[4] + 1) + 0
opt_pars[5] = (2 - 1) * opt_pars[5] / (opt_pars[5] + 1) + 1
opt_pars

# show the fitted mean function parameters
opt_pars[1]
opt_pars[2]
opt_pars_3
opt_pars_3 + max(newtraindata$x2) # this is the actual sigma_c
opt_pars[4]
opt_pars[5]

# now finally obtain the fitted mean function
mf_fit = mf(newtraindata,opt_pars[1],opt_pars[2],opt_pars_3,opt_pars[4],opt_pars[5])

# we can visualize the fitted mean function
vector_considered = newtraindata$x1
unique(vector_considered)
colors = c('hotpink','green','blue','orange','red','cyan','purple','pink','yellow','brown')

vector_for_color = c()
for (i in vector_considered) {
  for (j in 1:length(unique(vector_considered))) {
    if (i == unique(vector_considered)[j]) {
      vector_for_color = append(vector_for_color,colors[j])
    }
  }
}

plot(newtraindata$x2,newtraindata$y,col=vector_for_color,pch=20)
points(newtraindata$x2,mf_fit,col=vector_for_color,pch=8)

# ----------------------------- obtain the residual --------------------------------
# obtain the residual from the fitted mean function
r = newtraindata$y - mf_fit

# keep in mind that the three variables below will not be modified later
m_r = min(r) 
s_r = max(r) - min(r) 
r = (r - m_r)/s_r 

# we can visualize the residual
hist(r)
plot(traindata$x1,r)
plot(traindata$x2,r)
length(r)

# Note: we now have traindata, newtraindata (which is the same as traindata), and testdata. 

# keep in mind that the four variables below will not be modified later
max_minus_min_x1_new = max(newtraindata$x1) - min(newtraindata$x1)
max_minus_min_x2_new = max(newtraindata$x2) - min(newtraindata$x2)
min_x1_new = min(newtraindata$x1)
min_x2_new = min(newtraindata$x2)

newtraindata$x1 = (newtraindata$x1 - min_x1_new)/max_minus_min_x1_new
newtraindata$x2 = (newtraindata$x2 - min_x2_new)/max_minus_min_x2_new

traindata = newtraindata

testdata$x1 = (testdata$x1 - min_x1_new)/max_minus_min_x1_new
testdata$x2 = (testdata$x2 - min_x2_new)/max_minus_min_x2_new

# Reminder: we now have normalized traindata, newtraindata, and testdata, traindata and newtraindata are the same

# ------------------------------------ update newtraindata ----------------------------
# so that it represents the average residual instead

chosen_x1 = newtraindata$x1
chosen_x2 = newtraindata$x2
chosen_y = newtraindata$y

unique_x_df = unique(newtraindata[,1:2])
unique_x1_vector = unique_x_df$x1
unique_x2_vector = unique_x_df$x2

y_avg = c()
r_avg = c() # keep in mind that it will not be modified later
for (i in 1:length(unique_x1_vector)) {
  
  matching_elements = chosen_y[chosen_x1 == unique_x1_vector[i] & chosen_x2 == unique_x2_vector[i]]
  matching_elements2 = r[chosen_x1 == unique_x1_vector[i] & chosen_x2 == unique_x2_vector[i]]
  
  y_avg = append(y_avg, mean(matching_elements))
  r_avg = append(r_avg, mean(matching_elements2))
}

newtraindata = data.frame (x1 = unique_x1_vector, x2=unique_x2_vector, y = r_avg)

# -------------------------- GP: optimization for the covariance function parameters ------------------
dst_euclidean = function(a, b) sqrt(sum((a - b)^2)) # calculate the distance between any two vectors

# obtain the distance matrix, which will be needed for the covariance matrix later
mdist = matrix(nrow=nrow(newtraindata),ncol=nrow(newtraindata))
for (i in 1:nrow(newtraindata)) {
  for (j in 1:nrow(newtraindata)) {
    mdist[i,j] = dst_euclidean(newtraindata[i,1:2], newtraindata[j,1:2])
  }
}

mdist[1:5,1:5]
b_upper_bound = max(mdist)
b_lower_bound = max(10^(-3),min(mdist))

mle1=function(pars){
  
  # Calculate the log-likelihood (which is the optimization objective) 
  # Note: the parameters of the mean function have already been fitted
  # Input: 
  #  - the covariance function parameters
  # Output: 
  #  - negative of the log-likelihood (since we are doing minimization)
  
  pars=exp(pars)
  a=pars[1]
  b=(b_upper_bound - b_lower_bound) * pars[2] / (pars[2]+1) + b_lower_bound
  c=pars[3]
  
  # the corresponding covariance matrix
  S = a * exp(-0.5*(mdist/b)^2)
  diag(S) = diag(S)  + c + 10^-10
  
  mu = 0
  temp3 <- cholesky(S,parallel=T)
  temp4<-2*sum(log(abs(diag(temp3))))
  temp5<-forwardsolve(t(temp3),r_avg-mu) # the residual is being used here
  temp6<-t(r_avg-mu)%*%backsolve(temp3,temp5) # the residual is being used here
  temp<-(temp4+temp6)/2
  #cat(-temp,"\n")
  return(temp)
}

start.pt = c(0,0,0,0)
optimizing=nlm(mle1,start.pt,stepmax=3,steptol=10^-10,iterlim=300)
optimizing # check the convergence, code 1 and code 2 are more acceptable
log_likelihood = -optimizing$minimum

# obtain the fitted covariance function parameters (these will not be modified later)
a = exp(optimizing$estimate[1])
b=(b_upper_bound - b_lower_bound) * exp(optimizing$estimate[2]) / (exp(optimizing$estimate[2])+1) + b_lower_bound
c = exp(optimizing$estimate[3])

# obtain the final fitted covariance matrix
K = a * exp(-0.5*(mdist/b)^2)
diag(K)= diag(K) + c + 10^-10

K[1:5,1:5]
invcov = solve(K)

# ----------------------- GP: apply Bayes' rule to obtain the final predictive distribution ----------------

krig1D = function(train_mat,test_mat,y_in,cov_mat,invcov_mat){
  
  # Apply the Bayes' rule and predict for the test data
  # Input:
  #  - train_mat: the training X (we will ignore the third column)
  #  - test_mat: the testing X (we will ignore the third column)
  #  - y_in: the training y, which should be the residual
  #  - cov_mat: the fitted covariance matrix
  #  - invcov_mat: the inverse of the fitted covariance matrix
  # Output:
  #  - the predictive distribution

  num_train = dim(train_mat)[1] 
  num_test = dim(test_mat)[1]
  
  out.pred = out.var = vector(length = num_test)
  out = matrix(nrow=num_test,ncol=2)

  K = cov_mat
  Kinv = invcov_mat
  
  dst_euclidean = function(a, b) sqrt(sum((a - b)^2))
  
  for (i in 1:num_test){ # each test point
    test.coords = test_mat[i,1:2]
    
    dist_vec = matrix(nrow=num_train,ncol=1) # it is a column matrix
    for (j in 1:num_train) {
      dist_vec[j] = dst_euclidean(train_mat[j,1:2], test.coords)
    }
    
    k = a * exp(-0.5*(dist_vec/b)^2) # we use the fitted covariance function parameters here

    term1 = Kinv %*% k
    obs = as.matrix(y_in)
    pred = t(term1) %*% obs # this is the predictive mean
    
    # predvar = a - k %*% term1

    out.pred[i] = pred
    #out.var[i] = predvar
    
    out[i,] = c(out.pred[i],out.var[i])
  }
  return(out)
}


gp.out = krig1D(train_mat=as.matrix(newtraindata),
                y_in=r_avg, # need to use the residual
                test_mat=as.matrix(testdata), # this is for the actual test points
                cov_mat=K,invcov_mat=invcov)

gp_train.out = krig1D(train_mat=as.matrix(newtraindata),
                y_in=r_avg, # need to use the residual
                test_mat=as.matrix(traindata),
                cov_mat=K,invcov_mat=invcov)

# ------------------------------ to plot the interpolated training fit (part 1) ---------------------------
# Goal: 
#  * need to obtain a dataframe that is similar to traindata but interpolated
#  * also need to obtain the following three variables, which will be needed in part 2:
#   - interpolated_mf
#   - vector_for_color_interpolated
#   - fits_interpolated

interpolated = traindata[1,1:3] # for initialization, but we will cut it off later

for (i in 1:length(unique(traindata$x1))) { # we want to preserve the original order of temperature and interpolate only for stress
  subset = traindata[traindata$x1 == unique(traindata$x1)[i],1:3]
  subset = subset[order(subset$x2),1:3]
  
  for (j in 1:(nrow(subset)-1)) {
    left_value = subset$x2[j]
    right_value = subset$x2[j+1]
    middle_value = (left_value + right_value) / 2
    
    interpolated = rbind(interpolated,
                         c(unique(traindata$x1)[i], left_value, 0)) # fill the y column with zeros for now
    interpolated = rbind(interpolated,
                         c(unique(traindata$x1)[i], middle_value, 0))
  }
  interpolated = rbind(interpolated,
                       c(unique(traindata$x1)[i], right_value, 0))
  
}
interpolated = interpolated[2:nrow(interpolated),1:3] # cut off the first row

interpolated_mf = interpolated # interpolated but in the original scale
interpolated_mf$x1 = interpolated_mf$x1*max_minus_min_x1_new + min_x1_new
interpolated_mf$x2 = interpolated_mf$x2*max_minus_min_x2_new + min_x2_new

# colors for the interpolated training fit plotting
colors = c('blue','green') # 300, 350

vector_considered = interpolated_mf$x1
unique(vector_considered)

vector_for_color_interpolated = c()
for (i in vector_considered) {
  for (j in 1:length(unique(vector_considered))) {
    if (i == unique(vector_considered)[j]) {
      vector_for_color_interpolated = append(vector_for_color_interpolated,colors[j])
    }
  }
}

# now we obtain the interpolated training fit
gp_interpolated.out = krig1D(train_mat=as.matrix(newtraindata),
                             y_in=r_avg, 
                             test_mat=as.matrix(interpolated), # just change this part
                             cov_mat=K,invcov_mat=invcov) 

gp.interpolated = gp_interpolated.out[,1] # the predictive mean
gp.interpolated = gp.interpolated*s_r + m_r

fits_interpolated = mf(interpolated_mf,opt_pars[1],opt_pars[2],opt_pars_3,opt_pars[4],opt_pars[5]) + gp.interpolated

# -------------------------------- obtain the training fit and the (test) prediction ------------------------
# these are the predictive means
gp.preds = gp.out[,1]
gp.fits = gp_train.out[,1]

gp.preds = gp.preds*s_r + m_r
gp.fits = gp.fits*s_r + m_r

# note that r is no longer needed from now on

traindata_mf = traindata 
testdata_mf = testdata

traindata_mf$x1 = traindata$x1*max_minus_min_x1_new + min_x1_new
traindata_mf$x2 = traindata$x2*max_minus_min_x2_new + min_x2_new
testdata_mf$x1 = testdata$x1*max_minus_min_x1_new + min_x1_new
testdata_mf$x2 = testdata$x2*max_minus_min_x2_new + min_x2_new

after_LS_preds = mf(testdata_mf,opt_pars[1],opt_pars[2],opt_pars_3,opt_pars[4],opt_pars[5])
preds = after_LS_preds + gp.preds

after_LS_fits = mf(traindata_mf,opt_pars[1],opt_pars[2],opt_pars_3,opt_pars[4],opt_pars[5])
fits = after_LS_fits + gp.fits

# --------------------------- now we visualize the fit and the prediction ----------------------------------

# ------------------------ first plot the background (which is the original training data) ----------------------
x_to_plot1 = traindata_mf$x1
x_to_plot2 = traindata_mf$x2
y_to_plot = traindata_mf$y

vector_considered = x_to_plot1
unique(vector_considered) 

colors = c('blue','green') # 300, 350

vector_for_color = c()
for (i in vector_considered) {
  for (j in 1:length(unique(vector_considered))) {
    if (i == unique(vector_considered)[j]) {
      vector_for_color = append(vector_for_color,colors[j])
    }
  }
}

#png(file='2_D.png',width=1500,height=1500,unit='px',res=300)
par(cex.main=0.9)

xlim_range = c(5.8,6.8)

plot(x_to_plot2,y_to_plot,col=vector_for_color,
     xlim=xlim_range,ylim=c(0,600),pch=20,
     xlab=expression(bold(paste("Stress, ", sigma, " (GPa)"))),
     ylab=expression(paste(bold("Nucleation Time, "),bold(italic("t")["nuc"]),bold(" (ps)"))),
     cex.lab=1, cex.axis=1, xaxt='n', yaxt='n',panel.first = grid())
axis(side = 1, lwd = 1.5)
axis(side = 2, lwd = 1.5)

#------------------------ next we plot the main characters (which are the fit and the prediction) ---------- 

################## plot the mean function fit and the training fit (optional)
fits_real = fits

# new_x_to_plot1 = traindata$x1
# new_x_to_plot2 = traindata$x2
# new_y_to_plot = traindata$y
# 
# new_x_to_plot1 = traindata$x1*max_minus_min_x1_new + min_x1_new
# new_x_to_plot2 = traindata$x2*max_minus_min_x2_new + min_x2_new
# 
# vector_considered = new_x_to_plot1
# unique(vector_considered)
# 
# vector_for_color2 = c()
# for (i in vector_considered) {
#   for (j in 1:length(unique(vector_considered))) {
#     if (i == unique(vector_considered)[j]) {
#       vector_for_color2 = append(vector_for_color2,colors[j])
#     }
#   }
# }
# 
# # plot the mean function fit
# points(new_x_to_plot2,after_LS_fits,col=vector_for_color2,pch=8)
# for (j in 1:length(unique(vector_considered))) {
#   temp_indexes = new_x_to_plot1 == unique(vector_considered)[j]
#   lines(new_x_to_plot2[temp_indexes][order(new_x_to_plot2[temp_indexes])],after_LS_fits[temp_indexes][order(new_x_to_plot2[temp_indexes])],col=vector_for_color2[temp_indexes],pch=25)
# }
# 
# # plot the training fit
# points(new_x_to_plot2,fits_real,col=vector_for_color2,pch=25)
# for (j in 1:length(unique(vector_considered))) {
#   temp_indexes = new_x_to_plot1 == unique(vector_considered)[j]
#   lines(new_x_to_plot2[temp_indexes][order(new_x_to_plot2[temp_indexes])],fits_real[temp_indexes][order(new_x_to_plot2[temp_indexes])],col=vector_for_color2[temp_indexes],pch=25)
# }

################## plot the interpolated training fit (part 2)
new_x_to_plot1 = interpolated_mf$x1
new_x_to_plot2 = interpolated_mf$x2

vector_considered = new_x_to_plot1
unique(vector_considered)
# note that we already have vector_for_color_interpolated

train_lty = c(6,5) # 300, 350
train_lty = c(2,2)

fits_interpolated_real = fits_interpolated

#points(new_x_to_plot2,fits_interpolated_real,col=vector_for_color_interpolated,pch=25)
for (j in 1:length(unique(vector_considered))) { # each temperature
  temp_indexes = new_x_to_plot1 == unique(vector_considered)[j]
  lines(new_x_to_plot2[temp_indexes][order(new_x_to_plot2[temp_indexes])],
        fits_interpolated_real[temp_indexes][order(new_x_to_plot2[temp_indexes])],
        col=vector_for_color_interpolated[temp_indexes],lwd=2,lty=train_lty[j])
}
op__ = par(family = "sans")
par(op__)
#dev.off()

################## plot the prediction 
test_y_used = testdata_NRMSD$y # it should be same as testdata$y, but just in case

test_x_to_plot2 = testdata$x2

test_x_to_plot2 = testdata$x2 * max_minus_min_x2_new + min_x2_new

preds_real = preds

test_color = 'red' # 325
test_lty = 1

points(test_x_to_plot2,test_y_used,col=test_color,pch=17)
lines(test_x_to_plot2[order(test_x_to_plot2)],preds_real[order(test_x_to_plot2)],col=test_color,lwd=2,lty=test_lty)
#dev.off()

# --------------------------- finally calculate the errors --------------------
train_y_used = traindata_NRMSD$y # it should be same as traindata$y, but just in case

RMSD = sqrt(mean((test_y_used - preds_real)^2))
NRMSD = 100*RMSD/(max(test_y_used)-min(test_y_used))

print(paste('current RMSD:', RMSD))
print(paste('current NRMSD:', NRMSD))
