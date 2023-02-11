// priors, in a loop over RR variables (k in 1:Kd2)
Mvec[k] ~ normal(limKd2[k, 1] + halfRangeKd2[k], halfRangeKd2[k]);
