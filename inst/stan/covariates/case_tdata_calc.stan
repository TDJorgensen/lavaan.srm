  // calculate observed limits for case-level observations
  for (k in 1:Kp) {
    limKp[k,1] = min(Yp[ : , k] );
    limKp[k,2] = max(Yp[ : , k] );
    halfRangeKp[k] = (limKp[k,2] - limKp[k,1]) / 2;
  }
