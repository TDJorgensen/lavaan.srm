// appears in nested for-loops:

// for (k in 1:Kd2)
//   for (d in 1:Nd)
Mu_d[d, idx1] += (S_g[k] * GG[ IDg[d], k]); // add group effects to expected values
Mu_d[d, idx2] += (S_g[k] * GG[ IDg[d], k]);
