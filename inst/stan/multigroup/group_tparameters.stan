// appears in nested for-loops:

// for (k in 1:Kd2)
//   for (d in 1:Nd)
Yd2hat[d, idx1] += (S_g[k] * GG[ IDg[d], k]); // add group effects to expected values
Yd2hat[d, idx2] += (S_g[k] * GG[ IDg[d], k]);
