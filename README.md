# LMM-PROBE
This is a package to perform sparse high-dimensional linear mixed modeling based on a partitioned empirical Bayes ECM algorithm. 
Please refer to the package manual for more details on the lmmprobe function.  

To install the package, please follow the code snippet below: 
```
library(devtools)
install_github("anjazgodic/lmmprobe")
```

Here is an example for conducting analysis using LMM-PROBE: 
```
library(lmmprobe)
data(SLE)
ep <- 0.05
alpha <- 0.05
Y = SLE$Y
Z = SLE$Z
V = SLE$V
full_res <- lmmprobe(Y = Y, Z = Z, V = V, ep = ep, alpha = alpha)
```
