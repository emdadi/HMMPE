# HMMPE
HMM Parameter Estimation
Implementations of AntMarkov algorithm for parameter estimation hidden Markov model.

The implemenation of AntMarkov is presented in "mainAntMakov.m" file. Which runs the AntMarkov alorithm on a cell of sequences and returns the estimated matrices of emissions and transitions and the number of iterations.
In this package we gather the implementation of other methods that are developed by [1], [2]. Hsu method is proposed in [3] and its implementation is also available at [4] 
The implemenation of tabu search based algorithm is done as proposed by [5]. 
Files "learnHMM.m", "learnSNNMF.m" and "tabu.m" are the implementations of Hsu. et. al method, SNNMMf and tabu-based method.
The file "mainProgNormal.m" includes the code for executing all algorithms simultanously on a set of sequences and aggregate the results of all algorithms in output files. 
It calls the mentiond codes and also calls off-th-shelf matlab codes such as HMMviterbi and HMMtrain.
Thus, one can implement methods seperately by running their codes or simultanously by running "mainProgNormal.m".




[1] R. Mattila, \On identication of hidden markov models using spectral and non-negative matrix factorization methods," 2015

[2] C. Mattfeld, \Implementing spectral methods for hidden markov models with real-valued emissions," arXiv preprint arXiv:1404.7472, 2014.

[3] Hsu, Daniel, Sham M. Kakade, and Tong Zhang. "A spectral algorithm for learning hidden Markov models."  Journal of Computer and System Sciences 78.5 (2012): 1460-1480.

[4] https://github.com/cmgithub/spectral.

[5] T.-Y. Chen, X.-D. Mei, J.-s. Pan, and S.-H. Sun, \Optimization of hmm by the tabu search algorithm," J. Inf. Sci. Eng., vol. 20, no. 5, pp. 949{957, 2004.
