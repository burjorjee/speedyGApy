speedyGApy
==========

A numpy + matplotlib based port of [speedyGA](http://www.mathworks.com/matlabcentral/fileexchange/15164) to Python. 

SpeedyGApy is a fast, extensible, single-file, barebones, vectorized genetic algorithm that uses uniform crossover (UGA), sigma scaling, and stochastic universal sampling. 

+ The function `seapEvolve()` uses the advanced visualization capabilities of matplotlib to showcase a computational efficiency of the genetic algorithm discussed [here](http://blog.hackingevolution.net/2013/01/20/foga-2013-slides/). That is `seapEvolve()` empirically shows how a UGA can efficiently compute the effective attributes of a 4-bit stochastic effective attribute parity problem. Check for yourself that varying the length of the chromosomes does not affect the expected number of fitness evaluations required for the red dots, marking the locations of the 4 effective attributes, to diverge. Changing the effective attributes, i.e. varying the location of the red dots, also has no effect on the expected time-to-divergence. 

+ The function `staircaseFunctionEvolve()` provides proof of concept that a UGA is capable of using the computational power showcased by `seapEvolve()` to implement a global optimization heuristic that I call hyperclimbing

SpeedyGApy depends on numpy and matplotlib, which can typically be installed by executing

    easy_install numpy 
    easy_install matplotlib

On POSIX systems, you may be required to run `easy_install` as a superuser

The functions `seapEvolve()` and `staircaseFunctionEvolve()` will be called in sequence when speedyGA.py is run from the command line

    python speedyGA.py
    
Enjoy!
