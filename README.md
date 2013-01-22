speedyGApy
==========

A numpy + matplotlib based port of [speedyGA](http://www.mathworks.com/matlabcentral/fileexchange/15164) to Python. 

SpeedyGApy is a fast, extensible, single-file, barebones, vectorized genetic algorithm with uniform crossover (UX), sigma scaling, 
and stochastic universal sampling. 

By default speedyGApy showcases a computational efficiency of the genetic algorithm discussed 
[here](http://blog.hackingevolution.net/2013/01/20/foga-2013-slides/). It efficiently computes the 
effective attributes of a 4-bit stochastic effective attribute parity problem. Varying the length of the chromosomes 
does not effect the number of fitness evaluations required for the red dots 
(marking the locations of the 4 effective attributes) to diverge.

SpeedyGApy depends on numpy and matplotlib, which can typically be installed by executing

    easy_install numpy 
    easy_install matplotlib

Enjoy!
