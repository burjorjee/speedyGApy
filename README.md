speedyGApy
==========

A numpy + matplotlib based port of [speedyGA](http://www.mathworks.com/matlabcentral/fileexchange/15164) to Python. 

SpeedyGApy is a fast, extensible, single-file, barebones, vectorized genetic algorithm with uniform crossover, sigma scaling, and stochastic universal sampling. 

The script can be used to reproduce the experiments that form the basis for the [Hyperclimbing Hypothesis](http://s3.amazonaws.com/burjorjee/www/hyperclimbing_hypothesis_2013.pdf)---a new explanation for adaptation in genetic algorithms with uniform crossover.

+ The function `seapEvolve()` uses the the 4-Bit Stochastic Effective Attribute Parity problem discussed [here](http://blog.hackingevolution.net/2013/01/20/foga-2013-slides/) and [here](http://blog.hackingevolution.net/2009/06/29/red-dots-blue-dots/), and the advanced visualization capabilities of matplotlib to showcase a computational efficiency of the genetic algorithm. Check for yourself that varying the length of the chromosomes does not affect the expected number of fitness evaluations required for the red dots, marking the locations of the 4 effective attributes, to diverge. Changing the effective attributes, i.e. varying the location of the red dots, also has no effect on the expected time-to-divergence. 

+ The function `staircaseFunctionEvolve()` provides proof of concept that a UGA is capable of using the computational power showcased by `seapEvolve()` to implement a global optimization heuristic called [hyperclimbing](http://s3.amazonaws.com/burjorjee/www/hyperclimbing_hypothesis_2013.pdf)

SpeedyGApy depends on the Python packages numpy and matplotlib, which can typically be installed by executing 

    easy_install numpy 
    easy_install matplotlib

On POSIX systems, you may be required to run `easy_install` as a superuser


SpeedyGA Usage Instructions:

    usage: python speedyGA.py [-h] [--fitnessFunction {staircase,seap}]
                       [--probCrossover PROBCROSSOVER]
                       [--probMutation PROBMUTATION] [--popSize POPSIZE]
                       [--bitstringLength BITSTRINGLENGTH] [--gens GENS]
    
    Run SpeedyGA on `seap` or `staircase`, two fitness functions tailored to
    provide proof of concept for the Hyperclimbing Hypothesis. More details at
    http://blog.hackingevolution.net/2013/01/20/foga-2013-slides/
    
    optional arguments:
      -h, --help            show this help message and exit
      --fitnessFunction {staircase,seap}
                            The fitness function to use (default: staircase).
      --probCrossover PROBCROSSOVER
                            Number between 0 and 1 representing the fraction of
                            the population subject to crossover (default:1)
      --probMutation PROBMUTATION
                            The per bit mutation probability (default:0.005)
      --popSize POPSIZE     Size of the population (default:500)
      --bitstringLength BITSTRINGLENGTH
                            Length of a chromosome in the population (default:500)
      --gens GENS           The number of generations (default:500)
      
Enjoy!
