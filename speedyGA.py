import argparse
from functools import partial
from random import random
from matplotlib.pyplot import figure, plot, xlabel, ylabel, title, axis, hold
from numpy import *
from numpy.random import rand, randn

def visualizeGen(pop, gen, avgFitness, maxFitness, figNum=1):
    popSize,length = pop.shape
    f = figure(figNum)
    hold(False)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    plot(arange(length), bitFreqs,'b.', markersize=5)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def visualizeRun(avgFitnessHist, maxFitnessHist, figNum=2):
    f = figure(figNum)
    hold(False)
    plot(arange(len(avgFitnessHist)), avgFitnessHist, 'k-')
    hold(True)
    plot(arange(len(maxFitnessHist)), maxFitnessHist, 'c-')
    xlabel('Generation')
    ylabel('Fitness')
    f.show()


def evolve(fitnessFunction,
            length,
            popSize,
            maxGens,
            probMutation,
            probCrossover=1,
            sigmaScaling=True,
            sigmaScalingCoeff=1,
            SUS=True,
            visualizeGen=visualizeGen,
            visualizeRun=visualizeRun):

    maskReposFactor = 5
    uniformCrossoverMaskRepos = rand(popSize/2, (length+1)*maskReposFactor) < 0.5
    mutMaskRepos = rand(popSize, (length+1)*maskReposFactor) < probMutation

    avgFitnessHist = zeros(maxGens+1)
    maxFitnessHist = zeros(maxGens+1)

    pop = zeros((popSize, length), dtype='int8')
    pop[rand(popSize, length)<0.5] = 1

    for gen in xrange(maxGens):

        fitnessVals = fitnessFunction(pop)
        fitnessVals = transpose(fitnessVals)
        maxFitnessHist[gen] = fitnessVals.max()
        avgFitnessHist[gen] = fitnessVals.mean()

        print "gen = %.3d   avgFitness = %3.3f  maxfitness = %3.3f" % (gen, avgFitnessHist[gen], maxFitnessHist[gen])
        if visualizeGen:
            visualizeGen(pop, gen=gen, avgFitness=avgFitnessHist[gen], maxFitness=maxFitnessHist[gen])
        if sigmaScaling:
            sigma = std(fitnessVals)
            if sigma:
                fitnessVals = 1 + fitnessVals - fitnessVals.mean()/sigmaScalingCoeff
                fitnessVals[fitnessVals<0] = 0
            else:
                fitnessVals = ones(1,popSize)

        cumNormFitnessVals = cumsum(fitnessVals/fitnessVals.sum())
        if SUS:
            markers = random.random() + arange(popSize,dtype='float')/popSize
            markers[markers>1] = markers[markers >1] - 1
        else:
            markers = rand(1, popSize)
        markers = sort(markers)
        parentIndices = zeros(popSize, dtype='int16')
        ctr = 0
        for idx in xrange(popSize):
            while markers[idx]>cumNormFitnessVals[ctr]:
                ctr += 1
            parentIndices[idx] = ctr
        random.shuffle(parentIndices)

        # deterimine the first parents of each mating pair
        firstParents = pop[parentIndices[0:popSize/2],:]
        # determine the second parents of each mating pair
        secondParents = pop[parentIndices[popSize/2:],:]

        temp = floor(random.random() * length * maskReposFactor-1)
        masks = uniformCrossoverMaskRepos[:, temp:temp+length]
        reprodIndices = rand(popSize/2)<1-probCrossover
        masks[reprodIndices, :] = False
        firstKids = firstParents
        firstKids[masks] = secondParents[masks]
        secondKids = secondParents
        secondKids[masks] = firstParents[masks]
        pop = vstack((firstKids, secondKids))

        temp = floor(random.random()*length*(maskReposFactor-1))
        masks = mutMaskRepos[:, temp:temp+length]
        pop[masks] = pop[masks] + 1
        pop = remainder(pop, 2)

    visualizeRun(avgFitnessHist, maxFitnessHist)

########## Stochastic Effective Attribute Parity ###########
def stochasticEffectiveAttributeParity(pop, pivLoci):
    return remainder(pop[:, pivLoci].sum(axis=1),2)*0.5-.25+randn(len(pop))

def seapVisualizeGen(pop, gen, avgFitness, maxFitness, pivLoci, figNum=1):
    popSize,length = pop.shape
    f = figure(figNum)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    hold(False)
    plot(arange(length), bitFreqs,'b.', markersize=2)
    hold(True)
    plot(pivLoci, bitFreqs[pivLoci], 'r.', markersize=15)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def seapEvolve(length, probMutation, probCrossover, popSize, maxGens):
    pivLoci = floor(rand(4)*length).astype('int16')
    evolve(partial(stochasticEffectiveAttributeParity,pivLoci=pivLoci),
        length,
        popSize,
        maxGens,
        probMutation,
        probCrossover,
        visualizeGen=partial(seapVisualizeGen,pivLoci=pivLoci, figNum=1),
        visualizeRun=partial(visualizeRun,figNum=2))

####### Staircase Function ##########
def staircaseFunction(pop, L, V, delta, sigma):
    popSize, _ = pop.shape
    m, n = L.shape
    fitnessVals=randn(popSize)*sigma
    for i, chrom in enumerate(pop):
        for j in xrange(m):
            if all(chrom[L[j,:]] == V[j,:]):
                fitnessVals[i] += delta
            else:
                fitnessVals[i] -= delta/(2**n - 1)
                break
    return fitnessVals

def staircaseFunctionVisualize(pop, gen, avgFitness, maxFitness, L, figNum=1):
    m,n = L.shape
    colorMap = {0:'r', 1:'b', 2:'g', 3:'c', 4:'m', 5:'y', 6:'k'}
    popSize,length = pop.shape
    f = figure(figNum)
    hold(False)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    plot(arange(length), bitFreqs,'b.', markersize=2)
    hold(True)
    for i in xrange(L.shape[0]):
        plot(L[i,:], bitFreqs[L[i,:]], colorMap[i%len(colorMap)]+'.', markersize=30-2*i)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def staircaseFunctionEvolve(length, numSteps, order, delta, sigma, probMutation, probCrossover, popSize, maxGens):
    L = arange(length)
    random.shuffle(L)
    L=L[:order*numSteps]
    L.shape=(-1,order)
    V=ones(L.shape, dtype='int8')
    evolve(partial(staircaseFunction, L=L, V=V, delta=delta, sigma=sigma),
        length,
        popSize,
        maxGens,
        probMutation,
        probCrossover=probCrossover,
        visualizeGen=partial(staircaseFunctionVisualize,L=L, figNum=3),
        visualizeRun=partial(visualizeRun, figNum=4))

######################################
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Run SpeedyGA on `seap` or `staircase`, two fitness functions '+
                                                 'tailor made to provide proof of concept for the Hyperclimbing Hypothesis. '+
                                                 'More details at http://blog.hackingevolution.net/2013/01/20/foga-2013-slides/')
    parser.add_argument('--fitnessFunction', default="staircase", choices=['staircase','seap'], help="The fitness function to use (default: staircase).")
    parser.add_argument('--probCrossover', type =float, default=1, help="Number between 0 and 1 representing the fraction of the population subject to crossover (default:1)" )
    parser.add_argument('--probMutation', type =float, default=0.003, help="The per bit mutation probability (default:0.003)" )
    parser.add_argument('--popSize', type =int, default=1000, help="Size of the population (default:1000)")
    parser.add_argument('--bitstringLength', type =int, default=500, help="Length of a chromosome in the population (default:500)")
    parser.add_argument('--gens', type =int, default=500, help="The number of generations (default:500)")

    args = parser.parse_args()
    if args.fitnessFunction=="seap":
        seapEvolve(length=args.bitstringLength,
                    probMutation=args.probMutation,
                    probCrossover=args.probCrossover,
                    popSize=args.popSize,
                    maxGens=args.gens)
    else:
        staircaseFunctionEvolve(length=args.bitstringLength,
                                numSteps=10,
                                order=4,
                                delta=0.3,
                                sigma=1,
                                probCrossover=args.probCrossover,
                                probMutation=args.probMutation,
                                popSize=args.popSize,
                                maxGens=args.gens)
    raw_input('Hit Enter to end ...')
