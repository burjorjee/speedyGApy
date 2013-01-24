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

def visualizeRun(avgFitnessHist, maxFitnessHist, maxGens, figNum=2):
    f = figure(figNum)
    hold(False)
    plot(arange(maxGens + 1), avgFitnessHist, 'k-')
    hold(True)
    plot(arange(maxGens + 1), maxFitnessHist, 'c-')
    xlabel('Generation')
    ylabel('Fitness')
    f.show()


def evolve(fitnessFunction,
            length,
            popSize=500,
            maxGens=500,
            probCrossover=1,
            probMutation=0.003,
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
        masks = uniformCrossoverMaskRepos[:, temp:temp+length-1]
        reprodIndices = rand(popSize, 1)<1-probCrossover
        masks[reprodIndices, :] = False
        firstKids = firstParents
        firstKids[masks] = secondParents[masks]
        secondKids = secondParents
        secondKids[masks] = firstParents[masks]
        pop = vstack((firstKids, secondKids))

        temp = floor(random.random()*length*(maskReposFactor-1))
        masks = mutMaskRepos[:, temp:temp+length-1]
        pop[masks] = pop[masks] + 1
        pop = remainder(pop, 2)

    visualizeRun(avgFitnessHist, maxFitnessHist, maxGens)

########## Stochastic Effective Attribute Parity ###########
def stochasticEffectiveAttributeParity(pop, pivLoci):
    return remainder(pop[:, pivLoci].sum(axis=1),2)*0.5-.25+randn(len(pop))

def seapVisualizeGen(pop, gen, avgFitness, maxFitness, pivLoci, figNum=1):
    popSize,length = pop.shape
    f = figure(figNum)
    hold(False)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    plot(arange(length), bitFreqs,'b.', markersize=2)
    hold(True)
    plot(pivLoci, bitFreqs[pivLoci], 'r.', markersize=15)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def seapEvolve():
    length = 1000
    pivLoci = floor(rand(4)*length).astype('int16')
    evolve(partial(stochasticEffectiveAttributeParity,pivLoci=pivLoci), length,
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
    colorMap = {0:'r', 1:'b', 2:'g', 3:'c', 4:'m'}
    popSize,length = pop.shape
    f = figure(figNum)
    hold(False)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    plot(arange(length), bitFreqs,'b.', markersize=2)
    hold(True)
    for i in xrange(L.shape[0]):
        plot(L[i,:], bitFreqs[L[i,:]], colorMap[i%n]+'.', markersize=30-2*i)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def staircaseFunctionEvolve():
    length = 500
    n = 4
    L = arange(length)
    random.shuffle(L)
    L=L[:n*10]
    L.shape=(-1,n)
    V=ones(L.shape, dtype='int8')
    evolve(partial(staircaseFunction, L=L, V=V, delta=0.3, sigma=1),
        length,
        visualizeGen=partial(staircaseFunctionVisualize,L=L, figNum=3),
        visualizeRun=partial(visualizeRun, figNum=4))
######################################

if __name__=="__main__":
    seapEvolve()
    staircaseFunctionEvolve()
    raw_input('Hit Enter to end ...')
