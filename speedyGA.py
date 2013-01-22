from functools import partial
from random import random
from matplotlib.pyplot import figure, plot, xlabel, ylabel, title, axis, hold
from numpy import *
from numpy.random import rand, randn

def visualizeGen(pop, gen, avgFitness, maxFitness):
    popSize,length = pop.shape
    f = figure(1)
    hold(False)
    bitFreqs = pop.sum(axis=0).astype('float')/popSize
    plot(arange(length), bitFreqs,'b.', markersize=2)
    axis([0, length, 0, 1])
    title("Generation = %s, Average Fitness = %0.3f " % (gen,avgFitness))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def evolve(fitnessFunction,
            length,
            popSize=500,
            maxGens=500,
            probCrossover=1,
            probMutation=0.005,
            sigmaScaling=True,
            sigmaScalingCoeff=1,
            SUS=True,
            visualizeGen=visualizeGen):

    maskReposFactor = 5
    uniformCrossoverMaskRepos = rand(popSize/2, (length+1)*maskReposFactor) < 0.5
    mutMaskRepos = rand(popSize, (length+1)*maskReposFactor) < probMutation

    avgFitnessHist = zeros(maxGens+1)
    maxFitnessHist = zeros(maxGens+1)

    pop = zeros((popSize, length), dtype='int8')
    pop[rand(popSize, length)<0.5] = 1

    f = figure(1)
    for gen in xrange(maxGens):

        fitnessVals = transpose(fitnessFunction(pop))
        maxFitnessHist[gen] = max(fitnessVals)
        avgFitnessHist[gen] = mean(fitnessVals)

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

    f=figure(2)
    hold(False)
    plot(arange(maxGens+1), avgFitnessHist,'k-')
    hold(True)
    plot(arange(maxGens+1), maxFitnessHist, 'c-')
    xlabel('Generation')
    ylabel('Fitness')
    f.show()

def stochasticEffectiveAttributeParity(pop, pivLoci):
    return remainder(pop[:, pivLoci].sum(axis=1),2)*0.5-.25+randn(len(pop))

def seapEvolve():
    length = 1000
    pivLoci = floor(rand(4)*length).astype('int16')
    evolve(partial(stochasticEffectiveAttributeParity,pivLoci=pivLoci), length,
        visualizeGen=partial(seapVisualizeGen,pivLoci=pivLoci))

def seapVisualizeGen(pop, gen, avgFitness, maxFitness, pivLoci):
    popSize,length = pop.shape
    f = figure(1)
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

if __name__=="__main__":
    seapEvolve()