from contextlib import closing
from matplotlib.pyplot import plot, figure, hold, axis, ylabel, xlabel, savefig, title
from numpy import sort, logical_xor, transpose, logical_not
from numpy.numarray.functions import cumsum, zeros
from numpy.random import rand, shuffle
from numpy import mod, floor
import time
import cloud
from durus.file_storage import FileStorage
from durus.connection import Connection

def bitFreqVisualizer(effectiveAttrIndices, bitFreqs, gen):
    f = figure(1)
    n = len(bitFreqs)
    hold(False)
    plot(range(n), bitFreqs,'b.', markersize=10)
    hold(True)
    plot(effectiveAttrIndices, bitFreqs[effectiveAttrIndices],'r.', markersize=10)
    axis([0, n-1, 0, 1])
    title("Generation = %s" % (gen,))
    ylabel('Frequency of the Bit 1')
    xlabel('Locus')
    f.canvas.draw()
    f.show()

def showExperimentTimeStamps():
    with closing(FileStorage("soda_results.durus")) as durus:
        conn = Connection(durus)
        return conn.get_root().keys()

def neap_uga(m, n, gens, probMutation, effectiveAttrIndices, probMisclassification, bitFreqVisualizer=None):
    """ neap = "noisy effective attribute parity"
    """
    pop = rand(m,n)<0.5
    bitFreqHist= zeros((n,gens+1))
    for t in range(gens+1):

        print "Generation %s" % t

        bitFreqs = pop.astype('float').sum(axis=0)/m
        bitFreqHist[:,t] = transpose(bitFreqs)

        if bitFreqVisualizer:
            bitFreqVisualizer(bitFreqs,t)

        fitnessVals = mod(pop[:, effectiveAttrIndices].astype('byte').sum(axis=1) +
                          (rand(m) < probMisclassification).astype('byte'),2)
        totalFitness = sum (fitnessVals)
        cumNormFitnessVals = cumsum(fitnessVals).astype('float')/totalFitness

        parentIndices = zeros(2*m, dtype='int16')
        markers = sort(rand(2*m))
        ctr = 0
        for idx in xrange(2*m):
            while markers[idx]>cumNormFitnessVals[ctr]:
                ctr += 1
            parentIndices[idx] = ctr
        shuffle(parentIndices)

        crossoverMasks = rand(m, n) < 0.5
        newPop = zeros((m, n), dtype='bool')
        newPop[crossoverMasks] = pop[parentIndices[:m], :][crossoverMasks]
        newPop[logical_not(crossoverMasks)] = pop[parentIndices[m:], :][logical_not(crossoverMasks)]

        mutationMasks = rand(m, n)<probMutation
        pop = logical_xor(newPop,mutationMasks)

    return bitFreqHist[0, :], bitFreqHist[-1, :]

def f(gens):
        k = 7
        n= k + 1
        effectiveAttrIndices = range(k)
        probMutation = 0.004
        probMisclassification = 0.20
        popSize = 1500
        jid = cloud.call(neap_uga, **dict(m=popSize,
                       n=n,
                       gens=gens,
                       probMutation=probMutation,
                       effectiveAttrIndices=effectiveAttrIndices,
                       probMisclassification=probMisclassification))
        print "Kicked off trial %s" % jid
        return jid

def cloud_result(jid):
    result = cloud.result(jid)
    print "Retrieved results for trial %s" % jid
    return result

def run_trials():
    numTrials = 3000
    gens = 1000
    from multiprocessing.pool import ThreadPool as Pool
    pool = Pool(50)

    jids = pool.map(f,[gens]*numTrials)
    print "Done spawning trials. Retrieving results..."

    results = pool.map(cloud_result, jids)
    firstLocusFreqsHists = zeros((numTrials,gens+1), dtype='float')
    lastLocusFreqsHists = zeros((numTrials,gens+1), dtype='float')
    print "Done retrieving results. Press Enter to serialize..."

    raw_input()

    for i, result in enumerate(results):
        firstLocusFreqsHists[i, :], lastLocusFreqsHists[i, :] = result

    with closing(FileStorage("soda_results.durus")) as durus:
        conn = Connection(durus)
        conn.get_root()[str(int(floor(time.time())))] = (firstLocusFreqsHists, lastLocusFreqsHists)
        conn.commit()

    pool.close()
    pool.join()

def render_results(timestamp=None):

    with closing(FileStorage("soda_results.durus")) as durus:
        conn = Connection(durus)
        db = conn.get_root()
        if not timestamp:
            timestamp = sorted(db.keys())[-1]
        firstLocusFreqsHists, lastLocusFreqsHists = db[timestamp]
    print "Done deserializing results. Plotting..."

    x = [(2, 'First', firstLocusFreqsHists, "effective"),
         (3, 'Last', lastLocusFreqsHists, "non-effective")]

    for i, pos, freqsHists, filename in x :
        freqsHists = freqsHists[:,:801]
        f = figure(i)
        hold(False)
        plot(transpose(freqsHists), color='grey')
        hold(True)
        maxGens = freqsHists.shape[1]-1
        plot([0, maxGens], [.05,.05], 'k--')
        plot([0, maxGens], [.95,.95], 'k--')
        axis([0, maxGens, 0, 1])
        xlabel('Generation')
        ylabel('1-Frequency of the '+pos+' Locus')
        f.canvas.draw()
        f.show()
        savefig(filename+'.png', format='png', dpi=200)

if __name__ == "__main__":
    run_trials()
    render_results()