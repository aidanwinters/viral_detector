from itertools import product
import numpy as np
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def getKeys(k):
    return [''.join(x) for x in product(['A', 'T', 'C', 'G'], repeat = k)]

def genDictionary(k):
    return dict(zip(getKeys(k), np.arange(4**k)))

def getFreq(genome, k, window = None, make_continuous = True):

    dict = genDictionary(k)
    temp = np.zeros(len(dict), dtype="int32")

    # print out info:
    # number of k-mer combos
    # the size of the Genome
    #the window size and how many windows there are in genome

    fasta_seq = SeqIO.parse(open(genome),'fasta')
    final_seq = '';
    if make_continuous:
        #combine all the different sequences (important if there are chromosomes??)
        final_seq = [(final_seq + str(f.seq)) for f in fasta_seq][0]

    if window == None:
        print('true')
        window = len(final_seq)

    print('---------------------------')
    print('Number of K-mer Combos:', 4**k)
    print('Size of the Genome:', len(final_seq))
    print('Number of Windows', len(final_seq)/ window)
    print('---------------------------')

    freqs = []
    windows = np.arange(stop=len(final_seq), start=0, step=window)
    windows = np.append(windows, len(final_seq))

    for w in np.arange(len(windows) - 1):
        curr_freq = np.array(temp.copy())
        chunk = final_seq[windows[w]:(windows[w+1])]
        for j in np.arange(len(chunk) - k):
            curr_freq[dict[chunk[j:(j+k)]]] += 1
        freqs.append(curr_freq / window)

    return(freqs)



def main():
    k_val = 3
    freqs = getFreq('random.fna',k_val, window = 100000)
    df = pd.DataFrame(freqs, columns=getKeys(k_val))
    df.plot()
    plt.show()

if __name__ == "__main__":

	main()
