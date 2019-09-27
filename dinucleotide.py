import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import string
import tqdm as tq


"""
This module contains the dinucletotide class.

The dinculeotide class contains methods that calculate dinucleotide relative abudances for a given genome.
The class also methods for "simulating" genomes and plotting global dinucleotide frequencies.

The dinucleotide object takes 3 arguments: A fasta file, a window size and a skip parameter.

"""


class dinucleotide(object):

    di_RC = {'CG': 'CG', 'GC': 'GC', 'CC': 'GG', 'GG': 'CC', 'TT': 'AA', 'AA': 'TT',
             'TG': 'CA', 'CA': 'TG', 'AG': 'CT', 'GT': 'AC', 'GA': 'TC', 'TC': 'GA', 'TA': 'TA', 'AT': 'AT', 'CT': 'AG', 'AC': 'GT'}
    nuc_RC = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, fasta_file=None, window_size=50000, skip=2500, seq=None):

        self.window_size = window_size
        self.skip = skip
        self.header = None
        self.seq = None

        if fasta_file != None:
            self.header, self.seq = dinucleotide.parse_fasta(fasta_file)
        elif seq != None:
            self.header = None
            self.seq = seq
        else:
            raise Error('Enter a Fasta file or Sequence string.')

        self.length = len(self.seq)
        # added by Aidan
        self.k = 2

        self.deltas = []
        # These three variables will be filled when the following methods are called
        self.global_nuc_freqs = None
        self.global_di_freqs = None
        self.global_di_abundance = None

        # When you make an istance of this class the global counts and relative abundance will immediately
        # be calculated
        self.global_counts_freqs()
        self.global_relative_abundance()

        self.suspect_sequence = ""

        self.peaks = {}

    def global_counts_freqs(self):
        """Calculate GLOBAL nucleotide and dinucleotide frequencies.
        """

        self.global_nuc_freqs, self.global_di_freqs = dinucleotide.Calc_frequencies(
            self.seq, 2)
        return

    def global_relative_abundance(self):
        """
                Calculate GLOBAL nucleotide and dinucleotide relative abundances.
        """
        self.global_di_abundance = dinucleotide.calc_DN_relative_abund(
            self.global_nuc_freqs, self.global_di_freqs)
        return

    def calculate_local_freqs(self):
        """

                This method calculates local dinucleotide frequencies and compares them to their genomes global values by
                taking the difference between them.

                Specifically, this method will loop through the entire genome using a slidding window, skipping by some given amount (skip value (self.skip))

                Self.deltas will be a list of deltas that we later plot

                Returns:
                    -A list of deltas (distances) that is saved as the class attribute self.deltas
                    -Keeps track of sequence with the largest delta
        """

        max_delta = -999
        print("Calculating local dinucleotide frequencies")
        for i in tq.tqdm(range(0, len(self.seq) - self.window_size + 1, self.skip)):

            # print("Current window ", str(int(i / self.skip)))

            window_seq = self.seq[i:i + self.window_size]
            window_nuc_freqs, window_di_freqs = dinucleotide.Calc_frequencies(
                window_seq, 2)

            window_relative_abund = dinucleotide.calc_DN_relative_abund(
                window_nuc_freqs, window_di_freqs)

            sigma = 0

            # For each window, we want to sum the differences between each dinucleotide relative abundance
            # in the window AND the global relative abundance.
            # For example: AA relative abundance(in current window) - AA relative abundance(global)
            # We sum these

            for k, v in window_relative_abund.items():
                sigma += abs(v - self.global_di_abundance[k])

            # Multiply the sum by 1/16
            window_value = float(1 / 16) * sigma

            if window_value > max_delta:
                self.suspect_sequence = window_seq

            self.peaks[(i, i + self.window_size)] = window_value
            self.deltas.append(window_value)

        return

    def return_suspect_sequence(self):
        return self.suspect_sequence

    def grab_top_peaks(self, top_n):
        """
                        Sorts self.peaks and then reassigns it to only include top 5 deltas

            Params:
                -top_n: number of top deltas to return (int)

            returns: updates self.peaks to only include top peaks! This is sloppy but like its ok
            because we just need to breathe and take this one step at a time
        """
        newA = dict(sorted(self.peaks.items(),
                           key=operator.itemgetter(1), reverse=True)[:top_n])

        #self.peaks = newA
        return newA

    def write_suspect_to_fasta(self, outfile, top_n_peaks=0):
        """
                This method creates a fasta file that contains sequences associated with the top peaks.
                Default behavior writes sequence associated with largest peak.
                Params:
                        -outfile
                        -top_5_peaks: take on boolean value (Default is False)

                Returns: Fasta file is created
        """

        if top_n_peaks > 0:

            newA = self.grab_top_peaks(top_n_peaks)

            seq = []
            for k, v in newA.items():
                seq.append(self.seq[k[0]:k[1]])

        else:
            seq = [self.suspect_sequence]

        fh = open(outfile, "w")

        for temp_seq in seq:
            header = dinucleotide.random_generator(6) + "\n"
            temp = temp_seq + "\n"
            fh.write(header)
            fh.write(temp)

        fh.close()

        return

    def plot_deltas(self):
        """
                Plots all deltas for the genome
        """

        x_axis = [i * self.skip for i in range(len(self.deltas))]
        plt.plot(x_axis, self.deltas)
        plt.ylabel('Window delta relative abundance')
        plt.xlabel('Genome position')
        plt.show()

    def plot_deltas_wAnnotation(self, gff_file, annot="phage"):
        """
                Plots all deltas for the genome with annotated regions
        """

        roi = []

        with open(gff_file, 'r') as gff:
            for line in gff:
                if annot in line:
                    split = line.split('\t', 5)
                    roi.append((int(split[3]), int(split[4])))

        x_axis = [i * self.skip for i in range(len(self.deltas))]
        plt.plot(x_axis, self.deltas)

        for r in roi:
            plt.axvspan(r[0], r[1], color='red', alpha=0.25)

        plt.ylabel('Window delta relative abundance')
        plt.xlabel('Genome position')
        plt.legend()
        plt.show()

    def plot_distribution(self):
        """
                Plots the distribution of deltas
        """

        plt.hist(self.deltas, bins=100, alpha=0.5, label='x')
        # plt.hist(y, bins, alpha=0.5, label='y')
        plt.show()

    def getProportions(self):
        """
        Calculates nucleotides content. This method is used specifically for the simulateSequence method.
        Also this should just be part global_counts_freqs() because this is MOST DEF redundant
        """

        props = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

        for s in self.seq:
            props[s] += 1

        for k in props.keys():
            props[k] = props[k] / len(self.seq)

        return list(props.values())

    def simulateSequence(self):
        """
        Simulates DNA with specific nucleotide frequencies
        """

        draw = np.random.choice(['A', 'T', 'C', 'G'],
                                len(self.seq), p=self.getProportions())
        seq = "".join(draw)

        return seq

    @staticmethod
    def plot_compare_genome(genome_A, genome_B):
        """
                Plots deltas for two dinucleotide frequencies.
                params:
                    -genome_A: an instance of dinucleotide
                    -genome_B: another instance of dinucleotide

                returns: matplotliv figure
        """

        x_axis = [i * genome_A.skip for i in range(len(genome_A.deltas))]
        plt.plot(x_axis, genome_A.deltas)
        plt.plot(x_axis, genome_B.deltas)
        plt.ylabel('Window delta relative abundance')
        plt.xlabel('Genome position')
        plt.show()

    def plot_compare_distribution(genome_A, genome_B):
        """
            Plot
        """

        plt.hist(genome_A.deltas, bins=100, alpha=0.5)
        plt.hist(genome_B.deltas, bins=100, alpha=0.5)
        plt.show()

    @staticmethod
    def Calc_frequencies(sequence, size=2):
        """
            Given a sequence and a size (kmer length), calculate the frequencies of kmers of length "size" ()

            Params:
                -sequence: A string
                -size: An int
        """

        di_counts, nuc_counts = dict(), dict()

        # loop through sequence using a incredibly sophisticated sliding window approach
        for i in range(len(sequence) - size + 1):
            """
                We are accounting for the double stranded nature of DNA, although this shouldnt matter.
            """
            forward_mer = sequence[i:i+size]  # AA
            reverse_mer = dinucleotide.di_RC[forward_mer]  # TT

            rev_nuc = dinucleotide.nuc_RC[forward_mer[0]]

            nuc_counts[forward_mer[0]] = nuc_counts.get(forward_mer[0], 0) + 1
            nuc_counts[rev_nuc] = nuc_counts.get(rev_nuc, 0) + 1

            di_counts[forward_mer] = di_counts.get(forward_mer, 0) + 1
            di_counts[reverse_mer] = di_counts.get(reverse_mer, 0) + 1

        # ooff not cute yall but to get the count just right we need to add the last nucleotide to the count
        nuc_counts[forward_mer[0]] = nuc_counts.get(forward_mer[0], 0) + 1
        nuc_counts[rev_nuc] = nuc_counts.get(rev_nuc, 0) + 1

        length = len(sequence) * 2.0  # double stranded DNA my fuckin guyyyy

        nuc_freq = {k: v/float(length) for k, v in nuc_counts.items()}

        effective_length = 2 * (length - size + 1)
        di_freq = {k: v/effective_length for k, v, in di_counts.items()}

        return nuc_freq, di_freq

    @staticmethod
    def calc_DN_relative_abund(nuc_freq, di_freq):
        """
            Calculate dinucleotide relative abundance.

            Params:
                -nuc_freq: A dictionary containing frequencies for individual nucleotides
                -di_freq: A dictionary containing frequencies of all dinucleotides
        """

        return {k: v/(nuc_freq[k[0]] * nuc_freq[k[1]]) for k, v in di_freq.items()}

    @staticmethod
    def parse_fasta(file_name):
        """
                This function takes in a fasta file and returns its header
                and sequence. If fasta files contains several chromsomes than
                the sequenes will be concatenated.

                Prams

        """
        seq = ""
        headers = []

        with open(file_name) as fh:
            for line in fh:
                if(line[0] != '>'):
                    seq += line.strip()
                else:
                    headers.append(line.strip())
        return headers, seq

    @staticmethod
    def random_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for x in range(size))

    def __str__(self):
        """
        Special method that return the string representation of this class

        """
        return "This dinucleotide instance uses the file {} and contains a sequence of length {} ".format(self.fasta_file, str(self.length))


def main():

    bacteria1 = dinucleotide(fasta_file="random.fna")
    bacteria1.calculate_local_freqs()
    bacteria1.plot_deltas()


if __name__ == "__main__":
    main()
