import dinucleotide as dn



def main():
    if(len(sys.argv) <= 1):
        raise NameError('Did not provide a FASTA file argument.')
    # ps, length = getProportions(sys.argv[1])
    # newSeq = genSequence(ps, length)
    #
    # seq_record = SeqRecord(Seq(newSeq, SingleLetterAlphabet()), id="simulated_dna")
    #
    # file_name = 'sim.fna'
    #
    # out = SeqIO.write(seq_record, file_name, "fasta")
    #
    # sim = dn.dinucleotide(file_name)
    # sim.calculate_local_freqs()
    #
    # bacteria1 = dn.dinucleotide("random.fna")
    # bacteria1.calculate_local_freqs()
    #
    # bacteria1.plot_both(sim)

    h, seq = dn.dinucleotide.parse_fasta(sys.argv[1])


    return None

if __name__ == "__main__":
	main()
