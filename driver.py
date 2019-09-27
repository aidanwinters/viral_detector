import dinucleotide as dn
import blast
import sys


def main():
    if(len(sys.argv) <= 1):
        raise NameError('Did not provide a FASTA file argument.')

    bacteria1 = dn.dinucleotide(sys.argv[1])
    bacteria1.calculate_local_freqs()
    # bacteria1.plot_distribution()
    bacteria1.plot_deltas()
    bacteria1.plot_deltas_wAnnotation("ncbi-genomes-2019-09-27/GCF_000008865.2_ASM886v2_genomic.gff")


    #simulate a new sequence
    fake_seq = bacteria1.simulateSequence()
    sim_genome = dn.dinucleotide(seq=fake_seq)
    sim_genome.calculate_local_freqs()

    dn.dinucleotide.plot_compare_genome(bacteria1, sim_genome)
    dn.dinucleotide.plot_compare_distribution(bacteria1, sim_genome)

    #blast section

    #create fasta for uspect seq
    suspect = "suspect.fna"

    bacteria1.write_suspect_to_fasta(suspect)

    query = suspect #this will be the file from dinucleotie program
    subject = "./dbs/Actino_DB" #I hard coded this path so change it to make it work, wilson changed
    dbt = "nucl"
    outfile = "practice.out" #random ass outfile name
    useblast = blast.Blast(subject,dbt,query,outfile)
    useblast.makedb()
    useblast.blast()
    blast_results = blast.BlastReport(outfile)

    return None


if __name__ == "__main__":
    main()
    # Hasan was here
    # wilson hacl
