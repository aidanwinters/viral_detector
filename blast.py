import subprocess as sp
import sys,os,time,re
from collections import defaultdict


class Blast(object):

    """
        Document this class
    """

    def __init__(self,fasta_in,dbt,query,out = 'blastout.txt'):
        self.subject = fasta_in
        self.dbtype = dbt
        self.query = query
        self.outfile = out

    def makedb(self,force = False):
        """Create a custom database from a multi-FASTA file of sequences with this minimal command:
        makeblastdb -in mydb.fsa dbtype nucl -parse_seqids
        subprocess.Popen requires a list of arguments"""
        print ("Creating goddamn subject database...get stoked")
        print  (os.path.isfile(self.subject + ".nhr") and os.path.isfile(self.subject + ".nin") and os.path.isfile(self.subject + ".nsq"))

        if (os.path.isfile(self.subject + ".nhr") and os.path.isfile(self.subject + ".nin") and os.path.isfile(self.subject + ".nsq")) or force:
            return "Database already exists continuing"
        else:
           sp.Popen(['./ncbi-blast-2.9.0+/bin/makeblastdb','-in',self.subject,'-dbtype',self.dbtype]).wait()

        return "Finished creating blastdb my guy"

    def blast(self):
        """

        Run blastn or blastp

        Parameters for -db and -query are provided by the class instance

        """

        print("Running blast... OK dude, please chill")
        start = time.time()

        #use blastn
        if self.dbtype == "nucl":
        	sp.Popen(['./ncbi-blast-2.9.0+/bin/blastn','-db',self.subject,'-query',self.query,'-outfmt','7 stitle qstart sstart qend send qcovs qseqid','-out',self.outfile]).wait()
        #use blastp
        elif self.dbtype == "prot":
            sp.Popen(['./ncbi-blast-2.9.0+/bin/blastp','-db',self.subject,'-query',self.query,'-outfmt','7 stitle qstart sstart qend send qcovs qseqid','-out',self.outfile]).wait()

        end = time.time()

        return "Done running blast in %0.6f seconds" % (end - start)

class BlastReport(object):

    def __init__(self,blastresfile = 'blastout.txt'):
        self.file = blastresfile
        self.results = []
        self.parse() #call parse when new instance of the class is created

    def parse(self):
        find_query = re.compile("#\sQuery:\s.*").search #regex
        pos = -1

        with open(self.file) as fh:
            for line in fh:
                if find_query(line):
                    query = re.sub("#\sQuery:\s","",line)
                    result_dict = defaultdict(defaultdict)
                    self.results.append(result_dict)
                    pos += 1
                elif line.startswith("#"):
                    continue
                else:
                    line = line.split('\t')
                    self.results[pos][query][line[0]] = {'query_coords': line[1] + ":" + line[3],
                                                   'subject_coords' : line[2] + ":" + line[4],
                                                   'percent_query_coverage' : line[5] }
        print(self.results)

    def qcovs(self):
        #this takes the data that was parsed from the blast report and returns the percent query coverage
        qcov_list = []
        for result in self.results:
            for query in result:
                for subject in result[query]:
                    qcov_list.append((query,subject,result[query][subject]['percent_query_coverage']))
        return qcov_list

def main():
    query = "random.fna" #this will be the file from dinucleotie program
    subject = "./dbs/Actino_DB" #I hard coded this path so change it to make it work, wilson changed
    dbt = "nucl"
    outfile = "practice.out" #random ass outfile name
    useblast = Blast(subject,dbt,query,outfile)
    useblast.makedb()
    useblast.blast()
    blast_results = BlastReport()

if __name__ == '__main__':
    main()
