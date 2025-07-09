# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# ReadAnalysis.py
# class to analyze a given folder of fasta and fastq sequence reads

import gzip
from pathlib import Path
import numpy as np
from pprint import pprint
import time
from datetime import datetime

class ReadAnalysis:
    # initialize analysis lists and statistics dictionary
    def __init__(self, directorypath):  # full path from home needs to be given
        # initialize given path to directory of files that should be analyzed
        self.directorypath = Path.home() / directorypath

        if not self.directorypath.exists() or not self.directorypath.is_dir():
            raise ValueError(f"The path {self.directorypath} does not exist or is not a directory.")

        # lists that hold content from analysis
        self.readlengthvalues = []
        self.basecountvalues = {base: [] for base in ['A', 'C', 'T', 'G', 'N']}
        self.gccontents = []
        self.ncontents = []
        self.readsperfastafile = []
        self.readsperfastqfile = []
        self.qualityscores = []

        # statistics dictionary (0 bc no analysis performed yet)
        self.stats = {
            'readlengths': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0,
                'q25': 0,
                'q75': 0
                # figure this out:
                # 'distribution_type': 'normal',  # or 'gamma', 'lognormal', etc.
                # 'parameters': {'loc': 150.2, 'scale': 25.1}  # for sampling
            },
            'basecounts': {
                'A': {
                    'mean': 0,
                    'median': 0,
                    'std': 0,
                    'min': 0,
                    'max': 0
                },
                'C': {
                    'mean': 0,
                    'median': 0,
                    'std': 0,
                    'min': 0,
                    'max': 0,
                },
                'T': {
                    'mean': 0,
                    'median': 0,
                    'std': 0,
                    'min': 0,
                    'max': 0
                },
                'G': {
                    'mean': 0,
                    'median': 0,
                    'std': 0,
                    'min': 0,
                    'max': 0
                },
                'N': {
                    'mean': 0,
                    'median': 0,
                    'std': 0,
                    'min': 0,
                    'max': 0
                },
                'totalbases': 0
            },
            'gccontent': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            },
            'ncontent': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            },
            'readsperfasta': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            },
            'readsperfastq': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            },
            'qualityscore': {
                'mean': 0,
                'median': 0,
                'std': 0,
                'min': 0,
                'max': 0
            }
        }

    # func to calculate statistics after data added to lists
    def calculateStats(self):
        if self.readlengthvalues:
            self.stats['readlengths']['mean'] = np.mean(self.readlengthvalues)
            self.stats['readlengths']['median'] = np.median(self.readlengthvalues)
            self.stats['readlengths']['std'] = np.std(self.readlengthvalues)
            self.stats['readlengths']['min'] = np.min(self.readlengthvalues)
            self.stats['readlengths']['max'] = np.max(self.readlengthvalues)
            self.stats['readlengths']['q25'] = np.percentile(self.readlengthvalues, 25)
            self.stats['readlengths']['q75'] = np.percentile(self.readlengthvalues, 75)

        self.stats['basecounts'] = {base: {} for base in ['A', 'C', 'T', 'G', 'N']}
        for base in ['A', 'C', 'T', 'G', 'N']:
            if self.basecountvalues[base]:
                base_values = np.array(self.basecountvalues[base])
                self.stats['basecounts'][base]['mean'] = np.mean(base_values)
                self.stats['basecounts'][base]['median'] = np.median(base_values)
                self.stats['basecounts'][base]['std'] = np.std(base_values)
                self.stats['basecounts'][base]['min'] = np.min(base_values)
                self.stats['basecounts'][base]['max'] = np.max(base_values)

        if self.gccontents:
            self.stats['gccontent'] = {
                'mean': np.mean(self.gccontents),
                'median': np.median(self.gccontents),
                'std': np.std(self.gccontents),
                'min': np.min(self.gccontents),
                'max': np.max(self.gccontents)
            }

        if self.ncontents:
            self.stats['ncontent'] = {
                'mean': np.mean(self.ncontents),
                'median': np.median(self.ncontents),
                'std': np.std(self.ncontents),
                'min': np.min(self.ncontents),
                'max': np.max(self.ncontents)
            }

        if self.readsperfastafile:
            self.stats['readsperfasta'] = {
                'mean': np.mean(self.readsperfastafile),
                'median': np.median(self.readsperfastafile),
                'std': np.std(self.readsperfastafile),
                'min': np.min(self.readsperfastafile),
                'max': np.max(self.readsperfastafile)
            }

        if self.readsperfastqfile:
            self.stats['readsperfastq'] = {
                'mean': np.mean(self.readsperfastqfile),
                'median': np.median(self.readsperfastqfile),
                'std': np.std(self.readsperfastqfile),
                'min': np.min(self.readsperfastqfile),
                'max': np.max(self.readsperfastqfile)
            }

        if self.qualityscores:
            self.stats['qualityscore'] = {
                'mean': np.mean(self.qualityscores),
                'median': np.median(self.qualityscores),
                'std': np.std(self.qualityscores),
                'min': np.min(self.qualityscores),
                'max': np.max(self.qualityscores)
            }

    # function to print analysis once completed
    def printAnalysis(self):
        self.calculateStats()
        pprint(self.stats)

    # function to return the stat analysis 

    # function to analyze a single sequence of a file
    # updates statistics as part of analysis
    def analyzeSequence(self, sequence, quality=None):
        seqlength = len(sequence)
        # add to variables/lists
        self.readlengthvalues.append(seqlength)
        self.stats['basecounts']['totalbases'] += seqlength

        # base composition - count all bases in one pass
        basecountsdict = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}
        gccount = 0

        for base in sequence:
            if base in basecountsdict:
                basecountsdict[base] += 1
                if base in 'GC':
                    gccount += 1

        # update base count values (append counts to the lists)
        for base in ['A', 'C', 'T', 'G', 'N']:
            self.basecountvalues[base].append(basecountsdict[base])

        # calculate percentages
        thisncontent = (basecountsdict['N'] / seqlength) * 100 if seqlength > 0 else 0
        thisgccontent = (gccount / seqlength) * 100 if seqlength > 0 else 0

        self.ncontents.append(thisncontent)
        self.gccontents.append(thisgccontent)

        # quality scores
        if quality:
            # convert ASCII chars to numeric phred scores
            scores = [ord(char) - 33 for char in quality]
            self.qualityscores.extend(scores)

    # function to determine how to open file (gzipped or not)
    def openFile(self, filename):
        filepath = Path(self.directorypath) / filename
        if filepath.name.endswith('.gz'):
            return gzip.open(filepath, 'rt')
        else:
            return open(filepath, 'r')

    # function to analyze entire fasta file
    def analyzeFastaFile(self, filename):
        print(f'analyzing FASTA file: {filename}')
        amofreads = 0

        with self.openFile(filename) as file:
            sequence = ""

            for line in file:
                line = line.strip()

                if line.startswith('>'):
                    # then, if there is a sequence, analyze it
                    if sequence:
                        # test print
                        # print(f'length: {len(sequence)} \n first 10 letters: {sequence[:10]} \n')
                        self.analyzeSequence(sequence)
                        amofreads += 1
                    sequence = ""  # clear the sequence

                else:  # if line doesn't start with >, then add it to the sequence string
                    sequence += line.upper()  # convert all values to uppercase

            # after going through the file, process the last sequence
            if sequence:
                self.analyzeSequence(sequence)
                amofreads += 1

        # add amofreads to list of sequence reads per fasta file
        self.readsperfastafile.append(amofreads)

    # func to analyze one entire fastq file
    def analyzeFastqFile(self, filename):
        print(f'analyzing FASTQ file: {filename}')
        amofreads = 0

        with self.openFile(filename) as file:
            linecount = 0
            sequence = ""
            quality = ""

            for line in file:
                line = line.strip()
                linecount += 1

                # structured like this bc the file considers each sequence to have 4 lines, either:
                if linecount % 4 == 1:  # header line
                    continue
                elif linecount % 4 == 2:  # sequence line
                    sequence = line.upper()  # convert to upper and assign line to sequence
                elif linecount % 4 == 3:  # + line (representation of header)
                    continue
                elif linecount % 4 == 0:  # quality line
                    quality = line  # assign line to quality
                    self.analyzeSequence(sequence, quality)  # analyze sequence only when have the sequence and qualities
                    amofreads += 1
                    # clear/refresh variables
                    sequence = ""
                    quality = ""

        # add amofreads to list of sequence reads per fastq file
        self.readsperfastqfile.append(amofreads)

    # funct to log the run time of the analysis
    def logruntime(self, seconds):
        now = datetime.now()
        date = now.strftime("%Y-%m-%d")
        starttime = now.strftime("%H:%M:%S")

        with open("runtimelog.txt", "a") as runtimefile:
            runtimefile.write(f"\n{date}\t{starttime}\t{seconds:.2f}")

    # func to test fasta file analysis
    def fastaTest(self):
        # test fasta file analysis with run time
        startanalysis = time.time()
        self.directorypath = Path('references/sequencereadfiles')
        self.analyzeFastaFile('PB644_EB816.hifi_reads.fasta.gz')
        endanalysis = time.time()
        runseconds = round(endanalysis - startanalysis, 2)

        self.logruntime(runseconds)
        print(f"Analysis took {runseconds} seconds to run")

        startprint = time.time()
        self.printAnalysis()
        endprint = time.time()
        print(f"Print stats took {endprint - startprint:.2f} seconds to run")

    # func to test fastq file analysis
    def fastqTest(self):
        # test fastq file analysis & log run time
        startanalysis = time.time()
        self.directorypath = Path('references/sequencereadfiles')
        self.analyzeFastqFile('PB644_EB820.hifi_reads.fastq.gz')
        endanalysis = time.time()
        runseconds = round(endanalysis - startanalysis, 2)

        self.logruntime(runseconds)
        print(f"Analysis took {runseconds} seconds to run")

        startprint = time.time()
        self.printAnalysis()
        endprint = time.time()
        print(f"Print stats took {endprint - startprint:.2f} seconds to run")

    # function to analyze all fastq files in a directory (folder)
    def analyzeAllFastq(self):
        # get all fastq files in the directory (.fastq & .fasq.gz)
        fastqfiles = list(self.directorypath.glob('*.fastq')) + list(self.directorypath.glob('*.fastq.gz'))

        for filepath in fastqfiles:
            self.analyzeFastqFile(filepath.name)

    # function to analyze all fasta files in a directory (folder)
    def analyzeAllFasta(self):
        # get all fasta files in the directory (.fasta & .fasta.gz)
        fastafiles = list(self.directorypath.glob('*.fasta')) + list(self.directorypath.glob('*.fasta.gz'))

        for filepath in fastafiles:
            self.analyzeFastaFile(filepath.name)

    # function to analyze all files in the given directory
    def analyzeDirectory(self):
        print(f'Analyzing directory: {self.directorypath}')
        self.analyzeAllFasta()
        self.analyzeAllFastq()

def main():
    folder = 'SynologyDrive/CompSci/DouglassResearch/references/SequenceReadFiles'
    try:
        analyzer = ReadAnalysis(folder)

        analyzer.fastqTest()
        analyzer.fastaTest()

        startanalysis = time.time()
        analyzer.analyzeDirectory()
        endanalysis = time.time()

        runseconds = round(endanalysis - startanalysis, 2)
        analyzer.logruntime(runseconds)
        print(f"Analysis took {runseconds} seconds to run")

        startprint = time.time()
        analyzer.printAnalysis()
        endprint = time.time()

        print(f"Print stats took {endprint - startprint:.2f} seconds to run")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    main()
