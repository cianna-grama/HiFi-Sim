# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# SimReads.py
# Class to generate synethetic reads
'''
hbfhwb
'''

import gzip
import random
import math
import numpy as np
import sys
from GenerateQscores import GenerateQscores
import time

class SimReads:
    def __init__(self, errorrate=0, insertf=0, deletef=0, subf=0):
        self.errorrate = errorrate
        self.insertf = insertf
        self.deletef = deletef
        self.subf = subf
        self.numreads = {'min': 0, 'max': 0, 'mean': 0, 'median': 0}
        self.readlength = {'min': 0, 'max': 0, 'mean': 0, 'median': 0, 'sd': 0}

    def readFasta(self, filename):
        '''
        Returns a dictionary of sequences for ease of reading through
        FASTA files. 
        Dictionary keys: seqeunce IDs.
        Dictionary values: sequence strings.
        '''

        filepath = "references/sequencereadfiles/" + filename
        sequences = {}
        currentseqid = None
        currentseq = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if currentseqid:
                        sequences[currentseqid] = "".join(currentseq)
                    currentseqid = line[1:].split()[0] # get id
                    currentseq = []
                else:
                    currentseq.append(line.upper()) 
            if currentseqid: # add last sequence
                sequences[currentseqid] = "".join(currentseq)
        return sequences
    
    def reverseComplement(self, dnasequence):
        '''
        Returns the reverse complement of a DNA sequence.
        '''

        # align complements according to dna structure
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} 

        # reverse the sequence, complement reversed bases, and join into a strand
        return "".join(complement.get(base, base) for base in reversed(dnasequence))

    def ishomopolymer(self, seq):
        '''
        Checks if a set of bases is a homopolymer.
        Homopolymer is a sequence of consecutive, identical 
        nucleotide bases. Ex: 'AAAAA'. 
        Helper for introduceErrors function.
        '''
        return len(seq) >= 3 and len(set(seq)) == 1

    def isrepeat(self, seq):
        '''
        Checks simple base repetitions. Ex: ATAT
        Helper for introduceErrors function.
        '''
        return any(seq.startswith(rep * 2) for rep in ['AT', 'TA', 'CG', 'GC', 'AG', 'GA', 'TC', 'CT'])

    def introduceErrors(self, sequenceref, errorrate, insertf, deletef, subf):
        '''
        Introduce sequencing errors at specified error rate. 
        Errors: insertions, deletions, substitutions.
            sequenceref (str): Sequence string to introduce errors to.
            errorate (int): Error rate in the form of fraction (int between 0.0 and 1.0).
            insertf (int): Insertion error fractional split of the error rate.
            deletef (int): Deletion error fractional split of the error rate.
            subf (int): Substitution error fractional split of the error rate.
        '''
        if errorrate == 0:
            return sequenceref

        # get total fraction of indels and substitutions
        sumf = insertf + deletef + subf

        # if no fractions requested, default to pacbio hifi percentages
        if sumf <= 0 and errorrate > 0:
            print("Error fractions sum to zero or less with inputted error rate")
            print("Default to PacBio Hifi read percentages: 25% insertions, 20% deletions, 55% substitutions")
            # calculate fractions into probability percentages 
            insertp, deletep, subp = 0.25, 0.20, 0.55
            print(f'insert probability: {insertp}\tdelete probability: {deletep}\tsubstitute probability: {subp}')

        # if fractions given, transform to probability percentages
        elif sumf > 0:
            insertp = insertf/sumf
            deletep = deletef/sumf
            subp = subf/sumf

        # if no fractions and no error rate, return original seqeunce as is
        else:
            return sequenceref

        # initialize varibles for output
        newbases = []
        indexref = 0
        dnabases = ['A', 'C', 'G', 'T']

        # process each base in the given reference sequeunce
        # go through sequence ref until indexref exceeds length of sequence
        while indexref < len(sequenceref):
            # current base is item at current index in the sequence
            currentrefbase = sequenceref[indexref]

            # note: random.random() selects a number between [0.0 and 1.0] at random
            errorchance = random.random()

            # this line will be true with a probability of given errorrate
            # probability that random number is between 0.0 and the errorrate is the given errorrate
            if errorchance <= errorrate:  
                errortyperoll = random.random() 

                # probability ranges:
                    # between 0 & deletep = delete
                    # between deletep & (deletep + insertp) = insert
                    # between (deletep + insertp) & 1 - substitute
                
                # deletion
                if errortyperoll < deletep:
                    # simulate variable deletion size based on real hifi distribution
                    deletionroll = random.random()
                    if deletionroll < 0.95:
                        deletionsize = 1
                    elif deletionroll < 0.98:
                        deletionsize = random.choice([2, 3])
                    elif deletionroll < 0.995:
                        deletionsize = random.choice([4, 5])
                    else:
                        deletionsize = random.randint(6, 10)

                    # skip that many bases in the reference sequence
                    indexref += deletionsize
                    continue # move to next base after deletion

                # insertion - matches real insertion environments 
                elif errortyperoll < deletep + insertp:
                    # determine insertion context
                    context = sequenceref[max(0, indexref - 3): indexref + 3]

                    insertionroll = random.random()
                    
                    if self.ishomopolymer(context):
                        # in homopolymer context, insert 1 base ~95% of the time
                        insertlength = 1
                    elif self.isrepeat(context):
                        # in repeat context, allow 2–3 bp insertions occasionally
                        if insertionroll < 0.95:
                            insertlength = 1
                        else:
                            insertlength = random.choice([2, 3])
                    else:
                        # rare chance for longer insertions
                        if insertionroll < 0.98:
                            insertlength = 1
                        elif insertionroll < 0.995:
                            insertlength = random.choice([2, 3])
                        else:
                            insertlength = random.randint(4, 6)

                    # insert synthetic bases
                    for i in range(insertlength):
                        newbases.append(random.choice(dnabases))
                    
                    # insert current base after the inserted ones
                    newbases.append(currentrefbase)
                    indexref += 1

                # substitution
                elif errortyperoll <= 1:
                    # create list of possible sub bases - ensures A is not subbed with A
                    possiblesubs = [base for base in dnabases if base != currentrefbase]
                    # edge case - if no possible subs (so current base not ATGC), random sub
                    if not possiblesubs:
                        newbases.append(random.choice(dnabases))
                    # expected path
                    else: 
                        newbases.append(random.choice(possiblesubs))
                    indexref += 1
                    
            # no error occurs - add orignal base to new list
            else: 
                newbases.append(currentrefbase)
                indexref += 1
                
        # finally, return the new bases joined as a sequence string
        return "".join(newbases)
    
    def introduceBiologicalVariation(self, seq, snprate=0.001, indelrate=0.0001):
        """
        Introduce SNPs and indels into the read sequence 
        to simulate biological variation from the reference.
        Returns a new sequence string with variation.
        """
        
        bases = ['A', 'C', 'G', 'T']
        i = 0
        newseq = []

        while i < len(seq):
            
            # random number between 0 and 1, determines if variation occurs at that base position
            r = random.random()    

            # insertions
            if r < indelrate / 2:
                newseq.append(random.choice(bases))
                # no increment- insertion adds a base

            # deletions
            elif r < indelrate:
                i += 1  # skip  base to simulate deletion

            # SNPs
            elif r < indelrate + snprate:
                original = seq[i]
                newbase = random.choice([b for b in bases if b != original])
                newseq.append(newbase)
                i += 1

            else:
                newseq.append(seq[i])
                i += 1
        return ''.join(newseq)
    

    def mutateGenome(self, refgenomedict, snprate=0.001, indelrate=0.0001):
        """
        Takes in dictionary of reference genome and 
        returns a new genome dictionary with SNPs and 
        indels introduced into each chromosome/contig.
        """
        mutated = {}
        for name, seq in refgenomedict.items():
            mutated[name] = self.introduceBiologicalVariation(seq, snprate, indelrate)
        return mutated
    
    def loadReference(self, filepath):
        """
        Loads a reference genome (plain or gzipped FASTA) into a dictionary.
        Each key is a sequence (ex: chromosome or contig) name and value is the nucleotide string.
        Helper to simulateHiFiReads function.
        """
        
        # this dictionary can contain chromasomes or contigs loaded in as references
        sequences = {}
        currentname = None
        currentseq = []

        # handle gzipped or plaintext FASTA
        openfunc = gzip.open if filepath.endswith(".gz") else open

        with openfunc(filepath, 'rt') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if currentname:
                        sequences[currentname] = ''.join(currentseq)
                    currentname = line[1:].split()[0]
                    currentseq = []
                else:
                    currentseq.append(line)
            if currentname:
                sequences[currentname] = ''.join(currentseq)

        return sequences
    
    def err(self, message):
        '''
        Prints error message and introduces an error.
        Helper to simulateHiFiReads function.
        '''
        print(f"Error: {message}", file=sys.stderr)

    def simulateHiFiReads(self, genomedictref, numreads, errorrate, insertf, deletef, subf, outputformat, 
                      outputname, readlengthmode='lognormal', lengthmean=0, lengthsd=0, 
                      empiricallengths=None, qscoremode='normal', qscoreparams=None):
        """
        Simulates HiFi reads from a given reference genome dictionary.
        Outputs reads in FASTA or FASTQ format with errors and optional biological variation.
        Length distrubution and quality score distribution are modifyable.
        
        Parameters:
            genomedictref (dict): Dictionary of chromosome/contig name → sequence.
            numreads (int): Number of reads to simulate.
            errorrate (float): Overall error rate.
            insertf, deletef, subf (float): Fractions of insertions, deletions, substitutions.
            outputformat (str): 'fasta' or 'fastq'.
            outputname (str): Output filename prefix.
            readlengthmode (str): 'lognormal' or 'empirical'.
            lengthmean (float): Mean read length (for lognormal mode).
            lengthsd (float): Std dev of read length (for lognormal mode).
            empiricallengths (list): List of real read lengths for 'empirical' mode.
            qscoremode (str): 'Normal' or 'random' or 'weighted'.
            qscoreparams (dict): Dictionary of parameters for quality score generation.
                For qscoremode='normal':
                    - 'meanphred' (int): Mean Phred score. (default: 35)
                    - 'stddev' (float): Standard deviation (optional).
                For qscoremode='random':
                    - No parameters required.
                For qscoremode='weighted':
                    - 'phredrange' (list of int): List of possible Phred scores.
                    - 'weights' (list of float): Relative weights matching phredrange.
        """
        
        # error catching - no reference genome given
        if not genomedictref:
            self.err('No reference genome given')
            return

        # extract data from the dictionary of reference sequences or chromasomes
        sequencenames = list(genomedictref.keys())
        sequencelengths = {name: len(seq) for name, seq in genomedictref.items()}
        totalgenomelen = sum(sequencelengths.values())

        # error catching - if genome is empty
        if totalgenomelen == 0:
            self.err('Genome is empty- total genome is length is 0')
            return

        # normalize and check output format 
        outputformat = outputformat.lower()
        if outputformat not in ['fasta', 'fastq']:
            self.err(f"Unsupported output format '{outputformat}'. Must be 'fasta' or 'fastq'.")
            return

        # calculate weight for chance of selection for each sequence in the file
        # purpose: when generating reads from a reference - longer sequences
        # so, longer sequences in the file are more likely to be used/selected based on their length percentage in reference to the whole genome
        sequenceweights = [length / totalgenomelen for length in sequencelengths.values()]

        # output filename
        outputfilename = f'{outputname}.{outputformat}'

        # initialize quality score generator
        if outputformat == 'fastq':
            qscoregen = GenerateQscores()

            if qscoreparams is None:
                qscoreparams = {}
            if qscoremode == 'normal':
                if 'meanphred' not in qscoreparams:
                    qscoreparams['meanphred'] = 35  
            elif qscoremode == 'weighted':
                if 'phredrange' not in qscoreparams or 'weights' not in qscoreparams:
                    self.err("'weighted' mode requires both 'phredrange' and 'weights' in qscoreparams.")
                    return
            elif qscoremode == 'random':
                pass  # no parameters required
            else:
                self.err(f"Unsupported qscoremode '{qscoremode}'. Choose 'normal', 'random', or weighted'")
                return
        
        with open(outputfilename, 'w') as outputfile:   
            # start timer and progress bar interval
            starttime = time.time()
            updateinterval = max(1, numreads // 100)
            # loop through for each read that needs to be created
            writtenreads = 0
            while writtenreads < numreads:
                # create a unique id 
                readID = f'HiFiRead{writtenreads+1}'

                # determine the read length 
                # user choice of lognormal or empirical distribution
                valid = False   # use to remove reads that are too small 
                while not valid:
                    if readlengthmode == 'empirical':
                        if not empiricallengths or len(empiricallengths) == 0:
                            self.err('empiricallengths list is empty or not provided.')
                            return
                        simreadlength = random.choice(empiricallengths)

                    elif readlengthmode == 'lognormal':
                        simreadlength = max(1, int(random.lognormvariate(np.log(lengthmean), lengthsd / lengthmean)))
                    else:
                        self.err(f"Unsupported readlengthmode '{readlengthmode}'")
                        return

                    if simreadlength >= 1000:
                        valid = True
                
                # select reference sequence from reference chromasome/sequences
                refseqname = random.choices(sequencenames, weights=sequenceweights, k=1)[0]
                refseqread = genomedictref[refseqname] 
                refseqlen = len(refseqread)

                # select start position on the reference
                # if simulated read length is longer than reference sequence, reduce length
                if simreadlength >= refseqlen:
                    startpos = 0
                    simreadlength = refseqlen
                # if simulated not longer than reference sequence, random start position
                else: 
                    startpos = random.randint(0, refseqlen - simreadlength)
                
                # assign reference segment to the read ranging from the start until read length
                refsegment = refseqread[startpos: startpos + simreadlength]
                
                # decide if forward (True) or reverse complemented (False)
                simreadseq = refsegment if random.choice([True, False]) else self.reverseComplement(refsegment)

                # introduce errors 
                alteredsimseq = self.introduceErrors(simreadseq, errorrate, insertf, deletef, subf)

                # edge case - sequence deleted: regerate sequence
                if not alteredsimseq:
                    self.err("Read reduced to empty after error introduction.\nRegenerating read.")
                    continue # regerate read without incrimenting i 

                if outputformat == 'fasta':
                    outputfile.write(f'>{readID} simulated HiFi read length: {len(alteredsimseq)}\n')
                    # output in chunks
                    for j in range(0, len(alteredsimseq), 60):
                        outputfile.write(alteredsimseq[j:j+60] + '\n')

                elif outputformat == 'fastq':            
                    # normal distribution
                    if qscoremode == 'normal':
                        qualitystr = qscoregen.normaldist(len(alteredsimseq), **qscoreparams)
                    # uniform random distribution 
                    elif qscoremode == 'random':
                        qualitystr = qscoregen.randomdist(len(alteredsimseq), **qscoreparams)
                    # weighted distribution 
                    elif qscoremode == 'weighted':        
                        qualitystr = qscoregen.weightedchoice(len(alteredsimseq), **qscoreparams)
                    # unsupported mode
                    else:
                        self.err(f"Error: Unsupported qscoremode '{qscoremode}'. Select 'normal', 'random', or 'weighted'.")
                        return
                        
                    # write fasta output to file 
                    outputfile.write(f'@{readID} simulated HiFi read length: {len(alteredsimseq)}\n')
                    outputfile.write(alteredsimseq + "\n")
                    outputfile.write("+\n")
                    outputfile.write(qualitystr + "\n")

                writtenreads += 1 # increment loop after valid read is written to output file

                if writtenreads % updateinterval == 0 or writtenreads == numreads:
                    elapsed = time.time() - starttime
                    percentdone = (writtenreads / numreads) * 100
                    avtimeperread = elapsed / writtenreads
                    estimatedtotal = avtimeperread * numreads
                    timeleft = estimatedtotal - elapsed

                    # progress bar
                    barlength = 40
                    filledlength = int(barlength * percentdone // 100)
                    bar = '=' * filledlength + '-' * (barlength - filledlength)
                    print(f"\r[{bar}] {percentdone:5.1f}% | {writtenreads}/{numreads} reads | "
                    f"Elapsed: {elapsed:5.1f}s | ETA: {timeleft:5.1f}s", end='', flush=True)

            totaltime = time.time() - starttime
            print(f"\nFinished generating {numreads} reads in {totaltime:.2f} seconds.")

        print(f"Successfully generated {numreads} reads in {outputfilename}")
        
    def simulateFromReference(self, referencepath, outputname,
                          outputformat, readstats, lengthstats,
                          qscoremode, qscoreparams, errorrate=0.01,
                          insertf=0.25, deletef=0.20, subf=0.55,
                          snprate=0.001, indelrate=0.0001, 
                          readlengthmode='lognormal', empiricallengths=None):
        """
        Wrapper for simulating synthetic HiFi reads using a reference genome
        and statistics (assumed to come from ReadAnalysis).

        Parameters:
            referencepath (str): Path to reference FASTA/FASTA.gz file.
            outputname (str): Output filename prefix.
            outputformat (str): 'fasta' or 'fastq', passed via command line.
            readstats (dict): Read count stats from ReadAnalysis.
            lengthstats (dict): Read length stats from ReadAnalysis.
            qscoremode (str): 'normal', 'random', or 'weighted'.
            qscoreparams (dict): Parameters for quality score generation.
            errorrate (float): Overall sequencing error rate.
            insertf/deletef/subf (float): Error profile fractions.
            snprate/indelrate (float): Biological variation parameters.
            readlengthmode (str): 'lognormal' or 'empirical'.
            empiricallengths (list): List of real read lengths (optional).
        """

        print("[INFO] Loading reference genome...")
        refgenome = self.loadReference(referencepath)

        print("[INFO] Applying biological variation...")
        mutatedgenome = self.mutateGenome(refgenome, snprate=snprate, indelrate=indelrate)

        # Determine number of reads
        mu = readstats.get('mean', 10000)
        sigma = readstats.get('std', max(1.0, 0.1 * mu))  # avoid 0 std
        min_reads = readstats.get('min', int(0.5 * mu))
        max_reads = readstats.get('max', int(1.5 * mu))

        sampled_reads = int(random.gauss(mu, sigma))
        numreads = max(min_reads, min(sampled_reads, max_reads))

        print(f"[INFO] Number of reads to simulate: {numreads}")

        # Read length stats
        lengthmean = lengthstats.get('mean', 15000)
        lengthsd = lengthstats.get('sd', 2000)

        print(f"[INFO] Simulating {numreads} {outputformat.upper()} reads to '{outputname}.{outputformat}'")
        self.simulateHiFiReads(
            genomedictref=mutatedgenome,
            numreads=numreads,
            errorrate=errorrate,
            insertf=insertf,
            deletef=deletef,
            subf=subf,
            outputformat=outputformat,
            outputname=outputname,
            readlengthmode=readlengthmode,
            lengthmean=lengthmean,
            lengthsd=lengthsd,
            empiricallengths=empiricallengths,
            qscoremode=qscoremode,
            qscoreparams=qscoreparams
        )
