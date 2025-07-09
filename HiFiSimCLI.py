# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# HiFiSimli.py
# Command line interface for HiFiSim

import argparse
import os
import random

from FileGlance import FileGlance
from ReadAnalysis import ReadAnalysis
from SimReads import SimReads
from GenerateQscores import GenerateQScores

def main():
    parser = argparse.ArgumentParser(
        description="Simulate PacBio HiFi reads with error models, read length distributions, and quality scores."
    )

    # === Required arguments ===
    parser.add_argument('--input', required=True, help='Path to reference genome file or folder')
    parser.add_argument('--outputname', required=True, help='Prefix for output file name')
    parser.add_argument('--outputformat', choices=['fasta', 'fastq'], required=True, help='Output format')

    # === Simulation configuration ===
    parser.add_argument('--numreads', type=int, help='Number of reads to simulate (optional)')
    parser.add_argument('--errorrate', type=float, default=0.0, help='Overall error rate (0.0–1.0)')
    parser.add_argument('--insertf', type=float, default=0.25, help='Fraction of insertions')
    parser.add_argument('--deletef', type=float, default=0.20, help='Fraction of deletions')
    parser.add_argument('--subf', type=float, default=0.55, help='Fraction of substitutions')

    # === Read length distribution ===
    parser.add_argument('--readlengthmode', choices=['lognormal', 'empirical'], default='lognormal')
    parser.add_argument('--lengthmean', type=float, default=3000.0, help='Mean read length (lognormal only)')
    parser.add_argument('--lengthsd', type=float, default=500.0, help='Std dev of read length (lognormal only)')

    # === Quality score generation ===
    parser.add_argument('--qscoremode', choices=['normal', 'weighted'], default='normal')
    parser.add_argument('--qscoreparams', nargs='*', help='Parameters for Q-score generation. '
                                                           'For "normal": [mean sd], for "weighted": [score1 weight1 score2 weight2 ...]')

    args = parser.parse_args()

    # === Step 1: Load reference genome ===
    if os.path.isdir(args.input):
        genome_dict = FileGlance(args.input).load_folder()
    else:
        genome_dict = FileGlance(args.input).load_file()

    # === Step 2: Analyze reference data ===
    stats = ReadAnalysis(genome_dict).get_read_statistics()
    empiricallengths = stats.get('read_lengths', None)

    # === Step 3: Set number of reads if not given ===
    if args.numreads is None:
        read_stats = stats.get('readsperfasta', {})
        minr = int(read_stats.get('min', 1000))
        maxr = int(read_stats.get('max', 5000))
        median = int(read_stats.get('median', (minr + maxr) // 2))
        args.numreads = random.randint(minr, maxr)
        print(f"[INFO] Number of reads not specified. Randomly selected: {args.numreads} (range: {minr}-{maxr})")

    # === Step 4: Quality Score Generator (optional object if used within SimReads) ===
    qscore_gen = GenerateQScores(mode=args.qscoremode, params=args.qscoreparams)

    # === Step 5: Simulate HiFi Reads ===
    SimReads().simulateHiFiReads(
        genomedictref=genome_dict,
        numreads=args.numreads,
        errorrate=args.errorrate,
        insertf=args.insertf,
        deletef=args.deletef,
        subf=args.subf,
        outputformat=args.outputformat,
        outputname=args.outputname,
        readlengthmode=args.readlengthmode,
        lengthmean=args.lengthmean,
        lengthsd=args.lengthsd,
        empiricallengths=empiricallengths,
        qscoremode=args.qscoremode,
        qscoreparams=args.qscoreparams
    )


if __name__ == '__main__':
    main()
