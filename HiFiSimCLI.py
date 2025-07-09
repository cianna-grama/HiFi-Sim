# Cianna Grama
# SSRI - Synthetic Sequencing Simulation
# HiFiSimcli.py
# Command line interface for HiFiSim

import argparse
import os
import sys

from SimReads import SimReads
from FileGlance import FileGlance

def parse_qscore_params(args):
    if args.qscoremode == 'normal':
        return {
            'meanphred': int(args.qscoreparams[0]) if args.qscoreparams else 30,
            'stddev': float(args.qscoreparams[1]) if args.qscoreparams and len(args.qscoreparams) > 1 else 5
        }
    elif args.qscoremode == 'weighted':
        return {
            'phredrange': [30, 35, 40],
            'weights': [0.3, 0.4, 0.3]
        }
    else:
        return {}

def main():
    parser = argparse.ArgumentParser(description="Simulate HiFi reads with or without prior read analysis.")

    # Required
    parser.add_argument('--input', required=True, help='Path to reference genome file (.fa or .fa.gz)')
    parser.add_argument('--outputname', required=True, help='Prefix for output file')
    parser.add_argument('--outputformat', choices=['fasta', 'fastq'], required=True, help='Output format')

    # Optional: read analysis
    parser.add_argument('--readfolder', help='Path to folder of real reads to analyze before simulating')

    # Manual simulation inputs
    parser.add_argument('--numreads', type=int, help='Number of reads to simulate')
    parser.add_argument('--errorrate', type=float, default=0.01, help='Total error rate (0.0–1.0)')
    parser.add_argument('--insertf', type=float, default=0.25, help='Fraction of insertions')
    parser.add_argument('--deletef', type=float, default=0.20, help='Fraction of deletions')
    parser.add_argument('--subf', type=float, default=0.55, help='Fraction of substitutions')

    parser.add_argument('--readlengthmode', choices=['lognormal', 'empirical'], default='lognormal')
    parser.add_argument('--lengthmean', type=float, default=3000, help='Mean read length (lognormal mode)')
    parser.add_argument('--lengthsd', type=float, default=500, help='Standard deviation (lognormal mode)')

    parser.add_argument('--qscoremode', choices=['normal', 'random', 'weighted'], default='normal')
    parser.add_argument('--qscoreparams', nargs='*', help='Q-score params: normal=[mean stddev], weighted=[score1 weight1 ...]')

    args = parser.parse_args()
    sim = SimReads()

    # Load reference genome
    if os.path.isdir(args.input):
        print("Error: Input should be a reference genome file (.fa or .fa.gz), not a folder.")
        sys.exit(1)
    genome_dict = sim.loadReference(args.input)

    print("[INFO] Previewing reference genome:")
    FileGlance(args.input)

    # Read analysis path is optional
    if args.readfolder:
        print(f"[INFO] Running ReadAnalysis on folder: {args.readfolder}")
        from ReadAnalysis import ReadAnalysis
        analysis = ReadAnalysis(args.readfolder)
        stats = analysis.get_read_statistics()

        readstats = stats.get('readsperfasta', {})
        lengthstats = stats.get('lengthstats', {})
        empiricallengths = stats.get('read_lengths', None)

        sim.simulateFromReference(
            referencepath=args.input,
            outputname=args.outputname,
            outputformat=args.outputformat,
            readstats=readstats,
            lengthstats=lengthstats,
            qscoremode=args.qscoremode,
            qscoreparams=parse_qscore_params(args),
            errorrate=args.errorrate,
            insertf=args.insertf,
            deletef=args.deletef,
            subf=args.subf,
            snprate=0.001,
            indelrate=0.0001,
            readlengthmode=args.readlengthmode,
            empiricallengths=empiricallengths
        )

    else:
        sim.simulateHiFiReads(
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
            empiricallengths=None,
            qscoremode=args.qscoremode,
            qscoreparams=parse_qscore_params(args)
        )

if __name__ == '__main__':
    main()
