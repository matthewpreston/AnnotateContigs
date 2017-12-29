#!/usr/bin/env python2
#
# AnnotateContigs.py - Annotation consensus assembler. Takes aligned contigs
# from a BLAST XML file and reassembles them into a set of annotated consensus
# sequences. Queries NCBI for coding locations in order to create spliced
# transcripts containing solely the CDS. Can provide a reference FASTA file that
# acted as the BLAST database for faster, offline usage.
#
# Four modes (and combinations of the four):
#   a OR A      Outputs aligned contigs
#   c OR C      Outputs a consensus sequence derived from aligned contigs
#   s OR S      Outputs a spliced consensus sequence derived from aligned
#               contigs
#   d OR D      Provides statistical data (gene length, breadth of coverage,
#               etc.)
#
# Example usage:
# python AnnotateContigs.py -i 0.7 -s 0.7 -n 4 blast.xml csd
# 
# Written By: Matt Preston and Xiaotong Yang
# Created on: May 16, 2016      V1.0
# Revised on: May 17, 2016      V1.1 Added splicing and mapping locations
#             May 22, 2016      V1.2 Added protein records
#             May 26, 2016      V1.3 Added modes, gene assembly, and splicing
#             May 30, 2016      V1.4 Split into modules, added offline mode
#             Jun 22, 2016      V1.5 Added stats for breadth and depth of cov.
#             Jul 26, 2016      V1.6 Upgraded consensus building
#             Sep 21, 2016      V1.7 Added identity and similarity thresholds
#             Jun  1, 2017      V1.8 Supported multithreading, complete overhaul
#             Jun 21, 2017      V1.8.1 Bug fixes (key=alignment.title not gene)
#             Jul  6, 2017      V1.8.2 Added % UTR and intron in stats file
#             Aug 11, 2017      V1.8.3 Fixed race conditions with threads

import argparse
import logging
import sys
from Bio import Entrez, SeqIO
from AnnotateContigs.Alignments import Alignments
from AnnotateContigs.TerminalCommands import CreateDir

class CheckFileAction(argparse.Action):
    """Ensures that file exists and is OK to read/write"""
    
    def __call__(self, parser, namespace, values, option_string=None):
        # It is actually better to open a file than it is to check existance,
        # else it could lead to annoying bugs (Why can't it find my file?!?!?
        # It's clearly right there you stupid program!)
        try:
            handle = open(values)
        except IOError as e:
            parser.error(e)
        else:
            handle.close()
        setattr(namespace, self.dest, values)
        
class CheckModeAction(argparse.Action):
    """Ensures that mode is in proper format"""
    
    def __call__(self, parser, namespace, values, option_string=None):
        if not set(values) & set(['a','A','c','C','s','S','d','D']):
            parser.error("Invalid mode: {}; expected at least one character in"
                         "'a','A','c','C','s','S','d','D'".format(values))
        setattr(namespace, self.dest, values)

def BoundedActionFactory(lowerBound, upperBound, lowerCompChar='(',
                         upperCompChar=')'):
    """
    Creates BoundedAction classes with desired bounds and bound checking
    functions
    """
                              
    class BoundedAction(argparse.Action):
        """Ensures that an option that is to be a float is within bounds"""
        
        def __init__(self, option_strings, dest, **kwargs):
            super(BoundedAction, self).__init__(option_strings, dest,
                                                     **kwargs)
            self.lowerBound = lowerBound
            self.upperBound = upperBound
            self.lowerCompChar = lowerCompChar
            self.upperCompChar = upperCompChar
            
        def __call__(self, parser, namespace, values, option_string=None):
            # Get comparison operators (raises exception if invalid character)
            lowerComp = BoundedAction._GetCompOp(lowerCompChar)
            upperComp = BoundedAction._GetCompOp(upperCompChar)
            if not BoundedAction._WithinBounds(values, self.lowerBound,
                                               self.upperBound,
                                               lowerComp, upperComp):
                parser.error("--{flag} must be within "
                             "{lowerCompChar}{lowerBound}, "
                             "{upperBound}{upperCompChar}".format(
                                flag          = self.dest,
                                lowerCompChar = self.lowerCompChar,
                                lowerBound    = self.lowerBound,
                                upperBound    = self.upperBound,
                                upperCompChar = self.upperCompChar))
            setattr(namespace, self.dest, values)
            
        @staticmethod
        def _GetCompOp(char):
            """Sees if char is valid and returns appropriate function"""
            
            if char == '[':
                return lambda n,lb: n >= lb
            elif char == '(':
                return lambda n,lb: n > lb
            elif char == ')':
                return lambda n,ub: n < ub
            elif char == ']':
                return lambda n,ub: n <= ub
            else:
                raise Exception("Unexpected char {}, should be either '[', "
                                "'(', ')', or ']'")
            
        @staticmethod
        def _WithinBounds(num, lowerBound, upperBound, lowerComp=(lambda n,lb: n > lb),
                 upperComp=(lambda n,ub: n < ub)):
            """Returns true if num is within bounds"""
            
            return lowerComp(num, lowerBound) and upperComp(num, upperBound)
    return BoundedAction

def main(args):
    Entrez.email = args.email
    # Set global logging
    logging.basicConfig(level=getattr(logging, args.verbosity),
                        format="%(asctime)s (%(threadName)-10s) "
                               "(%(name)s)-%(levelname)s: %(message)s")
    logger = logging.getLogger(__name__)
    logger.info("Start")
                                        
    # If a query file was given, store them in a dictionary where the keys are
    # the gene names
    queryRecords = None
    if args.reference is not None:
        referenceRecords = {}
        logger.info("Fetching query records")
        with open(args.reference, 'r') as referenceHandle:
            for record in SeqIO.parse(referenceHandle, "fasta"):
                logger.debug("Adding record: %s", record.description)
                referenceRecords[record.description] = record
        logger.info("Done fetching")
                                        
    # If a reference was given, store them in a dictionary where the keys are
    # the gene names
    referenceRecords = None
    if args.reference is not None:
        referenceRecords = {}
        logger.info("Fetching reference records")
        with open(args.reference, 'r') as referenceHandle:
            for record in SeqIO.parse(referenceHandle, "fasta"):
                logger.debug("Adding record: %s", record.description)
                referenceRecords[record.description] = record
        logger.info("Done fetching")

    alignments = Alignments(queryRecords, referenceRecords, args.identity, 
                            args.similarity, args.threshold)
       
    # Parse through BLAST XML
    logger.info("Parsing through BLAST XML file")
    alignments.ParseBlastOutput(args.blast_file, args.threads)
    logger.info("Number of genes found: %d", len(alignments))
       
    # See how many contigs are associated with each gene (for testing)
    for key in alignments.iterkeys():
        for trans in alignments[key].alternativeTranscripts:
            logger.debug("Key: %s; NumContigs: %s", key, len(trans.contigs))
    
    # Output to crap
    CreateDir(args.output_dir, True)
    if 'a' in args.mode or 'A' in args.mode:
        logger.info("Outputting alignments to: %s", args.alignment_dir)
        CreateDir(args.alignment_dir, True)
        alignments.OutputAlignments(args.threads, args.alignment_dir,
                                    minCoverage=args.coverage,
                                    maxTranscripts=args.protein_num)
    if 'c' in args.mode or 'C' in args.mode:
        logger.info("Outputting contigs to: %s", args.consensus_dir)
        CreateDir(args.consensus_dir, True)
        alignments.OutputConsenses(args.threads, args.consensus_dir,
                                   minCoverage=args.coverage,
                                   maxTranscripts=args.protein_num)
    if 's' in args.mode or 'S' in args.mode:
        logger.info("Outputting spliced contigs to: %s", args.spliced_dir)
        CreateDir(args.spliced_dir, True)
        alignments.OutputSplicedConsenses(args.threads, args.spliced_dir,
                                          minCoverage=args.coverage,
                                          maxTranscripts=args.protein_num)
    if 'd' in args.mode or 'D' in args.mode:
        logger.info("Outputting statistics to: %s", args.stats)
        alignments.OutputStatistics(args.stats, args.coverage)

    logger.info("End")
    sys.exit(0)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Add positional arguments
    parser.add_argument("blast_file",
                        metavar="<blast.xml>",
                        action=CheckFileAction,
                        help="""XML file containing BLAST output""")
    parser.add_argument("mode",
                        metavar="<MODE>",
                        action=CheckModeAction,
                        help="""Use 'a'/'A' for outputting aligned contigs, 
                            'c'/'C' for consensus sequence derived from aligned 
                            contigs, 's'/'S' for a spliced consensus sequence 
                            derived from aligned contigs, and/or 'd'/'D' for
                            statistical data""")
    # Add optional arguments
    parser.add_argument("-q", "--query",
                        metavar="FILE",
                        action=CheckFileAction,
                        help="""If possible, provide the FASTA file containing 
                            the query sequences used for BLASTing to the 
                            reference. This will be used to generate more 
                            statistics (%% 5' and 3' UTR, %% intronic) 
                            [None]""")
    parser.add_argument("-r", "--reference",
                        metavar="FILE",
                        action=CheckFileAction,
                        help="""Instead of querying NCBI for references, 
                            provide a reference FASTA file for offline (and 
                            faster) use.\nNote: must have been used as the 
                            BLAST database when creating the BLAST XML file
                            [None]""")
    parser.add_argument("-e", "--email",
                        metavar="EMAIL",
                        default="matthew.preston@mail.utoronto.ca",
                        help="""For NCBI to contact if a problem arises with 
                            the spamming of their servers [Mine dammit]""")
    parser.add_argument("-n", "--threads",
                        metavar="INT",
                        type=int,
                        default=1,
                        help="""Number of threads to use [1]""")
    parser.add_argument("-o", "--output_dir",
                        metavar="DIR",
                        default="AnnotateContigs",
                        help="""Directory to store all output data 
                            [AnnotateContigs/]""")
    parser.add_argument("-A", "--alignment_dir",
                        metavar="DIR",
                        default="{output_dir}/Alignments",
                        help="""Will output alignments to DIR if MODE contains 
                            \"a\" or \"A\" [$output_dir/Alignments/]""")
    parser.add_argument("-C", "--consensus_dir",
                        metavar="DIR",
                        default="{output_dir}/Consenses",
                        help="""Will output consensus sequences to DIR if MODE 
                            contains \/"c\" or \"C\" [$output_dir/Consenses/]
                            """)
    parser.add_argument("-S", "--spliced_dir",
                        metavar="DIR",
                        default="{output_dir}/SplicedConsenses",
                        help="""Will output spliced consensus sequences to DIR 
                            if MODE contains \"s\" or \"S\" 
                            [$output_dir/SplicedConsenses/]""")
    parser.add_argument("-d", "--stats",
                        metavar="FILE",
                        default="{output_dir}/Summary_Results.csv",
                        help="""Name of file to dump statistical output (i.e. 
                            length, percent coverage, breadth of coverage) to 
                            FILE if MODE to contain \"d\" or \"D\" 
                            [$output_dir/Summary_Results.csv]""")
    parser.add_argument("-c", "--coverage",
                        metavar="FLOAT",
                        type=float,
                        default="0.0",
                        action=BoundedActionFactory(0.0, 1.0, '[', ']'),
                        help="""Minimum coverage threshold in order to output 
                            alignment, must be within [0.0, 1.0] [0.0]""")
    parser.add_argument("-i", "--identity",
                        metavar="FLOAT",
                        type=float,
                        default="0.0",
                        action=BoundedActionFactory(0.0, 1.0, '[', ']'),
                        help="""Lower threshold to keep contiguous sequences 
                            sufficient sequence identity (exact matches), 
                            must be within [0.0, 1.0] [0.0]""")
    parser.add_argument("-s", "--similarity",
                        metavar="FLOAT",
                        type=float,
                        default="0.0",
                        action=BoundedActionFactory(0.0, 1.0, '[', ']'),
                        help="""Lower threshold to keep contiguous sequences 
                            with sufficient sequence similarity (similar 
                            matches), must be within [0.0, 1.0] [0.0]""")
    parser.add_argument("-t", "--threshold",
                        metavar="FLOAT",
                        type=float,
                        default="0.7",
                        action=BoundedActionFactory(0.5, 1.0, '(', ']'),
                        help="""Threshold needed to exceed for creating a 
                            consensus nucleotide, must be within (0.5, 1.0] 
                            [0.7]""")
    parser.add_argument("-p", "--protein_num",
                        metavar="INT",
                        type=float,
                        default=float("inf"),
                        action=BoundedActionFactory(0.0, float("inf"), '(', ')'),
                        help="""Threshold needed to exceed for creating a 
                            consensus nucleotide, must be greater than zero 
                            [all]""")  
    parser.add_argument("-v", "--verbosity",
                        metavar="LEVEL",
                        nargs='?',
                        type=lambda l: l.upper(),
                        const="INFO",
                        default="ERROR",
                        choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],
                        help="""Verbosity of output. LEVEL can be either DEBUG,
                            INFO, WARNING, ERROR, or CRITICAL [ERROR] (if flag 
                            raised but no level given: INFO)""")
    args = parser.parse_args()
    # Augment directory namings
    args.alignment_dir = args.alignment_dir.format(output_dir=args.output_dir)
    args.consensus_dir   = args.consensus_dir.format(output_dir=args.output_dir)
    args.spliced_dir   = args.spliced_dir.format(output_dir=args.output_dir)
    args.stats         = args.stats.format(output_dir=args.output_dir)
    main(args)