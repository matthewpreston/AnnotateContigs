#!/usr/bin/env python2
#
# Alignments.py - An object that parses through BLAST XML output derived from
# NCBI's nt database and outputs the aligned segments, as well as creates a
# tentative consensus sequence and splices to create alternative transcripts.
# Since the nt database is used, this object will automatically attempt to pull
# data from NCBI in order to retrieve splicing information about genes using
# their accession numbers. Can also utilize a given reference FASTA file to
# avoid the use (and bandwidth costs) of querying NCBI.
#
# Written By: Matt Preston
# Written On: May 31, 2017
# Revised On: Never

from functools import wraps
import logging
import Queue
import threading
import time
import sys

from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

from AlignmentParser import ThreadSafeAlignmentParser
from Consensus import ConsensusSeq
from Exceptions import GeneNotFoundException, FailedToParseException
from RecordObjects import AlignmentRecord
from SequenceOperations import _fix_location, CodingLocationToList, \
                               MappedDifference, Splice
from TerminalCommands import CreateDir
from ThreadPool import ThreadPool

def identity(seqA, seqB):
    """Calculates the percent sequence identity between two sequences"""
    if len(seqA) != len(seqB):
        raise ValueError("Sequences are of different length")
    return float(sum(a == b for a,b in zip(seqA, seqB))) * 100 / len(seqA)

class Alignments:
    """Stores alignments common to one gene."""
    
    def __init__(self, queryRecords=None, referenceRecords=None, identity=0.0,
                 similarity=0.0, threshold=0.7):
        self._threadPool = ThreadPool()
        self._alignments = {}
        self.queryRecords = queryRecords
        self.referenceRecords = referenceRecords
        self.identity = identity
        self.similarity = similarity
        self.threshold = threshold
        self.logger = logging.getLogger(self.__class__.__name__)

    def __getitem__(self, key):
        return self._alignments[key]
        
    def __len__(self):
        return len(self._alignments)
        
    @staticmethod
    def _CreateFeatureLocation(locations):
        """
        Creates a CompoundLocation object, which is made of a list of 
        FeatureLocation objects for ease of output
        """
        return sum([FeatureLocation(l[0], l[1]) for l in locations])
        
    @staticmethod
    def _GetNonCodingBases(mapped, location):
        """
        Gets the 5' UTR, 3' UTR, and intronic bases based on the extent of
        contig mapping versus splicing locations
        """
        
        # Find regions that are not translated into protein by removing exons
        tempLocation = CodingLocationToList(_fix_location(location)) # Tuple-ize
        diff = MappedDifference(mapped, tempLocation)
        # The first tuple (if it exists) is either 5' UTR or an intron
        if len(diff) and diff[0][1] < tempLocation[0][0]:
            temp = diff.pop(0)
            UTR5 = temp[1] - temp[0] + 1
        else:
            UTR5 = 0
        # The last tuple (if it exists) is either 3' UTR or an intron
        if len(diff) and diff[-1][0] > tempLocation[-1][1]:
            temp = diff.pop(-1)
            UTR3 = temp[1] - temp[0] + 1
        else:
            UTR3 = 0
        # The remaining tuples are now mapping intronic bases
        intronic = sum(m[1]-m[0]+1 for m in diff)
        return (UTR5, UTR3, intronic)
        
    def _CheckCriteria(self, altTrans, minCoverage, name=None):
        """
        Ensures that this AlternativeTranscriptRecord is okay to output, i.e.
        has contigs associated with it, contigs that map to it, and has suitable
        coverage
        """
        
        if not altTrans.contigs:
            self.logger.debug("No contigs for %s",
                              altTrans.proteinID if name is None else name)
            return False
            
        if not altTrans.MappingToSpliced():
            self.logger.debug("No mapping for %s",
                              altTrans.proteinID if name is None else name)
            return False
            
        self.logger.debug("Coverage (%.6f) vs minimum (%.6f) for %s", 
                          altTrans.Coverage(), minCoverage,
                          altTrans.proteinID if name is None else name)
        return altTrans.Coverage() >= minCoverage
    
    def _SpawnThreads(self, threads):
        """Spawns the appropriate number of new threads in the thread pool"""
        
        # Using our new custom thread pool class, spawn the appropriate number
        # of threads
        threads -= 1 # Main thread will be utilized as well, don't be wasteful
        if threads > self._threadPool.numThreads:
            threadNames = []
            for i in range(self._threadPool.numThreads, threads):
                threadNames.append("Thread-%i" % (i+1))
            self._threadPool.AddThreads(threads - self._threadPool.numThreads,
                                        threadNames)
     
    # Basic methods
    def keys(self):
        return self._alignments.keys()
        
    def values(self):
        return self._alignments.values()
        
    def items(self):
        return self._alignments.items()
                  
    def iterkeys(self):
        for k in self._alignments.iterkeys():
            yield k

    def itervalues(self):
        for v in self._alignments.itervalues():
            yield v
            
    def iteritems(self):
        for k, v in self._alignments.iteritems():
            yield (k, v)

    def ParseBlastOutput(self, blastFile, threads=1):
        """Parses BLAST XML alignments & adds them to its internal dictonary"""

        # Spawn the appropriate number of threads
        self._SpawnThreads(threads)
        
        # Create an alignment parsing object which is going to be abused of its
        # method to parse through alignments (so descriptive!)
        ap = ThreadSafeAlignmentParser(name="%s.ThreadSafeAlignmentParser"
                                            % self.__class__.__name__)
        
        # Fill up queue
        self.logger.debug("Parsing %s", blastFile)
        with open(blastFile, 'r') as blastHandle:
            for record in NCBIXML.parse(blastHandle):
                for alignment in record.alignments:
                    self._threadPool.AddTask(ap.ParseAlignment,
                            alignment        = alignment,
                            targetDict       = self._alignments,
                            referenceRecords = self.referenceRecords,
                            identity         = self.identity,
                            similarity       = self.similarity)

        # Also have the main thread parse through the alignments
        self.logger.debug("Main thread parsing through alignments")
        while not self._threadPool.tasks.empty():
            try:
                func, args, kwargs = self._threadPool.tasks.get(False)
                func(*args, **kwargs)
                self._threadPool.tasks.task_done()
            except Queue.Empty as e:
                pass
            except Exception as e:
                self.logger.error(e, exc_info=True)
                self._threadPool.tasks.task_done()
                
        # Wait for threads to finish
        self._threadPool.Join()
        self.logger.debug("Finished parsing through alignments")
        
    def OutputAlignment(self, key, outputTitle=None, format="fasta",
                        outputDir=".", minCoverage=0.0,
                        maxTranscripts=float("inf")):
        """Writes out the multiple alignments aligned to a desired gene"""
        
        record = self._alignments.get(key)
        if record is None:
            raise GeneNotFoundException(
                "No alignments pertaining to: %s" % key)
        
        # Get the desired number of longest transcripts
        name = record.gene if self.referenceRecords is not None else None
        geneLength = len(record.referenceSequence)
        count = 0
        for trans in record.GetLongestTranscripts(maxTranscripts):
            if count > maxTranscripts:
                return
        
            # Ensure that transcript is suitable for output
            if not self._CheckCriteria(trans, minCoverage, name):
                continue
            
            # Formulate output file name
            if outputTitle is None:
                if self.referenceRecords is None:
                    tempTitle = "%s/GeneName_%s_Accession_%s_ProteinID_%s" \
                                % (outputDir, record.gene, record.accession,
                                   trans.proteinID)
                else:
                    tempTitle = "%s/%s" % (outputDir, record.gene)
            else:
                tempTitle = "%s/%s" % (outputDir, outputTitle)
                         
            # Coalesce all contigs associated with this protein
            ID = tempTitle.split("/")[-1]
            records = []
            for i, contig in enumerate(trans.contigs):
                contigSeq = contig.GetPaddedContig(geneLength)
                dscrpt = str(FeatureLocation(contig.start, contig.end, 1))
                temp = SeqRecord(Seq(contigSeq, generic_dna),
                                 id="%s_Contig_%i" % (ID, i+1),
                                 description=dscrpt)
                records.append(temp)
                
            # Output
            self.logger.debug("Outputting gene: %s", key)
            SeqIO.write(records, "%s.%s" % (tempTitle, format), format)
            count += 1
    
    def OutputAlignments(self, threads=1, outputDirectory=".", format="fasta",
                         minCoverage=0.0, maxTranscripts=float("inf")):
        """Writes all multiple alignments out into a desination folder"""
        
        # Spawn the appropriate number of threads
        self._SpawnThreads(threads)
        
        # Fill up queue
        self.logger.debug("Adding genes to output their respective alignments")
        for key in self._alignments.iterkeys():
            self._threadPool.AddTask(self.OutputAlignment,
                    key=key,
                    outputTitle=None,
                    format=format,
                    outputDir=outputDirectory,
                    minCoverage=minCoverage,
                    maxTranscripts=maxTranscripts)
                    
        # Also have the main thread output alignments
        self.logger.debug("Main thread outputting alignments")
        while not self._threadPool.tasks.empty():
            try:
                func, args, kwargs = self._threadPool.tasks.get(False)
                func(*args, **kwargs)
                self._threadPool.tasks.task_done()
            except Queue.Empty as e:
                pass
            except Exception as e:
                self.logger.error(e, exc_info=True)
                self._threadPool.tasks.task_done()
                
        # Wait for threads to finish
        self._threadPool.Join()
        self.logger.debug("Finished outputting alignments")
    
    def OutputConsensus(self, key, outputTitle=None, format="fasta",
                        outputDir=".", minCoverage=0.0,
                        maxTranscripts=float("inf")):
        """
        Writes out the consensus sequence created by merging the aligned contigs
        of a desired gene
        
        Note: pretty sure this is 'cargo cult science' because it doesn't even
        make sense. I would have to inspect how each transcriptome was assembled
        i.e. look at what de novo assembler was used, and then dicern how to
        create a consensus sequence. For example, Trinity uses a De Brujin graph
        which it later splits to create 'isoforms'. When BLAST aligns these
        'isoforms' and gives the alignment file to this program, I would
        essentially be double counting (triple counting, etc.) the aligned
        sections of contigs. So for me, TODO make the Consensus.py library
        decide how to create a consensus sequence based on de novo assembler
        instead of naively counting and averaging!
        """
        
        record = self._alignments.get(key)
        if record is None:
            raise GeneNotFoundException(
                "No alignments pertaining to: %s" % key)
        
        # Get the desired number of longest transcripts
        name = record.gene if self.referenceRecords is not None else None
        geneLength = len(record.referenceSequence)
        count = 0
        for trans in record.GetLongestTranscripts(maxTranscripts):
            if count > maxTranscripts:
                return
        
            # Ensure that transcript is suitable for output
            if not self._CheckCriteria(trans, minCoverage, name):
                continue
            
            # Formulate output file name
            if outputTitle is None:
                if self.referenceRecords is None:
                    tempTitle = "%s/GeneName_%s_Accession_%s_ProteinID_%s" \
                                % (outputDir, record.gene, record.accession,
                                   trans.proteinID)
                else:
                    tempTitle = "%s/%s" % (outputDir, record.gene)
            else:
                tempTitle = "%s/%s" % (outputDir, outputTitle)
                         
            # Create a consensus sequence
            consensus = trans.Consensus(geneLength, self.threshold)
            ID = tempTitle.split("/")[-1]
            dscrpt = str(self._CreateFeatureLocation(trans.mapped))
            records = [SeqRecord(Seq(consensus, generic_dna),
                                 id="%s_Consensus" % ID,
                                 description=dscrpt)]
                
            # Output
            self.logger.debug("Outputting key: %s", key)
            SeqIO.write(records, "%s.%s" % (tempTitle, format), format)
            count += 1
    
    def OutputConsenses(self, threads=1, outputDirectory=".", format="fasta",
                        minCoverage=0.0, maxTranscripts=float("inf")):
        """Write all consensus sequences out into a destination folder"""
        
        # Spawn the appropriate number of threads
        self._SpawnThreads(threads)
        
        # Fill up queue
        self.logger.debug("Adding genes to output their respective alignments")
        for key in self._alignments.iterkeys():
            self._threadPool.AddTask(self.OutputConsensus,
                    key=key,
                    outputTitle=None,
                    format=format,
                    outputDir=outputDirectory,
                    minCoverage=minCoverage,
                    maxTranscripts=maxTranscripts)
                    
        # Also have the main thread output alignments
        self.logger.debug("Main thread outputting consenses")
        while not self._threadPool.tasks.empty():
            try:
                func, args, kwargs = self._threadPool.tasks.get(False)
                func(*args, **kwargs)
                self._threadPool.tasks.task_done()
            except Queue.Empty as e:
                pass
            except Exception as e:
                self.logger.error(e, exc_info=True)
                self._threadPool.tasks.task_done()
                
        # Wait for threads to finish
        self._threadPool.Join()
        self.logger.debug("Finished outputting consenses")
    
    def OutputSplicedConsensus(self, key, outputTitle=None, format="fasta",
                               outputDir=".", minCoverage=0.0,
                               maxTranscripts=float("inf")):
        """Write out the spliced consensus transcript of a desired gene"""
        
        record = self._alignments.get(key)
        if record is None:
            raise GeneNotFoundException(
                "No alignments pertaining to: %s" % key)
        
        # Get the desired number of longest transcripts
        name = record.gene if self.referenceRecords is not None else None
        geneLength = len(record.referenceSequence)
        count = 0
        for trans in record.GetLongestTranscripts(maxTranscripts):
            if count > maxTranscripts:
                return
        
            # Ensure that transcript is suitable for output
            if not self._CheckCriteria(trans, minCoverage, name):
                continue
            
            # Formulate output file name
            if outputTitle is None:
                if self.referenceRecords is None:
                    tempTitle = "%s/GeneName_%s_Accession_%s_ProteinID_%s" \
                                % (outputDir, record.gene, record.accession,
                                   trans.proteinID)
                else:
                    tempTitle = "%s/%s" % (outputDir, record.gene)
            else:
                tempTitle = "%s/%s" % (outputDir, outputTitle)
                         
            # Create a spliced consensus sequence
            splicedConsensus = trans.SplicedConsensus(geneLength,self.threshold)
            ID = tempTitle.split("/")[-1]
            dscrpt = str(self._CreateFeatureLocation(trans.MappingToSpliced()))
            records = [SeqRecord(Seq(splicedConsensus, generic_dna),
                                 id="%s_Spliced_Consensus" % ID,
                                 description=dscrpt)]
                
            # Output
            self.logger.debug("Outputting key: %s", key)
            SeqIO.write(records, "%s.%s" % (tempTitle, format), format)
            count += 1
    
    def OutputSplicedConsenses(self, threads=1, outputDirectory=".",
                               format="fasta", minCoverage=0.0,
                               maxTranscripts=float("inf")):
        """
        Write all spliced consensus transcripts all into a destination folder
        """
        
        # Spawn the appropriate number of threads
        self._SpawnThreads(threads)
        
        # Fill up queue
        self.logger.debug("Adding genes to output their respective alignments")
        for key in self._alignments.iterkeys():
            self._threadPool.AddTask(self.OutputSplicedConsensus,
                    key=key,
                    outputTitle=None,
                    format=format,
                    outputDir=outputDirectory,
                    minCoverage=minCoverage,
                    maxTranscripts=maxTranscripts)
                    
        # Also have the main thread output alignments
        self.logger.debug("Main thread outputting spliced consenses")
        while not self._threadPool.tasks.empty():
            try:
                func, args, kwargs = self._threadPool.tasks.get(False)
                func(*args, **kwargs)
                self._threadPool.tasks.task_done()
            except Queue.Empty as e:
                pass
            except Exception as e:
                self.logger.error(e, exc_info=True)
                self._threadPool.tasks.task_done()
                
        # Wait for threads to finish
        self._threadPool.Join()
        self.logger.debug("Finished outputting spliced consenses")

    def OutputStatistics(self, outputFile, minCoverage=0.0,
                         maxTranscripts=float("inf")):
        """
        Outputs a variety of statistics to a file (breadth of coverage, number
        of aligned contigs, percent coding and non-coding, etc.)"""
        
        # Sort by gene name, then by accession number
        temp = [(k,v.gene,v.accession) for k,v in self._alignments.iteritems()]
        temp.sort(key=(lambda (k,g,a): (g.lower(),a,k)))
        sortedKeys = [k for k,g,a in temp]
        numGenes = 0
        geneLengths = 0.0
        totalContigs = 0.0
        transcriptLengths = 0.0
        totalBasesMapped = 0.0
        totalBreadthOfCoverage = 0.0
        totalUTR5 = 0.0
        totalUTR3 = 0.0
        totalIntronic = 0.0
        totalExonic = 0.0
        totalCDSBases = 0.0
        totalCDSCoverage = 0.0
        totalIdentities = 0.0
        transcripts = 0

        with open(outputFile, 'w') as o:
            if self.referenceRecords is None:
                o.write("Gene,Accession.Version,Gene Length,Number of Contigs,"
                        "Protein ID,Transcript Length,Bases Mapped,"
                        "Breadth of Coverage (%),5' UTR (%),3' UTR (%),"
                        "Intronic (%),Exonic (%),CDS Bases Mapped,"
                        "CDS Coverage (%),Identity (%)\n")
            else:
                o.write("Gene,Gene Length,Number of Contigs,Bases Mapped,"
                        "Breadth of Coverage (%),5' UTR (%),3' UTR (%),"
                        "Intronic (%),Exonic (%),CDS Bases Mapped,"
                        "CDS Coverage (%),Identity (%)\n")
                        
            if not sortedKeys:
                o.write("No genes found\n")
                return
                
            for key in sortedKeys:
                record = self._alignments[key]
                geneLength = len(record.referenceSequence)
                transcriptsProcessed = 0
                name =record.gene if self.referenceRecords is not None else None
                for trans in record.GetLongestTranscripts(maxTranscripts):
                    # Ensure that transcript is suitable for output
                    if not self._CheckCriteria(trans, minCoverage, name):
                        continue

                    numContigs = len(trans.contigs)
                    if self.referenceRecords is None:
                        # Write gene name and its length
                        o.write("%s,%s.%s,%i,%i," \
                                % (record.gene, record.accession,
                                   record.version, geneLength, numContigs))

                        # Write protein ID and its length
                        transcriptLength = len(trans.location)
                        transcriptLengths += transcriptLength
                        o.write("%s,%i," % (trans.proteinID, transcriptLength))
                    else:
                        # Write gene name and its length
                        o.write("%s,%i,%i," % (record.gene, geneLength,
                                               numContigs))
                        transcriptLength = geneLength

                    # Get the number of bases mapped to this gene
                    bm = trans.BasesMapped()
                    totalBasesMapped += bm
                    o.write("%i," % bm)

                    # Get the breadth of coverage of this gene
                    boc = float(bm) * 100 / transcriptLength
                    totalBreadthOfCoverage += boc
                    o.write("%f," % boc)
                    
                    # Get 5' UTR, 3' UTR, and intronic percentages
                    UTR5, UTR3, intronic = self._GetNonCodingBases(trans.mapped,
                            trans.location)
                    UTR5 = float(UTR5) * 100 / bm
                    UTR3 = float(UTR3) * 100 / bm
                    intronic = float(intronic) * 100 / bm
                    exonic = 100 - UTR5 - UTR3 - intronic
                    totalUTR5 += UTR5
                    totalUTR3 += UTR3
                    totalIntronic += intronic
                    totalExonic += exonic
                    o.write("%f,%f,%f,%f," % (UTR5,UTR3,intronic,exonic))
                    
                    # Get number of bases in CDS
                    cds = trans.Breadth()
                    totalCDSBases += cds
                    o.write("%i," % cds)
                    
                    # Get the CDS coverage of this gene
                    cdsCov = float(cds) * 100 / transcriptLength
                    totalCDSCoverage += cdsCov
                    o.write("%f," % cdsCov)
                    
                    # Get the percent identity between spliced reference
                    # and spliced consensus
                    ident = identity(Splice(record.referenceSequence,
                                            trans.location),
                                     trans.SplicedConsensus(geneLength,
                                                            self.threshold))
                    totalIdentities += ident
                    o.write("%f\n" % ident)
                    
                    transcriptsProcessed += 1
                if transcriptsProcessed:
                    numGenes += 1
                    geneLengths += geneLength
                    totalContigs += numContigs
                    transcripts += transcriptsProcessed

            # Calculate averages
            if numGenes:
                o.write("\nAverage,")
                if self.referenceRecords is None:
                    o.write(",")
                    o.write("%f," % (geneLengths / numGenes))
                    o.write("%f," % (totalContigs / numGenes))
                    o.write(",")
                    o.write("%f,"% (transcriptLengths / transcripts))
                else:
                    o.write("%f," % (geneLengths / numGenes))
                    o.write("%f," % (totalContigs / numGenes))
                o.write("%f," % (totalBasesMapped / transcripts))
                o.write("%f," % (totalBreadthOfCoverage / transcripts))
                o.write("%f," % (totalUTR5 / transcripts))
                o.write("%f," % (totalUTR3 / transcripts))
                o.write("%f," % (totalIntronic / transcripts))
                o.write("%f," % (totalExonic / transcripts))
                o.write("%f," % (totalCDSBases / transcripts))
                o.write("%f," % (totalCDSCoverage / transcripts))
                o.write("%f," % (totalIdentities / transcripts))
        
# Temporary main for testing
def main(argv):
    if len(argv) < 1:
        sys.stderr.write("USAGE: python blah.py <blast.xml> [numThreads=2] "
                         "[reference.fasta] [outDir]\n")
        sys.exit(-1)
    
    blastFile = argv[0]
    if len(argv) > 1:
        numThreads = int(argv[1])
    else:
        numThreads = 2
    if len(argv) > 2:
        referenceFile = argv[2]
    else:
        referenceFile = None
    if len(argv) > 3:
        outFolder = argv[3]
    else:
        outFolder = "."
    
    Entrez.email = "matthew.preston@mail.utoronto.ca"
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s (%(threadName)-10s) "
                               "(%(name)s)-%(levelname)s: %(message)s")
    logger = logging.getLogger(__name__)
    logger.info("Start")
    
    # If a reference was given, store them in a dictionary where the keys
    # the gene names
    referenceRecords = None
    referenceFile = None # Testing for internet, comment out for reference
    if referenceFile is not None:
        referenceRecords = {}
        logger.debug("Fetching reference records")
        with open(referenceFile, 'r') as referenceHandle:
            for record in SeqIO.parse(referenceHandle, "fasta"):
                logger.debug("Adding record: %s", record.description)
                referenceRecords[record.description] = record
    
    alignments = Alignments(referenceRecords)
    
    # Parse through BLAST XML
    alignments.ParseBlastOutput(blastFile, numThreads)
    logger.debug("Target dictionary length: %d", len(alignments))
    
    # See how many contigs are associated with each gene (for testing)
    for key in alignments.iterkeys():
        for trans in alignments[key].alternativeTranscripts:
            logger.debug("Key: %s; NumContigs: %s" % (key, len(trans.contigs)))
    
    # Output to crap
    logger.debug("Outputting alignments to: %s", outFolder)
    CreateDir(outFolder, True)
    CreateDir("%s/Alignments" % outFolder, True)
    CreateDir("%s/Consensus" % outFolder, True)
    CreateDir("%s/SplicedConsensus" % outFolder, True)
    alignments.OutputAlignments(numThreads, "%s/Alignments" % outFolder)
    alignments.OutputConsenses(numThreads, "%s/Consensus" % outFolder)
    alignments.OutputSplicedConsenses(numThreads, "%s/SplicedConsensus" % outFolder)
    logger.debug("Outputting stats")
    alignments.OutputStatistics("%s/Stats.csv" % outFolder)
    
    logger.info("End")
    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv[1:])