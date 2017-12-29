# Holds record structures used by AnnotateContigs and its dependencies

from Consensus import ConsensusSeq
import Exceptions
from MergeLocations import MergeLocations
from SequenceOperations import _fix_location, CodingLocationToList, \
                               MappedIntersection, Splice

class CodingLocationRecord:
    """
    Holds info on protein id and coding locations on a gene.

    Initialization:
        proteinID(string)
            Protein ID
        location(Bio.SeqFeature.SeqFeature OR Bio.SeqFeature.CompoundLocation)
            Location where mapping starts

            To get at each individual range of contiguous coding sequences, use:
            # Ex For NCBI accession: BK006917, 1st CDS is:
            # join(4332..6749,339984..340042,351193..351281,363963..364216)
            >>> record = CodingLocationRecord(id, location)
            >>> for part in record.location:
            ...     print part,
            ...
            [4331:6749](+) [339983:340042](+) [351192:351281](+) 
            [363962:364216](+)
    """

    def __init__(self, proteinID=None, location=[]):
        self.proteinID = proteinID
        self.location = location

    def __repr__(self):
        return "CodingLocationRecord(proteinID=%r, location=%r)" \
                % (self.proteinID, self.location)

class ContigRecord:
    """
    Holds info on a contig's sequence and it's mapping location on the reference
    gene.

    Initialization:
        seq(string)
            Sequence
        start(int)
            Location where mapping starts
        end(int)
            Location where mapping ends
    """

    def __init__(self, seq=None, start=None, end=None):
        self.seq = seq
        self.start = start
        self.end = end

    def __repr__(self):
        return "ContigRecord(seq=%r, start=%r, end=%r)" \
                % (self.seq, self.start, self.end)

    def GetPaddedContig(self, length):
        """Returns a '-' padded contig that is of a given length"""
        return "-"*(self.start-1) + self.seq + "-"*(length-self.end)
                
class AlternativeTranscriptRecord:
    """
    Holds NCBI protein ID, splicing instructions, locations of where the protein
    was mapped, and the contigs which map to it. 

    Initialization:
        proteinID(string)
            The protein ID provided by NCBI to what protein this alternative
            transcript translates to (used for identification purposes)
        location(list<CodingLocationRecord>)
            List of splicing instructions (start and stop)
        mapped(list<2-tuple>)
            Locations of where the protein was mapped by contigs in tuples of
            (start, end)
        contigs(list<ContigRecord>)
            List of contigs that map to this protein
    """

    def __init__(self, proteinID=None, location=[], mapped=[], contigs=[]):
        self._consensus = None
        self.proteinID = proteinID
        self.location = location
        self.mapped = mapped
        self.contigs = contigs

    def __repr__(self):
        return "AlternativeTranscriptRecord(proteinID=%r, location=%r, " \
                + "mapped=%r, contigs=%r)" \
                % (self.proteinID, self.location, self.mapped, self.contigs)
               
    def BasesMapped(self):
        """Returns the number of bases mapped to the entire gene sequence"""
        return sum(m[1]-m[0]+1 for m in self.mapped)
               
    def Breadth(self):
        """Returns the number of bases mapped to the final spliced transcript"""
        return sum(m[1]-m[0]+1 for m in self.MappingToSpliced())
                
    def Coverage(self):
        """
        Returns the coverage of the reference sequence based on mapping of the
        contigs
        """
        
        splicedLength = self.SplicedLength()
        if splicedLength == 0:
            raise Exceptions.ZeroLengthedAlternativeTranscriptException(
                "Zero lengthed alternative transcript coding protein ID: %s",
                proteinID)
        return float(self.Breadth()) / splicedLength
        
    def Consensus(self, length, threshold):
        """Creates a consensus sequence based on overlapping of contigs"""
        
        records = [c.GetPaddedContig(length) for c in self.contigs]
        return ConsensusSeq(threshold, records).GetConsensus()

    def MappingToSpliced(self):
        """Finds the extent of the contigs mapping to the final spliced mRNA"""
        return MappedIntersection(
                CodingLocationToList(_fix_location(self.location)), self.mapped)
        
    def SplicedConsensus(self, length, threshold):
        """
        Creates a spliced consensus sequence based on the overlapping of contigs
        and splicing locations
        """
        return Splice(self.Consensus(length, threshold), self.location)
            
    def SplicedLength(self):
        """Returns the length of the final spliced transcript"""
        return sum(p.end-p.start for p in self.location.parts)

class AlignmentRecord:
    """
    Stores all alignments for a certain gene. Accession, version, and ID are
    only needed if reference sequence is to be looked up online, else just
    provide one. The member 'proteins' will 

    Initialization:
        gene(string)
            Gene name
        accession(string)
            Accession number (GI numbers get phased out by Sept 2016)
        version(string)
            Version number
        ID(string)
            ID used to look up NCBI record
        referenceSequence(string)
            Reference sequence
        alternativeTranscripts(list<AlternativeTranscriptRecord>)
            List of alternative transcript records spliced from the gene
    """

    def __init__(self, gene=None, accession=None, version=None, ID=None, 
                 referenceSequence=None, alternativeTranscripts=[]):
        self.gene = gene
        self.accession = accession
        self.version = version
        self.ID = ID
        self.referenceSequence = referenceSequence
        self.alternativeTranscripts = alternativeTranscripts

    def __repr__(self):
        return "AlignmentRecord(gene=%r, accession=%r, version=%r, ID=%r, " \
                + "referenceSequence=%r, alternativeTranscripts=%r)" \
                % (self.gene, self.accession, self.version, self.ID,
                    self.referenceSequence, self.alternativeTranscripts)

    def Update(self, newRecord):
        """Adds additional alternative transcripts to record"""
        
        if self.alternativeTranscripts is None:
            raise Exceptions.EmptyAlternativeTranscriptsException(
                "AlignmentRecord.Update failed: No alternative transcripts")
        
        for newAltTrans in newRecord.alternativeTranscripts:
            for oldAltTrans in self.alternativeTranscripts:
                # Find the same transcript
                if newAltTrans.proteinID == oldAltTrans.proteinID:
                    # Update it
                    oldAltTrans.mapped = MergeLocations(
                            newAltTrans.mapped + oldAltTrans.mapped)
                    oldAltTrans.contigs += newAltTrans.contigs
                    break
        
    def ExomeLengths(self):
        """Returns a dictionary of {index: spliced transcript length}"""
        
        result = {}
        for i in range(len(self.alternativeTranscripts)):
            # Calculate the length of the spliced transcript
            length = 0
            for exon in self.alternativeTranscripts[i].location.parts:
                length += len(exon)
            result[i] = length
        return result
        
    def GetLongestTranscripts(self, n=float("inf")):
        """Returns a list of the nth longest alternative transcripts"""
        
        entries = [(k,v) for k,v in self.ExomeLengths().iteritems()]
        entries.sort(key=(lambda (k,v): v), reverse=True)
        result = []
        for entry, length in entries[:min(n, len(entries))]:
            result.append(self.alternativeTranscripts[entry])
        return result