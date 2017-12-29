from functools import wraps
import logging
import threading
import time
try:
    from urllib.error import HTTPError   # for Python 3
except ImportError as err:
    from urllib2 import HTTPError        # for Python 2

from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from Exceptions import *
from MergeLocations import MergeLocations
from RecordObjects import *
from TerminalCommands import CheckFileName

class AlignmentParser:
    """
    Contains methods to parse alignment into gene, accession, version, id,
    reference sequence, and contig.
    Initialization:
        loggerName(str)
            Name of logger for this class
    Members:
        logger(logging.Logger)
            Logger object with the given loggerName
    Methods:
        CheckIndel(ref, contig)
            Returns a fixed contig with indels removed
        GetGene(alignment)
            Returns gene abbreviation from alignment
        GetAccessionAndVersion(alignment)
            Returns accession and version numbers from alignment
        GetIDAndReferenceSequence(accession)
            Returns id and reference sequence from accession number
        GetLocations(ID)
            Returns a list of coding locations
        GetContigs(alignment, identity=1.0, similarity=1.0)
            Returns a list of contigs that align to the reference
        GetMappingandContigs(alignment, location, identity=1.0, similarity=1.0)
            Returns a list of contigs and the extent of mapping
        GetProteins(alignment, ID=None, identity=1.0, similarity=1.0)
            Returns a list of proteins/isoforms transcribed by the gene
    """

    def __init__(self, loggerName=__name__):
        self.logger = logging.getLogger(loggerName)
    
    def __repr__(self):
        return "AlignmentParser(loggerName=%r)" % self.logger.name

    def CheckIndel(self, ref, contig):
        """Rids of indels due to variation / sequencing error in contigs"""

        ref = list(ref)
        contig = list(contig)
        if len(ref) != len(contig):
            assert False, "AlignmentParser.CheckIndel Different sizes found"

        # Remove insertion (when aligned, reference will contain '-')
        for i in range(len(ref)):
            if ref[i] == '-':
                contig[i] = '!'
        ref = [n for n in ref if n != '-']
        contig = [n for n in contig if n != '!']

        # Remove deletion (fill with 'N')
        for i in range(len(contig)):
            if contig[i] == '-':
                contig[i] = 'N'
        return "".join(contig)
    
    def GetGene(self, alignment):
        """Returns gene abbreviation from alignment"""

        try:
            gene = alignment.hit_def

            # Usually, gene is surrounded by parentheses
            if gene.find("(") != -1:
                # Get the gene which is in the last set of parentheses
                # At least NCBI queries strongly tend to follow this
                gene = gene.split("(")[-1].split(")", 1)[0]

            # Augment genes with chars that may confuse file creation
            gene = CheckFileName(str(gene),
                    ['\\','/',':','*','?','\"','<','>','|',' ',','], '_')
        except AttributeError as e:
            self.logger.error("AlignmentParser.GetAccessionAndVersion failed")
            self.logger.error("Malformed alignment object")
            gene = None
        except Exception as e:
            self.logger.error("AlignmentParser.GetAccessionAndVersion failed")
            self.logger.error(e, exc_info=True)
            gene = None

        self.logger.debug("Gene: %s", gene)
        return gene
    
    def GetAccessionAndVersion(self, alignment):
        """Gets accession and version"""

        # NCBI gets rid of GI numbers in Sept 2016, so try statement
        # treats it as if GI numbers are still in effect, while except
        # statement believes they are phased out
        try:
            temp = alignment.title.split("|", 4)[3].split(".", 1)
            if len(temp) != 2:
                raise IndexError
            accession = temp[0]
            version = temp[1]
        except AttributeError:
            self.logger.error("AlignmentParser.GetAccessionAndVersion failed")
            self.logger.error("Malformed alignment object")
            accession = None
            version = None
        except IndexError as e:
            self.logger.warning("GI numbers may have been phased out, retrying")
            try:
                temp = alignment.title.split(" ", 1)[0].split(".", 1)
                if len(temp) != 2:
                    raise IndexError
                accession = temp[0]
                version = temp[1]
            except IndexError:
                self.logger.error("AlignmentParser.GetAccessionAndVersion failed")
                self.logger.error("Cannot find accession and version in: %s",
                                  alignment.title)
                accession = None
                version = None
        except Exception as e:
            self.logger.error("AlignmentParser.GetAccessionAndVersion failed")
            self.logger.error(e, exc_info=True)
            accession = None
            version = None

        self.logger.debug("Accession.version: %s.%s", accession, version)
        return (accession, version)
    
    def GetID(self, accession):
        """
        Gets the NCBI's internal ID for the accession number. To be used in
        conjunction with GetReferenceSequence() and GetLocations().
        """
        
        self.logger.debug("Attempting to get id")
        attempt = 1
        while attempt <= 3:
            try:
                netHandle = Entrez.esearch(db="nucleotide",
                                           term="%s[ACCN]" % accession)
                self.logger.debug("Got response, parsing id")
                record = Entrez.read(netHandle)
                if len(record["IdList"]):
                    ID = record["IdList"][0]
                else:
                    self.logger.warning("No IDs found from accession: %s",
                                        accession)
                    return None
                break
            except HTTPError as e:
                # NCBI/Internet problem
                self.logger.warning("Received error from server %s", e)
                self.logger.warning("Attempt %i of 3", attempt)
                attempt += 1
                time.sleep(15)
            except Exception as e:
                # BioPython throws a random exception, try again
                self.logger.warning("Received error from BioPython %s" % e)
                self.logger.warning("Attempt %i of 3", attempt)
                attempt += 1
                time.sleep(15)
            finally:
                try:
                    netHandle.close()
                except Exception as e:
                    pass
        if attempt > 3:
            self.logger.error("AlignmentParser.GetID failed")
            self.logger.error("Could not get ID from NCBI after 3 attempts, "
                              + "accession used: %s", accession)
            ID = None
        self.logger.debug("ID: %s", ID)
        return ID
    
    def GetReferenceSequence(self, ID):
        """Gets reference sequence from the internet using NCBI's internal ID"""

        self.logger.debug("Attempting to get reference sequence")
        attempt = 1
        while attempt <= 3:
            try:
                netHandle = Entrez.efetch(db="nucleotide", id=ID,
                                          rettype="fasta", retmode="text")
                self.logger.debug("Got response, parsing reference sequence")
                record = SeqIO.read(netHandle, "fasta")
                refSeq = str(record.seq)
                break
            except HTTPError as e:
                # NCBI/Internet problem
                self.logger.warning("Received error from server %s" % e)
                self.logger.warning("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(15)
            except Exception as e:
                # BioPython throws a random exception, try again
                self.logger.warning("Received error from BioPython %s" % e)
                self.logger.warning("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(15)
            finally:
                try:
                    netHandle.close()
                except Exception as e:
                    pass
        if attempt > 3:
            self.logger.error("AlignmentParser.GetReferenceSequence "
                              + "failed")
            self.logger.error("Could not get reference sequence from NCBI after"
                              + "3 attempts, accession used: %s", accession)
            refSeq = None
        self.logger.debug("ref[:10]: %s", refSeq[:10])
        return refSeq
    
    def GetLocations(self, ID):
        """Returns a list of coding locations"""

        # Don't try this at home kids! cd /; rm -rf *
        self.logger.debug("Attempting to get coding locations")
        attempt = 1
        while attempt <= 3:
            try:
                netHandle = Entrez.efetch(db="nucleotide", id=ID,
                                          rettype="gb", retmode="text")
                self.logger.debug("Got response, parsing coding locations")
                record = SeqIO.read(netHandle, "genbank")
                break
            except HTTPError as e:
                # NCBI/Internet problem
                self.logger.warning("Received error from server %s" % e)
                self.logger.warning("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(15)
            except Exception as e:
                # BioPython throws a random exception, try again
                self.logger.warning("Received error from BioPython %s" % e)
                self.logger.warning("Attempt %i of 3" % attempt)
                attempt += 1
                time.sleep(15)
            finally:
                try:
                    netHandle.close()
                except Exception as e:
                    pass
        if attempt > 3:
            self.logger.error("AlignmentParser.GetLocations failed")
            self.logger.error("Could not get coding locations from NCBI after"
                              + "3 attempts, ID used: %s", ID)
            return None

        # Parse these locations kthxbye
        locations = []
        for feature in record.features:
            if feature.type == "CDS":
                if "protein_id" in feature.qualifiers.keys():
                    locations.append(CodingLocationRecord(
                                     feature.qualifiers["protein_id"][0],
                                     feature.location))
                else:
                    locations.append(CodingLocationRecord("NA",
                                     feature.location))
        self.logger.debug("Splicing location count: %i", len(locations))
        return locations
        
    def GetContigs(self, alignment, identity=0.0, similarity=0.0):
        """Returns a list of contigs that align to the reference"""

        contigs = []
        try:
            for hsp in alignment.hsps:
                contig = self.CheckIndel(hsp.sbjct, hsp.query)

                # Check contig sequence identity (exact matches)
                if (float(hsp.identities) / hsp.align_length) < identity:
                    continue

                # Check contig sequence similarity (similar matches)
                if (float(hsp.positives) / hsp.align_length) < similarity:
                    continue

                # Collect start and end points of contig
                if hsp.sbjct_start < hsp.sbjct_end: # Foward
                    contigs.append(ContigRecord(contig, hsp.sbjct_start,
                                                hsp.sbjct_end))
                else: # Find reverse complement of contig for it to be forward
                    temp = Seq(contig, IUPAC.ambiguous_dna)
                    temp = temp.reverse_complement()
                    contigs.append(ContigRecord(str(temp),
                                                hsp.sbjct_end ,hsp.sbjct_start))
        except Exception as e:
            self.logger.error("AlignmentParser.GetContigs failed")
            self.logger.error(e, exc_info=True)
            return None

        self.logger.debug("Contig count: %i", len(contigs))
        return contigs
    
    def GetMappingandContigs(self, alignment, identity=0.0, similarity=0.0):
        """Returns a list of contigs and the extent of mapping"""
        # Try to get contigs
        contigs = self.GetContigs(alignment, identity, similarity)
        if contigs is None:
            self.logger.error("AlignmentParser.GetMappingandContigs failed")
            return (None, None)

        # Get mapping
        self.logger.debug("Calculating extent of mapping")
        mapped = [(contig.start, contig.end) for contig in contigs]
        mapped = MergeLocations(mapped)
        self.logger.debug("Mapped: %s", str(mapped))
        return (mapped, contigs)
    
    # DEPRECATED
    def GetMappingandContigs_Old(self, alignment, location, identity=0.0,
                             similarity=0.0):
        """Returns a list of contigs and the extent of mapping"""

        # Try to get contigs
        contigs = self.GetContigs(alignment, identity, similarity)
        if contigs is None:
            self.logger.error("AlignmentParser.GetMappingandContigs failed")
            return (None, None)

        self.logger.debug("Calculating extent of mapping")
        finalContigs = []
        temp = []
        for contig in contigs:
            for exon in location.parts:
                # BioPython's FeatureLocation class is stupid. When an exon
                # covers from [183, 295] (including nts at position 183 & 295),
                # the FeatureLocation class created from calling Entrez.efetch()
                # and SeqIO.parse() yields: [182:295](+). Testing whether a 
                # contig falls within this class yields:
                #
                # >>> 181 in exon
                # False
                # >>> 182 in exon
                # True (Should return False for our purposes)
                # >>> 183 in exon
                # True
                # >>> 294 in exon
                # True
                # >>> 295 in exon
                # False (Should return True for our purposes)
                #
                # As a solution, subtract 1 from end points for inclusion tests
                # Took some fiddling really to get it to work as intended
                if contig.start-1 in exon:
                    finalContigs.append(contig)
                    if contig.end-1 in exon:
                        temp.append((contig.start, contig.end))
                    else:
                        temp.append((contig.start, exon.end))
                elif contig.end-1 in exon:
                    finalContigs.append(contig)
                    temp.append((exon.start+1, contig.end))
        mapped = MergeLocations(temp)

        self.logger.debug("Mapped: %s", str(mapped))
        return (mapped, finalContigs)
    
    def GetAlternativeTranscripts(self, alignment, ID=None,
                                  identity=0.0, similarity=0.0):
        """
        Dangnabbit, get me them darned list o' alternatively spliced transcripts
        via them fancy alignments those new young folk be talking about. Aw hell
        gimme that ID too if ya couldja kindly?
        """

        if ID is None:
            ID = self.GetID(alignment)

        records = self.GetLocations(ID)
        if records is None:
            self.logger.error("AlignmentParser.GetProteins failed")
            return None

        # Find proteins based off of coding locations
        altRecs = []
        for record in records:
            altRec = AlternativeTranscriptRecord()
            altRec.proteinID = record.proteinID
            altRec.location = record.location
            altRec.mapped, altRec.contigs = self.GetMappingandContigs(
                                                alignment, identity, similarity)
            if altRec.mapped is None:
                self.logger.error("AlignmentParser.GetProteins failed")
                return None
            altRecs.append(altRec)
            self.logger.debug("Protein ID: %s", altRec.proteinID)
            self.logger.debug("Location: %s", altRec.location)
            self.logger.debug("Map: %s", altRec.mapped)
            self.logger.debug("Number of contigs: %i", len(altRec.contigs))
        return altRecs
    
    def GetReferenceRecord(self, alignment, referenceRecords):
        """Fetches desired reference record"""
        
        title = alignment.hit_def
        refRecord = referenceRecords.get(title)
        if refRecord is None:
            self.logger.error("AlignmentParser.GetReferenceRecord failed")
            self.logger.error("No reference sequence titled: %s", title)
            return None
        return refRecord
    
    def GetReferenceSequenceFromReference(self, alignment, referenceRecords):
        """
        Retrieves reference sequence based on the reference gene name used in
        the alignment
        """
        
        refRecord = self.GetReferenceRecord(alignment, referenceRecords)
        if refRecord is None:
            self.logger.error(
                    "AlignmentParser.GetReferenceSequenceFromReference failed")
            self.logger.error("Could not find reference record")
            return None
        return str(self.GetReferenceRecord(alignment, referenceRecords).seq)
    
    def GetAlternativeTranscriptsFromReference(self, alignment,
            referenceRecords, identity=1.0, similarity=1.0):
        """Retrieve contigs mapped to the reference"""
        
        # Get reference record
        refRecord = self.GetReferenceRecord(alignment, referenceRecords)
        if refRecord is None:
            return None
        
        # Get data
        location = FeatureLocation(0, len(refRecord), 1)
        contigs = self.GetContigs(alignment, identity, similarity)
        if contigs is None:
            return None
            
        temp = []
        for contig in contigs:
            temp.append((contig.start, contig.end))
        mapped = MergeLocations(temp)

        # Return that bad boy
        return [AlternativeTranscriptRecord(None, location, mapped, contigs)]
        
    def ParseAlignment(self, alignment, targetDict, referenceRecords=None,
                       identity=0.0, similarity=0.0):
        """
        Parses and BLAST alignment and adds it to the common target dictionary
        """
        
        # Check if we have already created an entry with the same gene
        key = alignment.title
        targetRecord = targetDict.get(key)
        
        if referenceRecords is None: # If no reference, use NCBI
            if targetRecord is None: # Uninitialized record
                record = AlignmentRecord()
                # Try to get the gene name
                gene = self.GetGene(alignment)
                if gene is None:
                    raise FailedToParseException("Failed to parse gene")
                record.gene = gene
                # Get accession.version
                acc, ver = self.GetAccessionAndVersion(alignment)
                if acc is None:
                    raise FailedToParseException("Failed to parse accession")
                record.accession = acc
                record.version = ver
                # Get NCBI internal ID
                ID = self.GetID(acc)
                if ID is None:
                    raise FailedToParseException("Failed to parse ID")
                record.ID = ID
                # Get reference sequence from NCBI
                refSeq = self.GetReferenceSequence(ID)
                if refSeq is None:
                    raise FailedToParseException("Failed to parse reference "
                                                 + "sequence via NCBI")
                record.referenceSequence = refSeq
                # Get alternative transcripts
                altTrans = self.GetAlternativeTranscripts(alignment, ID,
                                                          identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via NCBI")
                record.alternativeTranscripts = altTrans
            else: # Record has altready been initialized
                # Create another record regardless
                record = AlignmentRecord(gene=targetRecord.gene,
                        accession=targetRecord.accession,
                        version=targetRecord.version,
                        ID=targetRecord.ID,
                        referenceSequence=targetRecord.referenceSequence)
                # Get alternative transcripts
                altTrans = self.GetAlternativeTranscripts(alignment, record.ID,
                                                          identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via NCBI")
                record.alternativeTranscripts = altTrans
        else: # Reference available
            if targetRecord is None: # Uninitialized record
                record = AlignmentRecord()
                # Try to get the gene name
                gene = self.GetGene(alignment)
                if gene is None:
                    raise FailedToParseException("Failed to parse gene")
                record.gene = gene
                # Add given reference sequence
                refSeq = self.GetReferenceSequenceFromReference(alignment,
                        referenceRecords)
                if refSeq is None:
                    raise FailedToParseException("Failed to parse reference "
                                                 + "sequence via FASTA")
                record.referenceSequence = refSeq
                # Get alternative transcripts given reference
                altTrans = self.GetAlternativeTranscriptsFromReference(
                        alignment, referenceRecords, identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via FASTA")
                record.alternativeTranscripts = altTrans
            else: # Record has already been initialized
                # Create another record regardless
                record = AlignmentRecord(gene=targetRecord.gene,
                    accession=targetRecord.accession,
                    version=targetRecord.version,
                    ID=targetRecord.ID,
                    referenceSequence=targetRecord.referenceSequence)
                # Get alternative transcripts given reference
                altTrans = self.GetAlternativeTranscriptsFromReference(
                        alignment, referenceRecords, identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via FASTA")
                record.alternativeTranscripts = altTrans
            
        # Add to target dictionary
        targetDict[key].Update(record)
        
class ThreadSafeAlignmentParser(AlignmentParser):
    """
    Thread safe implementation of AlignmentParser, meant to be run within the
    target thread.
    """

    NCBI_DELAY = 0.333333334            # Delay between requests (across threads)
    lastRequest = 0.0                   # Static variable; last request time
    delayLock = threading.RLock()       # For making lastRequest thread-safe
    targetDictLock = threading.RLock()  # For making targetDict thread-safe
    
    def __init__(self, name=None):
        AlignmentParser.__init__(self, name)
    
    def delayNCBIRequest(function):
        """
        Delays the internet request to NCBI by a certain delay amongst all
        threads
        """
        
        @wraps(function)
        def wrapper(*args, **kwargs):
            with ThreadSafeAlignmentParser.delayLock:
                current = time.time()
                wait = ThreadSafeAlignmentParser.lastRequest \
                     + ThreadSafeAlignmentParser.NCBI_DELAY - current
                if wait > 0:
                    time.sleep(wait)
                    ThreadSafeAlignmentParser.lastRequest = current + wait
                else:
                    ThreadSafeAlignmentParser.lastRequest = current
            return function(*args, **kwargs)
        return wrapper
        
    # Decorate NCBI accessing inherited methods with delaying method
    GetID = delayNCBIRequest(
            AlignmentParser.GetID.__func__)
    GetReferenceSequence = delayNCBIRequest(
            AlignmentParser.GetReferenceSequence.__func__)
    GetLocations = delayNCBIRequest(
            AlignmentParser.GetLocations.__func__)
            
    def ParseAlignment(self, alignment, targetDict, referenceRecords=None,
                       identity=0.0, similarity=0.0):
        """
        Parses and BLAST alignment and adds it to the common target dictionary
        """
        
        # Check if we have already created an entry with the same gene
        key = alignment.title
        with ThreadSafeAlignmentParser.targetDictLock:
            targetRecord = targetDict.get(key)
        
        if referenceRecords is None: # If no reference, use NCBI
            if targetRecord is None: # Uninitialized record
                record = AlignmentRecord()
                # Try to get the gene name
                gene = self.GetGene(alignment)
                if gene is None:
                    raise FailedToParseException("Failed to parse gene")
                record.gene = gene
                # Get accession.version
                acc, ver = self.GetAccessionAndVersion(alignment)
                if acc is None:
                    raise FailedToParseException("Failed to parse accession")
                record.accession = acc
                record.version = ver
                # Get NCBI internal ID
                ID = self.GetID(acc)
                if ID is None:
                    raise FailedToParseException("Failed to parse ID")
                record.ID = ID
                # Get reference sequence from NCBI
                refSeq = self.GetReferenceSequence(ID)
                if refSeq is None:
                    raise FailedToParseException("Failed to parse reference "
                                                 + "sequence via NCBI")
                record.referenceSequence = refSeq
                # Get alternative transcripts
                altTrans = self.GetAlternativeTranscripts(alignment, ID,
                                                          identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via NCBI")
                record.alternativeTranscripts = altTrans
            else: # Record has altready been initialized
                # Create another record regardless (merging records needs to be 
                # atomic, thus we must gather all data first then merge in a
                # thread safe manner)
                record = AlignmentRecord(gene=targetRecord.gene,
                        accession=targetRecord.accession,
                        version=targetRecord.version,
                        ID=targetRecord.ID,
                        referenceSequence=targetRecord.referenceSequence)
                # Get alternative transcripts
                altTrans = self.GetAlternativeTranscripts(alignment, record.ID,
                                                          identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via NCBI")
                record.alternativeTranscripts = altTrans
        else: # Reference available
            if targetRecord is None: # Uninitialized record
                record = AlignmentRecord()
                # Try to get the gene name
                gene = self.GetGene(alignment)
                if gene is None:
                    raise FailedToParseException("Failed to parse gene")
                record.gene = gene
                # Add given reference sequence
                refSeq = self.GetReferenceSequenceFromReference(alignment,
                        referenceRecords)
                if refSeq is None:
                    raise FailedToParseException("Failed to parse reference "
                                                 + "sequence via FASTA")
                record.referenceSequence = refSeq
                # Get alternative transcripts given reference
                altTrans = self.GetAlternativeTranscriptsFromReference(
                        alignment, referenceRecords, identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via FASTA")
                record.alternativeTranscripts = altTrans
            else: # Record has already been initialized
                # Create another record regardless (merging records needs to be 
                # atomic, thus we must gather all data first then merge in a
                # thread safe manner)
                record = AlignmentRecord(gene=targetRecord.gene,
                    accession=targetRecord.accession,
                    version=targetRecord.version,
                    ID=targetRecord.ID,
                    referenceSequence=targetRecord.referenceSequence)
                # Get alternative transcripts given reference
                altTrans = self.GetAlternativeTranscriptsFromReference(
                        alignment, referenceRecords, identity, similarity)
                if altTrans is None:
                    raise FailedToParseException("Failed to parse alternative "
                                                 + "transcripts via FASTA")
                record.alternativeTranscripts = altTrans
            
        # Add to target dictionary, whilst keeping it thread safe
        with ThreadSafeAlignmentParser.targetDictLock:
            targetRecord = targetDict.get(key)
            if targetRecord is None:
                targetDict[key] = record
            else:
                targetDict[key].Update(record)