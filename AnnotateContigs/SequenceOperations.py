# Contains useful operations for dealing with sequences and location objects

from Bio.Alphabet import generic_dna
from Bio.GenBank import _FeatureConsumer
from Bio.GenBank.utils import FeatureValueCleaner
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

def _fix_location(location):
    """Increments start by one since parsing it subtracts one somehow"""
    
    fixed = []
    for part in location.parts:
        temp = FeatureLocation(part.start+1, part.end, part.strand)
        fixed.append(temp)
    return sum(fixed)

def Splice(sequence, location):
    """Creates a spliced sequence based on the splicing locations"""
    
    # Some class which contains a method to splice
    consumer = _FeatureConsumer(use_fuzziness=1,
                    feature_cleaner=FeatureValueCleaner())
    # Must format location into a string which the consumer can read and splice
    location = str(_fix_location(location)).replace(" ", "")
    location = location.replace("{","(").replace("}",")")
    location = location.replace("[","").replace("]","")
    location = location.replace("(+)","").replace("(-)","").replace(":","..")
    # Doesn't have to be "CDS", change it to whatever may seem suitable
    # However, you must input some string for the .location method to work
    consumer.feature_key("CDS")
    # Input splicing instructions
    consumer.location(location)
    
    # Construct a SeqRecord object for the splicing method to work
    record = SeqRecord(Seq(sequence, generic_dna), id="", description="")
    # Apply the splicing method
    record = consumer.data.features[0].location.extract(record)
    # Return just the sequence in string form
    return str(record.seq)
    
def MappedIntersection(loc1, loc2):
    """
    Find the intersection of these two mapping locations, assumes that both are
    sorted lists of 2-tuples (because they should be by my algorithms...)
    """
    
    result = []
    it1, it2 = 0, 0
    while it1 != len(loc1) and it2 != len(loc2):
        m1, m2 = loc1[it1], loc2[it2]
        if m1[1] < m2[0]:   # If m1 occurs before m2, skip and get next
            it1 += 1
        elif m2[1] < m1[0]: # If m2 occurs before m1, skip and get next
            it2 += 1
        else: # Some overlap
            result.append((max(m1[0],m2[0]),min(m1[1],m2[1])))
            if m1[1] <= m2[1]:   # If m1's end occurs before m2, get next
                it1 += 1
            elif m2[1] <= m1[1]: # If m2's end occurs before m1, get next
                it2 += 1
    return result
    
def MappedDifference(loc1, loc2):
    """
    Find the difference of these two mapping locations, assumes that both are
    sorted lists of 2-tuples (because they should be by my algorithms...)
    """
    
    temp = [l for l in loc1] # Make a copy since we'll be modifying contents
    result = []
    it1, it2 = 0, 0
    while it1 != len(temp) and it2 != len(loc2):
        m1, m2 = temp[it1], loc2[it2]
        if m1[1] < m2[0]:
            result.append(m1)
            it1 += 1
        elif m1[1] == m2[0]:
            result.append((m1[0],m2[0]-1))
            it1 += 1
        elif m1[0] < m2[0] and m1[1] < m2[1]:
            result.append((m1[0],m2[0]-1))
            it1 += 1
        elif m1[0] < m2[0] and m1[1] == m2[1]:
            result.append((m1[0],m2[0]-1))
            it1 += 1
            it2 += 1
        elif m1[0] < m2[0] and m1[1] > m2[1]:
            result.append((m1[0],m2[0]-1))
            temp[it1] = (m2[1]+1,m1[1])
            it2 += 1
        elif m1[0] == m2[0] and m1[1] < m2[1]:
            it1 += 1
        elif m1[0] == m2[0] and m1[1] == m2[1]:
            it1 += 1
            it2 += 2
        elif m1[0] == m2[0] and m1[1] > m2[1]:
            temp[it1] = (m2[1]+1,m1[1])
            it2 += 1
        # Now we know that m1[0] > m2[0], so stop testing that
        elif m1[1] < m2[1]:
            it1 += 1
        elif m1[1] == m2[1]:
            it1 += 1
            it2 += 2
        elif m1[0] < m2[1]:
            temp[it1] = (m2[1]+1,m1[1])
            it2 += 1
        elif m1[0] == m2[1]:
            temp[it1] = (m2[1]+1,m1[1])
            it2 += 1
        else: # m1[0] > m2[1]
            it2 += 1
    return result + temp[it1:] # Add the remainder
    
def CodingLocationToList(location):
    """Converts BioPython's CodingLocationRecord into a list of 2-tuples"""
    result = []
    for part in location.parts:
        result.append((int(part.start), int(part.end)))
    return result
