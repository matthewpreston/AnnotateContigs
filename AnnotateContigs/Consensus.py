from math import log

class ConsensusNT:
    def __init__(self, threshold):
        """
        Threshold is used to determine which nt to use as consensus
        Ex. (threshold = 0.7)
            ACGTAGTC
            ACSTMTWC
            TCSTCCAC
        yields
            Site    aveA    aveC    aveG    aveT    Consensus
            1       2/3     0       0       1/3     W
            2       0       1       0       0       C
            3       0       1/3     2/3     0       S
            4       0       0       0       1       T
            5       1/2     1/2     0       0       M
            6       0       1/3     1/3     1/3     B
            7       1/2     0       0       1/2     W
            8       0       1       0       0       C
        """
        self.A = 0
        self.C = 0
        self.G = 0
        self.T = 0
        if threshold <= 0.5:
            raise Exception("Threshold <= 0.5, doesn't make sense\n")
        self.threshold = threshold
        
    def GetAve(self):
        """Returns a tuple of nucleotide averages at this site"""
        total = self.A + self.C + self.G + self.T
        if total == 0:
            return (0,0,0,0)
        aveA = float(self.A) / total
        aveC = float(self.C) / total
        aveG = float(self.G) / total
        aveT = float(self.T) / total
        return (aveA, aveC, aveG, aveT)
        
    def GetConNT(self):
        """Returns nucleotide character (should probably do something about '-'"""
        nts = ['A','C','G','T']
        aves = self.GetAve()
        #Test if any nts in this position
        if aves == (0,0,0,0):
            return '-'
        #Set up dict
        result = {}
        for nt, ave in zip(nts, aves):
            result[nt] = ave
        #Try unambiguous singlet nts
        singlet = sorted([tuple([result[nt], nt]) for nt in nts],
                         key=lambda x:x[0],
                         reverse=True)
        #Check if one of the unambiguous singlet nts holds majority
        if singlet[0][0] >= self.threshold and singlet[0][0] != singlet[1][0]:
            return singlet[0][1]
        #Try ambiguous doublet nts
        doublet = sorted([tuple([(result['A'] + result['C']), 'M']),
                          tuple([(result['A'] + result['G']), 'R']),
                          tuple([(result['A'] + result['T']), 'W']),
                          tuple([(result['C'] + result['G']), 'S']),
                          tuple([(result['C'] + result['T']), 'Y']),
                          tuple([(result['G'] + result['T']), 'K'])],
                         key=lambda x:x[0],
                         reverse=True)
        #Check if one of the ambiguous doublet nts holds majority
        if doublet[0][0] >= self.threshold and doublet[0][0] != doublet[1][0]:
            return doublet[0][1]
        #Try ambiguous triplet nts
        triplet = sorted([tuple([(result['A']+result['C']+result['G']), 'V']),
                          tuple([(result['A']+result['C']+result['T']), 'H']),
                          tuple([(result['A']+result['G']+result['T']), 'D']),
                          tuple([(result['C']+result['G']+result['T']), 'B'])],
                         key=lambda x:x[0],
                         reverse=True)
        #Check if one of the ambiguous triplet nts holds majority
        if triplet[0][0] >= self.threshold and triplet[0][0] != triplet[1][0]:
            return triplet[0][1]
        #By process of elimination, it must be 'N'
        return 'N'

class ConsensusSeq:
    def __init__(self, threshold, records):
        """Takes fasta records and threshold (see ConsensusNT.__init__())"""
        self.seq = []
        for tup in zip(*records):
            temp = ConsensusNT(threshold)
            for i, nt in enumerate(tup):
                if nt == 'A':
                    temp.A += 1 
                elif nt == 'C':
                    temp.C += 1 
                elif nt == 'G':
                    temp.G += 1 
                elif nt == 'T':
                    temp.T += 1 
                elif nt == 'K':
                    temp.G += 0.5 
                    temp.T += 0.5 
                elif nt == 'M':
                    temp.A += 0.5 
                    temp.C += 0.5 
                elif nt == 'R':
                    temp.A += 0.5 
                    temp.G += 0.5 
                elif nt == 'Y':
                    temp.C += 0.5 
                    temp.T += 0.5 
                elif nt == 'S':
                    temp.C += 0.5 
                    temp.G += 0.5 
                elif nt == 'W':
                    temp.A += 0.5 
                    temp.T += 0.5 
                elif nt == 'B':
                    temp.C += float(1) / 3 
                    temp.G += float(1) / 3 
                    temp.T += float(1) / 3 
                elif nt == 'D':
                    temp.A += float(1) / 3 
                    temp.G += float(1) / 3 
                    temp.T += float(1) / 3 
                elif nt == 'H':
                    temp.A += float(1) / 3 
                    temp.C += float(1) / 3 
                    temp.T += float(1) / 3 
                elif nt == 'V':
                    temp.A += float(1) / 3 
                    temp.C += float(1) / 3 
                    temp.G += float(1) / 3 
                elif nt == 'N':
                    temp.A += 0.25 
                    temp.C += 0.25 
                    temp.G += 0.25 
                    temp.T += 0.25 
                else:
                    pass
            self.seq.append(temp)
            
    def GetConsensus(self):
        """Returns consensus sequence"""
        return "".join([n.GetConNT() for n in self.seq])
        
    def GetConStats(self):
        """Returns a list of tuples containing nt averages at each site"""
        return [n.GetAve() for n in self.seq]