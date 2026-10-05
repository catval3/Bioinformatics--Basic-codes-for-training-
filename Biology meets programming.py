
def main():
    sequence = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTGTGGATCATACTGATCTGTGCTTCCAATTACGGAGCGGGATTACGAAAGCATGGCTTCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGA"
    k = 5
    t = 10
    DNA = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

    print(f"The most frequent pattern is: {FrequentWords(sequence, k)}")
    print(f"The reverse complement of this sequence is:{ReverseComplement(sequence)}")
    print(f"The list of index number of ATA patterns: {PatternMatching("ATA", sequence)}")
    print(f"Skew Array:{SkewArray(sequence)}")
    Motifs = GreedyMotifSearch(DNA, k, t)
    print(Motifs)
    print(Score(Motifs))

#generating a frequence map, that counts how many times every possible k-mer (substring of length k) appears in a DNA sequence and stores the counts in a dictionary.
def FrequencyMap(Text, k):
    freq = {}       #start with a blank dictionary named “freq”
    n = len(Text)  #n is set to the length of the sequence
    for i in range(n-k+1): 
        Pattern = Text[i:i+k] 
        freq[Pattern] = 0  #adding a new key to the dictionary, the key is the pattern we got. 
    for i in range(n-k+1):  
        Pattern=Text[i:i+k] 
        freq[Pattern]+= 1 
    return freq

#after generating a frequency map, this code finds the most frequently occurring k-mer(s) in a DNA sequence and returns all k-mers tied for the highest frequency.
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
        # add each key to words whose corresponding frequency value is equal to m
    return words

#Produces the reverse-complementary DNA sequence of a given pattern. This is useful when looking for the same DNA pattern on the opposite strand.
def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

#Reverses the order of the nucleotides in a DNA sequence.
def Reverse(Pattern):
    Reverse = ""
    for char in Pattern:
        if char in ["A","T","C","G"]:
            Reverse = char + Reverse 
    return Reverse
#Converts each nucleotide into its complementary nucleotide: A ↔ T and C ↔ G.
def Complement(Pattern):
    complement = ""
    for char in Pattern:
        if char == "A":
            complement += char.replace("A", "T")
        elif char == "T":
            complement += char.replace("T", "A")
        elif char == "C":
            complement += char.replace("C", "G")
        elif char == "G":
            complement += char.replace("G", "C")
    return complement

#Finds the positions/index numbers where a specific DNA pattern occurs in a genome or DNA sequence.
def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome)-len(Pattern)+1): 
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

#Creates a skew array that tracks the difference between the number of G and C nucleotides as you move through the genome. It is commonly used to help identify regions such as the origin of replication.
def SkewArray(Genome):
    n = len(Genome)
    Skew = [0] 
    for i in range(n):
        if Genome[i] == "G":
            Skew.append(Skew[i] + 1) 
        elif Genome[i] == "C":
            Skew.append(Skew[i] - 1)
        elif Genome[i] == "A":
            Skew.append(Skew[i])
        elif Genome[i] == "T":
            Skew.append(Skew[i])
    return Skew 

#To find a common DNA pattern (motif) that appears, with some variation, across several DNA sequences.
def Count(Motifs): 
    count = {} 
    k = len(Motifs[0])
    for symbol in "ACGT": 
        count[symbol] = []
        for j in range(k): 
             count[symbol].append(0) 
    t = len(Motifs) 
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j] 
            count[symbol][j] += 1 
    return count
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    count = Count(Motifs) 
    for symbol in "ACGT":
        profile[symbol] = [] 
        for j in range(k): 
            profile[symbol].append(count[symbol][j] / t)  
    return profile
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    for i in range(len(Motifs)):
        for j in range(len(Motifs[i])):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score
def Consensus(Motifs):
    count = Count(Motifs)
    k = len(Motifs[0])
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
def ProfileMostProbableKmer(Text, k, Profile):
    max_pr = -1
    most_probable = Text[0:k]
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        p = Pr(kmer, Profile)
        if p > max_pr:
            max_pr = p
            most_probable = kmer
    return most_probable
def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p
def GreedyMotifSearch(Dna, k, t):
    n = len(Dna[0])
    best_motifs = [Dna[i][0:k] for i in range(t)]
    for i in range(n - k + 1):
        motifs = [Dna[0][i:i + k]]
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
    return best_motifs


main()