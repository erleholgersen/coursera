### A1.py #########################################################################################


### COURSE DEFINED FUNCTIONS ######################################################################

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33

### naive_with_rc #################################################################################
def naive_with_rc(p, t):

    occurrences = []

    # get matches for string itself
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record

    # get reverse complement
    rc_p = reverseComplement(p)

    # if reverse complement equals string itself, return
    if p == rc_p:
        return occurrences

    # get reverse complement matches
    # this implementation is definitely not optimal
    # could combine with loop above, but that  would require more finesse
    for i in range(len(t) - len(rc_p) + 1):  # loop over alignments
        match = True
        for j in range(len(rc_p)):  # loop over characters
            if t[i+j] != rc_p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences


lambdavirus_genome = readGenome('lambda_virus.fa')

### QUESTIONS #####################################################################################

# Question 1
q1_occurrences = naive_with_rc('AGGT', lambdavirus_genome)
q1_answer = len(q1_occurrences)

# Question 2
q2_occurrences = naive_with_rc('TTAA', lambdavirus_genome)
q2_answer = len(q2_occurrences)

# Question 3
q3_occurrences = naive_with_rc('ACTAAGT', lambdavirus_genome)
q3_answer = min(q3_occurrences)

# Question 4
q4_occurrences = naive_with_rc('AGTCGA', lambdavirus_genome)
q4_answer = min(q4_occurrences)

# Question 5

def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatch_count = 0 # keep track of number of mismatches
        match = True
        
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatch_count += 1
            if mismatch_count > 2:
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

q5_occurrences = naive_2mm('TTCAAGCC', lambdavirus_genome)
q5_answer = len(q5_occurrences)

# Question 6
q6_occurrences = naive_2mm('AGGAGGTT', lambdavirus_genome)
q6_answer = min(q6_occurrences)

# Question 7
reads, qualities = readFastq('ERR037900_1.first1000.fastq')

# get max read length
read_lengths = []
for read in reads:
    read_lengths.append(len(read))
    
max_read_length = max(read_lengths)


# examine quality scores for each offset
qualities_by_offset = []

for i in range(0, max_read_length):
    numeric_qualities = []
    
    # for each read, get numeric quality score at that position
    for read in reads:
        numeric_qualities.append( phred33ToQ(read[i]) )
    
    qualities_by_offset.append(numeric_qualities)

# print average quality score
for i in range(0, len(qualities_by_offset)):
    position_scores = qualities_by_offset[i]
    
    print i, sum(position_scores)/len(position_scores)

# offset 66 has higher quality than all the rest
# print reads at that position
for read in reads[1:20]:
    print read[66]

# many recorded as N.. suspicious



