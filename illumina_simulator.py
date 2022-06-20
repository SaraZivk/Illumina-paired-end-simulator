import random
import pysam
from scipy.stats import binom
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

genes = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}

def cigarise(aligned, l_src, l_read):
    # make CIGAR string from uples of overlapping segments
    cigar = []
    aligned_src  = aligned[0]
    aligned_read = aligned[1]
    index_src  = 0
    
    current_src  = 0
    current_read = 0    
    while current_src < l_src or current_read < l_read:        
        if index_src < len(aligned_read):
            start_src  = aligned_src[index_src][0]
            end_src    = aligned_src[index_src][1]
            
            start_read = aligned_read[index_src][0]
            end_read   = aligned_read[index_src][1]            
        
        if start_src - current_src == 0:
            if start_read - current_read == 0:
                cigar.append((0, end_src-start_src)) # match
                current_src  = end_src
                current_read = end_read
                index_src += 1
                continue
            else:
                cigar.append((1, start_read - current_read)) # insertion   
                current_read = start_read
                continue            
        if start_src - current_src > 0:
            delta_src  = start_src  - current_src
            delta_read = start_read - current_read
            current_src  = start_src
            current_read = start_read
            if delta_src == delta_read:
                cigar.append((0, delta_src)) # mismatch
                continue            
            if delta_read > delta_src:
                cigar.append((0, delta_src)) # mismatch
                cigar.append((1, delta_read - delta_src)) # insertion
                continue
            if delta_read > 0:
                cigar.append((0, delta_read)) # mismatch
                cigar.append((2, delta_src - delta_read)) # deletion
                continue
            cigar.append((2, delta_src)) # deletion
            continue
        if start_src - current_src < 0:
            delta_src  = l_src  - current_src
            delta_read = l_read - current_read
            current_src  = l_src
            current_read = l_read
            if delta_src == delta_read:
                cigar.append((0, delta_src)) # mismatch
                continue
            if delta_src > delta_read:
                cigar.append((2, delta_src - delta_read)) # deletion
                if delta_read > 0:
                    cigar.append((0, delta_read)) # mismatch
                continue
            if delta_src < delta_read:
                cigar.append((1, delta_read - delta_src)) # insertion
                if delta_src > 0:
                    cigar.append((0, delta_src)) # mismatch
                continue
    return cigar

def errorise(fragment, error_rate_snip, error_rate_del, error_rate_ins):
    # insert errors with given error rates
    i = 0
    while i < len(fragment):
        
        while random.random() < error_rate_del:
            if i == 0:
                continue
            fragment = fragment[:i-1] + fragment[i:]
            i -= 1
            
        if random.random() < error_rate_snip:
            snip = random.randint(0, 3)
            while genes[snip] == fragment[i]:
                snip = random.randint(0, 3)
            fragment = fragment[:i-1]+genes[snip]+fragment[i:]
                        
        while random.random() < error_rate_ins:
            ins = random.randint(0, 3)
            fragment = fragment[:i]+genes[ins]+fragment[i:]
            i += 1
        i += 1
        
    return fragment
    
def qualitize(fragment, max_diff):
    # compute nucleotide qualities based on normal distribution
    dist = [100 * binom.pmf(r, 2 * max_diff, 0.5) for r in range(max_diff, 2*max_diff+1)]
    qual_val = []
    
    for i in range(len(fragment)):
        start = max(i - max_diff, 0)
        end   = min(i + max_diff, len(fragment))
        
        q = 0
        for j in range(start, end):
            if fragment[i] == fragment[j]:
                q += dist[abs(i-j)]
        
        qual_val.append(q)   
    return qual_val
    
def illumina_simulator(reference_path, coverage, read_length, \
                       insert_size, avg_quality, max_diff, \
                       error_rate_snip, error_rate_ins, error_rate_del):
    
    output_R1 = "output_r1_" + reference_path[:-5] + "fastq"
    output_R2 = "output_r2_" + reference_path[:-5] + "fastq"
    output_sam = "output_" + reference_path[:-5] + "bam"
    
    genome_length = 0
    cumsum = [0]
    ids = []
    #%% compute number and positions of fragments
    record_dict = SeqIO.index(reference_path, "fasta")
    for seq_record in SeqIO.parse(reference_path, "fasta"):
        genome_length += len(seq_record)
        cumsum.append(cumsum[-1] + len(seq_record))
        ids.append(seq_record.id)

    N_reads = int(genome_length / read_length / 2) # N  ->>> * coverage, later 
    
    positions = []    
    for i in range(coverage):    
        positions.extend( \
            random.sample( \
                range(0, genome_length - insert_size + 1), N_reads ) )
    N_reads *= coverage
    #%% Fragmentation and error insertions
    all_source_sequences = []
    all_r1_sources       = []
    all_r2_sources       = []
    all_r1_qual_val      = []
    all_r2_qual_val      = []
    start_avg_quality    = 0
    n_qualities          = 0
    
    for pos in positions:
        i = 0
        while cumsum[i+1] < pos:
            i += 1
        j = i
        while pos + insert_size < cumsum[j]:
            j += 1
        start = pos - cumsum[i]
        if cumsum[i+1] - insert_size >= pos:        
            end   = start + insert_size
            source_sequence = MutableSeq(record_dict[ids[i]].seq[start : end])
        else:
            source_sequence = MutableSeq(record_dict[ids[i]].seq[start : ])
            end = insert_size - cumsum[i+1] + pos
            i += 1
            while i < j: # if the sequence extends over more lines in fasta file
                source_sequence.extend(record_dict[ids[i]].seq)
                end -= len(record_dict[ids[i]])
                i += 1
            source_sequence.extend(record_dict[ids[i]].seq[ : end])
        all_source_sequences.append(source_sequence)
        
        # errors are inserted into source fragment in 2 ways for 2 reads
        fragment1 = errorise(str(source_sequence), error_rate_snip, error_rate_del, error_rate_ins)
        fragment2 = errorise(str(source_sequence), error_rate_snip, error_rate_del, error_rate_ins)
                
        all_r1_sources.append(Seq(fragment1))
        all_r2_sources.append(Seq(fragment2))
        
        # nucleotide qualities for both fragments
        qual_val1 = qualitize(fragment1, max_diff)   
        qual_val2 = qualitize(fragment2, max_diff)        
        
        qual_val1 = qual_val1[:min(len(fragment1), read_length)]
        qual_val2 = qual_val2[:-min(len(fragment2), read_length)-1:-1]
            
        all_r1_qual_val.append(qual_val1)
        all_r2_qual_val.append(qual_val2)
        
        # for calculating average nucl qual 
        start_avg_quality += sum(qual_val1) + sum(qual_val2)
        n_qualities += len(qual_val1) + len(qual_val2)
        
    start_avg_quality /= n_qualities
    multiplier = avg_quality / start_avg_quality # scaling factor
    #%% scaling nucleotide qualities
    
    all_r1_qual = []
    all_r2_qual = []
    
    for qq in all_r1_qual_val:
        qual = ""
        for q in qq:
            q *= multiplier
            q = min(93, int(q))
            q = max(1, int(q))
            qual += chr(q+33)
        all_r1_qual.append(qual) 
    
    for qq in all_r2_qual_val:
        qual = ""
        for q in qq:
            q *= multiplier
            q = min(93, int(q))
            q = max(1, int(q))
            qual += chr(q+33)
        all_r2_qual.append(qual) 
    #%% Output files    
    header = { 'HD': {'VN': '1.0'},
                'SQ': [{'LN': 1575, 'SN': 'chr1'},
                       {'LN': 1584, 'SN': 'chr2'}] }
    
    f1 = open(output_R1, "w")
    f2 = open(output_R2, "w")
    outf = pysam.AlignmentFile(output_sam, "w", header=header)
    
    aligner = Align.PairwiseAligner()
    #%% Printing into output files
    for i in range(N_reads):    
        
        source_seq = all_source_sequences[i]
        
        fragment1  = all_r1_sources[i]
        read1 = fragment1[:min(len(fragment1), read_length)]
        qual1 = all_r1_qual[i]
        src1  = source_seq[:min(len(fragment1), read_length)]
        
        fragment2  = all_r2_sources[i]
        ll2 = min(len(fragment2), read_length)
        read2 = fragment2[:-ll2-1:-1]
        qual2 = all_r2_qual[i]
        src2  = source_seq[:-ll2-1:-1]
        
        f1.write("@%s\n%s\n+\n%s\n" % (str(i), read1, qual1))
        f2.write("@%s\n%s\n+\n%s\n" % (str(i), read2.complement(), qual2))
            
        # calculating alignments ang CIGAR strings
        alignments = aligner.align(src1, read1)
        cigar1 = cigarise(alignments[0].aligned, len(src1), len(read1))   
        
        first_base1 = alignments[0].aligned[0][0][0]
        templ_len1  = alignments[0].aligned[0][-1][1] - first_base1
        
        alignments = aligner.align(src2, read2)
        cigar2 = cigarise(alignments[0].aligned, len(src2), len(read2))        
        first_base2 = alignments[0].aligned[0][0][0]
        templ_len2  = alignments[0].aligned[0][-1][1] - first_base2
    
        a1 = pysam.AlignedSegment()
        a1.query_name = '*' # ‘*’ indicates the information is unavailable
        a1.query_sequence = str(read1) 
        a1.flag = 3
        a1.reference_id = 0
        a1.reference_start = positions[i] + first_base1
        a1.mapping_quality = 255 # A value 255 indicates that the mapping quality is not available.
        a1.cigartuples = tuple(cigar1)
        a1.next_reference_id = 0 
        a1.next_reference_start = positions[(i+1) % N_reads]
        a1.template_length = templ_len1        
        a1.query_qualities = pysam.qualitystring_to_array(qual1)
        outf.write(a1)
    
        a2 = pysam.AlignedSegment()
        a2.query_name = '*' # ‘*’ indicates the information is unavailable
        a2.query_sequence = str(read2.complement())
        a2.flag = 19 
        a2.reference_id = 0 
        a2.reference_start = positions[i] + insert_size - ll2 + first_base2
        a2.mapping_quality = 255 # A value 255 indicates that the mapping quality is not available.
        a2.cigartuples = tuple(cigar2) 
        a2.next_reference_id = 0 
        a2.next_reference_start = positions[(i+1) % N_reads] 
        a2.template_length = templ_len2
        a2.query_qualities = pysam.qualitystring_to_array(qual2)
        outf.write(a2)
        
    f1.close()
    f2.close()
    outf.close()
#############################################################
    
path = "example_human_reference.fasta" # input
path = "PhiX_genome.fasta"

coverage = 30        
read_length = 200    

insert_size = 350   
avg_quality = 50   
max_diff    = 3     
error_rate_snip = 1 / 500  
error_rate_ins  = 2.8 * 10e-6  
error_rate_del  = 5.1 * 10e-6  


illumina_simulator(path, coverage, read_length, insert_size, avg_quality, \
                   max_diff, error_rate_snip, error_rate_ins, error_rate_del)
    
