

def test_get_overlap_candidates():
    pass

def test_min_depth_with_enough_alignments_no_backtracking():
    pass

def test_min_depth_with_backtracking():
    pass

def test_get_raw_reads_whole_file():
    pass

def test_get_raw_reads_empty_region():
    pass

def test_get_raw_reads_nonempty_region():
    pass

def test_get_next_reads():
    pass

def test_get_depths_on_empty_reads():
    pass

def test_get_depths():
    pass


   bamfile = '~/projects/samtools_primer/tutorial/alignments/sim_reads_aligned.sorted.bam'

   ref_len_str = !samtools view -c $bamfile
   ref_length = int(ref_len_str[0])


        #candidate_sequences = [filter(lambda q: q.seq_length + i >= under_index, seqs) for i, seqs in enumerate(self.seq_matrix[:under_index]) ]  
        #candidate_sequences = get_candidate_sequences(under_index)
        #farthest_overlap_seq = max(candidate_sequences, key=lambda seq: seq.overlap)
        #farthest_overlap_seq = max(filter(lambda seq: seq.picked, candidate_sequences), key=lambda seq: seq.overlap)
'''
    def yield_next_candidate_row(self, under_index):
        rows = []
        for i, depth in self.depth_array[:under_index+1]:
            if depth >= self.min_depth:
                rows.append(self.seq_matrix[i])
            else:
                break
        return rows
'''
def new_subsample(bamfile, min_depth=20, ref=None, start=None, end=None):
   #ref_len_str = !wc -l $bamfile | grep -Eo ^[0-9]+ #samtools view -c $bamfile
   ref_len_str = !samtools view -c $bamfile
   ref_length = int(ref_len_str[0])
   matrix = DepthMatrix(ref_length, min_depth)
   all_alignments = map(parse_alignment, get_raw_reads(bamfile))
   for pos, depth in enumerate(running_depths):
       needed_depth = target_depth - int(depth)
       # for some reason there are position 0 alignments at the end of the samtools_primer file...  not sure why the graph is ignoring them
       if pos == 0 or needed_depth == 0: continue
       all_next_reads = get_next_reads(all_alignments, pos) 
       # random.sample expects N smaller than or equal size of colleciton
       sample_size = needed_depth if needed_depth <= len(all_next_reads) else len(all_next_reads)
       next_reads = random.sample(all_next_reads, sample_size)
       if next_reads:
          next_depths = np.array(get_depths(next_reads))
          #raise Exception
          running_depths[pos:pos+len(next_depths)] += next_depths[:ref_length-pos]
   return running_depths

def naive_subsample_bam_depths(bamfile, target_depth=20, ref=None, start=None, end=None):
   #ref_len_str = !wc -l $bamfile | grep -Eo ^[0-9]+ #samtools view -c $bamfile
   ref_len_str = !samtools view -c $bamfile
   ref_length = int(ref_len_str[0])
   running_depths = np.zeros(ref_length)
   all_alignments = map(parse_alignment, get_raw_reads(bamfile))
   start_reads = get_next_reads(all_next_reads)
   for pos, depth in enumerate(running_depths):
       needed_depth = target_depth - int(depth)
       # for some reason there are position 0 alignments at the end of the samtools_primer file...  not sure why the graph is ignoring them
       if pos == 0 or needed_depth == 0: continue
       # random.sample expects N smaller than or equal size of colleciton
       sample_size = needed_depth if needed_depth <= len(all_next_reads) else len(all_next_reads)
       next_reads = random.sample(all_next_reads, sample_size)
       if next_reads:
          next_depths = np.array(get_depths(next_reads))
          #raise Exception
          running_depths[pos:pos+len(next_depths)] += next_depths[:ref_length-pos]
   return running_depths


def get_bam_depths(bamfile, ref=None, start=None, end=None):
   #ref_len_str = !wc -l $bamfile | grep -Eo ^[0-9]+ #samtools view -c $bamfile
   ref_len_str = !samtools view -c $bamfile
   ref_length = int(ref_len_str[0])
   running_depths = np.zeros(ref_length)
   all_alignments = map(parse_alignment, get_raw_reads(bamfile))
   for pos, depth in enumerate(running_depths):
       if pos == 0: continue
       next_reads = get_next_reads(all_alignments, pos) 
       if next_reads:
          next_depths = np.array(get_depths(next_reads))
          #raise Exception
          running_depths[pos:pos+len(next_depths)] += next_depths[:ref_length-pos]
   return running_depths
