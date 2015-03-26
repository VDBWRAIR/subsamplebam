import re
import subprocess as sp
import multiprocessing
import subsamplebam

''' Not used '''
def get_num_alignments( bamfile):
   #ref_len_str = !samtools view -c $bamfile
   result = sp.check_output(['samtools view -c {0}'.format(bamfile)], shell=True)
   return int(result)

def get_ref_length(bamref): 
    ''' 
    :param str bamref: A line representing a ref. sequence form `samtools view -H`
    :returns the length of the reference sequence. 
    '''
    return int(re.findall(r'.*LN:(.+)', bamref)[0])

def get_ref_seq(bamref): 
    '''
    :param str bamref: A line representing a ref. sequence form `samtools view -H`
    :returns the ref idea of the sequence listed in the .sam header. 
    '''
    return re.findall(r'.*SN:(.*)\t', bamref)[0]

def parse_header(bamfile):
    '''
    :param str bamfile: name of bamfile
    :returns dict a dictionary of reference ids and their length
    ''' 
    view = sp.check_output(["samtools view -H {0}".format(bamfile)], shell=True)
    lines = [line for line in view.split('\n') if line.startswith('@SQ')]
    return {get_ref_seq(line) : get_ref_length(line) for line in lines}

def subsample_file(bamfile, depth=1000, allow_orphans=False):
    '''
    Given a .bam file, create (possibly multiple) subsampled (sorted, indexed) .bam files for each ref. sequence in parallel 
    :param str bamfile: name of .bam file
    :side-effect: create (subsampled) <seqname>.bam files for each reference sequence in bamfile.
    '''
    refseqs = parse_header(bamfile)
    args = [refinfo + (bamfile, depth, allow_orphans) for refinfo in refseqs.iteritems() ]
    print args
    subsample_reference(*args[0])
    return


    pool = multiprocessing.Pool() 
    # use pool.map instead of for loop
    pool.map(subsample_multi, args) #(refname, length, bamfile, depth, allow_orphans))
        

def subsample_multi(args):
    return subsample_reference(*args)


def subsample_reference(refname, reflength, bamfile, depth, allow_orphans=False):
    '''
    For running multiple references in parallel
    :param str seqname: name of  reference sequence as listed in .sam header
    :param int length: length of the reference sequence
    :returns str new subsampled .sam file as a string.
    ''' 
    import os; 
    folder = 'outputs'
    if not os.path.exists(folder):
        os.mkdir(folder)
    refname = '"{0}"'.format(refname)
    orphans = '-A' if allow_orphans else '' 
    bamfile = '{0}'.format(bamfile)
    subsam_cmd = "python subsample_mindepth.py {0} {1} --subsample {2} {3}".format(bamfile, refname, depth, reflength, orphans)
    filename = '''"{0}/{1}.min.{2}"'''.format(folder, refname, str(depth))
    index_cmd  = "samtools view -hSb - | samtools sort - {0}; samtools index {0}.bam;".format(filename)
    print "saving to {0}.bam".format(filename)
    with open('stderr.log', 'w') as errlog:
        print(subsam_cmd)
        subsam = sp.check_output(subsam_cmd)
        print(subsam)
        print(subsam)
#        subsample_proc = sp.Popen(
#                subsam_cmd, 
#                stdout=sp.PIPE, 
#                stderr=errlog,
#                shell=True)
        #subsam = subsample_proc.communicate()[0]
        #with open('{0}.bam'.format(refname), 'w') as out_bamfile:
        index_process = sp.Popen(
                index_cmd, 
                stdout=sp.PIPE,
                stderr=errlog,
                stdin=subsam, #subsample_proc.stdout, 
                shell=True) 
        index_process.wait()
    return

if __name__ == "__main__":
    args = subsamplebam.parse_args(wrapper=True)
    subsample_file(args.bamfile, args.subsample, args.count_orphans)

