#!/usr/bin/env python

# this script provides the functionality necessary to reverse
# engineer Illumina TSCA probes

def getSeq(ds):
    from subprocess import check_output, STDOUT
    temp = check_output('samtools faidx %s %s:%s-%s' % (ds['ref'], ds['chrom'], ds['low'], ds['high']), stderr=STDOUT, shell=True)

    finalSeq = ''
    for line in temp.split():
        line = line.decode('UTF-8')
        if '>' not in line:
            finalSeq += line

    #finalSeq = finalSeq.upper()
    return finalSeq

def getTargets(ID, seq):
    import primer3
    primers = primer3.bindings.designPrimers(
	{
	    'SEQUENCE_ID': ID,
	    'SEQUENCE_TEMPLATE': seq
	},
	{
	    'PRIMER_OPT_SIZE': 27,
	    'PRIMER_PICK_INTERNAL_OLIGO': 1,
	    'PRIMER_INTERNAL_MAX_SELF_END': 8,
	    'PRIMER_MIN_SIZE': 22,
	    'PRIMER_MAX_SIZE': 30,
	    'PRIMER_OPT_TM': 69.9,
	    'PRIMER_MIN_TM': 66.7,
	    'PRIMER_MAX_TM': 73.2,
	    'PRIMER_MIN_GC': 40.4,
	    'PRIMER_MAX_GC': 59.6,
	    'PRIMER_MAX_POLY_X': 100,
	    'PRIMER_INTERNAL_MAX_POLY_X': 100,
	    'PRIMER_SALT_MONOVALENT': 50.0,
	    'PRIMER_DNA_CONC': 50.0,
	    'PRIMER_MAX_NS_ACCEPTED': 0,
	    'PRIMER_MAX_SELF_ANY': 12,
	    'PRIMER_MAX_SELF_END': 8,
	    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
	    'PRIMER_PAIR_MAX_COMPL_END': 8,
	    'PRIMER_NUM_RETURN': 20
	})

    return primers

