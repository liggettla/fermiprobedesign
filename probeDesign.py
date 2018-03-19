#!/usr/bin/env python

# this script provides the functionality necessary to reverse
# engineer Illumina TSCA probes

def getSeq(ds):
    from subprocess import check_output, STDOUT

    allSeqs = {}
    for probe in ds['seqs']:
        low = probe['loc'] - 125
        high = probe['loc'] + 125
        temp = check_output('samtools faidx %s %s:%s-%s' % (ds['ref'], probe['chrom'], low, high), stderr=STDOUT, shell=True)

        finalSeq = ''
        for line in temp.split():
            line = line.decode('UTF-8')
            if '>' not in line:
                finalSeq += line

        allSeqs[probe['name']] = finalSeq

    return allSeqs

def getTargets(seq):
    import primer3
    allprimers = {}

    for name in seq:

        primers = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': name,
                'SEQUENCE_TEMPLATE': seq[name]
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

        allprimers[name] = primers
    return allprimers

# this will build a pandas DF with a bunch of possible oligos and thermo info about them
# this method will also sort the oligos based on a number of desireabl criteria
def possibleOligos(seq):
    import pandas as pd
    import primer3

    fullDF = pd.DataFrame(index=['Seq','Homo','Haripin','Stablity'])

    for probe in seq:
        template = seq[probe]
        upstream = seq[probe][:120]
        downstream = seq[probe][-120:]

        for oligolen in range(22,31):
            for start in range(len(upstream) - (oligolen + 5)):
                oligoUp = upstream[start+1:][:oligolen] # don't start at position 0
                oligoDown = downstream[start+6:][:oligolen] # this should start 5bp away from 5' end to avoid touching mutation

                # melting temperature
                tmUp = primer3.bindings.calcTm(oligoUp)
                tmDown = primer3.bindings.calcTm(oligoDown)

                # homodimerization
                homodimerUp = primer3.bindings.calcHomodimer(oligoUp)
                homodimerDown = primer3.bindings.calcHomodimer(oligoDown)
                if homodimerUp.structure_found:
                    homodimerUp = homodimerUp.dg
                else:
                    homodimerUp = 1000000. # if no structure set very high ΔG
                if homodimerDown.structure_found:
                    homodimerDown = homodimerDown.dg
                else:
                    homodimerDown = 1000000.

                # secondary structure
                hairpinUp = primer3.bindings.calcHairpin(oligoUp)
                hairpinDown = primer3.bindings.calcHairpin(oligoDown)
                if hairpinUp.structure_found:
                    hairpinUp = hairpinUp.dg
                else:
                    hairpinUp = 1000000. # if no structure set very high ΔG
                if hairpinDown.structure_found:
                    hairpinDown = hairpinDown.dg
                else:
                    hairpinDown = 1000000.

                # 3' end stability
                endstablilityUp = primer3.bindings.calcEndStability(oligoUp, template)
                endstablilityDown = primer3.bindings.calcEndStability(oligoDown, template)
                if endstablilityUp.structure_found:
                    endstablilityUp = endstablilityUp.dg
                else:
                    endstablilityUp = -1000000. # if no structure set very low ΔG
                if endstablilityDown.structure_found:
                    endstablilityDown = endstablilityDown.dg
                else:
                    endstablilityDown = -1000000.

                # GC clamp
                if oligoUp[-5:].count('G') + oligoUp[-5:].count('C') > 0 and oligoUp[-5:].count('G') + oligoUp[-5:].count('C') < 4:
                    gcClampUp = 'good'
                else:
                    gcClampUp = 'poor'
                if oligoDown[-5:].count('G') + oligoDown[-5:].count('C') > 0 and oligoDown[-5:].count('G') + oligoDown[-5:].count('C') < 4:
                    gcClampDown = 'good'
                else:
                    gcClampDown = 'poor'

                # GC content
                gcUp = 1
                gcDown = 1
                gcUp = (oligoUp.count('G') + oligoUp.count('C')) / len(oligoUp)
                gcDown = (oligoDown.count('G') + oligoDown.count('C')) / len(oligoDown)



                up = pd.DataFrame([oligoUp, tmUp, homodimerUp, hairpinUp, endstablilityUp, 'Up', oligolen, gcClampUp, gcUp, probe], index=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe'])
                down = pd.DataFrame([oligoDown, tmDown, homodimerDown, hairpinDown, endstablilityDown, 'Down', oligolen, gcClampDown, gcDown, probe], index=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe'])

                fullDF = fullDF.append(up.T, ignore_index=True)
                fullDF = fullDF.append(down.T, ignore_index=True)

    # eliminate undesireable oligos
    fullDF = fullDF[fullDF['Homo'] > -5000]
    fullDF = fullDF[fullDF['Hairpin'] > -2000]
    fullDF = fullDF[fullDF['GC'] > 0.404]
    fullDF = fullDF[fullDF['GC'] < 0.596]
    fullDF = fullDF[fullDF['gcClamp'] == 'good']
    #fullDF = fullDF[fullDF['Tm'] < 66.7]
    #fullDF = fullDF[fullDF['Tm'] > 73.2]
    fullDF = fullDF[fullDF['Tm'] > 60]



    fullDF = fullDF.reindex(columns=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe'])
    return fullDF

# this will return the first quarter of the DataFrame sorted by 3 prime stability
# this is just a crude but simpler method of optimizing 3 prime stability
def sampleStability(df):
    import pandas as pd

    oligosUp = up[up['Probe'] == probe]
    oligosUp = sort_values
    oligosDown = up[up['Probe'] == probe]

    #someting like this:

    # sample the top quarter of oligos sorted in
    # descending order of 3' end stability ΔG
    up = list(range(int(len(oligosUp) / 4)))
    oligosUp = oligosUp.sort_values(['Stability'], ascending=False).take(up)
    down = list(range(int(len(oligosDown) / 4)))
    oligosDown = oligosDown.sort_values(['Stability'], ascending=False).take(down)

    return oligosUp, oligosDown

# this will match up oligo pairs
def findPairs(up, down, seq):
    import pandas as pd
    import primer3

    for probe in seq:
        probeList = []
        oligosUp = up[up['Probe'] == probe]
        oligosDown = up[up['Probe'] == probe]

        for rowUp in up.itertuples():
            duplex = -100000
            currentPair = {}
            # row = Pandas(Index=1534, Seq='AGGGTGGAGGAAAGTTGGAATTCACAAACA', Tm=60.410845774884194, Homo=-4457.9524890635075, Hairpin=1000000.0,
            # Stability=878.8465328218736, Loc='Up', Length=30, gcClamp='good', GC=0.43333333333333335, Probe='TIIIN')
            oligoUp = rowUp[1]
            stabilityUp = rowUp[5]

            for rowDown in up.itertuples():
                # row = Pandas(Index=1534, Seq='AGGGTGGAGGAAAGTTGGAATTCACAAACA', Tm=60.410845774884194, Homo=-4457.9524890635075, Hairpin=1000000.0,
                # Stability=878.8465328218736, Loc='Up', Length=30, gcClamp='good', GC=0.43333333333333335, Probe='TIIIN')
                oligoDown = rowDown[1]
                stabilityDown = rowDown[5]

                dg = primer3.bindings.calcHeterodimer(oligoUp, oligoDown).dg
                # this finds the best pair by ΔG for each oligo
                if duplex < dg:
                    currentPair['Up'] = oligoUp
                    currentPair['Down'] = oligoDown
                    currentPair['Up3Prime'] = stabilityUp
                    currentPair['Down3Prime'] = stabilityDown


















































