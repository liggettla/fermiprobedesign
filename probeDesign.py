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

        allSeqs[probe['name']] = finalSeq.upper()

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

                # Oligo position
                # Position is given as the innermost position such that up - down yields capture length
                upPos = start + 1 + oligolen
                downPos = start + 6 + 131

                up = pd.DataFrame([oligoUp, tmUp, homodimerUp, hairpinUp, endstablilityUp, 'Up', oligolen, gcClampUp, gcUp, probe, upPos], index=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe','Pos'])
                down = pd.DataFrame([oligoDown, tmDown, homodimerDown, hairpinDown, endstablilityDown, 'Down', oligolen, gcClampDown, gcDown, probe, downPos], index=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe','Pos'])

                fullDF = fullDF.append(up.T, ignore_index=True)
                fullDF = fullDF.append(down.T, ignore_index=True)

    # filter oligos by specified parameters
    fullDF = eliminateOligos(fullDF)

    fullDF = fullDF.reindex(columns=['Seq','Tm','Homo','Hairpin','Stability','Loc','Length','gcClamp','GC','Probe','Pos'])
    return fullDF

# this will eliminate oligos based on desired filtering parameters
def eliminateOligos(df):
    import pandas as pd
    finalDF = pd.DataFrame()

    # this is useless at the moment but it may be advantageous in the future
    # to make the following parameters dynamic in order to yield similar numbers
    # of up and down oligos
    for loc in ['Up', 'Down']:
        temp = df[df['Loc'] == loc]

        temp = temp[temp['Homo'] > -5000]
        temp = temp[temp['Hairpin'] > -2000]
        temp = temp[temp['GC'] > 0.404]
        temp = temp[temp['GC'] < 0.596]
        temp = temp[temp['gcClamp'] == 'good']

        # this range was calculated by taking the melting temps of original probes
        temp = temp[temp['Tm'] > 52.3]
        temp = temp[temp['Tm'] < 67.1]

        finalDF = finalDF.append(temp, ignore_index=True)

    return finalDF

# this will return the first quarter of the DataFrame sorted by 3 prime stability
# this is just a crude but simpler method of optimizing 3 prime stability
def sampleStability(df, seq):
    import pandas as pd
    filteredDF = pd.DataFrame()

    for probe in seq:
        for loc in ['Up', 'Down']:
            temp = df[df['Probe'] == probe]
            temp = temp[temp['Loc'] == loc]

            # sample the top quarter of oligos sorted in

            # descending order of 3' end stability ΔG
            if len(temp) > 12:
                sample = list(range(int(len(temp) / 4)))
            else:
                sample = list(range(int(len(temp))))

            filteredDF = filteredDF.append(temp.sort_values(['Stability'], ascending=False).take(sample), ignore_index=True)

    return filteredDF

# this will match up oligo pairs
def findPairs(df, seq):
    import pandas as pd
    import primer3

    pairsList = []

    # append to dataframe like this
    #x.loc[index,column]=num

    # df.iloc[i] returns the ith row of df. i does not refer to the index label, i is a 0-based index

    up = df[df['Loc'] == 'Up']
    down = df[df['Loc'] == 'Down']

    duplex = -100000

    for probe in seq:
        probeList = []
        oligosUp = up[up['Probe'] == probe]
        oligosDown = down[down['Probe'] == probe]

        for rowUp in oligosUp.itertuples():
            downSeq = oligosDown.copy(deep=True) # don't modify original df
            downSeq['Hetero'] = downSeq['Seq'].apply(lambda x: primer3.bindings.calcHeterodimer(x, rowUp[1]).dg)
            downSeq = downSeq.sort_values(['Hetero'], ascending=False)

            # add to list of dicts which will be converted to pandas df
            # this is much faster than individually adding to a DataFrame
            pairsList.append({'Up':rowUp[1], 'Down':downSeq.iloc[0]['Seq'], 'Probe':probe, 'Hetero':downSeq.iloc[0]['Hetero'], 'CapLen':downSeq.iloc[0]['Pos'] - rowUp[11]})

    #convert final pairs to DataFrame
    finalDF = pd.DataFrame(pairsList)
    finalDF = finalDF.reindex(columns=['Probe','Up','Down','CapLen','Hetero'])

    return finalDF

# ensure probes do not cross react
def finalizeProbeSet(df, seq):
    from collections import defaultdict
    import pandas as pd
    import primer3

    finalProbes = []

    # create holding data structure containing the best pair from each probe
    for probe in seq:
        probeDF = df[df['Probe'] == probe]
        probeDF = probeDF.sort_values(['Hetero'], ascending=False)
        finalProbes.append({'Probe':probe, 'Up':probeDF.iloc[0]['Up'], 'Down':probeDF.iloc[0]['Down'], 'CapLen':1})


    unfinished = True
    while unfinished:

        pass



















































