from sequenceAnalyzer import FastAreader
import argparse 

parser = argparse.ArgumentParser(description='Summarize the distribution of ambiguity and missing data in a SARS-CoV-2 phylogeny') 
parser.add_argument('-msa', nargs='?', required=True,help='Path to msa file for SARS-CoV-2 phylogeny')
parser.add_argument('-min_MAF', nargs='?', required=True,
                    help='Minimum minor allele frequency for reported sites') 
args = vars(parser.parse_args())

ambDic = {}
NDic = {}
seqDic = {}
myReaderRef = FastAreader('NC_045512v2.fa')
for header, sequence in myReaderRef.readFasta():
    refHead = header
    refSeq = sequence.upper()
ambDicPerSite = {}
NDicPerSite = {}
myReaderFasta = FastAreader(args['msa'])
for header, sequence in myReaderFasta.readFasta():
    if sequence.replace('-','').upper() == refSeq:
        refSeq = sequence
    seqDic[header] = sequence.upper()
    NDic[header] = 0
    ambDic[header] = 0 
    ambDicPerSite[header] = {}
    NDicPerSite[header] = {}

NambDicPerSite = {}
altDic = {}
charDic = {'A':'','G':'','C':'','T':''}
refDic = {}
baseCount = 0
NCount = 0
ambCount = 0
NaltDic  ={}
siteAmbDic = {}
siteNDic = {}
for i in range(len(refSeq)-1):
    if refSeq[i] != '-':
        ambCount = 0
        altCount = 0
        NCount = 0
        baseCount += 1
        refDic[baseCount] = refSeq[i]
        for head in seqDic:
            alt = seqDic[head][i].upper()
            altDic[head] = {}
            NaltDic[head] = {}
#            ambDicPerSite[head] = {}
            if alt == 'N' or alt == '-':
                NDic[head] += 1
                NCount += 1
                NDicPerSite[head][baseCount] = 'N'
                ambDicPerSite[head][baseCount] = 'N'
            elif alt not in charDic:
                ambDic[head] += 1
                ambCount += 1
                ambDicPerSite[head][baseCount] = alt

            elif refSeq[i].upper() != alt:
                altCount += 1 

        Freq = float(altCount)/float(len(seqDic))
        if Freq >= float(args['min_MAF']) and Freq < (1.0-float(args['min_MAF'])):
            siteAmbDic[baseCount] = Freq
mafCount = {}
NmafCount = {}
for head in ambDicPerSite:
    mafCount[head] = 0
    NmafCount[head] = 0 
    for i in siteAmbDic:
        if i in NDicPerSite[head]:
            NmafCount[head] += 1
        elif i in ambDicPerSite[head]:
            mafCount[head] += 1

 
with open('ambiguity.tsv','w') as f:
    f.write('Sample\tAmb_Count\tN_Count\tAmb_Count_MAF>={0}\tN_Count_MAF>={0}\t'.format(str(args['min_MAF'])))
    for site in siteAmbDic:
        f.write('{1}_{0}_{2}\t'.format(site,refDic[site],siteAmbDic[site]))
    f.write('\n')

    for sample in ambDic:
        f.write('{0}\t{1}\t{2}\t{3}\t{4}\t'.format(sample,ambDic[sample],NDic[sample],mafCount[sample],NmafCount[sample]))
        for site in siteAmbDic:
            if site in ambDicPerSite[sample]:
                f.write('{0}\t'.format(ambDicPerSite[sample][site]))
            else:
                f.write('0\t')
        f.write('\n')


