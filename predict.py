import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import pandas as pd
import os, re, math, platform
from collections import Counter
import numpy as np
from Bio import Seq, SeqIO
from deepforest import CascadeForestClassifier
from deepforest import CascadeForestRegressor
import json

import argparse


parser = argparse.ArgumentParser(description='Run Model')
parser.add_argument('-i', type=str,required=True, help='Input Fasta')
parser.add_argument('-o', type=str,required=True, help='Output Result')
parser.add_argument('-predSpecie', type=bool, default=False)
parser.add_argument('--Ecoli', type=bool, default=False)
parser.add_argument('--Saureus', type=bool, default=False)
parser.add_argument('--Paeruginosa', type=bool, default=False)
parser.add_argument('--Bsubtilis', type=bool, default=False)
parser.add_argument('--Sepidermidis', type=bool, default=False)
parser.add_argument('--Mluteus', type=bool, default=False)
parser.add_argument('--Kpneumoniae', type=bool, default=False)
parser.add_argument('--Efaecalis', type=bool, default=False)
parser.add_argument('--Styphimurium', type=bool, default=False)
parser.add_argument('--Abaumannii', type=bool, default=False)
parser.add_argument('--predActivity', type=bool, default=False)
args = parser.parse_args()
inputfile = args.i
outfile = args.o 


# Miscellaneous
_AALetter = ['A', 'C', 'D', 'E', 'F', 'G', 'H',
             'I', 'K', 'L', 'M', 'N', 'P', 'Q',
             'R', 'S', 'T', 'V', 'W', 'Y']


"""
n_gram statistics
"""

_AALetter = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 
            'I', 'K', 'L', 'M', 'N', 'P', 'Q', 
            'R', 'S', 'T', 'V', 'W', 'Y']

def get_aan_corpus(n=2):
    '''
    Get AA corpus of n_gram (e.g. Di, Tri, etc.)
    Output: AA n_gram corpus ((e.g. Di:400, Tri:3000, etc.))
    '''
    n_corpus = []
    if n <= 2:
        for i in _AALetter:
            for j in _AALetter:
               n_corpus.append("{}{}".format(i, j))
        return n_corpus
    for i in get_aan_corpus(n - 1):
        for j in _AALetter:
            n_corpus.append("{}{}".format(i, j))
    return n_corpus


def get_ngram_counts(seq, n=2):
    '''
    Get n_gram statistics
    Input: peptide sequence and n
    Ouput: n_gram statistic (dictionary) {A.A Corp: Counts}
    '''
    # Get the name of ngram feature
    if n == 2:
        prefix = 'DPC'
    elif n == 3:
        prefix = 'TPC'
    else:
        prefix = '{}gram'.format(n)

    ngrams = [seq[i: i + n] for i in range(len(seq) - n + 1)]
    n_corpus = get_aan_corpus(n)
    ngram_stat = {}
    for aa_ng in n_corpus:
        ngram_stat['{}|{}'.format(prefix, aa_ng)] = ngrams.count(aa_ng) / len(ngrams) * 100
    return ngram_stat


def minSequenceLength(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(i[1]):
            minLen = len(i[1])
    return minLen


def minSequenceLengthWithNormalAA(fastas):
    minLen = 10000
    for i in fastas:
        if minLen > len(re.sub('-', '', i[1])):
            minLen = len(re.sub('-', '', i[1]))
    return minLen


"""
    input.fasta:      the input protein sequence file in fasta format.
    k_space:          the gap of two amino acids, integer, defaule: 5
    output:           the encoding file, default: 'encodings.tsv'
"""


def generateGroupPairs(groupKey):
    gPair = {}
    for key1 in groupKey:
        for key2 in groupKey:
            gPair['CKSAAGP|'+key1+'.'+key2] = 0
    return gPair


def cksaagp(fastas, gap = 5, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    if minSequenceLength(fastas) < gap+2:
        print('Error: all the sequence length should be greater than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
        return 0

    group = {'aliphatic': 'GAVLMI',
             'aromatic': 'FYW',
             'postivecharge': 'KRH',
             'negativecharge': 'DE',
             'uncharge': 'STCPNQ'}

    AA = 'ARNDCQEGHILKMFPSTWYV'

    groupKey = group.keys()

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    gPairIndex = []
    for key1 in groupKey:
        for key2 in groupKey:
            gPairIndex.append('CKSAAGP|'+key1+'.'+key2)

    encodings = []
    header = ['#']
    for g in range(gap + 1):
        for p in gPairIndex:
            header.append(p+'.gap'+str(g))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for g in range(gap + 1):
            gPair = generateGroupPairs(groupKey)
            sum = 0
            for p1 in range(len(sequence)):
                p2 = p1 + g + 1
                if p2 < len(sequence) and sequence[p1] in AA and sequence[p2] in AA:
                    gPair['CKSAAGP|'+index[sequence[p1]]+'.'+index[sequence[p2]]] = gPair['CKSAAGP|'+index[sequence[p1]]+'.'+index[sequence[p2]]] + 1
                    sum = sum + 1

            if sum == 0:
                for gp in gPairIndex:
                    code.append(0)
            else:
                for gp in gPairIndex:
                    code.append(gPair[gp] / sum)

        encodings.append(code)

    return encodings


"""
    input.fasta:      the input protein sequence file in fasta format.
    lambda:           the lambda value, integer, defaule: 30
    output:           the encoding file, default: 'encodings.tsv'
"""

def Rvalue(aa1, aa2, AADict, Matrix):
    return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)


def paac(fastas, lambdaValue=30, w=0.05, **kw):
    if minSequenceLengthWithNormalAA(fastas) < lambdaValue + 1:
        print('Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
        return 0

    dataFile = re.sub('codes$', '', os.path.split(os.path.realpath('__file__'))[0]) + r'\data\PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath('__file__'))[0]) + '/data/PAAC.txt'
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
        AAProperty1.append([(j-meanI)/fenmu for j in i])

    encodings = []
    header = ['#']
    for aa in AA:
        header.append('PAAC|' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('PAAC|lambda' + str(n))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(
                sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    return encodings


def GAAC(fastas, **kw):

group = {


'alphatic': 'GAVLMI',


'aromatic': 'FYW',


'postivecharge': 'KRH',


'negativecharge': 'DE',


'uncharge': 'STCPNQ'

}


groupKey = group.keys()


encodings = []

header = ['#']

for key in groupKey:


header.append("GAAC|"+key)

encodings.append(header)


for i in fastas:


name, sequence = i[0], re.sub('-', '', i[1])


code = [name]


count = Counter(sequence)


myDict = {}


for key in groupKey:



for aa in group[key]:




myDict[key] = myDict.get(key, 0) + count[aa]



for key in groupKey:



code.append(myDict[key]/len(sequence))


encodings.append(code)


return encodings


def GDPC(fastas, **kw):

group = {


'alphaticr': 'GAVLMI',


'aromatic': 'FYW',


'postivecharger': 'KRH',


'negativecharger': 'DE',


'uncharger': 'STCPNQ'

}


groupKey = group.keys()

baseNum = len(groupKey)

dipeptide = [g1+'.'+g2 for g1 in groupKey for g2 in groupKey]


index = {}

for key in groupKey:


for aa in group[key]:



index[aa] = key


encodings = []

header = ['#'] + ['GDPC|'+dipname for dipname in dipeptide]

encodings.append(header)


for i in fastas:


name, sequence = i[0], re.sub('-', '', i[1])



code = [name]


myDict = {}


for t in dipeptide:



myDict[t] = 0



sum = 0


for j in range(len(sequence) - 2 + 1):



myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] + 1



sum = sum +1



if sum == 0:



for t in dipeptide:




code.append(0)


else:



for t in dipeptide:




code.append(myDict[t]/sum)


encodings.append(code)


return encodings


def GTPC(fastas, **kw):

group = {


'alphaticr': 'GAVLMI',


'aromatic': 'FYW',


'postivecharger': 'KRH',


'negativecharger': 'DE',


'uncharger': 'STCPNQ'

}


groupKey = group.keys()

baseNum = len(groupKey)

triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]


index = {}

for key in groupKey:


for aa in group[key]:



index[aa] = key


encodings = []

header = ['#'] + ['GTPC|'+tname for tname in triple]

encodings.append(header)


for i in fastas:


name, sequence = i[0], re.sub('-', '', i[1])



code = [name]


myDict = {}


for t in triple:



myDict[t] = 0



sum = 0


for j in range(len(sequence) - 3 + 1):



myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1



sum = sum +1



if sum == 0:



for t in triple:




code.append(0)


else:



for t in triple:




code.append(myDict[t]/sum)


encodings.append(code)


return encodings

'''
Insert Iso_electric Point and net_charge(neutral) feature to the sequence data_frame
Input: sequence data_frame {IDx: Seq_x, ...}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., iep, net_charge}
'''


def insert_phycs(seq_df):
    seq_df = seq_df.copy()
    #  Function for compute Isoelectric Point or net_charge of peptide
    def get_ieq_nc(seq, is_iep=True):
        protparam = PA(seq)
        return protparam.isoelectric_point() if is_iep else protparam.charge_at_pH(7.0)

    # Calculating IsoElectricPoints and NeutralCharge
    data_size = seq_df.size
    seq_df['PHYC|IEP'] = list(map(get_ieq_nc, seq_df['Sequence'], [True] * data_size))  # IsoElectricPoints
    seq_df['PHYC|Net Charge'] = list(map(get_ieq_nc, seq_df['Sequence'], [False] * data_size))  # Charge(Neutral)

    # Calculating hydrophobic moment (My assume all peptides are alpha-helix)
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'eisenberg')
    descrpt.calculate_moment(window=1000, angle=100, modality='max')
    seq_df['PHYC|Hydrophobic Moment'] = descrpt.descriptor.reshape(-1)

    # Calculating "Hopp-Woods" hydrophobicity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'hopp-woods')
    descrpt.calculate_global()
    seq_df['PHYC|Hydrophobicity'] = descrpt.descriptor.reshape(-1)

    # Calculating Energy of Transmembrane Propensity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'tm_tend')
    descrpt.calculate_global()
    seq_df['PHYC|Transmembrane Propensity'] = descrpt.descriptor.reshape(-1)

    # Calculating Aromaticity
    descrpt = GlobalDescriptor(seq_df['Sequence'].values)
    descrpt.aromaticity()
    seq_df['PHYC|Aromacity'] = descrpt.descriptor.reshape(-1)

    # Calculating Levitt_alpha_helical Propensity
    descrpt = PeptideDescriptor(seq_df['Sequence'].values, 'levitt_alpha')
    descrpt.calculate_global()
    seq_df['PHYC|Alpha Helical Propensity'] = descrpt.descriptor.reshape(-1)

    # Calculating Aliphatic Index
    descrpt = GlobalDescriptor(seq_df['Sequence'].values)
    descrpt.aliphatic_index()
    seq_df['PHYC|Aliphatic Index'] = descrpt.descriptor.reshape(-1)

    # Calculating Boman Index
    descrpt = GlobalDescriptor(seq_df['Sequence'].values)
    descrpt.boman_index()
    seq_df['PHYC|Boman Index'] = descrpt.descriptor.reshape(-1)

    return seq_df


'''
Insert Amino acid composition to the sequence data_frame
Input: sequence data_frame {IDx: Seq_x}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., AAC_Ax ... AAC_Yx}
'''


def insert_aac(seq_df):
    seq_df = seq_df.copy()
    # Compute AAC for peptide in specific A.A
    def get_aac(seq, aa):
        return seq.count(aa) / len(seq) * 100

    # processing data_frame
    data_size = seq_df.size
    for ll in _AALetter:
        seq_df['AAC|{}'.format(ll)] = list(map(get_aac, seq_df['Sequence'], [ll] * data_size))
    return seq_df


'''
Insert n_grams Descriptor to the sequence data_frame
Input: sequence data_frame {IDx: Seq_x, ...}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., ngram_(1), .., ngram(20 ** n)}
'''


def insert_ngrams(seq_df, n=2):
    seq_df = seq_df.copy()
    data_size = seq_df.size

    ngrams_df = list(map(get_ngram_counts, seq_df['Sequence'], [n] * data_size))
    ngrams_df = pd.DataFrame(ngrams_df)  # Convert ngrams features to pd.DataFrame
    seq_df = pd.concat([seq_df, ngrams_df], axis=1)
    return seq_df


"""
Insert CKSAAGP Descriptor to the sequence data_frame
(Composition of k-spaced amino acid pairs)
Input: sequence data_frame {IDx: Seq_x, ...}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., CKSAAGP(0), ..., CKSAAGP(m)}
"""


def insert_cksaagp(seq_df, gap=2):
    seq_df.copy()
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = cksaagp(fastas, gap=gap)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]], axis=1)
    return seq_df


"""
Insert PAAC Descriptor to the sequence data_frame
(Pseudo Amino Acid Composition)
Input: sequence data_frame {IDx: Seq_x, ...}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., CKSAAGP(0), ..., CKSAAGP(m)}
"""


def insert_paac(seq_df, lamb=3, w=0.1):
    seq_df = seq_df.copy()
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    encoding = paac(fastas, lambdaValue=lamb, w=w)
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]], axis=1)
    return seq_df


"""
Insert Grouped n-gram Descriptor to the sequence data_frame
(Pseudo Amino Acid Composition)
Input: sequence data_frame {IDx: Seq_x, ...}
Output: data_frame of Peptide Seq {IDx: Seq_x, ..., CKSAAGP(0), ..., CKSAAGP(m)}
"""


def insert_gngram(seq_df, n=1):
    seq_df = seq_df.copy()
    fastas = [[idx, seq] for idx, seq in zip(seq_df['Id'], seq_df['Sequence'])]
    # encoding = paac(fastas, lambdaValue=lamb, w=w)
    if n is 1:
        encoding = GAAC(fastas)
    elif n is 2:
        encoding = GDPC(fastas)
    elif n is 3:
        encoding = GTPC(fastas)
    else:
        raise Warning("Invalid n-grams, no features added")
    encoding = pd.DataFrame(encoding[1:], columns=encoding[0])
    seq_df = pd.concat([seq_df, encoding.iloc[:, 1:]], axis=1)
    return seq_df


def construct_features(seq_df, paaclamb=4, paacw=0.5):
    """
    Construct Features for the AVPIden. We first investigated physiochemical feautres, AAC features, DiC features, 
    CKSAAGP features, PAAC features, and PHYC features.
    Parameters are pre-set according to the sequence identities.
    """
    seq_df = insert_aac(seq_df)
    seq_df = insert_ngrams(seq_df, n=2)
    seq_df = insert_cksaagp(seq_df, gap=3) # As the maximum motif length = 5.
    seq_df = insert_paac(seq_df, lamb=paaclamb, w=paacw)
    seq_df = insert_phycs(seq_df)
    return seq_df


def write_fasta(df, file_path, abbr_columns=None):
    """
    Save dataframe to a .fasta file, the df should contain at least columns named "Id" and "Sequence"
    
    df: dataframe for saving .fasta
    file_path: path(string) for saving the fasta file
    abbr_columns: string columns for adding abbreviations. Multiple abbr are splited by '|'.
    """
    Seqrecords = [SeqIO.SeqRecord(id=row['Id'], 
                              seq=Seq.Seq(row['Sequence']), 
                              description='|'.join(row[abbr_columns] if abbr_columns is not None else "")) \
             for idn, row in df.iterrows()]
    with open(file_path, 'w+') as fhandle:
        SeqIO.write(Seqrecords, fhandle, "fasta-2line")
        print("Saved {:d} sequences.".format(len(Seqrecords)))


def read_fasta(fname):
    '''
    Read fasta file to dictionary
    Input: path name of fasta
    Output: dataframe of Peptide Seq {ID1: Seq1, ID2: Seq2,...}
    '''
    with open(fname, "rU") as f:
        seq_dict = [(record.id, record.seq._data) for record in SeqIO.parse(f, "fasta")]
    seq_df = pd.DataFrame(data=seq_dict, columns=["Id", "Sequence"])
    return seq_df


def predictActivity(featureIndex,species, index, isPredActivity,threshold):

    if len(index) != 0:
        model_path_classification = "/home/AMPActiPred/model/DeepForest-multi-label/{}".format(species) #Please change this path to the path where you save the model.
        seq_df = featureIndex.iloc[index,:]
        seqFeature = featureIndex.iloc[index,2:].values

        model_multiLabel = CascadeForestClassifier()
        model_multiLabel.load(model_path_classification)
        y_prob = model_multiLabel.predict_proba(seqFeature)
        result_list = []

        if isPredActivity == True:
            model_path_regression = "/home/AMPActiPred/model/DeepForest-AMP-Regression/{}".format(species) #Please change this path to the path where you save the model.
            model_regression = CascadeForestRegressor()
            model_regression.load(model_path_regression)
            y_regression_pred = model_regression.predict(seqFeature).squeeze()
            print(y_regression_pred)
            print(y_regression_pred.shape==())
            print(len(y_prob))
        for i in range(0,len(y_prob)):
            if y_prob[i][1] >= threshold and isPredActivity == False:
                result_list.append(1)
            if y_prob[i][1] >= threshold and isPredActivity == True:
                if len(y_prob)==1:
                    result_list.append(float(y_regression_pred))
                else:
                    result_list.append(y_regression_pred[i])
            if y_prob[i][1] < threshold:
                result_list.append(0)
        print(result_list)
        return result_list
    else:
        return []





seqlist = []
idlist = []

for record in SeqIO.parse(inputfile,"fasta"):
    idlist.append(record.id)
    seqlist.append(str(record.seq))


df_bioseq = read_fasta(inputfile)
features_bioseq = construct_features(df_bioseq)
#features_bioseq.to_csv("/home/AMPActiPred/public_html/data/multiLabel.csv")

model_first_stage = CascadeForestClassifier()
model_first_stage.load("/home/AMPActiPred/model/first_stage_model_final") #Please change this path to the path where you save the model.

seq_feature = features_bioseq.iloc[:,2:].values
first_stage_pred = model_first_stage.predict(seq_feature)

isABPIndex = [i for i,x in enumerate(first_stage_pred) if x==1]
print(isABPIndex)

if args.predSpecie==True:
    if args.Ecoli == True:
        Ec_list = predictActivity(features_bioseq,"Ec",isABPIndex,args.predActivity,0.84)
    if args.Saureus == True:
        Sa_list = predictActivity(features_bioseq,"Sa",isABPIndex,args.predActivity,0.85)
    if args.Paeruginosa == True:
        Pa_list = predictActivity(features_bioseq,"Pa",isABPIndex,args.predActivity,0.48)    
    if args.Bsubtilis == True:
        Bs_list = predictActivity(features_bioseq,"Bs",isABPIndex,args.predActivity,0.35)           
    if args.Sepidermidis == True:
        Se_list = predictActivity(features_bioseq,"Se",isABPIndex,args.predActivity,0.2)         
    if args.Mluteus == True:
        Ml_list = predictActivity(features_bioseq,"Ml",isABPIndex,args.predActivity,0.15)     
    if args.Kpneumoniae == True:
        Kp_list = predictActivity(features_bioseq,"Kp",isABPIndex,args.predActivity,0.13)           
    if args.Efaecalis == True:
        Ef_list = predictActivity(features_bioseq,"Ef",isABPIndex,args.predActivity,0.12)          
    if args.Styphimurium == True:
        St_list = predictActivity(features_bioseq,"St",isABPIndex,args.predActivity,0.2)            
    if args.Abaumannii == True:
        Ab_list = predictActivity(features_bioseq,"Ab",isABPIndex,args.predActivity,0.05)     

finalList = []
EcIndex = SaIndex = PaIndex = BsIndex = StIndex = KpIndex = MlIndex = EfIndex = SeIndex = AbIndex = 0
for i in range(0, len(seq_feature)):
    result_obj = {}
    result_obj["Index"] = i
    result_obj["seqID"] = features_bioseq.iloc[i,:].Id
    result_obj["Sequence"] = features_bioseq.iloc[i,:].Sequence
    result_obj["isABP"] = True if first_stage_pred[i] == 1 else False
    result_obj["predSpecie"] = True if args.predSpecie == True else False
    result_obj["predActivity"] = True if args.predActivity == True else False
    result_obj["Ec"] = result_obj["Sa"] = result_obj["Pa"] = result_obj["Bs"] = result_obj["St"] = result_obj["Kp"] = result_obj["Ml"] = result_obj["Ef"] = result_obj["Se"] = result_obj["Ab"] = ""
    if first_stage_pred[i] == 1 and args.Ecoli == True:
        result_obj["Ec"] = Ec_list[EcIndex]
        EcIndex = EcIndex + 1

    if first_stage_pred[i] == 1 and args.Saureus == True:
        result_obj["Sa"] = Sa_list[SaIndex]
        SaIndex = SaIndex + 1
    
    if first_stage_pred[i] == 1 and args.Paeruginosa == True:
        result_obj["Pa"] = Pa_list[PaIndex]
        PaIndex = PaIndex + 1

    if first_stage_pred[i] == 1 and args.Bsubtilis == True:
        result_obj["Bs"] = Bs_list[BsIndex]
        BsIndex = BsIndex + 1

    if first_stage_pred[i] == 1 and args.Sepidermidis == True:
        result_obj["Se"] = Se_list[SeIndex]
        SeIndex = SeIndex + 1

    if first_stage_pred[i] == 1 and args.Mluteus == True:
        result_obj["Ml"] = Ml_list[MlIndex]
        MlIndex = MlIndex + 1

    if first_stage_pred[i] == 1 and args.Kpneumoniae == True:
        result_obj["Kp"] = Kp_list[KpIndex]
        KpIndex = KpIndex + 1

    if first_stage_pred[i] == 1 and args.Efaecalis == True:
        result_obj["Ef"] = Ef_list[EfIndex]
        EfIndex = EfIndex + 1

    if first_stage_pred[i] == 1 and args.Styphimurium == True:
        result_obj["St"] = St_list[StIndex]
        StIndex = StIndex + 1

    if first_stage_pred[i] == 1 and args.Abaumannii == True:
        result_obj["Ab"] = Ab_list[AbIndex]
        AbIndex = AbIndex + 1

    finalList.append(result_obj)

with open(outfile, 'w') as f:
    json.dump(finalList, f,sort_keys=True, indent=4, separators=(',', ': '))

(path, filename) = os.path.split(outfile)
features_bioseq.to_csv(path+"/extracted_feature.csv",index=False)
