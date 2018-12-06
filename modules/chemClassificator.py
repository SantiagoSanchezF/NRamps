''' Jose Santiago Sanchez Fragoso   -   2018

    chemClassificator   - AntiMicrobial propetide cleavage site predictorand classificator
    
    DEPPENDENCIES:
                    - Joblib
                    - Json
    EXAMPLES:
                    fasta = './pitd/proPitdApdAllDB.fasta'
                    model = './pitd/lun23/cleaveSVM.pkl'
                    with open("/scratch/santiago/cleavegeSite/pltcys/preforal/AAencode.json","r") as myfile:
                        AA = json.load(myfile)
                    cs = CleavageSite()
                    cs.predict_mature_seqs(fasta,model,AA) 
                    cs.print_mature_fasta()
'''
import json, re, getopt, sys
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
from rpy2.robjects.packages import importr

'''    FUNCTIONS    '''

class ChemClassificator:
    def __init__(self):
        self.df = {}
    def predict_mature_seqs(self,seqs,modelN,modelC):
        self.clfN = joblib.load(modelN)
        self.clfC = joblib.load(modelC)
        df = predict_Annot_matures(seqs,self.clfN,self.clfC)
        self.df['mature'] = df['mature']
        self.df['mature_annot'] = df['mature_annot']
        return self.df['mature'],self.df['mature_annot']
    def predict_amps(self,seqs, model, pep):
        self.dfchem = calc(seqs,pep)
        x = pd.DataFrame(self.dfchem)
        self.model = joblib.load(model)
        probs = self.model.predict_proba(x)
        abso = [[x for x in probs[n]].index(max(probs[n])) for n in range(len(probs))]
        probs = [round(max(probs[n]),2) for n in range(len(probs))]
        return abso,probs

## create a dataFrame with ids and sequences
def fasta_to_df(fasta):
    df = {}
    df['seqs'] = []
    df['ids'] = [] 
    ident = None
    seq = ""
    with open(fasta,"r") as myf:
        lines = myf.read().splitlines()
        while(not lines[-1].strip()):
            lines.pop() 
        last_line = lines[-1]
        for line in lines:
            i = line.split()[0]
            if re.search("^>",i):
                if ident and not len(seq) == 0:
                    df['ids'].append(ident)
                    df['seqs'].append(seq)
                ident = i.strip().split(">")[1]
                seq = ""
            if len(df['ids']) != len(set(df['ids'])):
                raise ValueError("In file " + fasta + ": \">" + ident + "\" identifier is used for more than one sequence. Identifiers must be unique. \nFirst non space characters are considered the id")
                return -1
            if not re.search("^>",i) and not re.search("[^ACEDGFIHKMLNQPSRTWVY]",i): 
                seq = seq + i.strip()
            if not re.search("^>",i) and re.search("[^ACEDGFIHKMLNQPSRTWVY]",i):
                raise ValueError("In file " + fasta + ": Only standar aminoacids one-letter code (ACEDGFIHKMLNQPSRTWVY) are allowed, Miss-use in sequence " + i + "\n")
                return -1
            if line == last_line and not len(seq) == 0:
                df['ids'].append(ident)
                df['seqs'].append(seq)
    return df 

## encode an octamer usig the given encoder
def code(octa):
     coder = {'A': [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'C': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'E': [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'D': [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'G': [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'F': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
 'I': [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'H': [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'K': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
 'M': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
 'L': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'N': [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'Q': [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'P': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
 'S': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
 'R': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
 'T': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
 'W': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
 'V': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
 'Y': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]}
     li = []
     for i in octa:
         li = li + coder[i]
     return li
##encode the 8 aa of termini given a seqeunce, NOT BELONGING TO CHEM CLASS BUT USES THE SAME FUNCTIONS
def encode_termini(seq):                                                      
    C = seq[len(seq)-8:len(seq)]                  
    N = seq[0:8]                                                              
    termini = code(N) + code(C)                                            
    return termini   

## get a list with the encoded octamers with-in a sequence
def seq_encode(seq):
     li = []
     for i in range(len(seq) - 7):
         li.append(code(seq[i:i+8]) + [i])
     return li
## predict a mature sequence gven the models of N/C-terminus Cleavage Site
def predict_mature(seq,modelN,modelC):
    if len(seq) < 8:
        return seq,"N-ter_CS negative","C-ter_CS negative"
    seq_encoded = seq_encode(seq)
    probsN = [x for x in modelN.predict_proba(seq_encoded)]
    probsC = [x for x in modelC.predict_proba(seq_encoded)]
    wtrmk = (-1,0,-1)        
    for pos in range(len(probsN)):
        for clus in range(len(probsN[0])):
            if probsN[pos][clus] > wtrmk[2] and clus !=0:
                wtrmk = (pos,clus,probsN[pos][clus])
    dum = [x for x in probsN[wtrmk[0]]]
    real = dum.index(max(dum))
    maxN = real
    indN = wtrmk[0] + 4
    probN = wtrmk[2]

    wtrmk = (-1,0,-1)        
    for pos in range(len(probsC)):
        for clus in range(len(probsC[0])):
            if probsC[pos][clus] > wtrmk[2] and clus !=0:
                wtrmk = (pos,clus,probsC[pos][clus])
    dum = [x for x in probsC[wtrmk[0]]]
    real = dum.index(max(dum))
    maxC = real
    indC = wtrmk[0] + 4
    probC = wtrmk[2]
    if maxN != 0:
        annotN = "N-ter_CS positive at position: " + str(indN) + ", proba:" + str(round(probN,2)) + " // "
    else:
        annotN = "N-ter_CS negative" + " // "
        indN = 0
    if maxC != 0:
        annotC = "C-ter_CS positive at position: " + str(indC) + ", CS-type:" + str(round(probC,2))
    else:
        annotC = "C-ter_CS negative" 
        indC = len(seq)
    if len(seq[indN:indC]) > 0:
        return seq[indN:indC],annotN,annotC
    else:
        return seq,annotN,annotC
## predict the mature sequences of a list, and annotate results
def predict_Annot_matures(seqs,modelN,modelC):                   
    dct = {}                                                            
    dct['mature'] = []                                                  
    dct['mature_annot'] = []
    for i in range(len(seqs)):
        a,b,c = predict_mature(seqs[i],modelN,modelC)
        dct['mature'].append(a)
        dct['mature_annot'].append( b  + c )
    return dct
## calulate the chemechical propierties of a list of sequences given the init Peptides library
def calc(seqs,pep):
    df = {}
    df['mw'] = []
    df['tiny'] = []
    df['small'] = []
    df['aro'] = []
    df['nonP'] = []
    df['chrgd'] = []
    df['basic'] = []
    df['acidic'] = []
    df['charge'] = []
    df['pI'] = []
    df['aindex'] = []
    df['instaindex'] = []
    df['boman'] = []
    df['hmoment100'] = []
    df['hmoment160'] = []
    for i in seqs:
        df['mw'].append(pep.mw(i)[0])
        df['tiny'].append(pep.aaComp(i)[0][-9])
        df['small'].append(pep.aaComp(i)[0][-8])
        df['aro'].append(pep.aaComp(i)[0][-6])
        df['nonP'].append(pep.aaComp(i)[0][-5])
        df['chrgd'].append(pep.aaComp(i)[0][-3])
        df['basic'].append(pep.aaComp(i)[0][-2])
        df['acidic'].append(pep.aaComp(i)[0][-1])
        df['charge'].append(pep.charge(i,pH=0,pKscale = "EMBOSS")[0])
        df['pI'].append(pep.pI(i,pKscale = "EMBOSS")[0])
        df['aindex'].append(pep.aIndex(i)[0])
        df['instaindex'].append(pep.instaIndex(i)[0])
        df['boman'].append(pep.boman(i)[0])
        df['hmoment100'].append(pep.hmoment(i, angle = 100, window = 11)[0])
        df['hmoment160'].append(pep.hmoment(i, angle = 160, window = 11)[0])
    return df

