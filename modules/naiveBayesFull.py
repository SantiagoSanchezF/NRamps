import os, json, re, math
'''
NaiveTrain is a software to create a model of Nave bayes classifier, given the df output (freq matrix) of patternFarming.py
USAGE:
    1. load this functions,
    2.load df.PattF.json as ptdf
    3. b = NaiveAmps()
    4.b.freqDict_to_model(ptdf.keys(),ptdf)
    5.b.classSeq('PATHtoFASTA.fa',CUT)
    6.pred sequences are in b.dataF['pred']. Alternativly you can save the file in a json
        b.save_Df(NameOfFile)
'''
class NaiveAmps():
    def __init__(self):
        self.pred = []
        self.model = {}
    def freqDict_to_model(self, subsetNames, freqDict):
        self.model = makeProbsDfModel(subsetNames, freqDict)
    def load_prev_model(self, filename):
        with open(filename,"r") as myfile:
            self.model = json.load(myfile)    
    def classSeq(self,model, seqs, cut):
        with open(model,"r") as myfile:
            self.model = json.load(myfile)
        self.pred = [0 if u < cut else 1 for u in [ seq_probs(x, self.model[1]) for x in seqs ] ]
        self.probs = [ seq_probs(x, self.model[1]) for x in seqs ] 
        self.annot = [self.model[0] for x in seqs]
    def save_Df(self,name):
        with open(name,"w") as myfile:
            json.dump(self.dataF,myfile)
        
# claculates probabiiti of being AMP
def seq_probs(sequence ,model):
    probAmpsLnSum = 0
    probNotAmpsLnSum = 0
    for i in range(len(model['pattRegex'])):
        if re.search(model['pattRegex'][i],sequence):
            probAmpsLnSum = probAmpsLnSum + math.log( model['probAmp'][i])
            probNotAmpsLnSum = probNotAmpsLnSum + math.log( model['probNotAmp'][i])
        else:
            probAmpsLnSum = probAmpsLnSum + math.log(1 - model['probAmp'][i])
            probNotAmpsLnSum = probNotAmpsLnSum + math.log(1 - model['probNotAmp'][i])
    probAmps = math.exp(probAmpsLnSum)
    probNotAmps = math.exp(probNotAmpsLnSum)
    probAmp = probAmps/(probAmps + probNotAmps)
    return probAmp

# add sequence and ids to df
def fasta_to_df(fasta):
    dataF = {}
    dataF['id'] = []
    dataF['seq'] = []
    dataF['pred'] = []
    with open(fasta,"r") as myfile:
        for i in myfile:
            if re.search("^>", i):
                dataF['id'].append(i.split()[0].split(">")[-1].split("\n")[0])
            elif re.search("^[A-Z]",i):
                dataF['seq'].append(i.split("\n")[0])
    if len(dataF['id']) != len(dataF['seq']):
        print("your fasta is missing a sequence or an id, migth give errors")
    return dataF

# gets size of each class
def get_class_freq(matrix):
    posN = sum(matrix['clase'])
    negN = len(matrix['clase']) - posN
    return posN,negN

# calculates binary event probability given a list
def pattProb(patList,classList,totAmps,totNotAmps, k = 0.5):
    probAmp = (sum([1 for x in range(len(classList)) if classList[x] == 1 and patList[x] == 1]) + k) / (totAmps + 2*k)
    probNotAmp = ( sum([patList[x] for x in range(len(classList)) if patList[x]== 1 and classList[x] == 0]) + k) / (totNotAmps + 2*k)
    return probAmp, probNotAmp

# creates a df with pattern and prob of Amp and Not Amp
def makeProbsDfModel(pattNames,freqDict):
    probDf = {}
    probDf['patt'] = [ x for x in pattNames if x not in ['seqs','clase']]
    probDf['probAmp'] = []
    probDf['probNotAmp'] = []
    probDf['pattRegex'] = []
    u,v = get_class_freq(freqDict)
    for i in probDf['patt']:
        a,b = pattProb(freqDict[i],freqDict['clase'],u,v)
        probDf['pattRegex'].append(patt_to_regex(i))
        probDf['probAmp'].append(a)
        probDf['probNotAmp'].append(b)
    return probDf

###transform pratt patts to python regex
def patt_to_regex(patt):
     single_exp = []
     splited_patt = patt.split("-")
     for i in splited_patt:
         if len(i) == 1:
             if i == 'x':
                 single_exp.append('[A-Z]')
             else:
                 single_exp.append(i)
             continue
         if re.search('^\[',i):
             single_exp.append(i)
             continue
         x = re.split(r'[(),]', i)
         if len(x) == 4 and x[0] == 'x':
             single_exp.append('[A-Z]' + '{' + x[1] + ',' + x[2] + '}')
         elif len(x) == 3 and x[0] == 'x':
             single_exp.append('[A-Z]' + '{' + x[1] + '}')
         elif len(x) == 4 and x[0] != 'x':
             single_exp.append(x[0] + '{' + x[1] + ',' + x[2] + '}')
         elif len(x) == 3 and x[0] != 'x':
             single_exp.append(x[0] + '{' + x[1] + '}')
     return ('').join(single_exp)
