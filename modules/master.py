'''
Software for AMP identification given a amino acids sequence
single-line multi-FASTA


        USAGE: 
            
            NRamps.py <FASTA> [OPTION]

        OPTIONS:
            -m predict mature sequences and classify with mature sequence

        DEPENDENCIES:
            sklearn library (Python)
            pandas  library (Python)
            importr library (Python)
            R       v > 3.0
            Peptides package (R)
      
        OTHERS
            All the files in the lib/ directory (including NRamps.py) must be in the same directory when executed.
            temporary files are written in tmpAmps directory, then deleted

'''
import warnings
def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
import sys,os,getopt,re
import naiveBayesFull as nb
#import chemClassificator as ch    
import chemClassificator as ch
from rpy2.robjects.packages import importr

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
        for line in range(len(lines)):
            i = lines[line].split()[0]
            if re.search("^>",i):
                if ident and not len(seq) == 0:
                    if len(seq) < 8:
                        None
                    else:
                        df['ids'].append(ident)
                        df['seqs'].append(seq)
                    
                ident = i.split()[0].split(">")[1]
                seq = ""
            if not re.search("^>",i) and not re.search("[^ACEDGFIHKMLNQPSRTWVY]",i): 
                
                seq = seq + i.strip()
            if line == len(lines) -1:
                df['ids'].append(ident)
                df['seqs'].append(seq)
    if len(df['seqs'][-1]) == 0:
        df['seqs'].pop()
        df['ids'].pop()
    return df['ids'],df['seqs']

cl_annot = {'0':'Negative',
            '1':'1, cluster rep. protein AC P0CAY5', 
            '2':'2, cluster rep. protein AC P85245', 
            '3':'3, cluster rep. protein AC Q2V318', 
            '4':'4, cluster rep. protein AC P83880', 
            '5':'5, cluster rep. protein AC P19656', 
            '6':'6, cluster rep. protein AC Q56XC2', 
            '7':'7, cluster rep. protein AC Q5USN8', 
            '8':'8, cluster rep. protein AC P82631', 
            '9':'9, cluster rep. protein AC P80214', 
            '10':'10, cluster rep. protein AC P69034', 
            '11':'11, cluster rep. protein AC Q718F4'}
### codee
libpath = os.path.dirname(os.path.realpath(sys.argv[0]))
try:
    fasta = sys.argv[1]
except:
    print( 'NRamps.py <fasta>')
    sys.exit(2)
try:
    opt = sys.argv[2]
except:
	opt = 0

ids,seqs = fasta_to_df(fasta)
#mature seqs"
chem = ch.ChemClassificator()

if opt == "-m":
    a,b = chem.predict_mature_seqs(seqs, libpath + "/../models/cleaveN_rf_model.pkl" , libpath + "/../models/cleaveC_rf_model.pkl") 

li = []
li2 = []
annot = []
pep = importr('Peptides')
modelChem = libpath + "/../models/chem_11_rf_model_UNSCALED.pkl"

if opt == "-m":
	cla,d = chem.predict_amps(a,modelChem,pep)
else:
    cla,d = chem.predict_amps(seqs,modelChem,pep)
for c in range(1,12):
    nv = nb.NaiveAmps()         
    if opt == "-m":                                  
        nv.classSeq(libpath + "/../models/new_NB_model1_" + str(c) + ".json", [x for x in a] ,0.5)            
    else:         
        nv.classSeq(libpath + "/../models/new_NB_model1_" + str(c) + ".json", [x for x in seqs] ,0.5)            
    li.append(nv.probs)
    li2.append(nv.annot)
clus = []
maxprobs = []              
for i in range(len(li[0])):
    clus.append([x[i] for x in li].index(max([x[i] for x in li])))  
    maxprobs.append(max([x[i] for x in li]))
for i in range(len(seqs)) :
    if maxprobs[i] > 0.5:
        annot.append("NaiveBayes class: " + cl_annot[str(clus[i])] +  ", prob:" + str(round(maxprobs[i],2)) ) 
    else:
        annot.append("NaiveBayes class: " + cl_annot['0'] + ", prob: " + str(round(maxprobs[i],2)))
for i in range(len(seqs)):
    if opt == "-m":
        print(">" + ids[i] + " // " + str(b[i]) + " // " + "RF class: " + cl_annot[str(cla[i])] + ", prob:" + str(d[i]) + " // " + annot[i]  + "\n" + a[i] )
    else:
        print(">" + ids[i] + " // " + "RF class: " + cl_annot[str(cla[i])] + ", prob:" + str(d[i]) + " // " + annot[i]  + "\n" + seqs[i] )

sys.exit()
