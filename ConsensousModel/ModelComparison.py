import cobra.io
from numpy import *
from os.path import join
from scipy.stats import *
from bisect import bisect_left
import cobra.manipulation
from copy import deepcopy
import scipy.stats
import DataCollection.EssentialGeneFinder
from pandas import DataFrame
import pandas
###DELRXN####
file = open('resistantReactionsKO.txt','r')
r = file.readlines()
r = [x[:-1] for x in r]
file.close()

file =open('ControlReactionsKO.txt','r')
c = file.readlines()
c = [x[:-1] for x in c]
file.close()

cu = list(set(c) - set(r))
ru = list(set(r) - set(c))

file = open('uniqueDelRxnC.txt','w')
file.write('Control\n')
[file.write(x+'\n') for x in cu]
file.close()
file = open('uniqueDelRxnR.txt','w')
file.write('Resistant\n')
[file.write(x+'\n') for x in ru]
file.close()


###GeneDeletionSimulation###
r = DataFrame.from_csv('resistantGeneDeletion.csv')
c = DataFrame.from_csv('ControlGeneDeletion.csv')
rG = set()
[rG.add((x,round(y,1))) for x,y in zip(r.index.values,r.flux.values)]
cG = set()
[cG.add((x,round(y,1))) for x,y in zip(c.index.values,c.flux.values)]
cu = cG - rG
ru = rG - cG
file = open('GDCompC.txt','w')
file.write('GeneName Flux\n')
[file.write(x[0] + ' ' + str(x[1])+'\n')  for x in cu]
file.close()
file = open('GDCompR.txt','w')
file.write('GeneName Flux\n')
[file.write(x[0] + ' ' + str(x[1])+'\n')  for x in ru]
file.close()

### RXN deletion ###
r = DataFrame.from_csv('resistantRxnDeletion.csv')
c = DataFrame.from_csv('ControlRxnDeletion.csv')
rR = set()
[rR.add((x,round(y,1))) for x,y in zip(r.index.values,r.flux.values)]
cR = set()
[cR.add((x,round(y,1))) for x,y in zip(c.index.values,c.flux.values)]
cu = cR - rR
ru = rR - cR

file = open('RDCompC.txt','w')
file.write('RxnName Flux\n')
[file.write(x[0] + ' ' + str(x[1])+'\n')  for x in cu]
file.close()
file = open('RDCompR.txt','w')
file.write('RxnName Flux\n')
[file.write(x[0] + ' ' + str(x[1])+'\n')  for x in ru]
file.close()

###FVA####
r = DataFrame.from_csv('resistantFVA.csv')
c = DataFrame.from_csv('ControlFVA.csv')
rR = set()
[rR.add((x,round(s,1),round(l,1))) for x,s,l in zip(r.index.values,r.minimum.values,r.maximum.values)]
cR = set()
[cR.add((x,round(s,1),round(l,1))) for x,s,l in zip(c.index.values,c.minimum.values,c.maximum.values)]
cu = cR - rR
ru = rR - cR
file = open('FVACompC.txt','w')
file.write('rxnName CMin CMax\n')
[file.write(x[0] + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n') for x in cu]
file.close()
file = open('FVACompR.txt','w')
file.write('rxnName RMin RMax\n')
[file.write(x[0] + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n') for x in ru]
file.close()

