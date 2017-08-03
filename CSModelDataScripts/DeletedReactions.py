from ABResistanceAnalysis import *
from pandas import DataFrame
import cobra
import cobra.manipulation

"""Script for finding list of deleted reactions from deleted genes given in textfile"""

data = creationOfDeletedReaction()
print data
file = open('removedReactionControl.txt','w')
file.write('Control\n')
[file.write(x+'\n') for x in data[0]]
file.close()
file = open('removedReactionExp.txt','w')
file.write('Experimental\n')
[file.write(x+'\n') for x in data[1]]
file.close()
