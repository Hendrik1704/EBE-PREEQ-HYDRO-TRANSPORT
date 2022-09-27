import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import os

PDG_SMASH_PATH = "./iSS/iSS_tables/pdg-SMASH.dat"

PDG_SMASH = pd.read_csv(PDG_SMASH_PATH,header=None, sep='\s\s+|,', engine='python')
PDG_SMASH.columns = ['PDG_ID','name','mass','width','gspin','baryon','strange','charm','bottom','gisospin','charge','decays']
PDG_SMASH = PDG_SMASH.dropna()

# generate arrays with the PDG_ID and the charge of the particles
PDG_ID = PDG_SMASH['PDG_ID'].to_numpy()
CHARGES = PDG_SMASH['charge'].to_numpy()

def get_PDG_ID_charge(id):
    charge = 0
    notinlist = 0
    for n in range(len(PDG_ID)):
        if PDG_ID[n] == id:
            charge = CHARGES[n]
            break
        if n == range(len(PDG_ID))[-1]:
            notinlist = 1
    if notinlist == 1:
        for n in range(len(PDG_ID)):
            if PDG_ID[n] == -id:
                charge = -CHARGES[n]
                break
    return int(charge)

# load the Oscar1999A file
parser = argparse.ArgumentParser()
parser.add_argument('inputFilePath', type=Path)
parser.add_argument('outputFilePath', type=Path)
p = parser.parse_args()

col_names = ['sample_idx','PDG_ID','px','py','pz','E','m','x','y','z','t']
OSCAR1999A_FILE = pd.read_csv(p.inputFilePath,names=col_names, sep='\s\s+|,', engine='python', skiprows=3)

is_NaN = OSCAR1999A_FILE.isnull()
row_has_NaN = is_NaN.any(axis=1)
rows_with_NaN = OSCAR1999A_FILE[row_has_NaN]

# count the number of events sampled from the hypersurface
NumberSamples = 0
event_idx = []
for i in range(len(OSCAR1999A_FILE)):
    if row_has_NaN.loc[i] == True:
        NumberSamples += 1
        event_idx.append(i)
event_idx.append(len(OSCAR1999A_FILE))

# generate new OSCAR2013 files
f = open(os.path.join(p.outputFilePath,'OSCAR%s'%(int(''.join(filter(str.isdigit, str(p.inputFilePath)))))), 'w')
f.write('#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n')
f.write('# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n')
for ev in range(NumberSamples):
    f.write('# event %d out %d\n'%(ev+1,OSCAR1999A_FILE.loc[event_idx[ev]].at['PDG_ID']))
    for line in range(event_idx[ev]+1,event_idx[ev+1]):
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['t']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['x']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['y']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['z']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['m']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['E']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['px']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['py']))
        f.write('%s '%(OSCAR1999A_FILE.loc[line].at['pz']))
        id = OSCAR1999A_FILE.loc[line].at['PDG_ID']
        f.write('%d '%(id))
        f.write('%d '%(OSCAR1999A_FILE.loc[line].at['sample_idx']))
        charge = get_PDG_ID_charge(int(id))
        f.write('%d'%(charge))
        f.write('\n')
    f.write('# event %d end 0 impact 0.000\n'%(ev+1))
f.close()
