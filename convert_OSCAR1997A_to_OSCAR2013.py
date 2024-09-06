import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import os

PDG_SMASH_PATH = "./iSS_tables/pdg-SMASH.dat"

PDG_SMASH = pd.read_csv(PDG_SMASH_PATH,header=None, sep='\s\s+|,', engine='python')
PDG_SMASH.columns = ['PDG_ID','name','mass','width','gspin','baryon','strange','charm','bottom','gisospin','charge','decays']
PDG_SMASH = PDG_SMASH.dropna()

# generate arrays with the PDG_ID and the charge of the particles
PDG_ID = PDG_SMASH['PDG_ID'].to_numpy()
CHARGES = PDG_SMASH['charge'].to_numpy()

# Cache for the charge of the particles which are already calculated
charge_cache = {}
def get_PDG_ID_charge(id):
    if id in charge_cache:
        return charge_cache[id]
    charge = 0
    notinlist = False
    for n in range(len(PDG_ID)):
        if PDG_ID[n] == id:
            charge = CHARGES[n]
            break
        if n == range(len(PDG_ID))[-1]:
            notinlist = True
    if notinlist:
        for n in range(len(PDG_ID)):
            if PDG_ID[n] == -id:
                charge = -CHARGES[n]
                break
    return int(charge)

def process_oscar1997_file(input_file, output_file, buffer_size=50000):
    num_skiplines = 3
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write('#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge\n')
        outfile.write('# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e\n')

        buffer = []  # Buffer to store lines
        event_count = 0
        for line in infile:
            if num_skiplines > 0:
                num_skiplines -= 1
                continue
            # If line has only four elements, it is a header
            if len(line.split()) == 4:
                # If event_count is >= 1 add a footer
                if event_count >= 1:
                    buffer.append(f'# event {event_count} end 0 impact 0.000\n')

                # The second number in the header is the number of hadrons in the event
                num_hadrons = int(line.split()[1])
                event_count += 1
                event_header = f'# event {event_count} out {num_hadrons}\n'
                buffer.append(event_header)

            parts = line.strip().split()
            if len(parts) == 11:
                sample_idx = int(parts[0])
                pdg_id = int(parts[1])
                px = float(parts[2])
                py = float(parts[3])
                pz = float(parts[4])
                E = float(parts[5])
                m = float(parts[6])
                x = float(parts[7])
                y = float(parts[8])
                z = float(parts[9])
                t = float(parts[10])
                charge = get_PDG_ID_charge(pdg_id)
                processed_line = f"{t:g} {x:g} {y:g} {z:g} {m:g} {E:g} {px:g} {py:g} {pz:g} {pdg_id} {sample_idx} {charge}\n"
                buffer.append(processed_line)

            # Write buffer to file once it reaches a certain size
            if len(buffer) >= buffer_size:
                outfile.writelines(buffer)
                buffer = []

        # Add footer for the last event
        buffer.append(f'# event {event_count} end 0 impact 0.000\n')
        # Write any remaining lines in the buffer to the file
        if len(buffer) > 0:
            outfile.writelines(buffer)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFilePath', type=Path)
    parser.add_argument('outputFilePath', type=Path)
    args = parser.parse_args()

    # Process the input file and write to the output file
    process_oscar1997_file(args.inputFilePath, args.outputFilePath)
