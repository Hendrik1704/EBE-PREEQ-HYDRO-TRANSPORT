import sys
import numpy as np
from scipy.interpolate import interp1d

def read_file(filename):
    data = []
    comment_line = ""
    
    try:
        with open(filename, 'r') as file:
            comment_line = next(file).strip()
            for line in file:
                # Split each line into a list of floating-point numbers
                row = [float(value) for value in line.split()]
                if len(row) > 0:
                    data.append(row)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e} in file '{filename}'. Make sure all values are numbers.")
        sys.exit(1)

    return np.array(data), comment_line

def new_e_and_bulk(e_kompost):
    e_k = np.array(e_kompost) # GeV/fm^3
    # s_kompost: sT = e + p_kompost = 4/3 e
    s_kompost = (4/3.)*e_k 
    
    # EOS
    eos_location = "./MUSIC/EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat"
    raw=np.fromfile(eos_location, dtype=(float,4))
    # 4 columns (energy density, local pressure, entropy density, local temperature)
    e=raw[:,0] # (GeV/fm^3)
    p=raw[:,1] # (GeV/fm^3)

    # entropy: sT = e+p
    entropy = [p+e for p, e in zip(p, e)]
    # energy as a function of entropy [e(sT)]
    energy_interp = interp1d(entropy, e, kind='linear', fill_value='extrapolate')
    # pressure as a function of energy [p(e)]
    p_interp = interp1d(e, p, kind='linear', fill_value='extrapolate')

    e_Music = energy_interp(s_kompost)
    
    ## same thing than do: s(e)*T(e)
    #sT = s*T
    #e_sT = interp1d(sT, e, kind='linear', fill_value='extrapolate')
    #e_Music = e_sT(s_kompost)
    
    Bulk_pressure = e_k/3 - p_interp(e_Music)
    
    return e_Music, Bulk_pressure # energy in GeV/fm^3
 
def main():
    # Check if all parameters are provided as command line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <filename_in> <filename_out>")
        sys.exit(1)

    # parse input parameters
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]

    # read the input file
    data, comment = read_file(filename_in)
    
    # define the energy and calculate the new energy
    energy_density = data[:, 3]
    energy_Music, bulk_pressure = new_e_and_bulk(energy_density)
    
    # add the bulk pressure as a column in data
    data = np.column_stack((data, bulk_pressure))

    # Replace the column 4: change the energy kompost to the new energy
    data[:, 3] = energy_Music

    # Write the final file: eta, x, y, new_energy, u^tau, u^x, u^y, u^eta,
    # pi^tautau, pi^taux, pi^tauy, pi^taueta, pi^xx, pi^xy, pi^xeta, pi^yy, pi^yeta, pi^etaeta
    with open(filename_out, 'w') as file:
        file.write(comment + '\n')
        for row in data:
            file.write(' '.join(map(str, row)) + '\n')

if __name__ == "__main__":
    main()
