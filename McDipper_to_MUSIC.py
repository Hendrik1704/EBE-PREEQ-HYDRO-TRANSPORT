import numpy as np
import sys
import re
from scipy.interpolate import interp1d
M_hbarc = 0.197326979

def read_energy_data(filename, Ns, NsLong):
    Energy = np.zeros((NsLong, Ns, Ns), dtype=np.float64)
    with open(filename, 'r') as f:
        # Extract tau
        comment_line = f.readline()
        match = re.search(r'#\s*tau0\s*=\s*([0-9.+-eE]+)', comment_line)
        if match:
            tau = float(match.group(1))
            
            # read file
            for line in f:
                ix, iy, ieta, energy = map(float, line.strip().split()[:4])
                ix, iy, ieta = int(ix), int(iy), int(ieta)
                Energy[ieta, ix, iy] = energy
        else:
            raise ValueError("Tau value not found in the comment line")
    return tau, Energy

def get_BulkPressure(p_interp, energy):
  return energy/3 - p_interp(energy)

def main():
    # Check if all parameters are provided as command line arguments
    if len(sys.argv) != 8:
        print("Usage: python script.py <eos_file> <Ns> <Nslong> <dx> <deta> <filename_in> <filename_out>")
        sys.exit(1)

    EOS_path = sys.argv[1]
    Ns = int(sys.argv[2])
    NsLong = int(sys.argv[3])
    dx = float(sys.argv[4])
    dy = dx
    deta = float(sys.argv[5])
    input_file = sys.argv[6]
    output_file = sys.argv[7]
    
    # EOS
    eos=np.fromfile(EOS_path, dtype=float).reshape(-1, 4)
    # 4 columns (energy density, local pressure, entropy density, local temperature)
    e=eos[:,0] # (GeV/fm^3)
    p=eos[:,1] # (GeV/fm^3)
    # pressure as a function of energy [p(e)]
    p_eos = interp1d(e, p, kind='linear', fill_value='extrapolate')

    energy_cut = 1e-15
    
    tau, Energy = read_energy_data(input_file, Ns, NsLong)

    comment_line = f"# tau_in_fm {tau} etamax= {NsLong} xmax= {Ns} ymax= {Ns} deta= {deta} dx= {dx} dy= {dy} \n"
    flow = [1.0, 0.0, 0.0, 0.0]

    with open(output_file, 'w') as file:
        file.write(comment_line)
        for ieta in range(NsLong):
            eta_pos = deta*(ieta-(NsLong-1)/2.0) 
            for ix in range(Ns):
                xpos = dx*(ix-(Ns-1)/2.0) ## fm
                for iy in range(Ns):
                    ypos = dy*(iy-(Ns-1)/2.0) ## fm
                    energy = Energy[ieta, ix, iy] ## GeV/fm^3
                    
                    if energy > energy_cut:
                        pitautau = pitaux = pitauy = pitaueta = pixy = pixeta = piyeta = 0
                        pixx = piyy = energy/(6*M_hbarc) # fm^-4
                        pietaeta = -energy/(3*tau*tau*M_hbarc) ## fm^-6

                        bulk = get_BulkPressure(p_eos, energy)  ## GeV/fm^3
                    else: 
                        pitautau = pitaux = pitauy = pitaueta = pixy = pixeta = piyeta = 0
                        pixx = piyy = pietaeta = bulk = 0.0

                    file.write(f"{eta_pos:.2f} {xpos:.2f} {ypos:.2f} {energy} {flow[0]} {flow[1]} {flow[2]} {flow[3]} {pitautau} {pitaux} {pitauy} {pitaueta} {pixx} {pixy} {pixeta} {piyy} {piyeta} {pietaeta} {bulk} \n")

if __name__ == "__main__":
    main()
