import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

HBARC = 0.197326979

def read_file(filename):
    data = []
    try:
        with open(filename, 'r') as file:
            for line in file:
                # Split each line into a list of floating-point numbers
                row = [float(value)*HBARC**3. for value in line.split()]
                if len(row) > 0:
                    data.append(row)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e} in file '{filename}'. Make sure all values are numbers.")
        sys.exit(1)

    return data

def get_Tmunu_Milne(data_read, tau_hydro):

    tau_inv_GeV = HBARC / tau_hydro

    Tmunu_matrix_Milne = []
    for point in range(len(data_read)):
        Tmunu = data_read[point][2:12]
    
        matrix = np.zeros((4,4))
        #KoMPoST version
        matrix[0][0] = Tmunu[0]
        matrix[1][1] = Tmunu[1]
        matrix[2][2] = Tmunu[2]
        matrix[3][3] = tau_inv_GeV**2.*Tmunu[3]
        matrix[0][1] = -Tmunu[4]
        matrix[0][2] = -Tmunu[5]
        matrix[0][3] = -tau_inv_GeV*Tmunu[6]
        matrix[1][2] = -Tmunu[7]
        matrix[2][3] = -tau_inv_GeV*Tmunu[8]
        matrix[3][1] = -tau_inv_GeV*Tmunu[9]

        matrix[1][0] = matrix[0][1]
        matrix[1][3] = matrix[3][1]
        matrix[2][0] = matrix[0][2]
        matrix[2][1] = matrix[1][2]
        matrix[3][0] = matrix[0][3]
        matrix[3][2] = matrix[2][3]

        Tmunu_matrix_Milne.append(matrix)

    return Tmunu_matrix_Milne

def get_energy_density_and_flow(Tmunu,tau_hydro):
    # T^mu_nu u^nu = energy_density u^mu
    # This function solves for the eigenvalues (energy) and the eigenvectors (flow)
    # of the T^{\mu\nu}

    tau_in_GeVinv = tau_hydro / HBARC

    g = np.zeros((4,4))
    g[0][0] = 1.0
    g[1][1] = -1.0
    g[2][2] = -1.0
    g[3][3] = -tau_in_GeVinv**2.

    imaginary_limit = 1e-30

    energy_density = []
    flow = []
    for point in range(len(Tmunu)):
        # T^mu_nu
        Tmunu_tmp = np.matmul(Tmunu[point],g)

        eigenvalues, eigenvectors = np.linalg.eig(Tmunu_tmp)

        found_a_candidate = False
        found_too_many_candidates = False
        eigenvalue_candidate = -1.0
        eigenvector_candidate = np.array([0.,0.,0.,0.])
        eigenvector_candidate_norm = 1.0
        # loop through the eigenvectors and eigenvalues
        for i in range(4):
            # get the eigenvalue
            eval_tmp = eigenvalues[i]
            eval_re = np.real(eval_tmp)
            eval_im = np.imag(eval_tmp)

            # continue if eigenvalue is positive and real
            if (eval_re > 0.) and (np.abs(eval_im/eval_re) < imaginary_limit):
                all_evec_elements_real = True
                norm = 0.0
                eigenvector = np.zeros(4)
                # loop over corresponding eigenvector elements 
                for j in range(4):
                    evec_elem = eigenvectors[j][i]
                    evec_elem_re = np.real(evec_elem)
                    evec_elem_im = np.imag(evec_elem)

                    # check if the element of the eigenvector is real
                    if (0.0 == evec_elem_im) or (np.abs(evec_elem_im/evec_elem_re) < imaginary_limit):
                        # g_mu_nu weight
                        norm += g[j][j]*evec_elem_re**2.
                        eigenvector[j] = evec_elem_re
                    else:
                        all_evec_elements_real = False
                        break
                
                # good if vector is timelike and components are real
                if (norm > 0.0) and all_evec_elements_real:
                    # if another real timelike candidate was already found we have a problem
                    if found_a_candidate:
                        similar_enough = 1e-3
                        if (np.abs(1.-eigenvalue_candidate/eval_re)>similar_enough) and\
                            (eigenvalue_candidate != eval_re) and\
                            (np.abs(1.-eigenvector_candidate[0]/eigenvector[0])>similar_enough) and\
                            (eigenvalue_candidate[0] != eigenvector[0]) and\
                            (np.abs(1.-eigenvector_candidate[1]/eigenvector[1])>similar_enough) and\
                            (eigenvalue_candidate[1] != eigenvector[1]) and\
                            (np.abs(1.-eigenvector_candidate[2]/eigenvector[2])>similar_enough) and\
                            (eigenvalue_candidate[2] != eigenvector[2]) and\
                            (np.abs(1.-eigenvector_candidate[3]/eigenvector[3])>similar_enough) and\
                            (eigenvalue_candidate[3] != eigenvector[3]):
                            found_too_many_candidates = True

                    eigenvalue_candidate = eval_re
                    eigenvector_candidate_norm = norm
                    eigenvector_candidate = eigenvector
                    found_a_candidate = True

        tmp_energy = eigenvalue_candidate

        time_component_sign = np.sign(eigenvector_candidate[0])
        tmp_flow = eigenvector_candidate * (time_component_sign/np.sqrt(eigenvector_candidate_norm))
        successful = not found_too_many_candidates and found_a_candidate

        if successful:
            energy_density.append(tmp_energy)
            flow.append(tmp_flow)
        else:
            raise ValueError('The eigensystem could not be parsed by the KoMPoST_to_MUSIC converter.')

    return energy_density, flow

def get_shear_and_bulk_conformal_EOS(Tmunu, tau_hydro, energy_density, flow, method=1):
    # return the bulk pressure and the shear tensor for all points on the grid
    # two different methods can be used, the default is the one used in KoMPoST
    tau_in_GeVinv = tau_hydro / HBARC
    
    g_inv = np.zeros((4,4))
    g_inv[0][0] = 1.0
    g_inv[1][1] = -1.0
    g_inv[2][2] = -1.0
    g_inv[3][3] = -1.0/tau_in_GeVinv**2.

    g = np.zeros((4,4))
    g[0][0] = 1.0
    g[1][1] = -1.0
    g[2][2] = -1.0
    g[3][3] = -tau_in_GeVinv**2.

    bulk_pressure = []
    pi_munu = []
    for point in range(len(Tmunu)):
        Ttautau = Tmunu[point][0][0] # GeV^4
        Txx = Tmunu[point][1][1] # GeV^4
        Tyy = Tmunu[point][2][2] # GeV^4
        Tetaeta = Tmunu[point][3][3] # GeV^6

        energy_density_point = energy_density[point] # GeV^4
        flow_point = flow[point]

        # conformal pressure
        pressure = energy_density_point / 3.0 # GeV^4
        # energy momentum tensor trace
        Tmunu_trace = Ttautau - Txx - Tyy - tau_in_GeVinv**2. * Tetaeta # GeV^4
        # bulk pressure
        p_bulk = -Tmunu_trace / 3.0 # GeV^4
        bulk_pressure.append(p_bulk)
        
        # compute u_mu
        flow_point_lower = [g[0][0]*flow_point[0],g[1][1]*flow_point[1],g[2][2]*flow_point[2],g[3][3]*flow_point[3]]
        # last component GeV^-2

        pi_munu_upper = np.zeros((4,4))
        if method == 1:
            # METHOD 1:
            # T^{\mu\nu} = -e*u^\mu*u^\nu + T^{\mu\alpha}*u_\alpha*u^\nu 
            #               + T^{\nu\alpha}*u_\alpha*u^\mu + 1/3*(g^{\mu\nu}-u^\mu*u^\nu)
            #               * (T^\alpha_\alpha - e) + \pi^{\mu\nu}
            for mu in range(4):
                for nu in range(4):
                    #\pi^{\mu\nu} = T^{\mu\nu} + e*u^\mu*u^\nu - T^{\mu\alpha}*u_\alpha*u^\nu 
                    #               - T^{\nu\alpha}*u_\alpha*u^\mu - 1/3*(g^{\mu\nu}-u^\mu*u^\nu)
                    #               * (T^\alpha_\alpha - e)
                    pi_munu_upper[mu][nu] = Tmunu[point][mu][nu]\
                        + energy_density_point*flow_point[mu]*flow_point[nu]\
                        - 1./3. * (g_inv[mu][nu] - flow_point[mu]*flow_point[nu]) * (Tmunu_trace - energy_density_point)
                    
                    for alpha in range(4):
                        pi_munu_upper[mu][nu] -= Tmunu[point][mu][alpha]*flow_point[nu]*flow_point_lower[alpha]
                        pi_munu_upper[mu][nu] -= Tmunu[point][nu][alpha]*flow_point[mu]*flow_point_lower[alpha]
        elif method == 2:
            # METHOD 2:
            # T^{\mu\nu} = e*u^\mu*u^\nu - (p + PI)*(g^{\mu\nu} - u^\mu*u^\nu) + \pi^{\mu\nu}
            for mu in range(4):
                for nu in range(4):
                    #\pi^{\mu\nu} = T^{\mu\nu} - e*u^\mu*u^\nu + (p + PI)*(g^{\mu\nu} - u^\mu*u^\nu)
                    pi_munu_upper[mu][nu] = Tmunu[point][mu][nu]\
                        - energy_density_point*flow_point[mu]*flow_point[nu]\
                        + (pressure + p_bulk) * (g_inv[mu][nu] - flow_point[mu]*flow_point[nu])
        else:
            raise ValueError("Method in 'get_shear_and_bulk_conformal_EOS' does not exist")

        pi_munu.append(pi_munu_upper)

    # pi^{t/x/y t/x/y} [GeV^4], pi^{t/x/y eta} [GeV^5], pi^{eta eta} [GeV^6]
    return bulk_pressure, pi_munu

def interpolate_e(e_kompost):
    e_k = np.array(e_kompost) # GeV^4
    # s_kompost: sT = e + p_kompost = 4/3 e
    s_kompost = (4./3.)*(e_k/HBARC**3.) # Needs to be in GeV/fm^3
    
    # EOS
    eos_location = "./MUSIC/EOS/hotQCD/hrg_hotqcd_eos_SMASH_binary.dat"
    eos = np.fromfile(eos_location, dtype=(float,4))
    # 4 columns (energy density, local pressure, entropy density, local temperature)
    e = eos[:,0] # (GeV/fm^3)
    p = eos[:,1] # (GeV/fm^3)

    # entropy: sT = e+p
    entropy = [p+e for p, e in zip(p, e)]
    # energy as a function of entropy [e(sT)]
    energy_interp = interp1d(entropy, e, kind='linear', fill_value='extrapolate')

    e_Music = energy_interp(s_kompost)
    
    ## same thing than do: s(e)*T(e)
    #sT = s*T
    #e_sT = interp1d(sT, e, kind='linear', fill_value='extrapolate')
    #e_Music = e_sT(s_kompost)
    
    return e_Music # energy in GeV/fm^3

def write_MUSIC_output(filename_out, energy_density, flow, shear_tensor, tau_hydro, Nx, Ny, dx, dy):
    # Generate x_pos, y_pos, and pos_tuple
    x_pos = []
    y_pos = []
    for ix in range(Nx):
        for iy in range(Ny):
            x_pos.append(dx*(ix - Nx/2.))
            y_pos.append(dy*(iy - Ny/2.))
    
    try:
        with open(filename_out, 'w') as file:
            file.write(f"# tau_in_fm {tau_hydro} etamax= 1 xmax= {Nx} ymax= {Ny} deta= 0 dx= {dx} dy= {dy}\n")
            for ix in range(Nx):
                for iy in range(Ny):
                    q = ix * Ny + iy # index for loop order in MUSIC
                    p = iy * Nx + ix # correct index for MUSIC loop order in KoMPoST list
                    pitautau = shear_tensor[p][0][0] / HBARC**4.
                    pitaux = shear_tensor[p][0][1] / HBARC**4.
                    pitauy = shear_tensor[p][0][2] / HBARC**4.
                    pitaueta = shear_tensor[p][0][3] / HBARC**5.
                    pixx = shear_tensor[p][1][1] / HBARC**4.
                    pixy = shear_tensor[p][1][2] / HBARC**4.
                    pixeta = shear_tensor[p][1][3] / HBARC**5.
                    piyy = shear_tensor[p][2][2] / HBARC**4.
                    piyeta = shear_tensor[p][2][3] / HBARC**5.
                    pietaeta = shear_tensor[p][3][3] / HBARC**6.

                    #Final units: x, y [fm], energy [GeV/fm^3], only u^tau [fm^-1], pi^{t/x/y} [fm^-4] pi^{eta t/x/y} [fm^-5] pi^{etaeta} [fm^-6]
                    row = [0, np.round(x_pos[q],2), np.round(y_pos[q],2), energy_density[p], flow[p][0], flow[p][1], flow[p][2], \
                        flow[p][3]/HBARC, pitautau, pitaux, pitauy, pitaueta, pixx, \
                        pixy, pixeta, piyy, piyeta, pietaeta]
                    file.write(" ".join(map(str, row)) + "\n")
    except IOError as e:
        print(f"Error writing to file '{filename_out}': {e}")

def generate_output_for_MUSIC(filename_out, data_read, tau_hydro, Nx, Ny, dx, dy):
    Tmunu_milne = get_Tmunu_Milne(data_read, tau_hydro)
    energy_density, flow = get_energy_density_and_flow(Tmunu_milne,tau_hydro)
    bulk_pressure, shear_tensor = get_shear_and_bulk_conformal_EOS(Tmunu_milne, tau_hydro, energy_density, flow, method=1)
    #bulk_pressure, shear_tensor = get_shear_and_bulk_conformal_EOS(Tmunu_milne, tau_hydro, energy_density, flow, method=2)
    
    energy_Music = interpolate_e(energy_density)
    
    write_MUSIC_output(filename_out, energy_Music, flow, shear_tensor, tau_hydro, Nx, Ny, dx, dy)

def main():
    # Check if all parameters are provided as command line arguments
    if len(sys.argv) != 8:
        print("Usage: python script.py <filename_in> <filename_out> <tau_hydro> <Nx> <Ny> <dx> <dy>")
        sys.exit(1)

    # parse input parameters
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    tau_hydro = float(sys.argv[3])
    Nx = int(sys.argv[4])
    Ny = int(sys.argv[5])
    dx = float(sys.argv[6])
    dy = float(sys.argv[7])

    # read the input file
    data = read_file(filename_in)

    # this creates similar output as the hydrodynamic output of KoMPoST
    generate_output_for_MUSIC(filename_out,data,tau_hydro,Nx,Ny,dx,dy)

if __name__ == "__main__":
    main()
