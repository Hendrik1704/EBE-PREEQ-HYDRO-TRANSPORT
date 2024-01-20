#!/bin/bash

echo "Execute the EBE evolution for the energy momentum tensors stored in 'input_energy_momentum_tensors'"
echo "The energy momentum tensors names should be of the format 'Tmunu_Event#_Ns#.dat'"
echo "The first number is the EventID and the second one the number of grid points in each direction Ns"
echo "The first line of the energy momentum tensors should contain a header!"

tau_EKT=0.2005
tau_hydro=0.8
eta_s=0.16
grid_spacing=0.1
hydro_oversampling=1

mkdir KoMPoST_output
mkdir KoMPoST_output_transformed
mkdir MUSIC_FOsurfaces
mkdir MUSIC_InputParameters
mkdir iSS_output
mkdir iSS_output_converted
mkdir smash_output

cd KoMPoST
INPUT_TMUNU_PATH="../input_energy_momentum_tensors/*"
for FILE in $INPUT_TMUNU_PATH
do

echo "Processing file: $FILE"
# Extract event number and NS from the file name
EVENTNUMBER=$(echo "$FILE" | grep -o -E '[0-9]+' | head -n1)
NS=$(echo "$FILE" | grep -o -E 'Ns([0-9]+)' | sed 's/Ns//')
echo "Create KoMPoST input file for event: $EVENTNUMBER"
echo "Create KoMPoST input file for grid with Ns: $NS"

# KoMPoST input file
# Change parameters here if needed!
cat <<EOF >> parameters_KoMPoST.ini
[KoMPoSTInputs]
tIn = $tau_EKT;
tOut = $tau_hydro;
inputfile = $FILE
outputfiletag = ../KoMPoST_output/$EVENTNUMBER.Tmunu

[KoMPoSTParameters]
EtaOverS = $eta_s
EtaOverSTemperatureScale = 0.0
# 0 for free-streaming, 1 for "KoMPoST" EKT evolution
EVOLUTION_MODE = 1
# 0 or 1
ENERGY_PERTURBATIONS = 1
# 0 or 1
MOMENTUM_PERTURBATIONS = 0
DECOMPOSITION_METHOD = 1

[EventInput]
afm = $grid_spacing
ns = $NS
xstart = 0
xend = `echo ${NS} 1 | awk '{print $1-$2}'`
ystart = 0
yend = `echo ${NS} 1 | awk '{print $1-$2}'`
EOF

echo "Execute KoMPoST"
./KoMPoST.exe parameters_KoMPoST.ini
# clear the KoMPoST directory of all files which are not needed any more
rm parameters_KoMPoST.ini
rm *.txt
cd ../KoMPoST_output
# delete all the unused output to save disc space
find . -type f ! -name '*music_init_flowNonLinear_pimunuTransverse.txt' -delete
cd ..

echo "Transforming $FILE into MUSIC input"
python3 KoMPoST_to_MUSIC.py ./KoMPoST_output/$EVENTNUMBER.Tmunu.music_init_flowNonLinear_pimunuTransverse.txt ./KoMPoST_output_transformed/$EVENTNUMBER.Tmunu.txt
cd KoMPoST
done

cd ../MUSIC

INPUT_TMUNU_PATH_MUSIC="../KoMPoST_output_transformed/*"
for FILE in $INPUT_TMUNU_PATH_MUSIC
do

echo "Processing file: $FILE"
echo "Create MUSIC input file MODE 2"

id=$(echo "$FILE" | grep -o -E '[0-9]+')
EVENTNUMBER=$(echo $id | cut -d' ' -f1)

# MUSIC input file
# Change parameters here if needed!
cat <<EOF >> ${EVENTNUMBER}.parameters_MUSIC.ini
###################################
# parameters list
###################################
#
echo_level 1    # controls the amount of messages in the terminal
mode 2          # this mode is evolution only
#
###################################
# parameters for initial conditions
###################################
#
Initial_profile 94               # Read in initial profile from a file, 94: full T^\mu\nu and read bulk
initialize_with_entropy 0       # 0: with energy density
#
Initial_Distribution_input_filename $FILE
#
s_factor  1.0   # normalization factor for initial profile
#
#######################################
# parameters for hydrodynamic evolution
#######################################
#
boost_invariant 1       # whether the simulation is boost-invariant
#
# grid information
Initial_time_tau_0 $tau_hydro       # starting time of the hydrodynamic evolution (fm/c)
Total_evolution_time_tau 50.    # the maximum allowed running evolution time (fm/c) (need to be set to some large number)
Delta_Tau 0.005                 # time step to use in the evolution [fm/c]
#
EOS_to_use 91                 # type of the equation of state
                              # 0: ideal gas
                              # 1: EOS-Q from azhydro
                              # 2: lattice EOS s95p-v1
                              #    (from Huovinen and Petreczky)
                              # 3: lattice EOS s95p with partial
                              #    chemical equilibrium (PCE) at 150 MeV
                              #    (see https://wiki.bnl.gov/TECHQM
                              #         /index.php/QCD_Equation_of_State)
                              # 4: lattice EOS s95p PCE at 155 MeV
                              # 5: lattice EOS s95p PCE at 160 MeV
                              # 6: lattice EOS s95p PCE at 165 MeV
                              # 7: lattice EOS s95p-v1.2 for UrQMD
                              # 9: lattice EOS hotQCD with UrQMD
                              # 91: lattice EOS hotQCD with SMASH
                              # 14: lattice EOS hotQCD at finite muB
#
# transport coefficients
Viscosity_Flag_Yes_1_No_0 1     # turn on viscosity in the evolution
Include_Shear_Visc_Yes_1_No_0 1 # include shear viscous effect
Shear_to_S_ratio $eta_s           # value of \eta/s
T_dependent_Shear_to_S_ratio  0 # flag to use temperature dep. \eta/s(T)
Include_Bulk_Visc_Yes_1_No_0 0  # include bulk viscous effect
T_dependent_Bulk_to_S_ratio 8   # include Temperature-dependent \zeta/s(T)
Include_second_order_terms 1    # include second order non-linear coupling terms
Include_Rhob_Yes_1_No_0 0
turn_on_baryon_diffusion 0
kappa_coefficient 0.0
#
# switches to output evolution information
output_evolution_data 0y                 # flag to output evolution history to file
output_movie_flag 0
output_evolution_T_cut 0.145
outputBinaryEvolution  1                # output evolution file in binary format
output_evolution_every_N_eta  1         # output evolution file every Neta steps
output_evolution_every_N_y  1           # output evolution file every Ny steps
output_evolution_every_N_x  1           # output evolution file every Nx steps
output_evolution_every_N_timesteps  200  # output evolution every Ntime steps
#
#
###########################################
# parameters for freeze out and Cooper-Frye
###########################################
Do_FreezeOut_Yes_1_No_0 1                       # flag to find freeze-out surface
Do_FreezeOut_lowtemp 1                          # flag to include cold corona
freeze_out_method 4                             # method for hyper-surface finder
                                                # 4: Cornelius
freeze_surface_in_binary 0                      # switch to output surface file in binary format
average_surface_over_this_many_time_steps 10    # the step skipped in the tau
#
freeze_Ncell_x_step 1
freeze_Ncell_y_step 1
freeze_Ncell_eta_step 1
freeze_eps_flag 0
N_freeze_out 1
use_eps_for_freeze_out 0        # flag to use energy density as criteria to
                                # find freeze-out surface
                                # 0: use temperature, 1: use energy density
T_freeze 0.155                  # freeze out temperature
#
#
EndOfData
EOF

./MUSIChydro ${EVENTNUMBER}.parameters_MUSIC.ini
mv surface_* ${EVENTNUMBER}.surface.dat
mv ${EVENTNUMBER}.parameters_MUSIC.ini ../MUSIC_InputParameters

mv *.surface.dat ../MUSIC_FOsurfaces
rm *.dat
done

cd ../iSS

echo "Create iSS input file"
# iSS input file
# Change parameters here if needed!
cat <<EOF >> parameters_iSS.ini
hydro_mode = 1           # mode for reading in freeze out information
                         # 1: reads outputs from MUSIC assuming
                         #    boost-invariant

afterburner_type = 2     # 0: PDG_Decay
                         # 1: UrQMD
                         # 2: SMASH

turn_on_bulk = 0         # read in bulk viscous pressure
turn_on_rhob = 0         # read in net baryon chemical potential
turn_on_diff = 0         # read in baryon diffusion current

include_deltaf_shear = 1      # include delta f contribution from shear
include_deltaf_bulk = 0       # include delta f contribution from bulk
include_deltaf_diffusion = 0  # include delta f contribution from diffusion

bulk_deltaf_kind = 1     # 0 : 14-moment approximation (s95p-v0-PCE)
                         # 1 : relaxation time approximation [default]
                         # 11: OSU 14-moment
                         # 20: WSU 22-moment (NEoS-BQS) shear and bulk
                         # 21: WSU Chapman-Enskog (NEoS-BQS) shear and bulk

restrict_deltaf = 1      # flag to apply restriction on the size of delta f
deltaf_max_ratio = 1.0   # the maximum allowed size of delta f w.r.t f0

quantum_statistics = 1   # include quantum statistics or not (1: yes, 0: no)

output_samples_into_files = 1  # output particle samples into individual files
                               # for individual particle species
store_samples_in_memory = 0  # flag to store particle samples in memory
use_OSCAR_format = 1         # output results in OSCAR format
use_gzip_format = 0          # output results in gzip format (only works with
                             # store_samples_in_memory = 1)
use_binary_format = 0        # output results in binary format
perform_decays = 0           # flag to perform resonance decay
perform_checks = 0           # perform tests for particle samples
include_spectators = 0       # include spectators (filename: spectators.dat)

local_charge_conservation = 0  # flag to impose local charge conservation

calculate_vn = 0         # 1/0: whether to calculate the
                         # dN/(pt dpt dphi dy) and v_n flows
                         # (they not required for MC-sampling)

MC_sampling = 2          # 0/1/2/3: whether to perform Monte-Carlo sampling
                         # (not required for spectra calculation).
                         # 0: No sampling.
                         # 1: use dN_dxtdetady to sample.
                         # 2: use dN_dxtdy to sample.
                         # 3: use dN_pTdpTdphidy to sample
                         #    (overwrites calculate_vn to be 1).
                         # 4: use FSSW (fast)
                         # Since this parameter controls whether to
                         # calculate certain arrays, it controls whether to
                         # perform other related calculations (see below).

dN_dy_sampling_model = 30    # Controls how an non-integer dN_dy is sampled to
                             # produce an integer that can be used in actual
                             # sampling
                             # -- all numbers below 100 are reserved for
                             #    "conventional" sampling where actual particle
                             #    yields are used to constrain the sampling.
                             #    101-200 are reserved for sampling using total
                             #    energy flux from freeze-out surface.
                             # -- 1: The fractional part of dN_dy is used as a
                             #    probability to determine whether there is 1
                             #    more particle
                             # -- 10: Use NBD to sample the fractional part
                             #    of dN_dy. The probability is calculated so
                             #    that p=nb/(nb+k) and k is assumed to be
                             #    proportional to n, and the proportionality
                             #    constant is given by the parameter
                             #    dN_dy_sampling_para1.
                             # -- 20: Use NBD to sample the whole dN_dy.
                             #    k is assumed to be proportional to n,
                             #    with coefficient dN_dy_sampling_para1.
                             # -- 30: Use Poisson distribution to sample the
                             #    whole dN_dy. The mean value is set to be
                             #    dN_dy
                             # -- 110: Total energy flux from freeze-out
                             #    surface is used to constrain dN_dy.
                             #    Whenever the total energy carried away from
                             #    all the sampled particles exceed the total
                             #    energy flux, the sampling procedure stops.


dN_dy_sampling_para1 = 0.16  # Additional parameters for dN/dy sampling.
                             # -- For dN_dy_sampling_model==10 or 20,
                             #    this parameter sets the ratio k/n for NBD,
                             #    see explanation for dN_dy_sampling_model.

y_LB = -5.0                  # lower bound for y-sampling;
                             # used in "conventional" sampling
y_RB = 5.0                   # upper bound for y-sampling; used in
                             # "conventional" sampling

eta_s_LB = -5.0              # lower bound for eta_s sampling; used only when
                             # sampling using total energy flux
eta_s_RB = 5.0               # upper bound for eta_s sampling.

use_dynamic_maximum = 0      # 0/1: Whether to automatically reduce the
                             # guessed maximum after some calculations.
                             # Work only when MC_sampling is set to 2.
adjust_maximum_after = 100000    # Used only when use_dynamic_maximum=1.
                                 # After the number of sampling given by
                                 # this parameter the guessed maximum is
                                 # adjusted.
adjust_maximum_to = 1.2      # [1,inf]: When guessed maximum is adjusted,
                             # it is adjusted to the "observed maximum"
                             # multiplied by this value. Note that the
                             # "observed maximum" is measured relative to
                             # the guessed maximum. See code for details.

grouping_particles = 1       # 0/1: Particles will be re-order according to
                             # their mass. This parameter combined with
                             # grouping_tolerance parameter can make particles
                             # with similar mass and chemical potentials to be
                             # sampled together.
grouping_tolerance = 0.01    # If two particles adjacent in the table have
                             # mass and chemical potentials close within this
                             # relative tolerance, they are considered to be
                             # identical and will be sampled successively
                             # without regenerating the dN / (dxt deta dy)
                             # matrix for efficiency.

use_historic_flow_output_format = 0    # 1/0: The "historical flow format"
                                       # means to output flows for all
                                       # particles in a single v2** file and
                                       # to add lines with particles names to
                                       # the v2** file. Turn this option off
                                       # to enbrace the new way of outputting
                                       # which allows the using of parameter
                                       # grouping_particles to speed up the
                                       # calculations.

calculate_vn_to_order = 9              # v_n's are calculated up to this order

sample_upto_desired_particle_number = 0  # flag to run sampling until desired
                                         # particle numbers is reached
number_of_particles_needed = 100000      # number of hadrons to sample
number_of_repeated_sampling = $hydro_oversampling     # How many times should the sampling be
                                       # repeated.
maximum_sampling_events = 10000

sample_pT_up_to = -1                   # Up to this value will pT be sampled;
                                       # if<0 then use the largest value in
                                       # the pT table.
sample_y_minus_eta_s_range = 3         # y_minus_eta_s will be sampled between
                                       # +- this value. It is used only when
                                       # sampling using
                                       # sample_using_dN_dxtdy_4all_particles
                                       # function.

use_pos_dN_only = 0                    # 1/0: When set to 1, all negative
                                       # emission functions will be skipped.
                                       # Effects the both dN_ptdptdphidy and
                                       # dN_dxtdetady calculations.

minimum_emission_function_val = 1e-30  # If dN/(dx_t deta dy) is evaluated to
                                       # be smaller than this value, then it
                                       # is replaced by this value.

calculate_dN_dtau = 0      # Output dN_dtau table. Only applicable
                           # if MC_sampling parameter is set to 1.
bin_tau0 = 0.6             # used to generate bins for
                           # calculate_dN_dtau_using_dN_dxtdeta function
bin_dtau = 0.2             # used to generate bins for
                           # calculate_dN_dtau_using_dN_dxtdeta function
bin_tau_max = 17.0         # used to generate bins for
                           # calculate_dN_dtau_using_dN_dxtdeta function

calculate_dN_dx = 0        # Output dN_dx table. Only applicable
                           # if MC_sampling parameter is set to 1.
bin_x_min = -10.0          # used to generate bins for
                           # calculate_dN_dx_using_dN_dxtdeta function
bin_dx = 0.5               # used to generate bins
                           # for calculate_dN_dx_using_dN_dxtdeta function
bin_x_max = 10.0           # used to generate bins for
                           # calculate_dN_dx_using_dN_dxtdeta function

calculate_dN_dphi = 0      # Output dN_dphi table. Only applicable
                           # if calculate_vn parameter is set to 1.
calculate_dN_deta = 1      # Output dN_deta table. Only applicable
                           # if MC_sampling parameter is set to 1.
calculate_dN_dxt = 1       # Output dN_dxt table. Only applicable
                           # if MC_sampling parameter is set to 1.

output_dN_dxtdy_4all = 0   # Output dN_dxtdy table. Only applicable
                           # if MC_sampling parameter is set to 2.

randomSeed = 0            # If <0, use system clock.
EOF

OUTPUT_PATH_MUSIC="../MUSIC_FOsurfaces/*"
for FILE in $OUTPUT_PATH_MUSIC
do

echo "Processing file: $FILE"
mkdir results
cp ${FILE} results

id=$(echo "$FILE" | grep -o -E '[0-9]+')
EVENTNUMBER=$(echo $id | cut -d' ' -f1)

cp ../MUSIC_FOsurfaces/${EVENTNUMBER}.surface.dat results
cp ../MUSIC_InputParameters/${EVENTNUMBER}.parameters_MUSIC.ini results
mv results/${EVENTNUMBER}.surface.dat results/surface.dat
mv results/${EVENTNUMBER}.parameters_MUSIC.ini results/music_input


./iSS.e parameters_iSS.ini
rm -rf results
mv OSCAR.DAT ../iSS_output/OSCAR${EVENTNUMBER}
done

cd ..

echo "Convert OSCAR format from 1997A version to 2013 version"
OUTPUT_PATH_ISS="./iSS_output/*"
for FILE in $OUTPUT_PATH_ISS
do
python3 convert_OSCAR1997A_to_OSCAR2013.py ${FILE} ./iSS_output_converted/
done

rm ./iSS/parameters_iSS.ini
cd smash/build

OSCAR_FILES_PATH="../../iSS_output_converted/*"
for FILE in $OSCAR_FILES_PATH
do

id=$(echo "$FILE" | grep -o -E '[0-9]+')

echo "Processing file: $FILE"
echo "Create SMASH input file"

# SMASH input file
# Change parameters here if needed!
cat <<EOF >> parameters_smash.yaml
Logging:
    default: INFO

General:
    Modus:         List
    Time_Step_Mode: None
    Delta_Time:    0.1
    End_Time:      1000.0
    Randomseed:    1
    Nevents:       $hydro_oversampling

Output:
    Output_Interval: 10.0
    Particles:
        Format:          ["Oscar2013"]

Modi:
    List:
        # If the build directory is not located in the smash directory anymore,
        # the absolute path specified below will not work anymore.
        # You can alternatively pass the path directly from the command line
        # with the "-c" command:
        # ./smash -i <path to config file> -c 'Modi: { List: { File_Directory: <path to file that is read in> } }
        File_Directory: "../../iSS_output_converted"
        File_Prefix: "OSCAR"
        Shift_Id: $id
EOF

mkdir ../../smash_output/Event${id}
./smash -i ./parameters_smash.yaml -o ../../smash_output/Event${id}
rm parameters_smash.yaml
done
