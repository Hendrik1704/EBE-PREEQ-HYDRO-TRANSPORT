# EBE-PREEQ-HYDRO-TRANSPORT

In this repository you find a numerical framework for (2+1)D event-by-event simulations of relativistic heavy-ion collisions.

## Modules

The code packages can be downloaded from their corresponding git repositories using the bash script `GetModulesFromGit.sh` and compiled with the `CompileFramework.sh` script. In the latter you can adjust the number of cores used to compile the codes.

At the moment the included modules are:

- Initial conditions:
    - Here you have to provide energy-momentum tensor files in the directory `input_energy_momentum_tensors`. The energy momentum tensor names should be of the format `Tmunu_Event#_Ns#.dat` with the event number and the number of grid sites in each direction. As format for the tensors please use the one from the KoMPoST code (link below).

- Pre-equilibrium evolution:
    - [KoMPoST](https://github.com/Hendrik1704/KoMPoST.git) (This is a slightly modified KoMPoST version, where $\nu_{\mathrm{eff}}=40$.)

- Relativistic viscous hydrodynamics:
    - [MUSIC](https://github.com/MUSIC-fluid/MUSIC)

- Cooper-Frye particlization:
    - [iSS](https://github.com/chunshen1987/iSS)

- Hadronic transport:
    - [smash](https://github.com/smash-transport/smash)

## Usage

Use the bash script `ExecuteEBE.sh` to change the parameters for the different codes and execute it for a production run. This file also contains comments for the different input parameters.
Be careful with parameter changes, as some parameters have to be changed consistently for multiple modules to the same values.

- Some parameters like `tau_EKT`, `tau_hydro`, `eta_s` and `grid_spacing=0.1` can be set at the top of the file.

If you have to abort a computation you can clear the result directories by running the `CleanResults.sh` script. Keep in mind that this script also deletes the outputs of SMASH at the end.
