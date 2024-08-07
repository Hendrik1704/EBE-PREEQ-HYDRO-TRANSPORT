# EBE-PREEQ-HYDRO-TRANSPORT

In this repository you find a numerical framework for (2+1)D event-by-event simulations of relativistic heavy-ion collisions.

## Modules and instalation

The code packages can be downloaded from their corresponding git repositories using the bash script `./GetModulesFromGit.sh` and compiled with the `./CompileFramework.sh` script. In the latter you can adjust the number of cores used to compile the codes.

At the moment the included modules are:

- Initial conditions:
    - Here you have to provide energy-momentum tensor files in the directory `input_energy_momentum_tensors`. The energy momentum tensor names should be of the format `Tmunu_Event#_Ns#.dat` with the event number and the number of grid sites in each direction. As format for the tensors please use the one from the KoMPoST code (link below).

- Pre-equilibrium evolution:
    - [KoMPoST](https://github.com/Hendrik1704/KoMPoST.git) This is a slightly modified KoMPoST version, where $\nu_{\mathrm{eff}}=40$.

- Code to do the pre-equilibrium/hydro matching:
    - KoMPoST_to_MUSIC.py: From a conformal pre-equilibrium to a non-conformal hydrodynamics the user can choose to do the matching using energy or entropy. The output will be a modified energy density that will be used in MUSIC.

- Relativistic viscous hydrodynamics:
    - [MUSIC](https://github.com/MUSIC-fluid/MUSIC) Modified version which can read in the bulk pressure and modified energy density from the `KoMPoST_to_MUSIC.py` script.

- Cooper-Frye particlization:
    - [iSS](https://github.com/chunshen1987/iSS)

- Converter:
    - `convert_OSCAR1997A_toOSCAR2013.py`: Convert the output from iSS to smash input.

- Hadronic transport:
    - [smash](https://github.com/smash-transport/smash) (version 3.1 with [PYTHIA](https://pythia.org/) 8.310)

## Usage

### Parameters

Use the bash script `ExecuteEBE.sh` or `ExecuteEBE_cluster.sh` (depending of the way you will run) to change the parameters for the different codes and execute it for a production run.

- Some parameters like `tau_EKT`, `tau_hydro`, `eta_s`, `grid_spacing`, `hydro_oversampling`
and `type_of_matching` can be set at the top of the file. In case of a specific analysis, you can go along the file and modify specific code parameters. 

### Running

- For running in your computer you can only use 
```
./ExecuteEBE.sh
``` 
which will do the analysis for all files in the folder `input_energy_momentum_tensors`.
The executable will run in a loop of the input files for each program.

For running a larger number of events, on a cluster with `slurm` 
(e.g. the noctua cluster in Paderborn, Germany), you can also use
```
sbatch noctua_script.sh
```
which essentially will run `ExecuteEBE_cluster.sh` with the configurations 
of parameters.
The cluster scripts are located in the `cluster_support` directory.
This code is essentially the same code previously but prepared to run a 
parallelized array of events.
In this case, all the modules run for each event in different jobs (job array).
Pay attention in this case, because you always need to adjust the size of the 
array `#SBATCH --array=0-9`, it needs to match with the total number of events 
inside of the input folder.

:exclamation: The execution script `noctua_script.sh` is specifically designed 
for the noctua cluster at the Paderborn University in Germany.
You might have to write your own script to run it on a different machine.
For the python scripts you might have to create your own virtual environment.
This script is only meant as an example.

:exclamation: To prepare the framework on the cluster, you can use 
the `noctua_PrepareFrameworkProductionRun.sh` job script, which downloads and 
installs all the modules in a separate slurm job.

If you have to abort a computation you can clear the result directories by 
running the `CleanResults.sh` script.
Keep in mind that this script also deletes the outputs of SMASH at the end.
