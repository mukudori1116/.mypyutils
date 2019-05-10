#!/usr/bin/python3

# -*- coding: utf-8 -*-

import os
import argparse
import requests


def make(file, text):
    with open(file, "w") as f:
        f.write(text)


def make_heat_files(project):

    HEAT_RUN_SH = """
#!/bin/sh

input="com"
prev_tpr="../minimize/min2.gro"
prev_cpt="../minimize/min2.trr"
for (( i=1;i<10; i++));do
    if [ $i -ne 1 ];then
        j=$((i-1))
        prev_tpr=md$j.tpr
        prev_cpt=md$j.cpt
    fi
    if [ ! -f md${i}.gro ];then
        rm -f md$i.out.mdp md$i.tpr
        gmx grompp \
            -f md$i.mdp \
            -c $prev_tpr \
            -r $prev_tpr \
            -p ../top/${input}.top \
            -n ../top/index.ndx \
            -t $prev_cpt \
            -o md$i.tpr \
            -po md$i.out.mdp || exit $?
        rm -f md$i.trr md$i.cpt md$i.gro md$i.edr md$i.log
        gmx mdrun \
            -v -s md$i -deffnm md$i -cpi md$i.cpt
    else
        echo "md${i}.gro present. Start next cycle."
    fi
done
"""
    MD1_MDP = """
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define       = -DPOSRES1000
; RUN CONTROL PARAMETERS
integrator   = md
dt           = 0.002
nsteps       = 100000
;
nstxout      = 500
nstlog       = 500
nstenergy    = 500
;
nstlist      = 5
ns_type      = grid
pbc          = xyz
rlist        = 0.8
;
coulombtype  = PME
rcoulomb     = 0.8
;
vdwtype      = cut-off
rvdw         = 0.8
dispcorr     = enerpres
;
optimize_fft = yes
;
tcoupl       = v-rescale
tc_grps      = Protein Non-Protein
tau_t        = 0.1 0.1
ref_t        = 300.0 300.0
;
pcoupl       = no
pcoupltype   = isotropic
tau_p        = 1.0
compressibility = 4.5e-5
ref_p        = 1.0
refcoord_scaling = all
;
gen_vel      = yes
gen_seed    = -1
;
continuation= no; Restarting after previous RUN
constraints  = hbonds
constraint_algorithm = LINCS

annealing    = single single
annealing_npoints = 2 2
annealing_time = 0 200 0 200
annealing_temp = 0 300 0 300
"""
    make(f'{project}/gromacs/heat/run.sh', HEAT_RUN_SH)
    make(f'{project}/gromacs/heat/md1.mdp', MD1_MDP)
    DEPOSRES = [1000, 500, 200, 100, 50, 20, 10, 0]
    for i, N in enumerate(DEPOSRES):
        MDN_MDP = f"""
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define    = -DEPOSRES{N}
; RUN CONTROL PARAMETERS
integrator   = md
dt           = 0.002
nsteps       = 50000
;
nstxout      = 500
nstlog       = 500
nstenergy    = 500
;
nstlist      = 5
ns_type      = grid
pbc          = xyz
rlist        = 0.8
;
coulombtype  = PME
rcoulomb     = 0.8
;
vdwtype      = cut-off
rvdw         = 0.8
dispcorr     = enerpres
;
optimize_fft = yes
;
tcoupl       = v-rescale
tc_grps      = Protein Non-Protein
tau_t        = 0.1 0.1
ref_t        = 300.0 300.0
;
pcoupl       = berendsen
pcoupltype   = isotropic
tau_p        = 1.0
compressibility = 4.5e-5
ref_p        = 1.0
refcoord_scaling = all
;
gen_vel      = no
;
continuation= yes; Restarting after previous RUN
constraints  = hbonds
constraint_algorithm = LINCS
"""
        make(f'{project}/gromacs/heat/md{i+2}.mdp', MDN_MDP)


def make_mini_files(project):
    RUN_SH = """
#!/bin/sh


input="../top/com"

rm -f min1.out.mdp min1.tpr

# Read configure file (mdp), initial coordinate file (gro) and topology file (top)
# Generate first tpr file (binary)
gmx grompp \
  -f min1.mdp \
  -po min1.out.mdp \
  -c ${input}.gro \
  -r ${input}.gro \
  -p ${input}.top \
  -o min1.tpr || exit 1
rm -f min1.log min1.gro min1.trr min1.edr

# Read initial minimize tpr file, run minimize, and output trajectory file (trr) + coordinate file (gro)
gmx mdrun -v \
  -s min1.tpr \
  -o min1.trr \
  -c min1.gro \
  -e min1.edr \
  -g min1.log || exit 1
rm -f min2.out.mdp min2.tpr

# Read configure file (mdp), initial coordinate file (gro) and topology file (top)
# Generate second tpr file (binary)
gmx grompp \
  -f min2.mdp \
  -po min2.out.mdp \
  -c min1.gro \
  -r ${input}.gro \
  -p ${input}.top \
  -o min2.tpr || exit 1
rm -f min2.log min2.gro min2.trr min2.edr

# Run second minimize
gmx mdrun -v \
  -s min2.tpr \
  -o min2.trr \
  -c min2.gro \
  -e min2.edr \
  -g min2.log
"""
    MIN1 = """
define       = -DPOSRES1000
;
cutoff-scheme = Verlet
integrator   = steep
nsteps       = 200
;
nstlog       = 1
;
nstlist      = 100
ns_type      = grid
pbc          = xyz
rlist        = 1.0
;
coulombtype  = PME
rcoulomb     = 1.0
;
vdwtype      = cut-off
rvdw         = 1.0
;
optimize_fft = yes
;
constraints  = none
"""
    MIN2 = """
define       =
;
cutoff-scheme = Verlet
integrator   = steep
nsteps       = 200
;
nstlog       = 1
;
nstlist      = 100
ns_type      = grid
pbc          = xyz
rlist        = 1.0
;
coulombtype  = PME
rcoulomb     = 1.0
;
vdwtype      = cut-off
rvdw         = 1.0
;
optimize_fft = yes
;
constraints  = none
"""
    make(f'{project}/gromacs/minimize/run.sh', RUN_SH)
    make(f'{project}/gromacs/minimize/min1.mdp', MIN1)
    make(f'{project}/gromacs/minimize/min2.mdp', MIN2)


def make_prod_files(project):
    RUN_SH = """
#!/bin/sh

input=com

# if md${input}.tpr doesn't exist, grompp will make the tpr file.
if [ ! -e md${input}.tpr ];then
    gmx grompp \
        -f mdrunset.mdp \
        -c ../heat/md9.tpr \
        -r ../heat/md9.tpr \
        -t ../heat/md9.cpt \
        -p ../top/${input}.top \
        -n ../top/index.ndx \
        -o md${input}.tpr \
        -po md${input}.mdp
fi

# gmx mdrun
gmx mdrun \
 -v -s md${input} -deffnm ${input} -cpi ${input}.cpt
"""
    MDRUN_MDP = """
define            =
; Run parameters
integrator      = md            ; leap-frog integrator
dt              = 0.002     ; 2 fs
nsteps          = 2500000  ; 2 fs * 2500000 = 5 ns
cutoff-scheme   = Verlet

; Output control
nstxout         = 5000          ; save coordinates every 10 ps
nstvout         = 5000            ; save velocities every 10 ps
nstenergy       = 2500              ; save energies every 5 ps
nstlog          = 5000        ; update log file every 10 ps

; Bond parameters
continuation    = yes           ; Restarting after previous RUN
constraint_algorithm = LINCS    ; holonomic constraints
constraints          = hbonds  ; Convert the bonds with H-atoms to constraints
lincs_iter           = 1            ; accuracy of LINCS
lincs_order          = 4              ; also related to accuracy

; Langevin Dynamics
bd_fric = 2.0           ; Brownian dynamics friction coefficient(correspond to gamma_ln in AMBER)
ld_seed = 1732  ; random seed

; Neighborsearching
ns_type         = grid          ; search neighboring grid cells
nstlist         = 20              ; 10 fs
rlist           = 1.0       ; short-range neighborlist cutoff (in nm)
rcoulomb        = 1.0         ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0           ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4               ; cubic interpolation
fourierspacing  = 0.12      ; grid spacing for FFT(nm)

; Temperature coupling is on
tcoupl    = V-rescale   ; modified Berendsen thermostat
tc-grps   = nowation Water_and_ions   ; two coupling groups - more accurate
tau_t     = 0.1  0.1               ; time constant, in ps
ref_t     = 300  300      ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl         = Parrinello-Rahman      ; Pressure coupling on in NPT
pcoupltype     = isotropic              ; uniform scaling of box vectors
tau_p          = 5.0              ; time constant, in ps
ref_p          = 5.0                ; reference pressure, in bar
compressibility = 4.5e-5  ; isothermal compressibility of water, bar^-1
refcoord_scaling = com

; Periodic boundary conditions
pbc                 = xyz               ; 3-D PBC

; Dispersion correction
DispCorr     = EnerPres ; account for cut-off vdW scheme

; Velocity generation
gen_vel         = no            ; Velocity generation
gen_temp        = 300.0         ; temperature for Maxwell distribution
gen_seed        = 1732          ; used to initialize random generator for random velocities
"""
    make(f'{project}/gromacs/pr/run.sh', RUN_SH)
    make(f'{project}/gromacs/pr/mdrunset.mdp', MDRUN_MDP)


def make_project_tree(project_name: str) -> None:
    DIRECTORIES = ['tmp', 'top', 'traj/mmpbsa',
                   'gromacs/minimize', 'gromacs/heat', 'gromacs/pr', 'gromacs/top']
    for d in DIRECTORIES:
        os.makedirs(f"{project_name}/{d}")
    make_mini_files(project_name)
    make_heat_files(project_name)
    make_prod_files(project_name)


def fetch_pdb(project, pdb_code):
    uri = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
    r = requests.get(uri)
    try:
        r.raise_for_status()
        with open(f"{project}/tmp/{pdb_code}_origin.pdb", "w") as f:
            f.write(r.text)
    except requests.exceptions.HTTPError as e:
        print(e)
        print("PDB code should not exists")


if __name__ == "__main__":

    os.chdir(os.getcwd())

    parser = argparse.ArgumentParser('Make Project tree for MD simulation')
    parser.add_argument('project', help='Project name')
    parser.add_argument('-f', '--fetch')
    args = parser.parse_args()

    make_project_tree(args.project)
    if args.fetch is not None:
        fetch_pdb(args.project, args.fetch)
