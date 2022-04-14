#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Ta2CO CT lat'
#SBATCH --output=StdOut.%j
#SBATCH --error=StdErr.%j

# Number of nodes (here 4 nodes with 16 CPUs each)
# The total number of nodes passed to mpirun will be nodes*ppn

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=spot-fsv2-64 
#SBATCH --qos=spot-fsv2-64

# Specify the account type and usage limits

#SBATCH --account=prj10_phase1
#SBATCH --time=05:00:00

#SBATCH --mail-user=as3359@bath.ac.uk
#SBATCH --mail-type=END

module purge
module use /apps/fsv2/modules/2021a/all/
module load QuantumESPRESSO/6.8-intel-2021a


# How to run QE codes:
PW_COMMAND='mpirun -np 2 pw.x - nk 1'
PH_COMMAND='mpirun -np 2 ph.x'
PH_IMAGE_COMMAND="mpirun -np 2 ph.x -ni 1"
DYNMAT_COMMAND='mpirun -np 1 dynmat.x'
PLOTBAND_COMMAND='mpirun -np 1 plotband.x'
BANDS_COMMAND='mpirun -np $NCORES bands.x'
Q2R_COMMAND='mpirun -np 1 q2r.x'
MATDYN_COMMAND='mpirun -np 1 matdyn.x'

export OMP_NUM_THREADS=2

TMP_DIR='/scratch/as3359/Mos2_ct_lat'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/FYP/Ta2CO/lat_const_CT'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR


#========================================================================

MX='Ta2CO'
atom1='Ta 180.94788  Ta.pbesol-spfn-rrkjus_psl.1.0.0.UPF'
atom2='C   12.0107   C.pbesol-n-rrkjus_psl.1.0.0.UPF'
atom3='O   15.9994  O.pbesol-n-rrkjus_psl.1.0.0.UPF'

atom1_position='C             0.0000000000        0.0000000000        0.0000000000  '
atom2_position='Ta            0.5237441000        0.2989654000       -0.4007024000  '
atom3_position='Ta            -0.0010323000        0.6027030000        0.4306046000  '
atom4_position='O             0.5210222000        0.3011869000        0.8013113000  '


cell_parameter1='1.023061336   0.000070588   0.000000000 '
cell_parameter2='-0.511469537   0.885657589   0.000000000'
cell_parameter3='0.000000000   0.000000000  10.000000000 '

for i in {5700..6000..100} ; do 
celldim=$(bc <<< "scale=3;$i/1000")
# self-consistent calculation
cat > $MX.scf.$celldim.in << EOF
 &control
    calculation='scf'
    restart_mode='restart',
    prefix='$MX',
    tstress = .true.
    tprnfor = .true.
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    max_seconds = 7200
 /
 &system
    ibrav= 0, 
    celldm(1)= $celldim,
    nat=  4, ntyp= 3,   
    ecutwfc = 150.00,
    ecutrho = 600.00
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-12
 /
 &ions
    ion_dynamics='bfgs'
 /
 &cell
    cell_dynamics='bfgs'
 /
CELL_PARAMETERS (alat)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (alat)
 $atom1_position
 $atom2_position
 $atom3_position
 $atom4_position
 K_POINTS {automatic}
 15 15 1 0 0 0 
EOF

$PW_COMMAND -input $MX.scf.$celldim.in > $MX.scf.$celldim.out

rm -rf $TMP_DIR

done
exit 0 