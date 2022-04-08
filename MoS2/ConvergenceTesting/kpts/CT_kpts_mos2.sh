#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='MoS2 CT Kpts'
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

TMP_DIR='/scratch/as3359/Ta2C_ct_kpts'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/FYP/MoS2/kpts_CT'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR


#========================================================================
# Variables required 

MX='MoS2'

iecutwfc=200

for k in {20..50..2} ; do 
# self-consistent calculation
cat > $MX.scf.$k.in << EOF
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
    celldm(1)= 6.072,
    nat=  3, ntyp= 2,   
    ecutwfc = $iecutwfc,
    ecutrho = $(( 5*$iecutwfc )),
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
   0.976383887  -0.000017274   0.000000000
  -0.487931186   0.846186578   0.000000000
   0.000000000   0.000000000   5.000000000 
ATOMIC_SPECIES
 Mo  95.94  Mo.pbesol-spn-kjpaw_psl.1.0.0.UPF
 S   32.065    S.pbesol-nl-kjpaw_psl.1.0.0.UPF 
ATOMIC_POSITIONS (alat)
Mo            0.0000000000        0.0000000000        1.6880765000
S             0.4883508000        0.2821722000        2.1737857000
S             0.4883508000        0.2821722000        1.2023670000
K_POINTS {automatic}
 $k $k 1 0 0 0 
EOF

$PW_COMMAND -input $MX.scf.$k.in > $MX.scf.$k.out

rm -rf $TMP_DIR

done
exit 0 