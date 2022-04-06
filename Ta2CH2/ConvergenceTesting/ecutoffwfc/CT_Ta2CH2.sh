#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Ta2CH2 CT_ecutwfc'
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
#SBATCH --time=06:00:00

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

TMP_DIR='/scratch/as3359/Ta2CH2-CT-ecutoff'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/FYP/Ta2CH2/ecutwfc_CT'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR
#========================================================================
# variables

MX='Ta2CH2'
atom1='Ta 180.94788  Ta.pbesol-spfn-rrkjus_psl.1.0.0.UPF'
atom2='C   12.0107   C.pbesol-n-rrkjus_psl.1.0.0.UPF'

atom1_position='C             0.0000000000        0.0000000000        0.0000000000 0 0 0'
atom2_position='Ta            0.6666666670        0.3333333330       -0.0410149592 1 1 1'
atom3_position='Ta            0.3333333330        0.6666666670        0.0410149592 1 1 1'
atom4_position='H            -0.0033389048       -0.0005507363        0.0924039106 1 1 1'
atom5_position='H             0.0005507363        0.0033389048        -0.0924039106 1 1 1'


cell_parameter1='1.040617955   0.000000000   0.000000000'
cell_parameter2='-0.520308977   0.901201584   0.000000000'
cell_parameter3='0.000000000   0.000000000  10.000000000'

#========================================================================

for iecutwfc in {60..300..5} ; do 
# self-consistent calculation
cat > $MX.scf.$iecutwfc.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='$MX',
    tstress = .true.
    tprnfor = .true.
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    max_seconds = 28800
 /
 &system
    ibrav= 0, 
    celldm(1)= 5.550,
    nat=  5, ntyp= 3,  
    ecutwfc = $iecutwfc,
    ecutrho = $(( 5*$iecutwfc )),
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-12
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
C             0.0000000000        0.0000000000        0.0000000000
Ta            0.5088008000        0.2948413000       -0.4364511000
Ta           -0.0007447000        0.5869407000        0.4364532000
H            -0.0017928000       -0.0015314000        0.8341930000
H             0.0002769000        0.0024412000       -0.8341930000
K_POINTS {automatic}
12 12 1 0 0 0
EOF

$PW_COMMAND -input $MX.scf.$iecutwfc.in > $MX.scf.$iecutwfc.out

rm -rf $TMP_DIR

done
exit 0 