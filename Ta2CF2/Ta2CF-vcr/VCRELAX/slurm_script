#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Ta2CF_vcr5'
#SBATCH --output=StdOut.o.%j
#SBATCH --error=StdErr.e.%j

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

TMP_DIR='/scratch/as3359/Ta2CF-VCrelax-5'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/Ta2CF/VCRELAX-5'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR

#========================================================================
#assign variables
fthresh=1.0d-3
ethresh=1.0d-4
ecut=65.0
rcut=600.0
elthresh=1.0d-7
pthresh=1.0
celldm=5.8997

nats=4
ntyp=3

atom1='Ta 180.94788  Ta.pbesol-spfn-rrkjus_psl.1.0.0.UPF'
atom2='C   12.0107   C.pbesol-n-rrkjus_psl.1.0.0.UPF'
atom3='F   18.9984  F.pbesol-n-rrkjus_psl.1.0.0.UPF'

atom1_position='Ta       0.3333333330        0.6666666670        0.0426707701    0   0   1'
atom2_position='C        0.0000000000        0.0000000000        0.0902850154    0   0   1'
atom3_position='Ta       0.6666666670        0.3333333330        0.1419119551   0   0   1'
atom4_position='F        0.0000000000        0.0000000000        0.2166719324   1   1   1'

nkx=12
nky=12

cell_parameter1='0.982215553  -0.000000000   0.000000000'
cell_parameter2='-0.491107776   0.850623622   0.000000000'
cell_parameter3='0.000000000   0.000000000   6.000000000'

prefix='Ta2CF'
# vc relax calculation
cat > $prefix.vcrelax.in << EOF
 &control
    calculation='vc-relax'
    restart_mode='from_scratch',
    prefix='Ta',
    tstress = .true.
    tprnfor = .true.
    forc_conv_thr=$fthresh,
    etot_conv_thr=$ethresh,
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    max_seconds = 56600
 /
 &system
    ibrav= 0, 
    celldm(1)= $celldm,
    nat=  $nats, ntyp= $ntyp,   
    ecutwfc =$ecut,
    ecutrho = $rcut,
    occupations='smearing',
    smearing='m-v', 
    degauss=0.005,
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  $elthresh
    diagonalization = 'cg'
 /
 &ions
    ion_dynamics='bfgs'
 /
 &cell
    cell_dynamics='bfgs', cell_dofree ='2Dxy'
 /
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (crystal)
 $atom1_position
 $atom2_position
 $atom3_position
 $atom4_position
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.vcrelax.in > $prefix.vcrelax.out

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep 'tau(' $prefix.vcrelax.out)
raw_cell_params=$(grep -A 3 'CELL_PARAMETERS' $prefix.vcrelax.out)
raw_f=$(grep 'Total force' $prefix.vcrelax.out | tail -1)
raw_p=$(grep 'bar' $prefix.vcrelax.out | tail -1)
current_time=$(grep 'time' $prefix.vcrelax.out | tail -1)

cat > Ta2C_H_raw_cartesian_atomic_positions << EOF
$raw_cart_at_pos
EOF

cat > Ta2C_H_raw_cell_parameters << EOF
$raw_cell_params
EOF

cat > Ta2C_H_raw_force << EOF
$raw_f
EOF

cat > Ta2C_H_raw_pressure << EOF
$raw_p
EOF

# Extract relevant data from these lines and store in the correct format for use as an input

cart_at_pos=$(awk '{print $2,$7,$8,$9}' Ta2C_H_raw_cartesian_atomic_positions)
current_cell_parameters=$(tail -3 Ta2C_H_raw_cell_parameters)
current_force=$(awk '{print $1,$2,$3,$4 }' Ta2C_H_raw_force)
current_pressure=$(awk '{print $5,$6 }' Ta2C_H_raw_pressure)

cat > Ta2C_H_cartesian_atomic_positions << EOF
$cart_at_pos
EOF

cat > tempfile_1 << EOF
$current_cell_parameters
EOF

# Split atomic positions file in two in order to separate initial and final sets

head -$nats Ta2C_H_cartesian_atomic_positions > Ta2C_H_init_cartesian_atomic_positions
tail -$nats Ta2C_H_cartesian_atomic_positions > testfile
mv testfile Ta2C_H_cartesian_atomic_positions

# Assign data to variables for future use

input_atomic_positions=$(cat Ta2C_H_init_cartesian_atomic_positions)
current_atomic_positions=$(cat Ta2C_H_cartesian_atomic_positions)

# Write key data to separate output file

initial_header="Initial positions:"
echo $initial_header >> Ta2C_H_atomic_positions.txt
cat Ta2C_H_init_cartesian_atomic_positions >> Ta2C_H_atomic_positions.txt

vcr1_cp_header="Cell parameters after vc-relax 1:"
echo $vcr1_cp_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

vcr1_ap_header="Atomic positions after vc-relax 1:"
echo $vcr1_ap_header >> Ta2C_H_atomic_positions.txt
cat Ta2C_H_cartesian_atomic_positions >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

# Initial relax calculation
cat > $prefix.relax.in << EOF
 &control
    calculation='relax'
    restart_mode='from_scratch',
    prefix='Ta',
    tstress = .true.
    tprnfor = .true.
    forc_conv_thr=$fthresh,
    etot_conv_thr=$ethresh,
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    max_seconds = 28800
 /
 &system
    ibrav= 0, 
    celldm(1)= $celldm,
    nat=  $nats, ntyp= $ntyp,   
    ecutwfc = $ecut,
    ecutrho = $rcut,
    occupations='smearing',
    smearing='m-v', 
    degauss=0.001,
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-12
    diagonalization = 'cg'
 /
 &ions
    ion_dynamics='bfgs'
 /
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (alat)
$current_atomic_positions
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS { alat }
$current_cell_parameters
EOF

$PW_COMMAND -input $prefix.relax.in > $prefix.relax.out

rm -rf $TMP_DIR