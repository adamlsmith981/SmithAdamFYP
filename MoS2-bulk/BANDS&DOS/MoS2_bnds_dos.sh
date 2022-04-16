#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='MoS2_bulk_BND_DOS'
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
#SBATCH --time=02:00:00

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
BANDS_COMMAND='mpirun -np 2 bands.x'
Q2R_COMMAND='mpirun -np 1 q2r.x'
MATDYN_COMMAND='mpirun -np 1 matdyn.x'
DOS_COMMAND='mpirun -np 2 dos.x'

export OMP_NUM_THREADS=2

TMP_DIR='/scratch/as3359/MoS2_bulk_bandsdos'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/FYP/MoS2-BULK/bands_dos_results'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR

#========================================================================
#assign variables
fthresh=1.0d-3
ethresh=1.0d-4
ecut=150.00
rcut=700.00
elthresh=1.0d-12
pthresh=1.0
celldm=6.072

nats=6
ntyp=2

nbnds=32

atom1='Mo  95.94  Mo.pbesol-spn-kjpaw_psl.1.0.0.UPF'
atom2='S   32.065    S.pbesol-nl-kjpaw_psl.1.0.0.UPF'

atom1_position='Mo            0.0000000000        0.0000000000        0.9844997000'
atom2_position='S             0.4887119000        0.2819702000        1.4701535000'
atom3_position='S             0.4887135000        0.2819635000        0.4988906000'
atom4_position='Mo            0.0000000000        0.0000000000        3.0155003000'
atom5_position='S            -0.4887135000       -0.2819635000        3.5011094000'
atom6_position='S            -0.4887119000       -0.2819702000        2.5298465000'


nkx=15
nky=15

cell_parameter1='0.977325997  -0.000176330   0.000000000'
cell_parameter2='-0.488539747   0.846361973   0.000000000'
cell_parameter3='0.000000000   0.000000000   8.000000000'

prefix='MoS2'
# scf calculation
cat > $prefix.scf.in << EOF
 &control
    calculation='scf'
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
    degauss=0.001,
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  $elthresh
 /
 &ions
    ion_dynamics='bfgs'
 /
ATOMIC_SPECIES
 $atom1
 $atom2
ATOMIC_POSITIONS (alat)
 $atom1_position
 $atom2_position
 $atom3_position
 $atom4_position
 $atom5_position
 $atom6_position
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.scf.in > $prefix.scf.out

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep 'tau(' $prefix.scf.out)
raw_cell_params=$(grep -A 3 'CELL_PARAMETERS'$prefix.scf.out)
raw_f=$(grep 'Total force' $prefix.scf.out | tail -1)
raw_p=$(grep 'bar' $prefix.scf.out | tail -1)
current_time=$(grep 'time' $prefix.scf.out | tail -1)

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

vcr1_cp_header="Cell parameters after scf:"
echo $vcr1_cp_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

vcr1_ap_header="Atomic positions after scf:"
echo $vcr1_ap_header >> Ta2C_H_atomic_positions.txt
cat Ta2C_H_cartesian_atomic_positions >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

# nscf calculation
cat > $prefix.nscf.in << EOF
 &control
    calculation='nscf'
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
 /
ATOMIC_SPECIES
 $atom1
 $atom2
ATOMIC_POSITIONS (alat)
 $atom1_position
 $atom2_position
 $atom3_position
 $atom4_position
 $atom5_position
 $atom6_position
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.nscf.in > $prefix.nscf.out

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep -A $nats 'ATOMIC_POSITIONS' $prefix.nscf.out)
raw_f=$(grep 'Total force' $prefix.nscf.out | tail -1)
raw_p=$(grep 'bar' $prefix.nscf.out| tail -1)
current_time=$(grep 'time' $prefix.nscf.out| tail -1)

cat > Ta2C_H_raw_cartesian_atomic_positions << EOF
$raw_cart_at_pos
EOF

cat > Ta2C_H_raw_force << EOF
$raw_f
EOF

cat > Ta2C_H_raw_pressure << EOF
$raw_p
EOF

# Extract relevant data from these lines and store in the correct format for use as an input

current_atomic_positions=$(tail -$nats Ta2C_H_raw_cartesian_atomic_positions)
current_force=$(awk '{print $1,$2,$3,$4 }' Ta2C_H_raw_force)
current_pressure=$(awk '{print $5,$6 }' Ta2C_H_raw_pressure)

cat > tempfile_1 << EOF
$current_atomic_positions
EOF

# Write key data to seperate output file

r1_ap_header="Atomic positions after nscf:"
echo $r1_ap_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

# band calculation
cat > $prefix.bands.in << EOF
 &control
    calculation='bands'
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
    nat=  $nats, ntyp= $ntyp, nbnd= $nbnds, 
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
 /
ATOMIC_SPECIES
 $atom1
 $atom2
ATOMIC_POSITIONS (alat)
 $atom1_position
 $atom2_position
 $atom3_position
 $atom4_position
 $atom5_position
 $atom6_position
K_POINTS crystal_b
4
0.0000000000     0.0000000000     0.0000000000 30.0
0.6666666667    -0.3333333333     0.0000000000 30.0
0.5000000000    -0.5000000000     0.0000000000  30.0
0.0000000000     0.0000000000     0.0000000000 1.0
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.bands.in &> $prefix.bands.out
$PW_COMMAND <$prefix.bands.in > $prefix.bands.out

cat > $prefix.bandx.in << EOF
&BANDS
    prefix="Ta"
    outdir="$TMP_DIR/"
    filband="$prefix.Bandx.dat"
/
EOF
$BANDS_COMMAND <$prefix.bandx.in > $prefix.bandx.out

# DOS calculations 
cat > $prefix.dos.in << EOF
&DOS
    prefix="Ta"
    outdir="$TMP_DIR/"
    fildos="$prefix.dos.dat"
    emin=-9.0
    emax=16.0
/
EOF
$DOS_COMMAND <$prefix.dos.in > $prefix.dos.out

rm -rf $TMP_DIR