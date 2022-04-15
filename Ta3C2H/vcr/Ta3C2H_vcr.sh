#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Ta3C2H_loop'
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
#SBATCH --time=12:00:00

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

TMP_DIR='/scratch/as3359/Ta3C2H_loop'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/FYP/Ta3C2H/vcr-loop'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR

#========================================================================
#assign variables
fthresh=1.0d-3
ethresh=1.0d-4
ecut=60.00
rcut=300.00
elthresh=1.0d-12
pthresh=1.0
celldm=5.550

nats=6
ntyp=3

atom1='Ta 180.94788  Ta.pbesol-spfn-rrkjus_psl.1.0.0.UPF'
atom2='C   12.0107   C.pbesol-n-rrkjus_psl.1.0.0.UPF'
atom3='H    1.000   H.pbesol-rrkjus_psl.0.1.UPF'

atom1_position='Ta            0.0000231000        0.6073201000        0.6349311000 1 1 1 '
atom2_position='C             0.0000093000        0.0000234000        1.0244802000 1 1 1 '
atom3_position='Ta            0.5259269000        0.3036382000        1.4781889000 1 1 1 '
atom4_position='C            -0.0000061000        0.6072475000        1.9318977000 1 1 1 '
atom5_position='Ta           -0.0000199000       -0.0000492000        2.3214464000 1 1 1 '
atom6_position='H             0.3333333333        0.3333333333        2.6014464000 1 1 1 '


nkx=15
nky=15

cell_parameter1='1.051850560   0.000005600   0.000000000 '
cell_parameter2='-0.525920430   0.910903436   0.000000000'
cell_parameter3='0.000000000   0.000000000  10.000000000 '

prefix='Ta3C2'

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
 &cell
    cell_dynamics='bfgs', cell_dofree ='2Dxy'
 /
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
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
ATOMIC_POSITIONS (alat)
$current_atomic_positions
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS { alat }
$current_cell_parameters
EOF

$PW_COMMAND -input $prefix.relax.in > $prefix.relax.out

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep -A $nats 'ATOMIC_POSITIONS' $prefix.relax.out)
raw_f=$(grep 'Total force' $prefix.relax.out | tail -1)
raw_p=$(grep 'bar' $prefix.relax.out| tail -1)
current_time=$(grep 'time' $prefix.relax.out | tail -1)

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

r1_ap_header="Atomic positions after relax 1:"
echo $r1_ap_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

loop_counter=0

# Find whether relax thresholds have been reached 

f=$(awk '{print $4 }' Ta2C_H_raw_force)
fdelta=`echo $f - $fthresh | bc -l`
bigfdelta=`echo $fdelta \* 1000000 | bc -l`
bigfdelta=${bigfdelta%.*}

p=$(awk '{print $6 }' Ta2C_H_raw_pressure)
pp=${p%.*}
if [ $pp -lt 0 ] 
then p=`echo 0 - $p | bc -l` 
fi
pdelta=`echo $p - $pthresh | bc -l`
bigpdelta=`echo $pdelta \* 1000000 | bc -l`
bigpdelta=${bigpdelta%.*}

# If they have not, loop through vc-relax/relax cycles until they have been

while [ $bigfdelta -gt 0 ]||[ $bigpdelta -gt 0 ]  

do

loop_counter=$(($loop_counter+1))

#
# vc-relax calculation
#

cat > $prefix.scf.vcrelax$loop_counter.in << EOF
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
    conv_thr =  $elthresh
 /
 &ions
    ion_dynamics='bfgs',
 /
&cell
    cell_dynamics='bfgs', cell_dofree='2Dxy'
 /
CELL_PARAMETERS { alat }
$current_cell_parameters
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (alat)
$current_atomic_positions
K_POINTS {automatic}
 $nkx $nky 1 0 0 0
EOF

$ECHO " running the vc-relax calculation ; \c"
$PW_COMMAND -input $prefix.scf.vcrelax$loop_counter.in &> $prefix.scf.vcrelax$loop_counter.out
$ECHO " done"

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep -A $nats 'ATOMIC_POSITIONS' $prefix.scf.vcrelax$loop_counter.out)
raw_cell_params=$(grep -A 3 'CELL' $prefix.scf.vcrelax$loop_counter.out)
raw_f=$(grep 'Total force' $prefix.scf.vcrelax$loop_counter.out | tail -1)
raw_p=$(grep 'bar' $prefix.scf.vcrelax$loop_counter.out | tail -1)
current_time=$(grep 'time' $prefix.scf.vcrelax$loop_counter.out | tail -1)

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

current_atomic_positions=$(tail - $nats Ta2C_H_raw_cartesian_atomic_positions)
current_cell_parameters=$(tail -3 Ta2C_H_raw_cell_parameters)
current_force=$(awk '{print $1,$2,$3,$4 }' Ta2C_H_raw_force)
current_pressure=$(awk '{print $5,$6 }' Ta2C_H_raw_pressure)

cat > tempfile_1 << EOF
$current_atomic_positions
EOF

cat > tempfile_2 << EOF
$current_cell_parameters
EOF

# Write key data to separate output file

vcr2_cp_header="Cell parameters after vc-relax $loop_counter:"
echo $vcr2_cp_header >> Ta2C_H_atomic_positions.txt
cat tempfile_2 >> Ta2C_H_atomic_positions.txt

vcr2_ap_header="Atomic positions after vc-relax $loop_counter:"
echo $vcr2_ap_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

#
# relax calculation
#

cat > $prefix.scf.relax$loop_counter.in << EOF
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
    ion_dynamics='bfgs',
 /
CELL_PARAMETERS { alat }
$current_cell_parameters
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (alat)
$current_atomic_positions
K_POINTS {automatic}
 $nkx $nky 1 0 0 0
EOF

$ECHO " running the relax calculation ; \c"
$PW_COMMAND -input $prefix.scf.relax$loop_counter.in &> $prefix.scf.relax$loop_counter.out
$ECHO " done"

# Get relevant unprocessed lines from output file that contain relevant data

raw_cart_at_pos=$(grep -A $nats 'ATOMIC_POSITIONS' $prefix.scf.relax$loop_counter.out)
raw_f=$(grep 'Total force' $prefix.scf.relax$loop_counter.out | tail -1)
raw_p=$(grep 'bar' $prefix.scf.relax$loop_counter.out | tail -1)
current_time=$(grep 'time' $prefix.scf.relax$loop_counter.out | tail -1)

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

# Write key data to separate output file

r2_ap_header="Atomic positions after relax $loop_counter:"
echo $r2_ap_header >> Ta2C_H_atomic_positions.txt
cat tempfile_1 >> Ta2C_H_atomic_positions.txt

echo $current_force >> Ta2C_H_atomic_positions.txt
echo $current_pressure >> Ta2C_H_atomic_positions.txt
echo $current_time >> Ta2C_H_atomic_positions.txt

# Determine if relax thresholds have now been met

f=$(awk '{print $4 }' Ta2C_H_raw_force)
fdelta=`echo $f - $fthresh | bc -l`
bigfdelta=`echo $fdelta \* 1000000 | bc -l`
bigfdelta=${bigfdelta%.*}

p=$(awk '{print $6 }' Ta2C_H_raw_pressure)
pp=${p%.*}
if [ $pp -lt 0 ] 
then p=`echo 0 - $p | bc -l` 
fi
pdelta=`echo $p - $pthresh | bc -l`
bigpdelta=`echo $pdelta \* 1000000 | bc -l`
bigpdelta=${bigpdelta%.*}

done

#
# final calculation using relaxed parameters
#

cat > $prefix.scf.final.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='Ta',
    tstress = .true.
    tprnfor = .true.
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
    conv_thr =  $elthresh
 /
CELL_PARAMETERS (alat)
$current_cell_parameters
ATOMIC_SPECIES
 $atom1
 $atom2
 $atom3
ATOMIC_POSITIONS (alat)
$current_atomic_positions
K_POINTS {automatic}
$nkx $nky 1 0 0 0
EOF

$PW_COMMAND -input $prefix.scf.final.in &> $prefix.scf.final.out
$PW_COMMAND <$prefix.scf.final.in > $prefix.scf.final.out

rm -rf $TMP_DIR