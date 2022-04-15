#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Ta2C Band Calc'
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
#SBATCH --time=04:00:00

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

export OMP_NUM_THREADS=2

TMP_DIR='/scratch/as3359/Ta2C-tmp'
PSEUDO_DIR='/shared/home/as3359/pseudo'
RESULTS_DIR='/shared/home/as3359/Ta2C/Bands'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR

#========================================================================
#assign variables
fthresh=1.0d-3
ethresh=1.0d-4
ecut=70.00
rcut=700.00
elthresh=1.0d-12
pthresh=1.0
celldm=5.550

nats=3
ntyp=2

nbnds=42

atom1='Ta 180.94788  Ta.pbesol-spfn-rrkjus_psl.1.0.0.UPF'
atom2='C   12.0107   C.pbesol-n-rrkjus_psl.1.0.0.UPF'

atom1_position='C             0.0000000000        0.0000000000        0.0000000000'
atom2_position='Ta            0.5202621000        0.3003735000       -0.4102757000'
atom3_position='Ta            0.0000000000        0.6007470000        0.4102757000'

nkx=20
nky=20

cell_parameter1='1.040524263   0.000000000   0.000000000'
cell_parameter2='-0.520262131   0.901120444   0.000000000'
cell_parameter3='0.000000000   0.000000000  10.000000000'

prefix='Ta2C'
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
    nbnd=$nbnds,
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
ATOMIC_SPECIES
 $atom1
 $atom2
ATOMIC_POSITIONS (alat)
 $atom1_position
 $atom2_position
 $atom3_position
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.scf.in > $prefix.scf.out


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
ATOMIC_POSITIONS (crystal)
 $atom1_position
 $atom2_position
 $atom3_position
K_POINTS {automatic}
 $nkx $nky 1 0 0 0 
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

$PW_COMMAND -input $prefix.nscf.in > $prefix.nscf.out


# bands calculation
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
    nat=  $nats, ntyp= $ntyp, nbnd= 32,  
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
ATOMIC_POSITIONS (crystal)
 $atom1_position
 $atom2_position
 $atom3_position
K_POINTS crystal_b
         4
        0.0000000000     0.0000000000     0.5000000000    30.0
        0.6666666667    -0.3333333333     0.5000000000    30.0
        0.5000000000    -0.5000000000     0.5000000000    30.0
        0.0000000000     0.0000000000     0.5000000000    1.0
CELL_PARAMETERS (alat= $celldm)
 $cell_parameter1
 $cell_parameter2
 $cell_parameter3
EOF

# $PW_COMMAND -input $prefix.bands.in &> $prefix.bands.out
$PW_COMMAND <$prefix.bands.in > $prefix.bands.out

# bands calculation
cat > $prefix.bandx.in << EOF
&BANDS
  prefix="Ta"
  outdir="$TMP_DIR/"
  filband="$prefix.Bandx.dat"
/
EOF
$BANDS_COMMAND <$prefix.bandx.in > $prefix.bandx.out
rm -rf $TMP_DIR