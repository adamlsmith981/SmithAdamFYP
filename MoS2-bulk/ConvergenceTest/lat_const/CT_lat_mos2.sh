#!/bin/bash
# lines starting with #SBATCH are options for the sbatch command
# version 26-10-2020

# Name of job and destinations for outputs

#SBATCH --job-name='Mos2-bulk CT lat'
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
#SBATCH --time=08:00:00

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
RESULTS_DIR='/shared/home/as3359/FYP/MoS2-BULK/lat_const_CT'

# not a restart

mkdir $TMP_DIR
mkdir $RESULTS_DIR
cp $0 $RESULTS_DIR
cd $RESULTS_DIR


#========================================================================

MX='MoS2'

for i in {6072..6080..1} ; do 
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
    nat=  6, ntyp= 2,   
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
   0.977325997  -0.000176330   0.000000000
  -0.488539747   0.846361973   0.000000000
   0.000000000   0.000000000   8.000000000  
ATOMIC_SPECIES
 Mo  95.94  Mo.pbesol-spn-kjpaw_psl.1.0.0.UPF
 S   32.065    S.pbesol-nl-kjpaw_psl.1.0.0.UPF 
ATOMIC_POSITIONS (alat)
Mo            0.0000000000        0.0000000000        0.9844997000
S             0.4887119000        0.2819702000        1.4701535000
S             0.4887135000        0.2819635000        0.4988906000
Mo            0.0000000000        0.0000000000        3.0155003000
S            -0.4887135000       -0.2819635000        3.5011094000
S            -0.4887119000       -0.2819702000        2.5298465000
K_POINTS {automatic}
 15 15 1 0 0 0 
EOF

$PW_COMMAND -input $MX.scf.$celldim.in > $MX.scf.$celldim.out

rm -rf $TMP_DIR

done
exit 0 