#!/bin/bash
#SBATCH --job-name=pigknee_sim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=interact
#SBATCH --account=cuuser_finou_biomedical_radiation
#SBATCH --output=output_%j.log

echo "[$(date)] Setting up simulation environment..."

module load spack
spack load gcc@12.3.0
spack load geant4@11.2.2

export CC=$(spack location -i gcc@12.3.0)/bin/gcc
export CXX=$(spack location -i gcc@12.3.0)/bin/g++
export GEANT4_DIR=$(spack location -i geant4@11.2.2)/lib/cmake/Geant4
export LD_LIBRARY_PATH=$(spack location -i geant4@11.2.2)/lib:$LD_LIBRARY_PATH

cd /scratch/kttucke/g4_pigknee/build || exit 1

PHOTON_COUNT=1000000
MACRO_FILE="temp_run.mac"
echo "/run/initialize" > ${MACRO_FILE}
echo "/run/beamOn ${PHOTON_COUNT}" >> ${MACRO_FILE}

echo "[$(date)] Starting Geant4 simulation..."
./pigknee ${MACRO_FILE} > run_${PHOTON_COUNT}.log 2>&1
RETURN_CODE=$?

echo "[$(date)] Simulation completed with exit code: $RETURN_CODE"
exit $RETURN_CODE
