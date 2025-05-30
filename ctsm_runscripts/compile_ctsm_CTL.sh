#! /bin/bash

# Script to compile CLM with case-specific settings
# For standalone CLM5.0 or CTSM development version

# based on compile script from Petra Sieber (https://github.com/pesieber/CTSM/blob/release-clm5.0/compile_cclm2.sh)
# adjusted by Inne Vanderkelen

set -e # failing commands will cause the shell script to exit


#==========================================
# Case settings
#==========================================

echo "*** Setting up case ***"

date=`date +'%Y%m%d-%H%M'` # get current date and time
startdate=`date +'%Y-%m-%d %H:%M:%S'`

COMPSET=IHistClm51Sp # I2000Clm50SpGs for release-clm5.0 
RES=hcru_hcru_mt13  #hcru_hcru_mt13 #f09_g17 for test glob
DOMAIN=glob # eur for CCLM2 (EURO-CORDEX), sa for South-America, glob for global 

CODE=CTSMdev # clm5.0 for official release, CTSMdev for latest master branch (requires ESMF installation)

COMPILER=gnu # gnu for gnu/gcc
COMPILERNAME=gcc # gcc for gnu/gcc

SCEN=CTL # scenario
EXP=test_${date}
FLAG_IRRIG=false
CASENAME=$COMPSET.$RES.$SCEN.spunup

DRIVER=nuopc # mct for clm5.0, mct or nuopc for CTSMdev, using nuopc requires ESMF installation (at least 8.4.1)
MACH=pizdaint
QUEUE=normal # USER_REQUESTED_QUEUE, overrides default JOB_QUEUE
WALLTIME="06:00:00" # USER_REQUESTED_WALLTIME, overrides default JOB_WALLCLOCK_TIME
PROJ=$(basename "$(dirname "${PROJECT}")") # extract project name (sm61/s1207)
NTASKS=24 # will be nr of NODES (was 24)
let "NCORES = $NTASKS * 12" # this will be nr of CPUS

NSUBMIT=29 #5  #24 # partition into smaller chunks, excludes the first submission
STARTYEAR=1985
ENDYEAR=2014
STARTDATE=$STARTYEAR"-01-01"
NYEARS=1 #5

# Set directories
export CLMROOT=$PWD # CLM code base directory  where this script is located
export CASEDIR=$SCRATCH/cases/$CASENAME # case directory on scratch
export CESMDATAROOT=$SCRATCH/CCLM2_inputdata # inputdata directory on scratch 
export CESMOUTPUTROOT=$SCRATCH/archive/$CASENAME # output directory on scratch

# Log output (use "tee" to send output to both screen and $outfile)
logfile=$HOME/case_logs/${CASENAME}_mylogfile.log
mkdir -p "$(dirname "$logfile")" && touch "$logfile" # create parent/child directories and logfile
print_log() {
    output="$1"
    echo -e "${output}" | tee -a $logfile
}

print_log "*** Case at: ${CASEDIR} ***"
print_log "*** Case settings: compset ${COMPSET}, resolution ${RES}, domain ${DOMAIN}, compiler ${COMPILER} ***"
print_log "*** Logfile at: ${logfile} ***"

# Sync inputdata on scratch because scratch will be cleaned every month (change inputfiles on $PROJECT!)
print_log "\n*** Syncing inputdata on scratch  ***"
#rsync -av /project/$PROJ/shared/CCLM2_inputdata/ $CESMDATAROOT/ | tee -a $logfile # also check for updates in file content
sbatch --account=$PROJ --export=ALL,PROJ=$PROJ transfer_cesm_inputdata.sh # xfer job to prevent overflowing the loginnode


#==========================================
# Load modules and find spack packages
#==========================================

# Find spack_esmf installation  (used in .cime/config_machines.xml and env_build.xml) and path of netcdf files
if [ $DRIVER == nuopc ]; then
    print_log "\n *** Finding spack_esmf ***"
    export ESMF_PATH=$(spack location -i esmf@8.4.1) # e.g. /project/s1207/ivanderk/spack-install/cray-cnl7-haswell/gcc-9.3.0/esmf-8.4.1-esftqomee2sllfsmjevw3f7cet6tbeb4/
    print_log "ESMF at: ${ESMF_PATH}"

    # direct to spack installation of ESMF (also in .cime/config_compilers.xml - but doesn't work yet for Inne)
    #export ESMFMKFILE=${ESMF_PATH}/lib/esmf.mk
    print_log "*** ESMFMKFILE: ${ESMFMKFILE} ***"

    print_log "*** LD_LIBRARY_PATH: ${LD_LIBRARY_PATH} ***"

fi


#==========================================
# Create case
#==========================================

print_log "\n*** Creating CASE: ${CASENAME} ***"

cd $CLMROOT/cime/scripts
./create_newcase --case $CASEDIR --compset $COMPSET --res $RES --mach $MACH --compiler $COMPILER --driver $DRIVER --project $PROJ --run-unsupported | tee -a $logfile



#==========================================
# Configure CLM
# Settings will appear in namelists and have precedence over user_nl_xxx
#==========================================

print_log "\n*** Modifying env_*.xml  ***"
cd $CASEDIR

# Set directory structure
./xmlchange RUNDIR="$CASEDIR/run" # by defaut, RUNDIR is $SCRATCH/$CASENAME/run
./xmlchange EXEROOT="$CASEDIR/bld"

# Change job settings (env_batch.xml or env_workflow.xml). Do this here to change for both case.run and case.st_archive
./xmlchange JOB_QUEUE=$QUEUE --force
./xmlchange JOB_WALLCLOCK_TIME=$WALLTIME
./xmlchange PROJECT=$PROJ

# Set run start/stop options and DATM forcing (env_run.xml)
./xmlchange RUN_TYPE=startup
./xmlchange RESUBMIT=$NSUBMIT
./xmlchange RUN_STARTDATE=$STARTDATE
./xmlchange STOP_OPTION=nyears,STOP_N=$NYEARS
# ./xmlchange NCPL_BASE_PERIOD="day",ATM_NCPL=48 # coupling freq default 30min = day,48

if [ $CODE == CTSMdev ] && [ $DRIVER == nuopc ]; then
	./xmlchange DATM_YR_START=$STARTYEAR,DATM_YR_END=$ENDYEAR,DATM_YR_ALIGN=$STARTYEAR # new variable names in CTSMdev with nuopc drive
else
    ./xmlchange DATM_CLMNCEP_YR_START=2004,DATM_CLMNCEP_YR_END=2004,DATM_CLMNCEP_YR_ALIGN=2004 # in clm5.0 and CLM_features, with any driver
fi

# Set the number of cores and nodes (env_mach_pes.xml)
./xmlchange COST_PES=$NCORES
./xmlchange NTASKS_CPL=-$NTASKS
./xmlchange NTASKS_ATM=-$NTASKS
./xmlchange NTASKS_OCN=-$NTASKS
./xmlchange NTASKS_WAV=-$NTASKS
./xmlchange NTASKS_GLC=-$NTASKS
./xmlchange NTASKS_ICE=-$NTASKS
./xmlchange NTASKS_ROF=-$NTASKS
./xmlchange NTASKS_LND=-$NTASKS

# ESMF interface and time manager (env_build.xml)
./xmlchange --file env_build.xml --id COMP_INTERFACE --val $DRIVER # mct is default in clm5.0, nuopc is default in CTSMdev (requires ESMF installation); adding --driver mct to create_newcase adds everything needed

if [ $DRIVER == mct ]; then
    ./xmlchange --file env_build.xml --id USE_ESMF_LIB --val "FALSE" # FALSE is default in clm5.0; since cesm1_2 ESMF is no longer necessary to run with calendar=gregorian
elif [ $DRIVER == nuopc ]; then
    ./xmlchange --file env_build.xml --id USE_ESMF_LIB --val "TRUE" # using the ESMF library specified by env var ESMFMKFILE (config_machines.xml), or ESMF_LIBDIR (not found in env_build.xml)
fi


#==========================================
# Set up the case (creates user_nl_xxx)
#==========================================

print_log "\n*** Running case.setup ***"
./case.setup -r | tee -a $logfile


#=========================================SOILLI=
# User namelists (use cat >> to append)
# Surface data: domain-specific
# Paramfile: can be exchanged for newer versions
# Domainfile: has to be provided to DATM
#==========================================
print_log "\n*** Modifying user_nl_*.xml  ***"


# For global domain keep the defaults (downloaded from svn trunc)
if [ $DOMAIN == glob ]; then
cat >> user_nl_clm << EOF

finidat='/scratch/snx3000/ivanderk/processing_4p1000/initdata_4p1000/IHistClm51Sp.hcru_hcru_mt13.CTL.clm2.r.2035-01-01-00000.nc'

fsurdat = "/scratch/snx3000/ivanderk/processing_4p1000/surfdata_4p1000/surfdata_360x720cru_16pfts_Irrig_CMIP6_simyr2000_c170824_$SCEN.nc"
do_transient_crops=false
do_transient_pfts=false
flanduse_timeseries=''

irrigate=$FLAG_IRRIG

hist_fincl1 = "QROOTSINK", "watsat", "watfc", "EFF_POROSITY", "SOILPSI", "SOILLIQ"
hist_fincl2 = "QROOTSINK", "TOTSOILLIQ","TOTSOILICE", "SOILWATER_10CM",  "QINFL" , "QOVER", "QDRAI", "ZWT", "FH2OSFC",'EFLX_LH_TOT', 'QIRRIG_FROM_SURFACE', 'H2OSOI',"watsat", "watfc","EFF_POROSITY", "SOILPSI","QSOIL","QVEGE","QVEGT","TLAI","QHR","SMP", "TSA", "TG", "TSOI", "SOILLIQ"

! h0 with monthly output, h1 subgrid monthly output.
hist_nhtfrq = 0, 0
hist_mfilt = 1, 1
hist_dov2xy = .true., .false.
! indicate if averaging over plant functional types, columns or land units is wanted
hist_type1d_pertape = ' ', ' '
EOF
fi


#==========================================
# Build
#==========================================

print_log "\n*** Building case ***"
./case.build --clean-all | tee -a $logfile
./case.build | tee -a $logfile


print_log "\n*** Finished building new case in ${CASEDIR} ***"


#==========================================
# Check and download input data
#==========================================

print_log "\n*** Downloading missing inputdata (if needed) ***"
print_log "Consider transferring new data to PROJECT, e.g. rsync -av ${SCRATCH}/CCLM2_inputdata /project/${PROJ}/shared/CCLM2_inputdata"
./check_input_data --download


#==========================================
# Preview and submit job
#==========================================

print_log "\n*** Preview the run ***"
./preview_run | tee -a $logfile

print_log "\n*** Submitting job ***"
./case.submit -a "-C gpu" | tee -a $logfile

squeue --user=$USER | tee -a $logfile
#less CaseStatus

enddate=`date +'%Y-%m-%d %H:%M:%S'`
duration=$SECONDS
print_log "Started at: $startdate"
print_log "Finished at: $enddate"
print_log "Duration to create, setup, build, submit: $(($duration / 60)) min $(($duration % 60)) sec"

print_log "\n*** Check the job: squeue --user=${USER} ***"
print_log "*** Check the case: in ${CASEDIR}, run less CaseStatus ***"


#==========================================
# Copy final CaseStatus to logs
#==========================================

# Notes:
#env_case = model version, components, resolution, machine, compiler [do not modify]
#env_mach_pes = NTASKS, number of MPI tasks (or nodes if neg. values) [modify before setup]
#env_mach_specific = controls machine specific environment [modify before setup]
#env_build = component settings [modify before build]
#env_batch = batch job settings [modify any time]
#env_run = run settings incl runtype, coupling, pyhsics/sp/bgc and output [modify any time]
#env_workflow = wallt
