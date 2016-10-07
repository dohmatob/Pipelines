#!/bin/bash

set -e
. Examples/Scripts/SetUpHCPPipeline.sh
. global/scripts/log.shlib  # Logging related functions


## Sanitize command-line args
if [ $# -lt 5 ]; then
    echo "Insufficient arguments!"
    echo "Usage:"
    echo "bash $0 <StudyFolder> <OutputFolder> <SubjectId> <TaskName> <UnwarpDir>"
    echo""
    echo "Example:"
    echo "$0 /home/elvis/nilearn_data/drago/storage/data/HCP/S500-1 ~/mnt/32-bit-system/home/elvis/hcp_preproc 100307 MOTOR RL"
    exit 1
fi


## Misc
: ${T1BASED:=0}

# grab input args
StudyFolder=$1
OutputFolder=$2
SubjectId=$3
TaskName=$4
Dir=$5
if [[ ${Dir} == "RL" ]]; then
    UnwarpDir="x"
else
    UnwarpDir="-x"
fi

SubjectFolder=${StudyFolder}/${SubjectId}
fMRIName=tfMRI_${TaskName}_${Dir}
fMRIOutputFolder=${OutputFolder}/${SubjectId}/${fMRIName}
T1wOutputFolder=${OutputFolder}/${SubjectId}/T1w
fMRITimeSeries=${SubjectFolder}/unprocessed/3T/${fMRIName}/${SubjectId}_3T_${fMRIName}
fMRISBRef=${SubjectFolder}/unprocessed/3T/${fMRIName}/${SubjectId}_3T_${fMRIName}_SBRef
SpinEchoPhaseEncodePositive=${SubjectFolder}/unprocessed/3T/${fMRIName}/${SubjectId}_3T_SpinEchoFieldMap_RL
SpinEchoPhaseEncodeNegative=${SubjectFolder}/unprocessed/3T/${fMRIName}/${SubjectId}_3T_SpinEchoFieldMap_LR

# nasty file organization to please HCP logic
ScoutInputName=${fMRISBRef}
ScoutInputBaseName=`basename ${ScoutInputName}`
T1wFolder=${SubjectFolder}/T1w # Location of T1w images
T1wImage=${T1wFolder}/T1w_acpc_dc
T1wRestoreImage=${T1wFolder}/T1w_acpc_dc_restore_brain

# make sure all input files are present
for input_file in ${fMRITimeSeries} ${fMRISBRef} ${T1wImage} ${T1wRestoreImage}; do
    if [ ! -e ${input_file}.nii.gz ]; then
	echo "Input file ${input_file}.nii.gz present. Qutting..."
	exit 1
    fi
done    

log_Msg "mkdir ${fMRIOutputFolder}"
mkdir -p "$fMRIOutputFolder"

## Motion Correction
MovementRegressor="Movement_Regressors" #No extension, .txt appended
MotionMatrixFolder="MotionMatrices"
MotionMatrixPrefix="MAT_"
mkdir -p ${fMRIOutputFolder}/MotionCorrection
if [ ! -e ${fMRIOutputFolder}/MotionCorrection/${fMRIName}_mc.par ]; then
    fMRIVolume/scripts/MotionCorrection.sh \
	${fMRIOutputFolder}/MotionCorrection \
	${fMRITimeSeries} \
	${fMRISBRef} \
	${fMRIOutputFolder}/${fMRIName}_mc \
	${fMRIOutputFolder}/${MovementRegressor} \
	${fMRIOutputFolder}/${MotionMatrixFolder} \
	${MotionMatrixPrefix} \
	"MCFLIRT"
fi


## Topup distortion correction using a blip-reversed SE pair "fieldmap" sequence
TopupConfig="global/config/b02b0.cnf"
DwellTime="0.00058"
UseJacobian="true"
GradientDistortionCoeffs="NONE" # Set to NONE to skip gradient distcorr
DistCorrWD=${fMRIOutputFolder}/DistortionCorrection
mkdir -p ${DistCorrWD}
if [ ! -e ${DistCorrWD}/WarpField.nii.gz ]; then
    global/scripts/TopupPreprocessingAll.sh \
        --workingdir=${DistCorrWD}/FieldMap \
        --phaseone=${SpinEchoPhaseEncodeNegative} \
        --phasetwo=${SpinEchoPhaseEncodePositive} \
        --scoutin=${ScoutInputName} \
        --echospacing=${DwellTime} \
        --unwarpdir=${UnwarpDir} \
        --owarp=${DistCorrWD}/WarpField \
        --ojacobian=${DistCorrWD}/Jacobian \
        --gdcoeffs=${GradientDistortionCoeffs} \
        --topupconfig=${TopupConfig} \
        --usejacobian=${UseJacobian}
fi

# create a spline interpolated image of scout (distortion corrected in same space)
log_Msg "create a spline interpolated image of scout (distortion corrected in same space)"
${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${ScoutInputName} -w ${DistCorrWD}/WarpField -o ${DistCorrWD}/${ScoutInputBaseName}_undistorted


## Coregistration: Use epi_ref (freesurfer-less BBR) to register EPI to T1w space
dof=6

# register undistorted scout image to T1w. This is just an initial registration, refined later in
# this script, but it is actually pretty good
log_Msg "register undistorted scout image to T1w"
if [ ! -e ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.mat ]; then
    global/scripts/epi_reg_dof --dof=${dof} --epi=${DistCorrWD}/${ScoutInputBaseName}_undistorted --t1=${T1wImage} --t1brain=${T1wRestoreImage} --out=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init
fi

# copy the initial registration into the final affine's filename, as it is pretty good
# we need something to get between the spaces to compute an initial bias field
cp ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.mat ${DistCorrWD}/fMRI2str.mat

# generate combined warpfields and spline interpolated images + apply bias field correction
log_Msg "generate combined warpfields and spline interpolated images and apply bias field correction"
if [ ! -e ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_warp.nii.gz ]; then
    ${FSLDIR}/bin/convertwarp --relout --rel -r ${T1wImage} --warp1=${DistCorrWD}/WarpField.nii.gz --postmat=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.mat -o ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_warp
fi


## My 1-step resampling procedure
# Save TR for later
TR_vol=`${FSLDIR}/bin/fslval ${fMRITimeSeries} pixdim4 | cut -d " " -f 1`
NumFrames=`${FSLDIR}/bin/fslval ${fMRITimeSeries} dim4`

# Apply combined transformations to fMRI (combines gradient non-linearity distortion, motion correction, and registration to T1w space, but keeping fMRI resolution)
for T1BASED in 1 0; do
    OutputfMRI=${fMRIOutputFolder}/tfMRI_${TaskName}_${Dir}_undistorted_mc
    if [ ${T1BASED} == "1" ]; then
	OutputfMRI=${OutputfMRI}2T1
    fi
    if [ ! -e ${OutputfMRI}.nii.gz ]; then
	mkdir -p ${DistCorrWD}/prevols
	mkdir -p ${DistCorrWD}/postvols    
	${FSLDIR}/bin/fslsplit ${fMRITimeSeries} ${DistCorrWD}/prevols/vol -t
	FrameMergeSTRING="" 
	k=0
	while [ $k -lt $NumFrames ] ; do
	    vnum=`${FSLDIR}/bin/zeropad $k 4`
	    
	    if [ ${T1BASED} == "0" ]; then
		flirt -in ${DistCorrWD}/prevols/vol${vnum} -ref ${T1wImage} \
		      -out ${DistCorrWD}/postvols/vol${k} \
		      -init ${fMRIOutputFolder}/${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} \
		      -applyxfm
	    else
		# combine
		${FSLDIR}/bin/convertwarp --relout --rel --ref=${DistCorrWD}/prevols/vol${vnum} \
			 --warp1=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_warp \
			 --postmat=${fMRIOutputFolder}/${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum} \
			 --out=${fMRIOutputFolder}/${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp
		
		# apply
		${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DistCorrWD}/prevols/vol${vnum} \
			 --warp=${fMRIOutputFolder}/${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_gdc_warp \
			 --ref=${T1wImage} --out=${DistCorrWD}/postvols/vol${k}
	    fi
	    echo ${DistCorrWD}/postvols/vol${k}.nii.gz
	    FrameMergeSTRING="${FrameMergeSTRING}${DistCorrWD}/postvols/vol${k}.nii.gz " 
	    k=`echo "$k + 1" | bc`
	done

	# Merge together results and restore the TR (saved beforehand)
	${FSLDIR}/bin/fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
    fi
    echo "${OutputfMRI}.nii.gz"
done
