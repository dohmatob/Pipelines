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
    echo 'for task in MOTOR LANGUAGE EMOTION RELATIONAL SOCIAL GAMBLING WM; do for x in /storage/data/HCP/S500-1/*/unprocessed/3T/tfMRI_MOTOR_LR/*_3T_tfMRI_MOTOR_LR.nii.gz; do sid=$(basename $(dirname $(dirname $(dirname $(dirname $x))))); echo $sid $task; done; done | xargs -n 2 -P 32 -i bash -c './intra.sh /storage/data/HCP/S500-1 /storage/workspace/elvis/hcp_new_preproc $@ LR' _ {}'
    echo "$0 /home/elvis/nilearn_data/drago/storage/data/HCP/S500-1 ~/mnt/32-bit-system/home/elvis/hcp_preproc 100307 MOTOR RL"
    exit 1
fi

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
FinalResolution=1.5  # downsample t1 images to this res; set to "orig" if you prefer original res

# make sure all input files are present
for input_file in ${fMRITimeSeries} ${fMRISBRef} ${T1wImage} ${T1wRestoreImage}; do
    if [ ! -e ${input_file}.nii.gz ]; then
	echo "Input file ${input_file}.nii.gz present. Qutting..."
	exit 1
    fi
done    

log_Msg "make output dirs"
mkdir -p "$fMRIOutputFolder"
mkdir -p "$T1wOutputFolder"

## Motion Correction
MovementRegressor=${fMRIOutputFolder}/"Movement_Regressors" #No extension, .txt appended
MotionMatrixFolder=${fMRIOutputFolder}/"MotionMatrices"
MotionMatrixPrefix="MAT_"
mkdir -p ${fMRIOutputFolder}/MotionCorrection
if [ ! -e ${fMRIOutputFolder}/MotionCorrection/${fMRIName}_mc.par ]; then
    fMRIVolume/scripts/MotionCorrection.sh \
	${fMRIOutputFolder}/MotionCorrection \
	${fMRITimeSeries} \
	${fMRISBRef} \
	${fMRIOutputFolder}/${fMRIName}_mc \
	${MovementRegressor} \
	${MotionMatrixFolder} \
	${MotionMatrixPrefix} \
	"MCFLIRT"
else
    log_Msg "motion correction"
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
else
    log_Msg "topup fieldmap generation and gradient unwarping"
fi

# create a spline interpolated image of scout (distortion corrected in same space)
log_Msg "create a spline interpolated image of scout (distortion corrected in same space)"
if [ ! -e ${DistCorrWD}/${ScoutInputBaseName}_undistorted.nii.gz ]; then
    ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${ScoutInputName} \
	-w ${DistCorrWD}/WarpField -o ${DistCorrWD}/${ScoutInputBaseName}_undistorted
fi
${FSLDIR}/bin/imcp ${DistCorrWD}/${ScoutInputBaseName}_undistorted ${fMRIOutputFolder}


## Coregistration: Use epi_ref (freesurfer-less BBR) to register EPI to T1w space
dof=6

# downsample t1 images to for speed
if [ ! ${FinalResolution} == "orig" ]; then
    for t1 in ${T1wImage} ${T1wRestoreImage}; do
	log_Msg "resample ${t1}.nii.gz to ${FinalResolution}mm isotropic"
	resampled=${T1wOutputFolder}/`basename ${t1}`${FinalResolution}mm
	if [ ! -e ${resampled}.nii.gz ]; then
	    ${FSLDIR}/bin/flirt -in ${t1} -ref ${t1} -applyisoxfm ${FinalResolution} -noresampblur -o ${resampled}
	fi
    done
    T1wImageResampled=${T1wOutputFolder}/`basename ${T1wImage}`${FinalResolution}mm
    T1wRestoreImageResampled=${T1wOutputFolder}/`basename ${T1wRestoreImage}`${FinalResolution}mm
else
    T1wImageResampled=${T1wImage}
    T1wRestoreImageResampled=${T1wRestoreImage}
fi

# register undistorted scout image to T1w. This is just an initial registration, 
# refined later in this script, but it is actually pretty good
CoregSBRef=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w${FinalResolution}mm_init
SBRef2T1wTransformMat=${CoregSBRef}.mat
log_Msg "Coregister undistorted scout image to (resampled) t1 space"
if [ ! -e ${SBRef2T1wTransformMat} ]; then
    global/scripts/epi_reg_dof --dof=${dof} --epi=${DistCorrWD}/${ScoutInputBaseName}_undistorted --t1=${T1wImageResampled} \
	--t1brain=${T1wRestoreImageResampled} --out=${CoregSBRef}
fi


## My 1-step resampling procedure
# Save TR for later
TR_vol=`${FSLDIR}/bin/fslval ${fMRITimeSeries} pixdim4 | cut -d " " -f 1`
NumFrames=`${FSLDIR}/bin/fslval ${fMRITimeSeries} dim4`

# Apply combined transformations to fMRI (combines gradient non-linearity distortion,
# motion correction, and registration to T1w space, but keeping fMRI resolution)
log_Msg "1-step resampling of fMRI"
for t1_based in 0 1; do
    OutputfMRI=${fMRIOutputFolder}/tfMRI_${TaskName}_${Dir}_undistorted_mc
    if [ "${t1_based}" == "1" ]; then
	OutputfMRI=${OutputfMRI}2T1w${FinalResolution}mm
    fi
    if [ ! -e ${OutputfMRI}.nii.gz ]; then
	mkdir -p ${DistCorrWD}/prevols
	mkdir -p ${DistCorrWD}/postvols    
	${FSLDIR}/bin/fslsplit ${fMRITimeSeries} ${DistCorrWD}/prevols/vol -t
	FrameMergeSTRING="" 
	k=0
	while [ $k -lt $NumFrames ] ; do
	    vnum=`${FSLDIR}/bin/zeropad $k 4`

	    # combine transformations from motion correction and BBR
	    if [ "${t1_based}" == "1" ]; then
		Affine=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}_all_affine.mat
		Ref=${T1wImageResampled}
		${FSLDIR}/bin/convert_xfm -omat ${Affine} \
		    -concat ${SBRef2T1wTransformMat} ${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}
	    else
		Affine=${MotionMatrixFolder}/${MotionMatrixPrefix}${vnum}
		Ref=${DistCorrWD}/${ScoutInputBaseName}_undistorted
	    fi

	    # apply combined warp fields
	    ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DistCorrWD}/prevols/vol${vnum} \
		--warp=${DistCorrWD}/WarpField --ref=${Ref} --out=${DistCorrWD}/postvols/vol${k} \
		--postmat=${Affine}
		
	    echo ${DistCorrWD}/postvols/vol${k}.nii.gz
	    FrameMergeSTRING="${FrameMergeSTRING}${DistCorrWD}/postvols/vol${k}.nii.gz " 
	    k=`echo "$k + 1" | bc`
	done
	
	# Merge together results and restore the TR (saved beforehand)
	${FSLDIR}/bin/fslmerge -tr ${OutputfMRI} $FrameMergeSTRING $TR_vol
    fi
    echo "${OutputfMRI}.nii.gz"
done
