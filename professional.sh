#!/bin/bash

# Example usage:
# bash professional.sh /home/elvis/nilearn_data/drago/storage/data/HCP/S500-1 ~/mnt/32-bit-system/home/elvis/hcp_preproc 100307 MOTOR RL
#

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
    echo "bash professional.sh /home/elvis/nilearn_data/drago/storage/data/HCP/S500-1 ~/mnt/32-bit-system/home/elvis/hcp_preproc 100307 MOTOR RL"
    exit 1
fi


## Misc
# grav input args
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

T1wImage=${StudyFolder}/${SubjectId}/T1w/T1w_acpc_dc
T1wRestoreImage=${StudyFolder}/${SubjectId}/T1w/T1w_acpc_dc_restore_brain

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
ScoutName="Scout"
OrigScoutName="${ScoutName}_orig"

ResultsFolder="Results"
T1wFolder=${SubjectFolder}/T1w #Location of T1w images

# files produced by HCP's T1 ==> std step. we can replace this with another method
AtlasSpaceFolder=${SubjectFolder}/MNINonLinear
AtlasTransform="acpc_dc2standard"
BiasField="BiasField_acpc_dc"
BiasFieldMNI="BiasField"

# files produced by HCP's freesurfer recon-all. we could replace this with an output from FSL's BET, e.g
FreeSurferBrainMask="brainmask_fs"

ResultsFolder=${AtlasSpaceFolder}/${ResultsFolder}/${fMRIName}
OutputfMRI2StandardTransform=${fMRIName}2standard
Standard2OutputfMRITransform=standard2${fMRIName}
fMRI2strOutputTransform=${fMRIName}2str
T1wAtlasName="T1w_restore"
JacobianOut="Jacobian"
FinalfMRIResolution=2
# note, this file doesn't exist yet, gets created by ComputeSpinEchoBiasField.sh during
# DistortionCorrectionAnd...
sebasedBiasFieldMNI=$AtlasSpaceFolder/Results/$fMRIName/${fMRIName}_sebased_bias.nii.gz
UseBiasFieldMNI="$sebasedBiasFieldMNI"

## T1 segmentation, BET, bias correction, and registration to MNI space
log_Msg "Affine T1w ==> MNI"
T1w2MNIMat=${T1wOutputFolder}/`basename ${T1wRestoreImage}`2std.mat
if [ ! -e ${T1w2MNIMat} ]; then
    ${FSLDIR}/bin/flirt -in ${T1wRestoreImage} -ref ${FSLDIR}/data/standard/MNI152_T1_1mm_brain -omat ${T1w2MNIMat} -v
fi

log_Msg "Nonlinear T1w ==> MNI"
mkdir -p ${T1wOutputFolder}/xfms
T1w2MNIWarp=${T1wOutputFolder}/xfms/${AtlasTransform}
T1w2MNIJac=${T1w2MNIWarp}_jac
if [ ! -e ${T1w2MNIWarp}.nii.gz ] || [ ! -e ${T1w2MNIJac}.nii.gz ]; then
    ${FSLDIR}/bin/fnirt --ref=${FSLDIR}/data/standard/MNI152_T1_1mm_brain \
	     --in=${T1wRestoreImage} --aff=${T1w2MNIMat} \
	     --cout=${T1w2MNIWarp} --iout=titi --interp=spline \
	     --jout=${T1w2MNIJac} -v
fi

log_Msg "Warping bias field"
if [ ! -e ${T1wOutputFolder}/${BiasFieldMNI}.nii.gz ]; then
    applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_1mm_brain \
	      --in=${T1wFolder}/${BiasField} --interp=spline \
	      --out=${T1wOutputFolder}/${BiasFieldMNI} -w ${T1w2MNIWarp} -v
fi
UseBiasFieldMNI=${T1wOutputFolder}/${BiasFieldMNI}

log_Msg "Warping brainmask"
if [ ! -e ${T1wOutputFolder}/${FreeSurferBrainMask}.nii.gz ]; then
    applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_2mm_brain \
	      --in=${T1wFolder}/${FreeSurferBrainMask} --interp=spline \
	      --out=${T1wOutputFolder}/${FreeSurferBrainMask} -w ${T1w2MNIWarp} -v
fi
AtlasSpaceFolder=${T1wOutputFolder}

log_Msg "mkdir ${fMRIOutputFolder}"
mkdir -p "$fMRIOutputFolder"


## Fake gradient distortion correction
log_Msg "NOT PERFORMING GRADIENT DISTORTION CORRECTION"
${FSLDIR}/bin/imcp ${ScoutInputName}.nii.gz "$fMRIOutputFolder"/"$OrigScoutName".nii.gz
${FSLDIR}/bin/imcp ${fMRITimeSeries} "$fMRIOutputFolder"/"$fMRIName"_gdc
${FSLDIR}/bin/fslroi "$fMRIOutputFolder"/"$fMRIName"_gdc "$fMRIOutputFolder"/"$fMRIName"_gdc_warp 0 3
${FSLDIR}/bin/fslmaths "$fMRIOutputFolder"/"$fMRIName"_gdc_warp -mul 0 "$fMRIOutputFolder"/"$fMRIName"_gdc_warp
${FSLDIR}/bin/imcp "$fMRIOutputFolder"/"$OrigScoutName" "$fMRIOutputFolder"/"$ScoutName"_gdc

# make fake jacobians of all 1s, for completeness
${FSLDIR}/bin/fslmaths "$fMRIOutputFolder"/"$OrigScoutName" -mul 0 -add 1 "$fMRIOutputFolder"/"$ScoutName"_gdc_warp_jacobian
${FSLDIR}/bin/fslroi "$fMRIOutputFolder"/"$fMRIName"_gdc_warp "$fMRIOutputFolder"/"$fMRIName"_gdc_warp_jacobian 0 1
${FSLDIR}/bin/fslmaths "$fMRIOutputFolder"/"$fMRIName"_gdc_warp_jacobian -mul 0 -add 1 "$fMRIOutputFolder"/"$fMRIName"_gdc_warp_jacobian


## Motion Correction
MotionCorrectionType="MCFLIRT"
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
	${MotionCorrectionType}
fi


## Topup distortion correction using a blip-reversed SE pair "fieldmap"
# sequence
TopupConfig="global/config/b02b0.cnf"
DwellTime="0.00058"
GradientDistortionCoeffs="NONE" # Set to NONE to skip gradient distcorr
UseJacobian="false"
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

# apply Jacobian correction to scout image (optional)
# gdc jacobian is already applied in main script, where the gdc call for the scout is
if [[ $UseJacobian == "true" ]]; then
    log_Msg "apply Jacobian correction to scout image"
    ${FSLDIR}/bin/fslmaths ${DistCorrWD}/${ScoutInputBaseName}_undistorted -mul ${DistCorrWD}/Jacobian ${DistCorrWD}/${ScoutInputBaseName}_undistorted
fi


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
if [ ! -e ${DistCorrWD}/Jacobian2T1w.nii.gz ]; then
    ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${DistCorrWD}/Jacobian.nii.gz -r ${T1wImage} --premat=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.mat -o ${DistCorrWD}/Jacobian2T1w.nii.gz
fi

# 1-step resample from input (gdc) scout - NOTE: no longer includes jacobian correction, if
# specified
if [ ! -e  ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.nii.gz ]; then
    ${FSLDIR}/bin/applywarp --rel --interp=spline -i ${ScoutInputName} -r ${T1wImage} -w ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_warp -o ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init
fi

# # resample phase images to T1w space these files were obtained by the import script from the
# # FieldMap directory, save them into the package and resample them
# # we don't have the final transform to actual T1w space yet, that occurs later in this script
# # but, we need the T1w segmentation to make the bias field, so use the initial registration above,
# # then compute the bias field again at the end
# Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc SBRef_dc"
# for File in ${Files}
# do
#     # NOTE: this relies on TopupPreprocessingAll generating _jac versions of the files
#     if [[ $UseJacobian == "true" ]]
#     then
#         ${FSLDIR}/bin/applywarp --interp=spline -i "${DistCorrWD}/FieldMap/${File}_jac" -r ${SubjectFolder}/T1w/T2w_acpc_dc.nii.gz --premat=${DistCorrWD}/fMRI2str.mat -o ${DistCorrWD}/${File}
#     else
#         ${FSLDIR}/bin/applywarp --interp=spline -i "${DistCorrWD}/FieldMap/${File}" -r ${SubjectFolder}/T1w/T2w_acpc_dc.nii.gz --premat=${DistCorrWD}/fMRI2str.mat -o ${DistCorrWD}/${File}
#     fi
# done

# bias correction
UseBiasField=${DistCorrWD}/ComputeSpinEchoBiasField/${fMRIName}_sebased_bias.nii.gz
log_Msg "compute spin-echo biasfield"
mkdir -p $DistCorrWD/ComputeSpinEchoBiasField
if [ ! -e ${UseBiasField} ]; then
    fMRIVolume/scripts/ComputeSpinEchoBiasField.sh \
	--workingdir=$DistCorrWD/ComputeSpinEchoBiasField \
	--subjectfolder=$SubjectFolder \
	--fmriname=$fMRIName \
	--corticallut=global/config/FreeSurferCorticalLabelTableLut.txt \
	--subcorticallut=global/config/FreeSurferSubcorticalLabelTableLut.txt \
	--smoothingfwhm=2 \
	--inputdir=$DistCorrWD
fi

# apply bias correction to scout images
Files="PhaseOne_gdc_dc PhaseTwo_gdc_dc"
for File in ${Files}
do
    # we need to apply the new bias field to them for output
    ${FSLDIR}/bin/fslmaths ${DistCorrWD}/${File} -div "$UseBiasField" ${DistCorrWD}/${File}_unbias
    
    #don't need the T1w versions
    #${FSLDIR}/bin/imcp ${DistCorrWD}/${File}_unbias ${SubjectFolder}/T1w/Results/${fMRIName}/${fMRIName}_${File}
done

# copy bias field and dropouts, etc to results dir
# ${FSLDIR}/bin/imcp $DistCorrWD/ComputeSpinEchoBiasField/${fMRIName}_sebased_bias ${UseBiasFieldMNI}

# apply Jacobian correction and bias correction options to scout image
if [[ $UseJacobian == "true" ]] ; then
    log_Msg "apply Jacobian correction to scout image"
    if [[ "$UseBiasField" != "" ]]
    then
        ${FSLDIR}/bin/fslmaths ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init -div ${UseBiasField} -mul ${DistCorrWD}/Jacobian2T1w.nii.gz ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init.nii_unbias.gz
    else
        ${FSLDIR}/bin/fslmaths ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init -mul ${DistCorrWD}/Jacobian2T1w.nii.gz ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_unbias.nii.gz
    fi
else
    log_Msg "do not apply Jacobian correction to scout image"
    if [[ "$UseBiasField" != "" ]]
    then
        ${FSLDIR}/bin/fslmaths ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init -div ${UseBiasField} ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_unbias.nii.gz
    fi
    # these all overwrite the input, no 'else' needed for "do nothing"
fi
${FSLDIR}/bin/imcp ${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_warp.nii.gz ${DistCorrWD}/fMRI2str.nii.gz

# move some files around
OutputTransformDir=${fMRIOutputFolder}/xfms
log_Msg "mkdir -p ${OutputTransformDir}"
mkdir -p ${OutputTransformDir}
log_Msg "Copying EPI ==> T1 output transforms"
${FSLDIR}/bin/imcp ${DistCorrWD}/Jacobian2T1w.nii.gz ${fMRIOutputFolder}/${JacobianOut}.nii.gz

# QA image (sqrt of EPI * T1w)
QAImage=qa
log_Msg 'generating QA image (sqrt of EPI * T1w)'
RegOutput=${DistCorrWD}/${ScoutInputBaseName}_undistorted2T1w_init_unbias
${FSLDIR}/bin/fslmaths ${T1wRestoreImage}.nii.gz -mul ${RegOutput}.nii.gz -sqrt ${QAImage}.nii.gz
log_Msg "END"
echo " END: `date`" >> ${DistCorrWD}/log.txt


## One Step Resampling
log_Msg "One Step Resampling"
log_Msg "mkdir -p ${fMRIOutputFolder}/OneStepResampling"
mkdir -p ${fMRIOutputFolder}/OneStepResampling
if [ ! -e ${fMRIOutputFolder}/${fMRIName}_nonlin.nii.gz ]; then
    fMRIVolume/scripts/OneStepResampling.sh \
	--workingdir=${fMRIOutputFolder}/OneStepResampling \
	--infmri=${fMRITimeSeries}.nii.gz \
	--t1=titi \ #${AtlasSpaceFolder}/${T1wAtlasName} \
	--fmriresout=${FinalfMRIResolution} \
	--fmrifolder=${fMRIOutputFolder} \
	--fmri2structin=${fMRIOutputFolder}/xfms/${fMRI2strOutputTransform} \
	--struct2std=${AtlasSpaceFolder}/xfms/${AtlasTransform} \
	--owarp=${fMRIOutputFolder}/xfms/${OutputfMRI2StandardTransform} \
	--oiwarp=${fMRIOutputFolder}/xfms/${Standard2OutputfMRITransform} \
	--motionmatdir=${fMRIOutputFolder}/${MotionMatrixFolder} \
	--motionmatprefix=${MotionMatrixPrefix} \
	--ofmri=${fMRIOutputFolder}/${fMRIName}_nonlin \
	--freesurferbrainmask=${AtlasSpaceFolder}/${FreeSurferBrainMask} \
	--biasfield=${AtlasSpaceFolder}/${BiasFieldMNI} \
	--gdfield=${fMRIOutputFolder}/${fMRIName}_gdc_warp \
	--scoutin=${fMRIOutputFolder}/${OrigScoutName} \
	--scoutgdcin=${fMRIOutputFolder}/${ScoutName}_gdc \
	--oscout=${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin \
	--ojacobian=${fMRIOutputFolder}/${JacobianOut}_MNI.${FinalfMRIResolution}
fi

# # create MNI space corrected fieldmap images
# log_Msg "create MNI space corrected fieldmap images"
# ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DistCorrWD}/PhaseOne_gdc_dc_unbias -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -r ${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin -o ${ResultsFolder}/${fMRIName}_PhaseOne_gdc_dc
# ${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DistCorrWD}/PhaseTwo_gdc_dc_unbias -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -r ${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin -o ${ResultsFolder}/${fMRIName}_PhaseTwo_gdc_dc

# # create MNINonLinear final fMRI resolution bias field outputs
# log_Msg "create MNINonLinear final fMRI resolution bias field outputs"
# if [[ ${BiasCorrection} == "SEBASED" ]]
# then
#     ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DistCorrWD}/ComputeSpinEchoBiasField/sebased_bias_dil.nii.gz -r ${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin -w ${SubjectFolder}/MNINonLinear/xfms/acpc_dc2standard.nii.gz -o ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_bias.nii.gz
#     ${FSLDIR}/bin/fslmaths ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_bias.nii.gz -mas ${fMRIOutputFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_bias.nii.gz
    
#     ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DistCorrWD}/ComputeSpinEchoBiasField/sebased_reference_dil.nii.gz -r ${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin -w ${SubjectFolder}/MNINonLinear/xfms/acpc_dc2standard.nii.gz -o ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_reference.nii.gz
#     ${FSLDIR}/bin/fslmaths ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_reference.nii.gz -mas ${fMRIOutputFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_sebased_reference.nii.gz
        
#     ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DistCorrWD}/ComputeSpinEchoBiasField/${fMRIName}_dropouts.nii.gz -r ${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin -w ${SubjectFolder}/MNINonLinear/xfms/acpc_dc2standard.nii.gz -o ${SubjectFolder}/MNINonLinear/Results/${fMRIName}/${fMRIName}_dropouts.nii.gz
# fi


## Intensity Normalization and Bias Removal
log_Msg "Intensity Normalization and Bias Removal"
echo "fMRIVolume/scripts/IntensityNormalization.sh \
    --infmri=${fMRIOutputFolder}/${fMRIName}_nonlin \
    --biasfield=${UseBiasFieldMNI} \
    --jacobian=${T1w2MNIJac} \ # ${fMRIOutputFolder}/${JacobianOut}_MNI.${FinalfMRIResolution} \
    --brainmask=${T1wOutputFolder}/${FreeSurferBrainMask} \ #${fMRIOutputFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution} \
    --ofmri=${fMRIOutputFolder}/${fMRIName}_nonlin_norm \
    --inscout=${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin \
    --oscout=${fMRIOutputFolder}/${fMRIName}_SBRef_nonlin_norm \
    --usejacobian=${UseJacobian}"
