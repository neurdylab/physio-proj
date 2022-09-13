#!/bin/bash
# Code for pre-processing NKI fMRI data
# on Vanderbilt cluster (ACCRE)
#
#

# "main" directory with raw data
maindir_raw=/data1/neurdylab/datasets/nki_rockland/raw

# "main" directory that will contain our processed data
maindir_proc=/data1/neurdylab/datasets/nki_rockland/proc

# paths on platypus
MNI_T1_2mm_template="/data1/neurdylab/MNI152_T1_2mm_brain.nii.gz"
scripts_path="/data1/neurdylab/scripts/vu_meica_pipeline"
afni_init="singularity exec --bind /data1:/data1 ${scripts_path}/afni_cmake_build_AFNI_21.1.03.sif" 

# module load for FSL
module load GCC/8.2.0  OpenMPI/3.1.4
module load FSL/6.0.1-Python-3.7.2

# for ants:
pyscripts_path="/data1/neurdylab/eegfmri_vu_pipeline/scripts/ants_python"

# start loop over subjects

for sub_id in `cat /data1/neurdylab/datasets/nki_rockland/scripts/batch_rerun.txt`

do
echo $sub_id

# dirs where subject's processed fmri files will be written:
proc_func_dir=${maindir_proc}/${sub_id}/ses-BAS1/func
proc_anat_dir=${maindir_proc}/${sub_id}/ses-BAS1/anat
mkdir -p $proc_func_dir $proc_anat_dir

# subject's raw functional & anatomic dirs:
raw_func_dir=${maindir_raw}/${sub_id}/ses-BAS1/func
raw_anat_dir=${maindir_raw}/${sub_id}/ses-BAS1/anat
anat_img=${raw_anat_dir}/${sub_id}_ses-BAS1_T1w.nii.gz

func_prefix=${sub_id}_ses-BAS1_task-rest_acq-1400_bold

log_path="${proc_func_dir}/logs"
mkdir -p $log_path
log_file="${log_path}/${sub_id}.log"
printf "Starting processing for ${sub_id}\n\n" > ${log_file}
echo "Log file created at ${log_file}"


# anatomic (T1) image processing
# ------------------------------

if ! [ -e ${proc_anat_dir}/T1_bet.nii.gz ]; then

        echo "Bet and unifize T1 anatomical, done once per subject..."

        # Brain Extraction (bet)
        cmd="bet ${anat_img} ${proc_anat_dir}/T1_bet.nii.gz -f 0.4 -o -m"
        printf "Running ${cmd}\n" >> $log_file
        $cmd
fi

T1_bet_img=${proc_anat_dir}/T1_bet.nii.gz

if ! [ -e ${proc_anat_dir}/T1_unifize.nii.gz ]; then

        # unifize
        ${afni_init} 3dcopy ${T1_bet_img} ${proc_anat_dir}/T1_bet
        ${afni_init} 3dUnifize -prefix ${proc_anat_dir}/T1_unifize ${proc_anat_dir}/T1_bet+orig
        ${afni_init} 3dcopy ${proc_anat_dir}/T1_unifize+orig ${proc_anat_dir}/T1_unifize.nii.gz

        rm -f ${proc_anat_dir}/T1_bet+orig*
        rm -f ${proc_anat_dir}/T1_unifize+orig*

fi

# Motion correction (re-alignment)
#----------------------------------

if ! [ -e ${proc_func_dir}/${func_prefix}_mo.nii.gz ]; then
        echo "mcflirt done once per subject"

        cmd="mcflirt -plots -in ${raw_func_dir}/${func_prefix}.nii.gz \
        -out ${proc_func_dir}/${func_prefix}_mo.nii.gz"
        printf "Running ${cmd}\n" >> $log_file
        echo $cmd
        $cmd

fi

# define T1_mcflirt.nii.gz filename, for convenience
mcflirt_img=${proc_func_dir}/${func_prefix}_mo


# Extract single BOLD volume (reference for registration)
if ! [ -e ${proc_func_dir}/oneVol.nii.gz ]; then

        echo "fslroi done once per subject"

        cmd="fslroi ${mcflirt_img}.nii.gz ${proc_func_dir}/oneVol.nii.gz 0 1"

        printf "Running ${cmd}\n" >> $log_file
        echo $cmd

        $cmd

fi

# register BOLD volume to same subject's T1

if ! [ -e ${proc_func_dir}/oneVol_EPI2T1.nii.gz ]; then

 echo "epi_reg done once per subject"

        cmd="epi_reg --epi=${proc_func_dir}/oneVol.nii.gz \
        --t1=${anat_img} \
        --t1brain=${T1_bet_img} \
        --out=${proc_func_dir}/oneVol_EPI2T1.nii.gz"
 echo $cmd
	printf "Running ${cmd}\n" >> $log_file
        $cmd

fi


# ANTS VERSION
# ------------------------------
if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz ]; then

echo "ants reg..."
mkdir -p $proc_func_dir/ants_out

applywarp --ref=${proc_anat_dir}/T1_unifize.nii.gz --in=${proc_func_dir}/${func_prefix}_mo.nii.gz --out=${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz \
--postmat=${proc_func_dir}/oneVol_EPI2T1.mat

fi

# module change for ANTs
module load GCC/6.4.0-2.28 OpenMPI/2.1.1
module load ANTs/2.3.0-Python-2.7.14


if ! [ -e ${proc_anat_dir}/ants_reg_tforms/tform1_0GenericAffine.mat ]; then

    echo "running ANTS!!"
    mkdir -p ${proc_anat_dir}/ants_reg_tforms

    # T1 to MNI
${pyscripts_path}/ants_doReg_accre.py \
                 --in_nii ${proc_anat_dir}/T1_unifize.nii.gz \
                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
                 --out_dir ${proc_anat_dir}/ants_reg_tforms/


# EPI to MNI one Vol
${pyscripts_path}/ants_applyReg_accre.py \
                 --in_nii  ${proc_func_dir}/oneVol_EPI2T1.nii.gz \
                 --out_nii ${proc_func_dir}/ants_out/oneVol_EPI2MNI.nii.gz  \
                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
                 --tform_dir ${proc_anat_dir}/ants_reg_tforms/
fi

# EPI to MNI for whole fMRI 4D
${pyscripts_path}/ants_applyReg_accre.py \
		--in_nii ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz \
		--out_nii ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI.nii.gz \
		--ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
		--tform_dir ${proc_anat_dir}/ants_reg_tforms/

# original
# ${pyscripts_path}/ants_applyReg_accre.py \
#                 --in_nii ${current_path}/${sub_id}-${scan_num}_EPI2T1.nii.gz \
#                 --out_nii ${current_path}/${sub_id}-${scan_num}_EPI2MNI.nii.gz \
#                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
#                 --tform_dir ${base_dir}/reg_tforms/


# module change back to FSL
module load GCC/8.2.0  OpenMPI/3.1.4
module load FSL/6.0.1-Python-3.7.2

# ------------------------------

#Post-registration smoothing and nuisance regression
# --------------------------

scripts_path="/data1/neurdylab/scripts/vu_meica_pipeline"
afni_init="singularity exec --bind /data1:/data1 ${scripts_path}/afni_cmake_build_AFNI_21.1.03.sif"

#Spatial Blurring
if ! [ -e ${proc_func_dir}/${func_prefix}_mo_EPI2MNI_sm.nii.gz ]; then

        echo "spatial blurring for subject"
        cmd="${afni_init} 3dmerge -1blur_fwhm 3.0 -doall -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz" "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI.nii.gz""
        echo $cmd
        printf "Running ${cmd}\n" >> $log_file
        $cmd

fi

#Calculate Mean

if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz ]; then

	echo "calculating mean for subject"

	cmd="${afni_init} 3dTstat -mean -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz" "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz""
	echo $cmd
	printf "Running ${cmd}\n" >> $log_file
	$cmd

fi

#Nuisance Regression

if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz ]; then
	echo "completing nuisance regression for subject"
	cmd="${afni_init} 3dDetrend -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz" -polort 4 "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz""
	echo $cmd
	printf "Running ${cmd}\n" >> $log_file
	$cmd

fi


#Add Back Mean

echo "adding back the mean"
${afni_init} 3dcalc -a "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz" \
             -b "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz" \
             -expr 'a+b' -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm_nr.nii.gz"


# QA for ants:
mkdir -p ${proc_func_dir}/ants_out/QA/imgs

slicer "${proc_func_dir}/ants_out/oneVol_EPI2MNI.nii.gz" "/data1/neurdylab/MNI152_T1_2mm_brain.nii.gz"  -s 2 \
        -x 0.35 "${proc_func_dir}/ants_out/QA/imgs/sla.png" -x 0.45 "${proc_func_dir}/ants_out/QA/imgs/slb.png" -x 0.55 "${proc_func_dir}/ants_out/QA/imgs/slc.png" -x 0.65 "${proc_func_dir}/ants_out/QA/imgs/sld.png" \
        -y 0.35 "${proc_func_dir}/ants_out/QA/imgs/sle.png" -y 0.45 "${proc_func_dir}/ants_out/QA/imgs/slf.png" -y 0.55 "${proc_func_dir}/ants_out/QA/imgs/slg.png" -y 0.65 "${proc_func_dir}/ants_out/QA/imgs/slh.png" \
        -z 0.35 "${proc_func_dir}/ants_out/QA/imgs/sli.png" -z 0.45 "${proc_func_dir}/ants_out/QA/imgs/slj.png" -z 0.55 "${proc_func_dir}/ants_out/QA/imgs/slk.png" -z 0.65 "${proc_func_dir}/ants_out/QA/imgs/sll.png"

pngappend "${proc_func_dir}/ants_out/QA/imgs/sla.png" + "${proc_func_dir}/ants_out/QA/imgs/slb.png" + "${proc_func_dir}/ants_out/QA/imgs/slc.png" + "${proc_func_dir}/ants_out/QA/imgs/sld.png" + "${proc_func_dir}/ants_out/QA/imgs/sle.png" \
	+ "${proc_func_dir}/ants_out/QA/imgs/slf.png" + "${proc_func_dir}/ants_out/QA/imgs/slg.png" + "${proc_func_dir}/ants_out/QA/imgs/slh.png" + "${proc_func_dir}/ants_out/QA/imgs/sli.png" + "${proc_func_dir}/ants_out/QA/imgs/slj.png" \
	+ "${proc_func_dir}/ants_out/QA/imgs/slk.png" + "${proc_func_dir}/ants_out/QA/imgs/sll.png" "${proc_func_dir}/ants_out/QA/${func_prefix}_EPI2MNI_ants.png"


# create write permissions
chmod g+w -R $proc_func_dir $proc_anat_dir

done

