function [data_clean] = clean_image(subj, mask, fmripath)
% Read fMRI data within the brain mask,
... clean the data (remove Nan,Inf, columns with constant values), save the nii file.
data_path = [fmripath 'sherlock_recall_s' num2str(subj) '.nii'];

if isfile([fmripath 'fmri_s' num2str(subj) '.nii']) && isfile([fmripath 'fmri_s' num2str(subj) '.mat'])
   fprintf(['all files exist for the subject ' num2str(subj) '.\n' 'skip the operation.\n']);
elseif isfile([fmripath 'fmri_s' num2str(subj) '.nii'])
    data = cosmo_fmri_dataset(data_path,'mask',mask);
    data_clean=cosmo_remove_useless_data(data);
    save([fmripath 'fmri_s' num2str(subj) '.mat'], 'data_clean');
    fprintf(['nii file exist for the subject ' num2str(subj) '.\n' 'saveing mat file.\n']);
elseif isfile([fmripath 'fmri_s' num2str(subj) '.mat'])
    data = cosmo_fmri_dataset(data_path,'mask',mask);
    data_clean=cosmo_remove_useless_data(data);
    cosmo_map2fmri(data_clean,[fmripath 'fmri_s' num2str(subj) '.nii']);
    fprintf(['mat file exist for the subject ' num2str(subj) '.\n' 'saveing nii file.\n']);
else
    data = cosmo_fmri_dataset(data_path,'mask',mask);
    data_clean=cosmo_remove_useless_data(data);
    save([fmripath 'fmri_s' num2str(subj) '.mat'], 'data_clean');
    cosmo_map2fmri(data_clean,[fmripath 'fmri_s' num2str(subj) '.nii']);
    fprintf(['saveing nii and mat files for the subject ' num2str(subj) '.\n']);
end
load([fmripath 'fmri_s' num2str(subj) '.mat'])
end