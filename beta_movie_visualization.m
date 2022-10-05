clearvars; clc; close all;
addpath(fullfile(pwd,'helpers'));
cd ../;
basepath = pwd;
fmripath = fullfile(basepath, 'recallFMRI/');
groupdir = fullfile(fmripath, 'results/group/');
name_feature={'SI'};
tel=1
load('group_weights_all.mat') % load beta 11-SI

for subj=[1:17]
    disp(['subject: ' num2str(subj)])
    clear data_clean;
    group_weights{1,subj}.samples=group_weights_all{1,subj}(11,:);
end
[idxs,group_intersect_cell]=cosmo_mask_dim_intersect(group_weights);
nsubj=numel(group_intersect_cell);

for subject_i=1:nsubj 
    stacked_group=group_intersect_cell{subject_i};
    stacked_group.samples=group_intersect_cell{subject_i}.samples;
    stacked_group.sa.chunks=ones(1,1)*subject_i; %number of subj
    stacked_group.sa.targets=ones(1,1); %predictor
    group_intersect{subject_i}=stacked_group;
end
group_clean=cosmo_stack(group_intersect,1,'drop_nonunique');
group_r=group_clean.samples;

group=mean(group_r);
group_clean.sa.chunks=ones(1,1)*1; %number of subj
group_clean.sa.targets=ones(1,1); %predictor
group_clean.samples=group;
cosmo_map2fmri(group_clean,[groupdir 'group_beta_' name_feature '_movie.nii']);




