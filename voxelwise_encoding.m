clearvars; clc; close all;
addpath(fullfile(pwd,'helpers'));
cd ../;
basepath = pwd;
fmripath = fullfile(basepath, 'recallFMRI/');
featurepath=fullfile(basepath, 'onset/');
mask = fullfile(basepath, 'mask/','isc_bin.nii');
groupdir = fullfile(fmripath, 'results/group/');
name_feature={'SI'};
tel=1
load('group_weights_all.mat') % load beta 11-SI
for subj=[1:4 6:17]
    cd(featurepath)
    load(['social_onset' num2str(subj) '.mat']); % load onset
    cd(basepath)
    tic
    disp(['subject: ' num2str(subj)])
    clear data_clean;
    savedir = fullfile(fmripath, ['results/' num2str(subj) '/']);
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    [data_clean] = clean_image(subj, mask, fmripath); %load image, mask it, and clean the data (Nan and Inf)
    numvoxel=size(data_clean.samples,2);  %how many voxels?
    group_weights{1,subj}.samples=group_weights_all{1,subj}(11,:); %beta for social interaction=11
    group{1}=group_weights{1,subj};
    group{2}=data_clean;
    [idxs,group_intersect_cell]=cosmo_mask_dim_intersect(group);
    for i=1:2 
        stacked_group=group_intersect_cell{i};
        stacked_group.samples=group_intersect_cell{i}.samples;
        stacked_group.sa.chunks=ones(1,1); %number of subj
        stacked_group.sa.targets=ones(1,1); %predictor
        group_intersect{i}=stacked_group;
    end
    
    %social
    s=onset_fmri_prepro(:,4)==1;
    ind_s=onset_fmri_prepro(s,1:3);
    social_ind=[];
    for i=1:size(ind_s,1)
        tem_ind=ind_s(i,1):ind_s(i,2);
        social_ind = [social_ind tem_ind];
    end
    social_ind=unique(social_ind)';
    %nonsocial
    ns=onset_fmri_prepro(:,4)==-1;
    ind_ns=onset_fmri_prepro(ns,1:3);
    nonsocial_ind=[];
    for i=1:size(ind_ns,1)
        tem_ind=ind_ns(i,1):ind_ns(i,2);
        nonsocial_ind = [nonsocial_ind tem_ind];
    end
    nonsocial_ind=unique(nonsocial_ind)';
    %remove overlapping
    remov_social=ismember(social_ind,nonsocial_ind);
    remov_nonsocial=ismember(nonsocial_ind,social_ind);
    social_ind(remov_social==1)=[];
    nonsocial_ind(remov_nonsocial==1)=[];
    socialX=[social_ind ones(size(social_ind))];
    nonsocialX=[nonsocial_ind (ones(size(nonsocial_ind))*-1)];
    X=sortrows([socialX; nonsocialX]); %make X (social vs. nonsocial)
    normX=normalize(X(:,2));
    % correlation between estimated y and true y.
    numvoxel=size(group_intersect{1}.samples,2);
    y_hat =  normX * group_intersect{1}.samples; %predicted bold signal with beta from MovieWatching
    y_true=normalize(group_intersect{2}.samples);
    
    R=zeros(numvoxel,1,'single');
    for v=1:numvoxel
        R(v) = corr(y_true(X(:,1),v), y_hat(:,v)); %prediction accuracy
    end
    
    group_intersect{2}.samples=R';
    performance=group_intersect{2};
    cosmo_map2fmri(performance,[savedir 'r_' name_feature{f} '_sub' num2str(subj) '.nii']);
    group_performance{:,tel}=performance;
    tel=tel+1;
end


save([groupdir 'group_' name_feature{f} '.mat'], 'group_performance');
[idxs,group_intersect_cell]=cosmo_mask_dim_intersect(group_performance); % remove un-shared voxels across subj
nsubj=numel(group_intersect_cell);
for subject_i=1:nsubj %tedious job to do.. we need to put chunks and target info in "Searchlight_sociality_intersect_cell"
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
cosmo_map2fmri(group_clean,[groupdir 'group_r_' name_feature{f} '.nii']);
min(group_clean.samples)
max(group_clean.samples)
mean(group)

[p, obs_stat, rand_stat, pvalue_corr] = randomize_r(group_r);
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,.05);

group=mean(group_r);
group(find(h==0))=0;
group_clean.sa.chunks=ones(1,1)*1; %number of subj
group_clean.sa.targets=ones(1,1); %predictor
group_clean.samples=group;
cosmo_map2fmri(group_clean,[groupdir 'group_permu_' name_feature '.nii']);


