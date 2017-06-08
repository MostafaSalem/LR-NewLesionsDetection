clear all;
clc
DataPath='Path to data folder';
VoxelSelectionFolderName='voxelselection1.5';
VoxelSelectionFileName = 't2_wmmasked_voxel_selection_0sigma';
Num_Patients=6; % the number of patients you want to train with

%The first 4 features (Basal images intensities) (Have to be normalized)
T1_Basal=[];T2_Basal=[];PD_Basal=[];FLAIR_Basal=[];
%The first 4 features (Basal images intensities) (Have to be normalized)
T1_FollowUp=[];T2_FollowUp=[];PD_FollowUp=[];FLAIR_FollowUp=[];
%The second 4 features (Diff between Basal &12M)
T1_Diff=[];T2_Diff=[];PD_Diff=[];FLAIR_Diff=[];
%Jacobian features
T1_Jacobian=[];T2_Jacobian=[];PD_Jacobian=[];FLAIR_Jacobian=[];
%Divergence features
T1_Divergence=[];T2_Divergence=[];PD_Divergence=[];FLAIR_Divergence=[];
%NormDiv features
T1_NormDiv=[];T2_NormDiv=[];PD_NormDiv=[];FLAIR_NormDiv=[];
%Response
Response=[];

display(['Voxel Selection Step: Working with Voxel Selectoin Image ', VoxelSelectionFileName]);

%Large lesoins only
Patients=ones(Num_Patients,1);
 % you can add any number of patients as you want
Patient_Nums = {'622','624','637','641','642','645'};
for patient=1:size(Patient_Nums,2)
    
    display(['Working with Patient No. ', Patient_Nums{patient}]);
    %Read the selectoin image
    voxelSelection = load_nifti([DataPath,Patient_Nums{patient},'/12M/',VoxelSelectionFolderName,'/',VoxelSelectionFileName]);
    voxelSelectionImage = logical(voxelSelection.img);
    %The lesion image
    gt = load_nifti([DataPath,Patient_Nums{patient},'/12M/lesionMask']);
    %adding all the lesion voxels to the training samples
    voxelSelectionImage(logical(gt.img))=1;
    
    lesoinVoxelsSelected = gt.img(voxelSelectionImage);
    lesoinVoxelsSelectedLabels=categorical(lesoinVoxelsSelected);
    % ###########################################################################################
    % ###################################Basal###################################################
    t1Basal = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/t1_moved']);
    intenisties = t1Basal.img(voxelSelectionImage);
    t1Basal.img = (t1Basal.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2Basal = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/t2_moved']);
    intenisties = t2Basal.img(voxelSelectionImage);
    t2Basal.img = (t2Basal.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    pdBasal = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/pd_moved']);
    intenisties = pdBasal.img(voxelSelectionImage);
    pdBasal.img = (pdBasal.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairBasal = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/flair_moved']);
    intenisties = flairBasal.img(voxelSelectionImage);
    flairBasal.img = (flairBasal.img - mean(intenisties))/std(intenisties);
    % ###########################################################################################
    % ###################################FollowUp################################################    
    t1FollowUp = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/t1_registered']);
    intenisties = t1FollowUp.img(voxelSelectionImage);
    t1FollowUp.img = (t1FollowUp.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2FollowUp = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/t2_corrected']);
    intenisties = t2FollowUp.img(voxelSelectionImage);
    t2FollowUp.img = (t2FollowUp.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    pdFollowUp = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/pd_corrected']);
    intenisties = pdFollowUp.img(voxelSelectionImage);
    pdFollowUp.img = (pdFollowUp.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairFollowUp = load_nifti([DataPath,Patient_Nums{patient},'/12M/preprocessed/flair_registered']);
    intenisties = flairFollowUp.img(voxelSelectionImage);
    flairFollowUp.img = (flairFollowUp.img - mean(intenisties))/std(intenisties);
    % #############################################################################################
    % ###################################Subtraction################################################
    t1Subtraction = load_nifti([DataPath,Patient_Nums{patient},'/12M/subtraction/t1_subtraction']);
    intenisties = t1Subtraction.img(voxelSelectionImage);
    t1Subtraction.img = (t1Subtraction.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2Subtraction = load_nifti([DataPath,Patient_Nums{patient},'/12M/subtraction/t2_subtraction']);
    intenisties = t2Subtraction.img(voxelSelectionImage);
    t2Subtraction.img = (t2Subtraction.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % --------------------------------------------------------------Patients_IndicesTesting-----------------------------
    pdSubtraction = load_nifti([DataPath,Patient_Nums{patient},'/12M/subtraction/pd_subtraction']);
    intenisties = pdSubtraction.img(voxelSelectionImage);
    pdSubtraction.img = (pdSubtraction.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairSubtraction = load_nifti([DataPath,Patient_Nums{patient},'/12M/subtraction/flair_subtraction']);
    intenisties = flairSubtraction.img(voxelSelectionImage);
    flairSubtraction.img = (flairSubtraction.img - mean(intenisties))/std(intenisties);
    % #############################################################################################
    % ###################################Jacobian##################################################
    t1Jacobian = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t1_multidemons_jacobian']);
    intenisties = t1Jacobian.img(voxelSelectionImage);
    t1Jacobian.img = (t1Jacobian.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2Jacobian = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t2_multidemons_jacobian']);
    intenisties = t2Jacobian.img(voxelSelectionImage);
    t2Jacobian.img = (t2Jacobian.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    pdJacobian = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/pd_multidemons_jacobian']);
    intenisties = pdJacobian.img(voxelSelectionImage);
    pdJacobian.img = (pdJacobian.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairJacobian = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/flair_multidemons_jacobian']);
    intenisties = flairJacobian.img(voxelSelectionImage);
    flairJacobian.img = (flairJacobian.img - mean(intenisties))/std(intenisties);
    
    % #############################################################################################
    % ###################################Divergence################################################
    t1Divergence = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t1_multidemons_divergence']);
    intenisties = t1Divergence.img(voxelSelectionImage);
    t1Divergence.img = (t1Divergence.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2Divergence = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t2_multidemons_divergence']);
    intenisties = t2Divergence.img(voxelSelectionImage);
    t2Divergence.img = (t2Divergence.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    pdDivergence = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/pd_multidemons_divergence']);
    intenisties = pdDivergence.img(voxelSelectionImage);
    pdDivergence.img = (pdDivergence.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairDivergence = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/flair_multidemons_divergence']);
    intenisties = flairDivergence.img(voxelSelectionImage);
    flairDivergence.img = (flairDivergence.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    % #############################################################################################
    % ###################################Norm*Div##################################################
    t1NormDiv = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t1_multidemons_norm_mul_divergence']);
    intenisties = t1NormDiv.img(voxelSelectionImage);
    t1NormDiv.img = (t1NormDiv.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    t2NormDiv = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/t2_multidemons_norm_mul_divergence']);
    intenisties = t2NormDiv.img(voxelSelectionImage);
    t2NormDiv.img = (t2NormDiv.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    pdNormDiv = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/pd_multidemons_norm_mul_divergence']);
    intenisties = pdNormDiv.img(voxelSelectionImage);
    pdNormDiv.img = (pdNormDiv.img - mean(intenisties))/std(intenisties);
    % -------------------------------------------------------------------------------------------
    % -------------------------------------------------------------------------------------------
    flairNormDiv = load_nifti([DataPath,Patient_Nums{patient},'/12M/deformation/flair_multidemons_norm_mul_divergence']);
    intenisties = flairNormDiv.img(voxelSelectionImage);
    flairNormDiv.img = (flairNormDiv.img - mean(intenisties))/std(intenisties);
    % #############################################################################################
    % #############################################################################################
    % #############################################################################################
    % #############################################################################################
    
    %Basal Features
    T1_Basal=[ T1_Basal; t1Basal.img(voxelSelectionImage)];
    T2_Basal=[ T2_Basal; t2Basal.img(voxelSelectionImage)];
    PD_Basal=[ PD_Basal; pdBasal.img(voxelSelectionImage)];
    FLAIR_Basal=[ FLAIR_Basal; flairBasal.img(voxelSelectionImage)];
    
    % #############################################################################################
    %FollowUp Features
    T1_FollowUp=[ T1_FollowUp; t1FollowUp.img(voxelSelectionImage)];
    T2_FollowUp=[ T2_FollowUp; t2FollowUp.img(voxelSelectionImage)];
    PD_FollowUp=[ PD_FollowUp; pdFollowUp.img(voxelSelectionImage)];
    FLAIR_FollowUp=[ FLAIR_FollowUp; flairFollowUp.img(voxelSelectionImage)];
    
    % #############################################################################################
    %Diff Features
    T1_Diff=[ T1_Diff; t1Subtraction.img(voxelSelectionImage)];
    T2_Diff=[ T2_Diff; t2Subtraction.img(voxelSelectionImage)];
    PD_Diff=[ PD_Diff; pdSubtraction.img(voxelSelectionImage)];
    FLAIR_Diff=[ FLAIR_Diff; flairSubtraction.img(voxelSelectionImage)];
    
    % #############################################################################################
    %Jacobian Features
    T1_Jacobian=[ T1_Jacobian; t1Jacobian.img(voxelSelectionImage)];
    T2_Jacobian=[ T2_Jacobian; t2Jacobian.img(voxelSelectionImage)];
    PD_Jacobian=[ PD_Jacobian; pdJacobian.img(voxelSelectionImage)];
    FLAIR_Jacobian=[ FLAIR_Jacobian; flairJacobian.img(voxelSelectionImage)];
    
    % #############################################################################################
    %Divergence Features
    T1_Divergence=[ T1_Divergence; t1Divergence.img(voxelSelectionImage)];
    T2_Divergence=[ T2_Divergence; t2Divergence.img(voxelSelectionImage)];
    PD_Divergence=[ PD_Divergence; pdDivergence.img(voxelSelectionImage)];
    FLAIR_Divergence=[ FLAIR_Divergence; flairDivergence.img(voxelSelectionImage)];
    
    % #############################################################################################
    %NormDiv Features
    T1_NormDiv=[ T1_NormDiv; t1NormDiv.img(voxelSelectionImage)];
    T2_NormDiv=[ T2_NormDiv; t2NormDiv.img(voxelSelectionImage)];
    PD_NormDiv=[ PD_NormDiv; pdNormDiv.img(voxelSelectionImage)];
    FLAIR_NormDiv=[ FLAIR_NormDiv; flairNormDiv.img(voxelSelectionImage)];
    % #############################################################################################
    %Response
    Response = [Response;lesoinVoxelsSelectedLabels];
    
end
% #############################################################################################
TrainingTable = table(T1_Basal,T2_Basal,PD_Basal,FLAIR_Basal,...
    T1_FollowUp,T2_FollowUp,PD_FollowUp,FLAIR_FollowUp,...
    T1_Diff,T2_Diff,PD_Diff,FLAIR_Diff,...
    T1_Jacobian,T2_Jacobian,PD_Jacobian,FLAIR_Jacobian,...
    T1_Divergence,T2_Divergence,PD_Divergence,FLAIR_Divergence,...
    T1_NormDiv,T2_NormDiv,PD_NormDiv,FLAIR_NormDiv,...
    Response);
% #############################################################################################
display('Training');
[LogisticRegressionModel, x]=trainLRClassifier(TrainingTable);
display('Saving The LR-Model');
save('Path to code folder/LogisticRegressionModel.m','LogisticRegressionModel');

