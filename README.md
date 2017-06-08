# LR-NewLesionsDetection: Logistic Regression Model for New Multiple Sclerosis Lesions Detection

__BrainTools_v2__ C++ toolbox based on ITK for brain MRI processing. It is a C++ repository with several fully automated tools to process brain MRI data from MS patients. It is an updated version of __BrainTools__ (https://github.com/marianocabezas/braintools) that supports version **3.20** of ITK. Currently it implements _bias correction_ using the __N4__ algorithm; _atlas registration_ using affine and b-splines transformations; _tissue segmentation_ based on <a name="deref1">[[1](#ref1)]</a>; _lesion segmentation and detection_ also based on [[1](#ref1)].

# Dependencies:
+ This code supports version **4.9.0** of ITK. Also, the CMakelists might need some tweaking to define the ITK path depending on your OS.
+ CMake is needed to create a Makefile or project for your desired IDE. We recommend using CMake version 3+ for better compatibility.
+ In order to complement the preprocessing steps with skull stripping we recommend installing a separate toolbox. For instance, [BET](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET) is a widely used toolbox in medical imaging. In the example Bash script we provide in this README we decided to use [ROBEX](https://sites.google.com/site/jeiglesias/ROBEX) due to their improved results and better integration with our pipeline. You can download [ROBEX](https://www.nitrc.org/projects/robex).
+ The **ICBM 452** atlas is not provided with this software. You can dowload it from the [LONI website](http://www.loni.usc.edu/atlases/Atlas_Detail.php?atlas_id=6). Once downloaded, the template file and the probabilistic maps should be located in ICBM452 folder in the directory where you have your compiled files.

# Defaults
To use these tools you will need to acquire T1-w, T2-w, PD-w and T2-FLAIR-w images for each timepoint and patient. In order to automatically detect each image we assume that the image name will contain any of the substrings defined in the following table:

Brain mask |  T1-w  | T2-w | PD-w | T2-FLAIR-w
---------- | ------ | ---- | ---- | ----------
brain_mask | MPRAGE |  T2  |  DP  | flair
brainmask  | mprage |  t2  |  PD  | FLAIR
BrainMask  |  MPR   |      |  dp  | dark_fluid
           |  mpr   |      |  pd  | darkfluid
           |   T1   |      |      |
           |   t1   |      |      |

We have also defined default names for the following variables used when executing the tools:
+ **images_folder**: preprocessed/
+ **transforms_folder**: transforms/
+ **atlas_folder**: atlas/
+ **segmentation_folder**: segmentation/
+ **subtraction_folder**: subtraction/
+ **deformation_folder**: deformation/

# Bias correction
The name of the tool that implements this step is PreTool. The syntax for this tool is as follows:

```bash
PreTool baseline_folder followup_folder
```

This preprocessing tool is adapted to longitudinal analysis, thus, both followup and baseline folder must be specified. Also, we assume that a brainmask is located inside the patient's folder.
The output will be written to the _preprocessed/_ subdirectory inside the timepoint's folder.

To adapt that tool check the file check the file __premain.cpp__ and make the necessary adjustments.

# Coregistration
The name of the tool that implements this step is CoregTool. The syntax for this tool is as follows:

```bash
CoregTool base_folder [images_folder] [transformations_folder]
```

This tool coregisters the preprocessed FLAIR and T1 images inside *base_folder/images_folder*. The affine transformations are then stored inside *base_folder/transforms_folder*.

To adapt that tool check the file __coregmain.cpp__ and make the necessary adjustments.

# Atlas registration
The name of the tool that implements this step is AtlasTool. The syntax for this tool is as follows:

```bash
AtlasTool base_folder [images_folder] [transforms_folder] [atlas_folder]
```

This tool registers the **ICBM 452** atlas to a preprocessed T1 image inside *base_folder/images_folder* with its coregistration transformations inside *base_folder/transforms_folder*. As such, you will need to run PreTool and CoregTool first. Also, we assume that a brainmask is located inside the patient's folder.

The registered probabilistic maps will be written to the *base_folder/atlas_folder* subdirectory inside the timepoint's folder.

To adapt that tool check the file __atlasmain.cpp__ and make the necessary adjustments.

# Tissue and lesion segmentation
The name of the tool that implements this step is TissueTool. The syntax for this tool is as follows:

```bash
TissueTool base_folder [images_folder] [transforms_folder] [atlas_folder] [segmentation_folder]
```

This tool segments tissues using T1-w, T2-w, PD-w images stored in *base_folder/images_folder* with its coregistration transformations inside *base_folder/transforms_folder*. The segmentation based on [[1](#ref1)] uses an atlas located inside *base_folder/atlas_folder* to create an initial tissue segmentation with lesion voxels missclassified as tissue. This segmentation also includes a *partial volume* class that should be reclassified as either *grey matter* or *cerebro-spinal fluid*. After this initial segmentation, lesions are segmented using the T2-FLAIR-w image.

Finally, the tissue and lesion masks are stored inside *base_folder/segmentation_folder*.

To adapt that tool check the file check the file __tissuemain.cpp__ and make the necessary adjustments.

# Subtraction
The name of the tool that implements this step is SubtractionTool. The syntax for this tool is as follows:

```bash
SubtractionTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder]
```

This tool applies a subtraction operation between the PD-w, T2-w and T2-FLAIR-w images inside the *baseline_folder/images_folder* and *followup_folder/images_folder* after applying and affine registration. Afterwards, a gaussian smoothing is applied and an initial new lesions mask is computed applying a threshold. All results are then saved inside *followup_folder/subtraction_folder*.

To adapt that tool check the file check the file __submain.cpp__ and make the necessary adjustments.

# Demons registration and deformation analysis

The name of the tool that implements this step is DeformableTool. The syntax for this tool is as follows:

```bash
DeformableTool baseline_folder followup_folder [images_folder] [transforms_folder] [deformation_folder]
```

This tool applies a multiresolution Demons algorithm from ITK between the PD-w, T2-w and T2-FLAIR-w images inside the *baseline_folder/images_folder* and *followup_folder/images_folder* after obtaining and initial affine transformation between them. For each deformation field, divergence, jacobian, Norm, and Norm*Divergence images are computed. All these outputs are stored in *followup_folder/deformation_folder*.

To adapt that tool check the file check the file __defomain.cpp__ and make the necessary adjustments.

__Logistic Regression Pipeline__
# Final new lesion segmentation

__LR-Model__: The name of the tool that implements this step is __PredictLRModel__. It is a Matlab .m file that (you have to specify the data path inside) obtains the final mask for the new and enhancing lesions by loading a trained logistic regression model [[2](#ref2)] and use it to test your data accordinglly (You can also train the model with your own data using the Matlab .m file __TrainLRModel__).

To adapt that PredictLRModel tool check the file check the file __PredictLRModel.m__ and make the necessary adjustments.
To train the model with your data, check the file __TrainLRModel.m__ and make the necessary adjustments.

# Example Bash Script (Preprocessing)
This script fully integrates the pipeline described in [[2](#ref2)] using ROBEX with an Ubuntu 14.04 machine. You have to write InputDataDir which is the path to your data. We assume that each patient has a folder with a numerical code inside InputDataDir and for each patient's folder two  timepoints are labeled Basal and 12M. We assume that the code folder contains a folder _ROBEX_ which contains the ROBEX Skullstripping binaries and Longitudinal_Lesion_Analysis_Tools which contains the BrainTools2 binaries.

```Bash Script
#!/bin/sh
#-------------------------------------------------------------------------------------------------------------------------------
InputDataDir='full path to the data folder' #The Data folder: the folder that contains all the patients data
CodeSource='full path to the code folder'                       #The Code folder: the folder that contains all the required binarires
#-------------------------------------------------------------------------------------------------------------------------------
# Loop througth the patients you want to preprocess
for PatientNum in 622 624 637 641 642 645 661 668 677 690 709 713 714 722 727 728 730 738 739 741 744 747 758 762 765 769 772 
do
#-------------------------------------------------------------------------------------------------------------------------------
echo "                                              Starting Working With Patient Num:  $PatientNum "

################################################################################################################
echo "                                              ---------------------------------"
echo "                                              ROBEX SkullStripping Basal&12M  Images    $PatientNum"
echo "                                              ---------------------------------"
   	$CodeSource/ROBEX/runROBEX.sh $InputDataDir/$PatientNum/Basal/PD.nii $InputDataDir/$PatientNum/Basal/brain.nii $InputDataDir/$PatientNum/Basal/brainmask.nii
	rm $InputDataDir/$PatientNum/Basal/brain.nii
	
   	$CodeSource/ROBEX/runROBEX.sh $InputDataDir/$PatientNum/12M/PD.nii $InputDataDir/$PatientNum/12M/brain.nii $InputDataDir/$PatientNum/12M/brainmask.nii
	rm $InputDataDir/$PatientNum/12M/brain.nii
################################################################################################################
echo "                                              ---------------------------------"
echo "                                              PreProcessing Basal&12M  Images    $PatientNum"
echo "                                              ---------------------------------"
        $CodeSource/Longitudinal_Lesion_Analysis_Tools/PreTool $InputDataDir/$PatientNum/Basal/ $InputDataDir/$PatientNum/12M/        
################################################################################################################
echo "                                              ---------------------------------"
echo "                                              Coregistration Basal&12M Images   $PatientNum"
echo "                                              ---------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/CoregTool $InputDataDir/$PatientNum/Basal/
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/CoregTool $InputDataDir/$PatientNum/12M/
################################################################################################################
echo "                      ---------------------------------------------------------------------------------------------------------------------"
echo "                      AtlasTool on Basal&12M to get the register the ICBM atlas to Basal&12M PD Space and generate the probability maps    $PatientNum"
echo "                      ---------------------------------------------------------------------------------------------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/AtlasTool $InputDataDir/$PatientNum/Basal/
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/AtlasTool $InputDataDir/$PatientNum/12M/       
################################################################################################################
echo "                                              --------------------------------------------------"
echo "                                              TissueTool on Basal&12M to get the WM &WML Masks   $PatientNum"
echo "                                              --------------------------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/TissueTool $InputDataDir/$PatientNum/Basal/
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/TissueTool $InputDataDir/$PatientNum/12M/       
################################################################################################################
echo "                                              ------------------------------------------------------------"
echo "                                              SubtractionTool on tim1&12M to get the positive activities   $PatientNum"
echo "                                              ------------------------------------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/SubtractionTool $InputDataDir/$PatientNum/Basal/ $InputDataDir/$PatientNum/12M/
       
################################################################################################################
echo "                                              ------------------------------------------------------------"
echo "                                            VoxelSelectoinTool on tim1&12M to get the voel selectoins images  $PatientNum"
echo "                                              ------------------------------------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/VoxelSelectionTool $InputDataDir/$PatientNum/12M/

################################################################################################################
echo "                                               -----------------------------------------------------------"
echo "                                               DeformableTool on tim1&12M to get the positive activities   $PatientNum"
echo "                                               -----------------------------------------------------------"
       $CodeSource/Longitudinal_Lesion_Analysis_Tools/DeformableTool $InputDataDir/$PatientNum/Basal/  $InputDataDir/$PatientNum/12M/
       
################################################################################################################

echo ""
echo "Done Patient Num:  $PatientNum  ...................................."

done
echo "########################################################"
echo ""
echo ""
echo "Done All Patients ...................................."
echo ""
echo ""
echo "########################################################"
```
# Example Matlab (Training) 

```Matlab Code
%The TrainingTable table contains the features and Responses (True labels)
TrainingTable = table(T1_Basal,T2_Basal,PD_Basal,FLAIR_Basal,...
    T1_FollowUp,T2_FollowUp,PD_FollowUp,FLAIR_FollowUp,...
    T1_Diff,T2_Diff,PD_Diff,FLAIR_Diff,...
    T1_Jacobian,T2_Jacobian,PD_Jacobian,FLAIR_Jacobian,...
    T1_Divergence,T2_Divergence,PD_Divergence,FLAIR_Divergence,...
    T1_NormDiv,T2_NormDiv,PD_NormDiv,FLAIR_NormDiv,...
    Response);

display('Training');
[LogisticRegressionModel, x]=trainLRClassifier(TrainingTable);
display('Saving The LR-Model');
%save the LogisticRegressionModel
save('Path to code folder/LogisticRegressionModel.m','LogisticRegressionModel');
```

# Example Matlab (Testing)  

```Matlab Code
%The TestingTable table contains the features
TestingTable = table(T1_Basal,T2_Basal,PD_Basal,FLAIR_Basal,...
        T1_FollowUp,T2_FollowUp,PD_FollowUp,FLAIR_FollowUp,...
        T1_Diff,T2_Diff,PD_Diff,FLAIR_Diff,...
        T1_Jacobian,T2_Jacobian,PD_Jacobian,FLAIR_Jacobian,...
        T1_Divergence,T2_Divergence,PD_Divergence,FLAIR_Divergence,...
        T1_NormDiv,T2_NormDiv,PD_NormDiv,FLAIR_NormDiv);
        
    display('Testing');
    
    load('Path to code folder/LogisticRegressionModel.m','LogisticRegressionModel');
    %The labels contains the generated labels (threshold the probability by 0.5) and Probability contains the probability map
    [labels,Probability]=LogisticRegressionModel.predictFcn(TestingTable);
```

<a name="ref1">[[1](#deref1)]</a> M. Cabezas, A. Oliver, E. Roura, J. Freixenet, J.C Vilanova, Ll. Ramió-Torrentà, A. Rovira, X. Lladó. [_Automatic multiple sclerosis lesion detection in brain MRI by FLAIR thresholding_](https://doi.org/10.1016/j.cmpb.2014.04.006). __Computer Methods and Programs in Biomedicine__, 115(3), pp. 147-161. 2014

<a name="ref2">[[2](#deref2)]</a> M. Salem, M. Cabezasa, S. Valverdea, D. Paretoc, A. Rovira, J. Salvia, A. Olivera, X. Lladó. [_Automatic Detection of New T2 Lesions in Multiple Sclerosis Using Supervised Learning and Deformation Fields_]. __under revision__, 2017.

