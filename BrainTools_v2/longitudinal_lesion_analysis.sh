#!/bin/sh
#-------------------------------------------------------------------------------------------------------------------------------
InputDataDir='Path to Data Folder' #The Data folder: the folder that contains all the patients data
CodeSource='Path to Code Folder'   #The Code folder: the folder that contains all the required binarires
#-------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------
# Loop througth the patients you want to preprocess
for PatientNum in 622 624 637 641 642 645 661 668 677 690 709 713 714 722 727 728 730 738 739 741 744 747 758 762 765 769 772 774 779 786 795 815 821 829 834 843
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

	
	
	
