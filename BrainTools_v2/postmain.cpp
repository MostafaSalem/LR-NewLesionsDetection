#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainregistration.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Posttool           |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;

    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName, baselineSegmentationFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName, followupSegmentationFolderName;
	FileNameType subtractionFolderName, deformationFolderName;
	switch (argc) {
		case 3:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
            baselineImagesFolderName = baselineFolderName + "preprocessed/";
            followupImagesFolderName = followupFolderName + "preprocessed/";
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 4:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 5:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
            baselineSegmentationFolderName = baselineFolderName + "segmentation/";
            followupSegmentationFolderName = followupFolderName + "segmentation/";
            subtractionFolderName = followupFolderName + "subtraction/";
            deformationFolderName = followupFolderName + "deformation/";
		case 6:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
            subtractionFolderName = followupFolderName + "subtraction/";
            deformationFolderName = followupFolderName + "deformation/";
		case 7:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
			subtractionFolderName = followupFolderName + argv[6];
            deformationFolderName = followupFolderName + "deformation/";
		case 8:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			baselineSegmentationFolderName = baselineFolderName + argv[5];
			followupSegmentationFolderName = followupFolderName + argv[5];
			subtractionFolderName = followupFolderName + argv[6];
			deformationFolderName = followupFolderName + argv[7];
		default:
            std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: PostTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder] [deformation_folder]" << std::endl;
			return EXIT_FAILURE;
			break;
	}

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------> Deformation
	//-- Jacobian
	FileNameType jacobian1Name = deformationFolderName + "avg3_multidemons_jacobian.nii.gz";
	ProbabilityImage avg3jacobian = brainio->ReadProbabilityImage(jacobian1Name);
	FileNameType jacobian2Name = deformationFolderName + "avg4_multidemons_jacobian.nii.gz";
	ProbabilityImage avg4jacobian = brainio->ReadProbabilityImage(jacobian2Name);

	FileNameType pdJacobianName = deformationFolderName + "pd_multidemons_jacobian.nii.gz";
	ProbabilityImage pdJacobian = brainio->ReadProbabilityImage(pdJacobianName);
	FileNameType t2JacobianName = deformationFolderName + "t2_multidemons_jacobian.nii.gz";
	ProbabilityImage t2Jacobian = brainio->ReadProbabilityImage(t2JacobianName);
	FileNameType t1JacobianName = deformationFolderName + "t1_multidemons_jacobian.nii.gz";
	ProbabilityImage t1Jacobian = brainio->ReadProbabilityImage(t1JacobianName);
	FileNameType flairJacobianName = deformationFolderName + "flair_multidemons_jacobian.nii.gz";
	ProbabilityImage flairJacobian = brainio->ReadProbabilityImage(flairJacobianName);
	//-- Divergence
	FileNameType divergence1Name = deformationFolderName + "avg3_multidemons_divergence.nii.gz";
	ProbabilityImage avg3divergence = brainio->ReadProbabilityImage(divergence1Name);
	FileNameType divergence2Name = deformationFolderName + "avg4_multidemons_divergence.nii.gz";
	ProbabilityImage avg4divergence = brainio->ReadProbabilityImage(divergence2Name);

	FileNameType pdDivergenceName = deformationFolderName + "pd_multidemons_divergence.nii.gz";
	ProbabilityImage pdDivergence = brainio->ReadProbabilityImage(pdDivergenceName);
	FileNameType t2DivergenceName = deformationFolderName + "t2_multidemons_divergence.nii.gz";
	ProbabilityImage t2Divergence = brainio->ReadProbabilityImage(t2DivergenceName);
	FileNameType t1DivergenceName = deformationFolderName + "t1_multidemons_divergence.nii.gz";
	ProbabilityImage t1Divergence = brainio->ReadProbabilityImage(t1DivergenceName);
	FileNameType flairDivergenceName = deformationFolderName + "flair_multidemons_divergence.nii.gz";
	ProbabilityImage flairDivergence = brainio->ReadProbabilityImage(flairDivergenceName);

	//------> PD 
	//-- Baseline
	//FileNameType pdBaselineName = followupImagesFolderName + "pd_moved.nii.gz";
	//ProbabilityImage pdBaseline = brainio->ReadProbabilityImage(pdBaselineName);
	//-- FollowUp
	//FileNameType pdFollowupName = followupImagesFolderName + "pd_corrected.nii.gz";
	//ProbabilityImage pdFollowup = brainio->ReadProbabilityImage(pdFollowupName);
	//-- Positive activity
	FileNameType pdBrainMaskedPositiveActivityName = subtractionFolderName + "pd_brainmasked_positive_activity.nii.gz";
	MaskImage pdBrainMaskedPositiveActivity = brainio->ReadMaskImage(pdBrainMaskedPositiveActivityName);
	FileNameType pdWMMaskedPositiveActivityName = subtractionFolderName + "pd_wmmasked_positive_activity.nii.gz";
	MaskImage pdWMMaskedPositiveActivity = brainio->ReadMaskImage(pdWMMaskedPositiveActivityName);

	//------> T2
	//-- Baseline
	//FileNameType t2BaselineName = followupImagesFolderName + "t2_moved.nii.gz";
	//ProbabilityImage t2Baseline = brainio->ReadProbabilityImage(t2BaselineName);
	//-- FollowUp
	//FileNameType t2FollowupName = followupImagesFolderName + "t2_corrected.nii.gz";
	//ProbabilityImage t2Followup = brainio->ReadProbabilityImage(t2FollowupName);
	//-- Positive activity
	FileNameType t2BrainMaskedPositiveActivityName = subtractionFolderName + "t2_brainmasked_positive_activity.nii.gz";
	MaskImage t2BrainMaskedPositiveActivity = brainio->ReadMaskImage(t2BrainMaskedPositiveActivityName);
	FileNameType t2WMMaskedPositiveActivityName = subtractionFolderName + "t2_wmmasked_positive_activity.nii.gz";
	MaskImage t2WMMaskedPositiveActivity = brainio->ReadMaskImage(t2WMMaskedPositiveActivityName);

	//------> FLAIR
	//-- Baseline
	//FileNameType flairBaselineName = followupImagesFolderName + "flair_moved.nii.gz";
	//ProbabilityImage flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	//-- FollowUp
	//FileNameType flairFollowupName = followupImagesFolderName + "flair_registered.nii.gz";
	//ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	//-- Positive activity
	FileNameType flairBrainMaskedPositiveActivityName = subtractionFolderName + "flair_brainmasked_positive_activity.nii.gz";
	MaskImage flairBrainMaskedPositiveActivity = brainio->ReadMaskImage(flairBrainMaskedPositiveActivityName);
	FileNameType flairWMMaskedPositiveActivityName = subtractionFolderName + "flair_wmmasked_positive_activity.nii.gz";
	MaskImage flairWMMaskedPositiveActivity = brainio->ReadMaskImage(flairWMMaskedPositiveActivityName);

	//------> T1
	FileNameType t1BrainMaskedPositiveActivityName = subtractionFolderName + "t1_brainmasked_positive_activity.nii.gz";
	MaskImage t1BrainMaskedPositiveActivity = brainio->ReadMaskImage(t1BrainMaskedPositiveActivityName);
	FileNameType t1WMMaskedPositiveActivityName = subtractionFolderName + "t1_wmmasked_positive_activity.nii.gz";
	MaskImage t1WMMaskedPositiveActivity = brainio->ReadMaskImage(t1WMMaskedPositiveActivityName);


	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t          Joint mask           " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	std::vector<MaskImage> masksBrain;
	if (pdBrainMaskedPositiveActivity != (MaskImage)NULL)
		masksBrain.push_back(pdBrainMaskedPositiveActivity);
	if (t2BrainMaskedPositiveActivity != (MaskImage)NULL)
		masksBrain.push_back(t2BrainMaskedPositiveActivity);
	if (flairBrainMaskedPositiveActivity != (MaskImage)NULL)
		masksBrain.push_back(flairBrainMaskedPositiveActivity);

	MaskImage intersection3Brain = BrainSegmentation::Intersection(masksBrain);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint3_brainmasked_positive_activity.nii.gz",
		intersection3Brain
	);
	if (t1BrainMaskedPositiveActivity != (MaskImage)NULL)
		masksBrain.push_back(t1BrainMaskedPositiveActivity);
	
	MaskImage intersection4Brain = BrainSegmentation::Intersection(masksBrain);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint4_brainmasked_positive_activity.nii.gz",
		intersection4Brain
	);

	ConnectComponentFilterType::Pointer connectedFilterBrain = ConnectComponentFilterType::New();
	connectedFilterBrain->SetInput(intersection3Brain);
    connectedFilterBrain->Update();
    ConnectedImage labelsBrain = connectedFilterBrain->GetOutput();	
	brainio->WriteConnectedImage(
		subtractionFolderName + "joint3_labels_brainmasked_positive_activity.nii.gz",
		labelsBrain
	);
	connectedFilterBrain->SetInput(intersection4Brain);
    connectedFilterBrain->Update();
    	 labelsBrain = connectedFilterBrain->GetOutput();	
	brainio->WriteConnectedImage(
		subtractionFolderName + "joint4_labels_brainmasked_positive_activity.nii.gz",
		labelsBrain
	);

	std::vector<MaskImage> masksWM;
	if (pdWMMaskedPositiveActivity != (MaskImage)NULL)
		masksWM.push_back(pdWMMaskedPositiveActivity);
	if (t2WMMaskedPositiveActivity != (MaskImage)NULL)
		masksWM.push_back(t2WMMaskedPositiveActivity);
	if (flairWMMaskedPositiveActivity != (MaskImage)NULL)
		masksWM.push_back(flairWMMaskedPositiveActivity);
	MaskImage intersection3WM = BrainSegmentation::Intersection(masksWM);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint3_wmmasked_positive_activity.nii.gz",
		intersection3WM
	);
	if (t1WMMaskedPositiveActivity != (MaskImage)NULL)
		masksWM.push_back(t1WMMaskedPositiveActivity);
	MaskImage intersection4WM = BrainSegmentation::Intersection(masksWM);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint4_wmmasked_positive_activity.nii.gz",
		intersection4WM
	);

	ConnectComponentFilterType::Pointer connectedFilterWM = ConnectComponentFilterType::New();
	connectedFilterWM->SetInput(intersection3WM);
    connectedFilterWM->Update();
    ConnectedImage labelsWM = connectedFilterWM->GetOutput();	
	brainio->WriteConnectedImage(
		subtractionFolderName + "joint3_labels_wmmasked_positive_activity.nii.gz",
		labelsWM
	);
	connectedFilterWM->SetInput(intersection4WM);
    connectedFilterWM->Update();
    	labelsWM = connectedFilterWM->GetOutput();	
	brainio->WriteConnectedImage(
		subtractionFolderName + "joint4_labels_wmmasked_positive_activity.nii.gz",
		labelsWM
	);

	
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t     Demons postprocessing     " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    //Brain Masked Positive Activities
	if (pdBrainMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage pdBrainMaskedLesions = BrainSegmentation::DeformationPostProcessing(pdBrainMaskedPositiveActivity, pdJacobian, pdDivergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "pd_demons_brainmasked_positive_activity.nii.gz",
			pdBrainMaskedLesions
		);
	}
	if (t2BrainMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage t2BrainMaskedLesions = BrainSegmentation::DeformationPostProcessing(t2BrainMaskedPositiveActivity, t2Jacobian, t2Divergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "t2_demons_brainmasked_positive_activity.nii.gz",
			t2BrainMaskedLesions
		);
	}
	if (flairBrainMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage flairBrainMaskedLesions = BrainSegmentation::DeformationPostProcessing(flairBrainMaskedPositiveActivity, flairJacobian, flairDivergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "flair_demons_brainmasked_positive_activity.nii.gz",
			flairBrainMaskedLesions
		);
	}
	if (t1BrainMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage t1BrainMaskedLesions = BrainSegmentation::DeformationPostProcessing(t1BrainMaskedPositiveActivity, t1Jacobian, t1Divergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "t1_demons_brainmasked_positive_activity.nii.gz",
			t1BrainMaskedLesions
		);
	}
	MaskImage lesionsBrain = BrainSegmentation::DeformationPostProcessing(intersection3Brain, avg3jacobian, avg3divergence);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint3_demons_brainmasked_positive_activity.nii.gz",
		lesionsBrain
	);
			 lesionsBrain = BrainSegmentation::DeformationPostProcessing(intersection4Brain, avg4jacobian, avg4divergence);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint4_demons_brainmasked_positive_activity.nii.gz",
		lesionsBrain
	);
	//WM Masked Positive Activities
		if (pdWMMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage pdWMMaskedLesions = BrainSegmentation::DeformationPostProcessing(pdWMMaskedPositiveActivity, pdJacobian, pdDivergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "pd_demons_wmmasked_positive_activity.nii.gz",
			pdWMMaskedLesions
		);
	}
	if (t2WMMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage t2WMMaskedLesions = BrainSegmentation::DeformationPostProcessing(t2WMMaskedPositiveActivity, t2Jacobian, t2Divergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "t2_demons_wmmasked_positive_activity.nii.gz",
			t2WMMaskedLesions
		);
	}
	if (flairWMMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage flairWMMaskedLesions = BrainSegmentation::DeformationPostProcessing(flairWMMaskedPositiveActivity, flairJacobian, flairDivergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "flair_demons_wmmasked_positive_activity.nii.gz",
			flairWMMaskedLesions
		);
	}
	if (t1WMMaskedPositiveActivity != (MaskImage)NULL) {
		MaskImage t1WMMaskedLesions = BrainSegmentation::DeformationPostProcessing(t1WMMaskedPositiveActivity, t1Jacobian, t1Divergence);
		brainio->WriteMaskImage(
			subtractionFolderName + "t1_demons_wmmasked_positive_activity.nii.gz",
			t1WMMaskedLesions
		);
	}
	MaskImage lesionsWM = BrainSegmentation::DeformationPostProcessing(intersection3WM, avg3jacobian, avg3divergence);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint3_demons_wmmasked_positive_activity.nii.gz",
		lesionsWM
	);
			 lesionsWM = BrainSegmentation::DeformationPostProcessing(intersection4WM, avg4jacobian, avg4divergence);
	brainio->WriteMaskImage(
		subtractionFolderName + "joint4_demons_wmmasked_positive_activity.nii.gz",
		lesionsWM
	);

	/*std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t     Ratio postprocessing      " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	masks.clear();
	if (pdPositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- PD" << std::endl;
		MaskImage pdLesions = BrainSegmentation::OnurPostProcessing(pdPositiveActivity, pdBaseline, pdFollowup);
		brainio->WriteMaskImage(
			subtractionFolderName + "pd_onur_positive_activity.nii.gz",
			pdLesions
		);
		masks.push_back(pdLesions);
	}
	if (t2PositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- T2" << std::endl;
		MaskImage t2Lesions = BrainSegmentation::OnurPostProcessing(t2PositiveActivity, t2Baseline, t2Followup);
		brainio->WriteMaskImage(
			subtractionFolderName + "t2_onur_positive_activity.nii.gz",
			t2Lesions
		);
		masks.push_back(t2Lesions);
	}
	if  (flairPositiveActivity != (MaskImage)NULL) {
        std::cout << "\t/-- FLAIR" << std::endl;
		MaskImage flairLesions = BrainSegmentation::OnurPostProcessing(flairPositiveActivity, flairBaseline, flairFollowup);
		brainio->WriteMaskImage(
			subtractionFolderName + "flair_onur_positive_activity.nii.gz",
			flairLesions
		);
		masks.push_back(flairLesions);
	}
    std::cout << "\t/-- JOINT" << std::endl;
	MaskImage onurLesions = BrainSegmentation::Intersection(masks);
	brainio->WriteMaskImage(
		subtractionFolderName + "pdt2_joint_onur_positive_activity.nii.gz",
		onurLesions
	);
*/
	delete(brainio);

    return EXIT_SUCCESS;

}
