#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainregistration.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Subtool            |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();


	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName, baselineSegmentationFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName, followupSegmentationFolderName;
	FileNameType subtractionFolderName;
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
		default:
			std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: SubtractionTool baseline_folder followup_folder [images_folder] [transforms_folder] [segmentation_folder] [subtraction_folder]" << std::endl; 
			return EXIT_FAILURE;
			break;
	}
    #ifdef WIN32
    CreateDirectory(subtractionFolderName.c_str(),NULL);
    #else
    mkdir(subtractionFolderName.c_str(), 0777);
    #endif


	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------> PD 
	//-- Baseline
	FileNameType pdBaselineName = baselineImagesFolderName + "pd_corrected_matched.nii.gz";
	ProbabilityImage pdBaseline = brainio->ReadProbabilityImage(pdBaselineName);
	//-- FollowUp
	FileNameType pdFollowupName = followupImagesFolderName + "pd_corrected.nii.gz";
	ProbabilityImage pdFollowup = brainio->ReadProbabilityImage(pdFollowupName);

	//------> T2
	//-- Baseline
	FileNameType t2BaselineName = baselineImagesFolderName + "t2_corrected_matched.nii.gz";
	ProbabilityImage t2Baseline = brainio->ReadProbabilityImage(t2BaselineName);
	//-- FollowUp
	FileNameType t2FollowupName = followupImagesFolderName + "t2_corrected.nii.gz";
	ProbabilityImage t2Followup = brainio->ReadProbabilityImage(t2FollowupName);

	//------> FLAIR
	//-- Baseline
	FileNameType flairBaselineName = baselineImagesFolderName + "flair_registered.nii.gz";
	ProbabilityImage flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	FileNameType flairBaselineTransformName = baselineTransformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform flairBaselineAffine;
	if (pdBaseline != (ProbabilityImage)NULL)
		flairBaselineAffine = brainio->ReadAffineTransform(flairBaselineTransformName);
	//-- FollowUp
	FileNameType flairFollowupName = followupImagesFolderName + "flair_corrected.nii.gz";
	ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	FileNameType flairFollowupTransformName = followupTransformsFolderName + "affineFLAIRtoSpace.tfm";
	if (pdFollowup != (ProbabilityImage)NULL) {
		AffineTransform flairFollowupAffine = brainio->ReadAffineTransform(flairFollowupTransformName);
		flairFollowup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairFollowup, flairFollowupAffine);
	}

	//------> T1
	//-- Baseline
	FileNameType t1BaselineName = baselineImagesFolderName + "t1_registered.nii.gz";
	ProbabilityImage t1Baseline = brainio->ReadProbabilityImage(t1BaselineName);
	FileNameType t1BaselineTransformName = baselineTransformsFolderName + "affineT1toSpace.tfm";
	AffineTransform t1BaselineAffine;
	if (pdBaseline != (ProbabilityImage)NULL)
		t1BaselineAffine = brainio->ReadAffineTransform(t1BaselineTransformName);
	//-- FollowUp
	FileNameType t1FollowupName = followupImagesFolderName + "t1_corrected.nii.gz";
	ProbabilityImage t1Followup = brainio->ReadProbabilityImage(t1FollowupName);
	FileNameType t1FollowupTransformName = followupTransformsFolderName + "affineT1toSpace.tfm";
	if (pdFollowup != (ProbabilityImage)NULL) {
		AffineTransform t1FollowupAffine = brainio->ReadAffineTransform(t1FollowupTransformName);
		t1Followup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t1Followup, t1FollowupAffine);
	}

	//------> Brain Masks files
	//-- Baseline
	MaskImage baselineBrainMask = brainio->ReadMaskImage(baselineFolderName + "brainmask.nii");
	//-- FollowUp
	MaskImage followupBrainMask = brainio->ReadMaskImage(followupFolderName + "brainmask.nii");
	
	//------> Segmentation files
	//-- Baseline
	MaskImage baselineWMMask = brainio->ReadMaskImage(baselineSegmentationFolderName + "wm_mask.nii.gz");
	//-- FollowUp
	MaskImage followupWMMask = brainio->ReadMaskImage(followupSegmentationFolderName + "wm_mask.nii.gz");


	/*--- Registration + subtraction ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Registration          " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We set the fixed and moving images we just preprocessed
    std::cout << "\t/-- Affine registration" << std::endl;
	FileNameType affineName = baselineTransformsFolderName + "baseline2followup.tfm";
	AffineTransform affine = brainio->ReadAffineTransform(affineName);
	if (affine == (AffineTransform)NULL) {
		if ((pdFollowup != (ProbabilityImage)NULL) && (pdBaseline != (ProbabilityImage)NULL)) {
			affine = BrainRegistration::MultiAffineRegistration(pdFollowup, pdBaseline,3000);
			brainio->WriteAffineTransform(
				baselineTransformsFolderName + "baseline2followup.tfm",
				affine
			);
		} else if ((pdFollowup != (ProbabilityImage)NULL) && (flairBaseline != (ProbabilityImage)NULL)) {
			affine = BrainRegistration::MultiAffineRegistration(pdFollowup, flairBaseline,3000);
			brainio->WriteAffineTransform(
				baselineTransformsFolderName + "baseline2followup.tfm",
				affine
			);
		}
		else {
			std::cerr << "The aren't any images to perform subtraction" << std::endl; 
			return EXIT_FAILURE;
		}
	}

    std::cout << "\t/-- Image resampling" << std::endl;

    std::cout << "\t -- Brain mask" << std::endl;
	MaskImage baselineBrainMaskMoved = BrainRegistration::ResampleImage<AffineTransformType>(followupBrainMask, baselineBrainMask, affine);

	std::cout << "\t -- WM mask" << std::endl;
	MaskImage baselineWMMaskMoved = BrainRegistration::ResampleImage<AffineTransformType>(followupWMMask, baselineWMMask, affine);
	
	std::cout << "\t -- PD" << std::endl;
	ProbabilityImage pdMoved;
	if ((pdFollowup != (ProbabilityImage)NULL) && (pdBaseline != (ProbabilityImage)NULL)) {
		pdMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, pdBaseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "pd_moved.nii.gz",
			pdMoved
		);
	}


	std::cout << "\t -- T2" << std::endl;
	ProbabilityImage t2Moved;
	if (t2Baseline != (ProbabilityImage)NULL) {
		t2Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t2Baseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "t2_moved.nii.gz",
			t2Moved
		);
	}

	std::cout << "\t -- FLAIR" << std::endl;
	ProbabilityImage flairMoved;
	if (flairBaseline != (ProbabilityImage)NULL) {
		flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairBaseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "flair_moved.nii.gz",
			flairMoved
		);
	}

std::cout << "\t -- T1" << std::endl;
	ProbabilityImage t1Moved;
	if (t1Baseline != (ProbabilityImage)NULL) {
		t1Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t1Baseline, affine);
		brainio->WriteProbabilityImage(
			followupImagesFolderName + "t1_moved.nii.gz",
			t1Moved
		);
	}

	// Final ROI
	MaskImage brainMask = BrainSegmentation::Union(baselineBrainMaskMoved, followupBrainMask);
	brainio->WriteMaskImage(
		followupSegmentationFolderName + "union_brain_mask.nii.gz",
		brainMask
	);

	MaskImage wmMask = BrainSegmentation::Union(baselineWMMaskMoved, followupWMMask);
	brainio->WriteMaskImage(
		followupSegmentationFolderName + "union_wm_mask.nii.gz",
		wmMask
	);

	// The last step is to apply the subtraction
    std::cout << "\t/-- Subtractions" << std::endl;
	ProbabilityImage subtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_subtraction.nii.gz");
	if ((subtractionPD == (ProbabilityImage)NULL) && (pdFollowup != (ProbabilityImage)NULL) && (pdMoved != (ProbabilityImage)NULL)) {
		subtractionPD = BrainSegmentation::Subtraction(pdFollowup, pdMoved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_subtraction.nii.gz",
			subtractionPD
		);
	}

	ProbabilityImage subtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_subtraction.nii.gz");
	if ((subtractionT2 == (ProbabilityImage)NULL) && (t2Followup != (ProbabilityImage)NULL) && (t2Moved != (ProbabilityImage)NULL)) {
		subtractionT2 = BrainSegmentation::Subtraction(t2Followup, t2Moved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_subtraction.nii.gz",
			subtractionT2
		);
	}

	ProbabilityImage subtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_subtraction.nii.gz");
	if ((subtractionT1 == (ProbabilityImage)NULL) && (t1Followup != (ProbabilityImage)NULL) && (t1Moved != (ProbabilityImage)NULL)) {
		subtractionT1 = BrainSegmentation::Subtraction(t1Followup, t1Moved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t1_subtraction.nii.gz",
			subtractionT1
		);
	}

	ProbabilityImage subtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_subtraction.nii.gz");
	if ((subtractionFLAIR == (ProbabilityImage)NULL) && (flairFollowup != (ProbabilityImage)NULL) && (flairMoved != (ProbabilityImage)NULL)) {
		subtractionFLAIR = BrainSegmentation::Subtraction(flairFollowup, flairMoved);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_subtraction.nii.gz",
			subtractionFLAIR
		);
	}

	// Masking with WMMask
	ProbabilityImage maskedWMSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_subtraction.nii.gz");
	if ((maskedWMSubtractionPD == (ProbabilityImage)NULL) && (subtractionPD != (ProbabilityImage)NULL)) {
		maskedWMSubtractionPD = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionPD, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_wmmasked_subtraction.nii.gz",
			maskedWMSubtractionPD
		);
	}

	ProbabilityImage maskedWMSubtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_wmmasked_subtraction.nii.gz");
	if ((maskedWMSubtractionT1 == (ProbabilityImage)NULL) && (subtractionT1 != (ProbabilityImage)NULL)) {
		maskedWMSubtractionT1 = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionT1, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t1_wmmasked_subtraction.nii.gz",
			maskedWMSubtractionT1
		);
	}
	ProbabilityImage maskedWMSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_subtraction.nii.gz");
	if ((maskedWMSubtractionT2 == (ProbabilityImage)NULL) && (subtractionT2 != (ProbabilityImage)NULL)) {
		maskedWMSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionT2, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_wmmasked_subtraction.nii.gz",
			maskedWMSubtractionT2
		);
	}

	ProbabilityImage maskedWMSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_subtraction.nii.gz");
	if ((maskedWMSubtractionFLAIR == (ProbabilityImage)NULL) && (subtractionFLAIR != (ProbabilityImage)NULL)) {
		maskedWMSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::MaskImage(subtractionFLAIR, wmMask);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_wmmasked_subtraction.nii.gz",
			maskedWMSubtractionFLAIR
		);
	}
	
	//WM Masked Smoothing With SIgma 0.5,1,1.5
	//PD
	ProbabilityImage gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionPD == (ProbabilityImage)NULL) && (maskedWMSubtractionPD != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionPD = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionPD, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_wmmasked_smoothed_0.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionPD
		);
	}
	gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionPD == (ProbabilityImage)NULL) && (maskedWMSubtractionPD != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionPD = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionPD, 1);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_wmmasked_smoothed_1sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionPD
		);
	}
	gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionPD == (ProbabilityImage)NULL) && (maskedWMSubtractionPD != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionPD = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionPD, 1.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "pd_wmmasked_smoothed_1.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionPD
		);
	}	
	//T2
	ProbabilityImage gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionT2 == (ProbabilityImage)NULL) && (maskedWMSubtractionT2 != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionT2, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_wmmasked_smoothed_0.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionT2
		);
	}
	gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionT2 == (ProbabilityImage)NULL) && (maskedWMSubtractionT2 != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionT2, 1);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_wmmasked_smoothed_1sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionT2
		);
	}
	gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionT2 == (ProbabilityImage)NULL) && (maskedWMSubtractionT2 != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionT2 = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionT2, 1.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "t2_wmmasked_smoothed_1.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionT2
		);
	}
	//FLAIR
	ProbabilityImage gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionFLAIR == (ProbabilityImage)NULL) && (maskedWMSubtractionFLAIR != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionFLAIR, 0.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_wmmasked_smoothed_0.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionFLAIR
		);
	}
	gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionFLAIR == (ProbabilityImage)NULL) && (maskedWMSubtractionFLAIR != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionFLAIR, 1);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_wmmasked_smoothed_1sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionFLAIR
		);
	}
	gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionFLAIR == (ProbabilityImage)NULL) && (maskedWMSubtractionFLAIR != (ProbabilityImage)NULL)) {
		gaussWMMaskedSubtractionFLAIR = BrainPreprocessing<ProbabilityImageType>::GaussianFiltering(maskedWMSubtractionFLAIR, 1.5);
		brainio->WriteProbabilityImage(
			subtractionFolderName + "flair_wmmasked_smoothed_1.5sigma_subtraction.nii.gz",
			gaussWMMaskedSubtractionFLAIR
		);
	}

	delete(brainio);

    return EXIT_SUCCESS;

}
