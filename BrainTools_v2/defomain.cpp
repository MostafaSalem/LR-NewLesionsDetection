#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainsegmentation.h>
#include <brainregistration.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|            Defotool           |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;

    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baselineFolderName, baselineImagesFolderName, baselineTransformsFolderName;
	FileNameType followupFolderName, followupImagesFolderName, followupTransformsFolderName;
	FileNameType deformationFolderName;

	switch (argc) {
		case 3:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
            baselineImagesFolderName = baselineFolderName + "preprocessed/";
            followupImagesFolderName = followupFolderName + "preprocessed/";
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 4:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
            baselineTransformsFolderName = baselineFolderName + "transforms/";
            followupTransformsFolderName = followupFolderName + "transforms/";
            deformationFolderName = followupFolderName + "deformation/";
			break;
		case 5:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
            deformationFolderName = followupFolderName + "deformation/";
		case 6:
			baselineFolderName = argv[1];
			followupFolderName = argv[2];
			baselineImagesFolderName = baselineFolderName + argv[3];
			followupImagesFolderName = followupFolderName + argv[3];
			baselineTransformsFolderName = baselineFolderName + argv[4];
			followupTransformsFolderName = followupFolderName + argv[4];
			deformationFolderName = followupFolderName + argv[5];
		default:
            std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: DeformableTool baseline_folder followup_folder [images_folder] [transforms_folder] [deformation_folder]" << std::endl;
			return EXIT_FAILURE;
			break;
	}
    #ifdef WIN32
    CreateDirectory(deformationFolderName.c_str(),NULL);
    #else
    mkdir(deformationFolderName.c_str(), 0777);
    #endif
	FileNameType baselineSearchName = baselineFolderName + "*.nii";
	FileNameType followupSearchName = followupFolderName + "*.nii";
	std::vector<FileNameType> toFindVector;

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------> BRAINMASK
	toFindVector.clear();
	toFindVector.push_back("brain_mask");
	toFindVector.push_back("brainmask");
	toFindVector.push_back("BrainMask");
	//-- Baseline
	FileNameType maskBaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	if (!maskBaselineName.empty())
		maskBaselineName = baselineFolderName + maskBaselineName;
	else
	{
		std::cerr << "ERROR: There is no brain mask with any of these names in folder "
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
		return EXIT_FAILURE;
	}
	MaskImage maskBaseline = brainio->ReadMaskImage(maskBaselineName);
	//-- FollowUp
	FileNameType maskFollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	if (!maskFollowupName.empty())
		maskFollowupName = followupFolderName + maskFollowupName;
	else
	{
		std::cerr << "ERROR: There is no brain mask with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
		return -1;
	}
	MaskImage maskFollowup = brainio->ReadMaskImage(maskFollowupName);

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
	
	//------> T1
	//-- Baseline
	FileNameType t1BaselineName = baselineImagesFolderName + "t1_corrected_matched.nii.gz";
	ProbabilityImage t1Baseline = brainio->ReadProbabilityImage(t1BaselineName);
	FileNameType t1BaselineTransformName = baselineTransformsFolderName + "affineT1toSpace.tfm";
	AffineTransform t1BaselineAffine = brainio->ReadAffineTransform(t1BaselineTransformName);
	//-- FollowUp
	FileNameType t1FollowupName = followupImagesFolderName + "t1_corrected.nii.gz";
	ProbabilityImage t1Followup = brainio->ReadProbabilityImage(t1FollowupName);
	FileNameType t1FollowupTransformName = followupTransformsFolderName + "affineT1toSpace.tfm";
	AffineTransform t1FollowupAffine = brainio->ReadAffineTransform(t1FollowupTransformName);
	if (t1FollowupAffine != (AffineTransform)NULL)
		t1Followup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t1Followup, t1FollowupAffine);
	//------> Affine from baseline to follow-up
	FileNameType affineName = baselineTransformsFolderName + "baseline2followup.tfm";
	AffineTransform affine2followup = brainio->ReadAffineTransform(affineName);


	//------> FLAIR
	//-- Baseline
	FileNameType flairBaselineName = baselineImagesFolderName + "flair_corrected_matched.nii.gz";
	ProbabilityImage flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	FileNameType flairBaselineTransformName = baselineTransformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform flairBaselineAffine = brainio->ReadAffineTransform(flairBaselineTransformName);
	//-- FollowUp
	FileNameType flairFollowupName = followupImagesFolderName + "flair_corrected.nii.gz";
	ProbabilityImage flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	FileNameType flairFollowupTransformName = followupTransformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform flairFollowupAffine = brainio->ReadAffineTransform(flairFollowupTransformName);
	if (flairFollowupAffine != (AffineTransform)NULL)
		flairFollowup = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairFollowup, flairFollowupAffine);
	
	
	/*--- Moving the baseline images ---*/
	MaskImage maskMoved = BrainRegistration::ResampleImage<AffineTransformType>(maskFollowup, maskBaseline, affine2followup);
	ProbabilityImage pdMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, pdBaseline, affine2followup);
	ProbabilityImage t2Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t2Baseline, affine2followup);
	

	AffineTransform t1Affine = AffineTransformType::New();
	if (t1BaselineAffine != (AffineTransform)NULL) {
		AffineTransformType::MatrixType t1BaselineTransformMatrix = t1BaselineAffine->GetMatrix();
		AffineTransformType::MatrixType baseline2FolloupMatrix = affine2followup->GetMatrix();
		AffineTransformType::MatrixType t1AffineTransformMatrix = baseline2FolloupMatrix * t1BaselineTransformMatrix;		
		t1Affine->SetMatrix(t1AffineTransformMatrix);
	}
	ProbabilityImage t1Moved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, t1Baseline, t1Affine);


	AffineTransform flairAffine = AffineTransformType::New();
	if (flairBaselineAffine != (AffineTransform)NULL) {
		AffineTransformType::MatrixType flairBaselineTransformMatrix = flairBaselineAffine->GetMatrix();
		AffineTransformType::MatrixType baseline2FolloupMatrix = affine2followup->GetMatrix();
		AffineTransformType::MatrixType flairAffineTransformMatrix = baseline2FolloupMatrix * flairBaselineTransformMatrix;		
		flairAffine->SetMatrix(flairAffineTransformMatrix);
	}
	ProbabilityImage flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pdFollowup, flairBaseline, flairAffine);
	
	/*--- Masking all images ---*/
	MaskImage mask = BrainSegmentation::Intersection(maskFollowup, maskMoved);
	ProbabilityImage pdMaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(pdMoved, mask);
	ProbabilityImage t1MaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(t1Moved, mask);
	ProbabilityImage t2MaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(t2Moved, mask);
	ProbabilityImage flairMaskedBaseline = BrainPreprocessing<ProbabilityImageType>::MaskImage(flairMoved, mask);

	ProbabilityImage pdMaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(pdFollowup, mask);
	ProbabilityImage t1MaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(t1Followup, mask);
	ProbabilityImage t2MaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(t2Followup, mask);
	ProbabilityImage flairMaskedFollowup = BrainPreprocessing<ProbabilityImageType>::MaskImage(flairFollowup, mask);

	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t   Demons multi Registration   " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;

	std::cout << "\tRegistration" << std::endl;
	DeformationField pdMRDeformation = brainio->ReadDeformationField(deformationFolderName + "pd_multidemons_deformation.nii.gz");
	if (pdMRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing PD" << std::endl;
		pdMRDeformation = BrainRegistration::MultiDemonsRegistration(pdMaskedFollowup, pdMaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "pd_multidemons_deformation.nii.gz",
			pdMRDeformation
		);
	}
	
	DeformationField t1MRDeformation = brainio->ReadDeformationField(deformationFolderName + "t1_multidemons_deformation.nii.gz");
	if (t1MRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing T1" << std::endl;
		t1MRDeformation = BrainRegistration::MultiDemonsRegistration(t1MaskedFollowup, t1MaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "t1_multidemons_deformation.nii.gz",
			t1MRDeformation
		);
	}

	DeformationField t2MRDeformation = brainio->ReadDeformationField(deformationFolderName + "t2_multidemons_deformation.nii.gz");
	if (t2MRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing T2" << std::endl;
		t2MRDeformation = BrainRegistration::MultiDemonsRegistration(t2MaskedFollowup, t2MaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "t2_multidemons_deformation.nii.gz",
			t2MRDeformation
		);
	}

	DeformationField flairMRDeformation = brainio->ReadDeformationField(deformationFolderName + "flair_multidemons_deformation.nii.gz");
	if (flairMRDeformation == (DeformationField)NULL) {
        std::cout << "\t/-- Processing FLAIR" << std::endl;
		flairMRDeformation = BrainRegistration::MultiDemonsRegistration(flairMaskedFollowup, flairMaskedBaseline);
        std::cout << "\t/---- Writing the results" << std::endl;
		brainio->WriteDeformationField(
			deformationFolderName + "flair_multidemons_deformation.nii.gz",
			flairMRDeformation
		);
	}

	std::vector<DeformationField> deformations;
	if (pdMRDeformation != (DeformationField)NULL)
		deformations.push_back(pdMRDeformation);
	if (t2MRDeformation != (DeformationField)NULL)
		deformations.push_back(t2MRDeformation);
	if (flairMRDeformation != (DeformationField)NULL)
		deformations.push_back(flairMRDeformation);
	DeformationField avg3MRDeformation = BrainRegistration::Sum(deformations);
	brainio->WriteDeformationField(
		deformationFolderName + "avg3_multidemons_deformation.nii.gz",
		avg3MRDeformation
	);
	if (t1MRDeformation != (DeformationField)NULL)
		deformations.push_back(t1MRDeformation);
	DeformationField avg4MRDeformation = BrainRegistration::Sum(deformations);
	brainio->WriteDeformationField(
		deformationFolderName + "avg4_multidemons_deformation.nii.gz",
		avg4MRDeformation
	);


	std::cout << "\tJacobian" << std::endl;
	if (pdMRDeformation != (DeformationField)NULL) {
		ProbabilityImage pdJacobian = BrainRegistration::Jacobian(pdMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_jacobian.nii.gz",
			pdJacobian
		);
	}

	if (t2MRDeformation != (DeformationField)NULL) {
		ProbabilityImage t2Jacobian = BrainRegistration::Jacobian(t2MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_jacobian.nii.gz",
			t2Jacobian
		);
	}

	if (t1MRDeformation != (DeformationField)NULL) {
		ProbabilityImage t1Jacobian = BrainRegistration::Jacobian(t1MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t1_multidemons_jacobian.nii.gz",
			t1Jacobian
		);
	}

	if (flairMRDeformation != (DeformationField)NULL) {
		ProbabilityImage flairJacobian = BrainRegistration::Jacobian(flairMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_jacobian.nii.gz",
			flairJacobian
		);
	}

	ProbabilityImage avg3Jacobian = BrainRegistration::Jacobian(avg3MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg3_multidemons_jacobian.nii.gz",
		avg3Jacobian
	);
	ProbabilityImage avg4Jacobian = BrainRegistration::Jacobian(avg4MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg4_multidemons_jacobian.nii.gz",
		avg4Jacobian
	);

	std::cout << "\tNorm" << std::endl;
	ProbabilityImage pdNorm;
	if (pdMRDeformation != (DeformationField)NULL) {
		pdNorm = BrainRegistration::Norm(pdMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_norm.nii.gz",
			pdNorm
		);
	}
	
	ProbabilityImage t2Norm;
	if (t2MRDeformation != (DeformationField)NULL) {
		t2Norm = BrainRegistration::Norm(t2MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_norm.nii.gz",
			t2Norm
		);
	}

	ProbabilityImage t1Norm;
	if (t1MRDeformation != (DeformationField)NULL) {
		t1Norm = BrainRegistration::Norm(t1MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t1_multidemons_norm.nii.gz",
			t1Norm
		);
	}

	ProbabilityImage flairNorm;
	if (flairMRDeformation != (DeformationField)NULL) {
		flairNorm = BrainRegistration::Norm(flairMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_norm.nii.gz",
			flairNorm
		);
	}

	ProbabilityImage avg3Norm = BrainRegistration::Norm(avg3MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg3_multidemons_norm.nii.gz",
		avg3Norm
	);

	ProbabilityImage avg4Norm = BrainRegistration::Norm(avg4MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg4_multidemons_norm.nii.gz",
		avg4Norm
	);
	
	std::cout << "\tDivergence" << std::endl;
	ProbabilityImage pdDivergence;
	if (pdMRDeformation != (DeformationField)NULL) {
		pdDivergence = BrainRegistration::Divergence(pdMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_divergence.nii.gz",
			pdDivergence
		);
	}

	ProbabilityImage t1Divergence;
	if (t1MRDeformation != (DeformationField)NULL) {
		t1Divergence = BrainRegistration::Divergence(t1MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t1_multidemons_divergence.nii.gz",
			t1Divergence
		);
	}
	ProbabilityImage t2Divergence;
	if (t2MRDeformation != (DeformationField)NULL) {
		t2Divergence = BrainRegistration::Divergence(t2MRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_divergence.nii.gz",
			t2Divergence
		);
	}

	ProbabilityImage flairDivergence;
	if (flairMRDeformation != (DeformationField)NULL) {
		flairDivergence = BrainRegistration::Divergence(flairMRDeformation);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_divergence.nii.gz",
			flairDivergence
		);
	}

	ProbabilityImage avg3Divergence = BrainRegistration::Divergence(avg3MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg3_multidemons_divergence.nii.gz",
		avg3Divergence
	);

	ProbabilityImage avg4Divergence = BrainRegistration::Divergence(avg4MRDeformation);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg4_multidemons_divergence.nii.gz",
		avg4Divergence
	);
	
	std::cout << "\tNormMulDivergence" << std::endl;
	if ( (pdNorm != (ProbabilityImage)NULL) && (pdDivergence != (ProbabilityImage)NULL)) {
		ProbabilityImage pdNormMulDivergence = BrainRegistration::NormMulDivergence(pdNorm, pdDivergence);
		brainio->WriteProbabilityImage(
			deformationFolderName + "pd_multidemons_norm_mul_divergence.nii.gz",
			pdNormMulDivergence
		);
	}

	if ( (t1Norm != (ProbabilityImage)NULL) && (t1Divergence != (ProbabilityImage)NULL)) {
		ProbabilityImage t1NormMulDivergence = BrainRegistration::NormMulDivergence(t1Norm, t1Divergence);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t1_multidemons_norm_mul_divergence.nii.gz",
			t1NormMulDivergence
		);
	}

	if ( (t2Norm != (ProbabilityImage)NULL) && (t2Divergence != (ProbabilityImage)NULL)) {
		ProbabilityImage t2NormMulDivergence = BrainRegistration::NormMulDivergence(t2Norm, t2Divergence);
		brainio->WriteProbabilityImage(
			deformationFolderName + "t2_multidemons_norm_mul_divergence.nii.gz",
			t2NormMulDivergence
		);
	}

	if ( (flairNorm != (ProbabilityImage)NULL) && (flairDivergence != (ProbabilityImage)NULL)) {
		ProbabilityImage flairNormMulDivergence = BrainRegistration::NormMulDivergence(flairNorm, flairDivergence);
		brainio->WriteProbabilityImage(
			deformationFolderName + "flair_multidemons_norm_mul_divergence.nii.gz",
			flairNormMulDivergence
		);
	}

	ProbabilityImage avg3NormMulDivergence = BrainRegistration::NormMulDivergence(avg3Norm, avg3Divergence);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg3_multidemons_norm_mul_divergence.nii.gz",
		avg3NormMulDivergence
	);

	ProbabilityImage avg4NormMulDivergence = BrainRegistration::NormMulDivergence(avg4Norm, avg4Divergence);
	brainio->WriteProbabilityImage(
		deformationFolderName + "avg4_multidemons_norm_mul_divergence.nii.gz",
		avg4NormMulDivergence
	);

	delete(brainio);

    return EXIT_SUCCESS;

}
