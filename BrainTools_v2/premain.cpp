#include <cstdlib>
#include <cstdio>
#ifdef WIN32
#include <windows.h>
#else

#endif
#include <brainio.h>
#include <brainpreprocessing.h>

using namespace std;

int main(int argc, char **argv)
{
	std::cout << "/-------------------------------\\" << std::endl;
    std::cout << "|            Pretool            |" << std::endl;
    std::cout << "\\-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	if (argc<3) {
		std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: PreTool baseline_folder followup_folder" << std::endl; 
		return EXIT_FAILURE;
	}
	FileNameType baselineFolderName = argv[1];
	FileNameType followupFolderName = argv[2];
	FileNameType baselineSearchName = baselineFolderName + "*.nii";
	FileNameType followupSearchName = followupFolderName + "*.nii";
    FileNameType baselineProcessedFolder = baselineFolderName + "preprocessed/";
    FileNameType followupProcessedFolder = followupFolderName + "preprocessed/";
	FileNameType correctedName, matchedName;
	FileNameType C3DN4Command;
	std::vector<FileNameType> toFindVector;

    #ifdef WIN32
	CreateDirectory(baselineProcessedFolder.c_str(),NULL);
	CreateDirectory(followupProcessedFolder.c_str(),NULL);
    #else
    mkdir(baselineProcessedFolder.c_str(), 0777);
    mkdir(followupProcessedFolder.c_str(), 0777);
    #endif

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// Get the names from Baseline and read the files (mask + image for baseline image and mask + image for follow-up)
	// First we search for the files in this order: brainmask, pd, t1, t2 and flair

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
	}
	MaskImage maskFollowup = brainio->ReadMaskImage(maskFollowupName);


	//------> PD
	toFindVector.clear();
	toFindVector.push_back("DP");
	toFindVector.push_back("PD");
	toFindVector.push_back("pd");
	toFindVector.push_back("dp");
	//-- Baseline
	FileNameType pdBaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	if (!pdBaselineName.empty())
		pdBaselineName = baselineFolderName + pdBaselineName;
	else
	{
		std::cerr << "ERROR: There is no PD image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage pdBaseline = brainio->ReadProbabilityImage(pdBaselineName);
	//-- FollowUp
	FileNameType pdFollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	if (!pdFollowupName.empty())
		pdFollowupName = followupFolderName + pdFollowupName;
	else
	{
		std::cerr << "ERROR: There is no PD image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage pdFollowup = brainio->ReadProbabilityImage(pdFollowupName);


	//------> T2
	toFindVector.clear();
	toFindVector.push_back("T2");
	toFindVector.push_back("t2");
	//-- Baseline
	FileNameType t2BaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	if (!t2BaselineName.empty())
		t2BaselineName = baselineFolderName + t2BaselineName;
	else
	{
		std::cerr << "ERROR: There is no T2 image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage t2Baseline = brainio->ReadProbabilityImage(t2BaselineName);
	//-- FollowUp
	FileNameType t2FollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	if (!t2FollowupName.empty())
		t2FollowupName = followupFolderName + t2FollowupName;
	else
	{
		std::cerr << "ERROR: There is no T2 image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage t2Followup = brainio->ReadProbabilityImage(t2FollowupName);


	//------> T1
	toFindVector.clear();
	toFindVector.push_back("MPRAGE");
	toFindVector.push_back("mprage");
	toFindVector.push_back("MPR");
	toFindVector.push_back("mpr");
	toFindVector.push_back("T1");
	toFindVector.push_back("t1");
	//-- Baseline
	FileNameType t1BaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	if (!t1BaselineName.empty())
		t1BaselineName = baselineFolderName + t1BaselineName;
	else
	{
		std::cerr << "ERROR: There is no T1/MPRAGE image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage t1Baseline = brainio->ReadProbabilityImage(t1BaselineName);
	//-- FollowUp
	FileNameType t1FollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	if (!t1FollowupName.empty())
		t1FollowupName = followupFolderName + t1FollowupName;
	else
	{
		std::cerr << "ERROR: There is no T1/MPRAGE image with any of these names in folder " 
			<< baselineFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
    ProbabilityImage t1Followup = brainio->ReadProbabilityImage(t1FollowupName);
	

	//------> T1
	toFindVector.clear();
	toFindVector.push_back("flair");
	toFindVector.push_back("FLAIR");
	toFindVector.push_back("dark_fluid");
	toFindVector.push_back("darkfluid");
	//-- Baseline
	FileNameType flairBaselineName = BrainIO::SearchForFile(baselineSearchName, toFindVector);
	ProbabilityImage flairBaseline;
	if (!flairBaselineName.empty()) {
		flairBaselineName = baselineFolderName + flairBaselineName;
		flairBaseline = brainio->ReadProbabilityImage(flairBaselineName);
	}
	//-- FollowUp
	FileNameType flairFollowupName = BrainIO::SearchForFile(followupSearchName, toFindVector);
	ProbabilityImage flairFollowup;
	if (!flairFollowupName.empty()) {
		flairFollowupName = followupFolderName + flairFollowupName;
		flairFollowup = brainio->ReadProbabilityImage(flairFollowupName);
	}
    
	/*--- Bias correction ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t        Bias correction        " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// First we normalise using N4, then we mask the images and finally we match the histograms
	// If the images exist we do nothing

	//------> PD
	//-- Baseline
	correctedName = baselineProcessedFolder + "pd_corrected.nii.gz";
	ProbabilityImage correctedPDBaseline = brainio->ReadProbabilityImage(correctedName);
	if ((correctedPDBaseline == (ProbabilityImage)NULL) && (pdBaseline != (ProbabilityImage)NULL)) {
		correctedPDBaseline = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(pdBaseline, maskBaseline);
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "pd_corrected.nii.gz",
			correctedPDBaseline
		);
	}
	//-- FollowUp
	correctedName = followupProcessedFolder + "pd_corrected.nii.gz";
	ProbabilityImage correctedPDFollowup = brainio->ReadProbabilityImage(correctedName);
	if ((correctedPDFollowup == (ProbabilityImage)NULL) && (pdFollowup != (ProbabilityImage)NULL)) {
		correctedPDFollowup = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(pdFollowup, maskFollowup);
		brainio->WriteProbabilityImage(
			followupProcessedFolder + "pd_corrected.nii.gz",
			correctedPDFollowup
		);
	}


	//------> T2
	//-- Baseline
	correctedName = baselineProcessedFolder + "t2_corrected.nii.gz";
	ProbabilityImage correctedT2Baseline = brainio->ReadProbabilityImage(correctedName);
	if ((correctedT2Baseline == (ProbabilityImage)NULL) && (t2Baseline != (ProbabilityImage)NULL)) {
		correctedT2Baseline = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(t2Baseline, maskBaseline);
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "t2_corrected.nii.gz",
			correctedT2Baseline
		);
	}
	//-- FollowUp
	correctedName = followupProcessedFolder + "t2_corrected.nii.gz";
	ProbabilityImage correctedT2Followup = brainio->ReadProbabilityImage(correctedName);
	if ((correctedT2Followup == (ProbabilityImage)NULL) && (t2Followup != (ProbabilityImage)NULL)) {
		correctedT2Followup = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(t2Followup, maskFollowup);
		brainio->WriteProbabilityImage(
			followupProcessedFolder + "t2_corrected.nii.gz",
			correctedT2Followup
		);
	}


	//------> T1
	//-- Baseline
	correctedName = baselineProcessedFolder + "t1_corrected.nii.gz";
	ProbabilityImage correctedT1Baseline = brainio->ReadProbabilityImage(correctedName);
	if ((correctedT1Baseline == (ProbabilityImage)NULL) && (t1Baseline != (ProbabilityImage)NULL)) {
		correctedT1Baseline = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(t1Baseline, maskBaseline);
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "t1_corrected.nii.gz",
			correctedT1Baseline
		);
	}
	//-- FollowUp
	correctedName = followupProcessedFolder + "t1_corrected.nii.gz";
	ProbabilityImage correctedT1Followup = brainio->ReadProbabilityImage(correctedName);
	if ((correctedT1Followup == (ProbabilityImage)NULL) && (t1Followup != (ProbabilityImage)NULL)) {
		correctedT1Followup = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(t1Followup, maskFollowup);
		brainio->WriteProbabilityImage(
			followupProcessedFolder + "t1_corrected.nii.gz",
			correctedT1Followup
		);
	}


	//------> FLAIR
	//-- Baseline
	correctedName = baselineProcessedFolder + "flair_corrected.nii.gz";
	ProbabilityImage correctedFLAIRBaseline = brainio->ReadProbabilityImage(correctedName);
	if ((correctedFLAIRBaseline == (ProbabilityImage)NULL) && (flairBaseline != (ProbabilityImage)NULL)) {
		correctedFLAIRBaseline = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(flairBaseline, maskBaseline);
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "flair_corrected.nii.gz",
			correctedFLAIRBaseline
		);
	}
	//-- FollowUp
	correctedName = followupProcessedFolder + "flair_corrected.nii.gz";
	ProbabilityImage correctedFLAIRFollowup = brainio->ReadProbabilityImage(correctedName);
	if ((correctedFLAIRFollowup == (ProbabilityImage)NULL) && (flairFollowup != (ProbabilityImage)NULL)) {
		correctedFLAIRFollowup = BrainPreprocessing<ProbabilityImageType>::N4BiasCorrection(flairFollowup, maskFollowup);
		brainio->WriteProbabilityImage(
			followupProcessedFolder + "flair_corrected.nii.gz",
			correctedFLAIRFollowup
		);
	}


	/*--- Histogram matching ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t       Histogram matching      " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// Histogram matching of all the images
	// If the images exist we do nothing
	matchedName = baselineProcessedFolder + "pd_corrected_matched.nii.gz";
	ProbabilityImage matchedPD = brainio->ReadProbabilityImage(matchedName);
	if ((matchedPD == (ProbabilityImage)NULL) && (correctedPDBaseline != (ProbabilityImage)NULL) && (correctedPDFollowup != (ProbabilityImage)NULL)) {
		matchedPD = BrainPreprocessing<ProbabilityImageType>::MatchHistogram(correctedPDBaseline, correctedPDFollowup);	
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "pd_corrected_matched.nii.gz",
			matchedPD
		);
	}

	matchedName = baselineProcessedFolder + "t2_corrected_matched.nii.gz";
	ProbabilityImage matchedT2 = brainio->ReadProbabilityImage(matchedName);
	if ((matchedT2 == (ProbabilityImage)NULL) && (correctedT2Baseline != (ProbabilityImage)NULL) && (correctedT2Followup != (ProbabilityImage)NULL)) {
		matchedT2 = BrainPreprocessing<ProbabilityImageType>::MatchHistogram(correctedT2Baseline, correctedT2Followup);	
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "t2_corrected_matched.nii.gz",
			matchedT2
		);
	}

	matchedName = baselineProcessedFolder + "t1_corrected_matched.nii.gz";
	ProbabilityImage matchedT1 = brainio->ReadProbabilityImage(matchedName);
	if ((matchedT1 == (ProbabilityImage)NULL) && (correctedT1Baseline != (ProbabilityImage)NULL) && (correctedT1Followup != (ProbabilityImage)NULL)) {
		matchedT1 = BrainPreprocessing<ProbabilityImageType>::MatchHistogram(correctedT1Baseline, correctedT1Followup);	
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "t1_corrected_matched.nii.gz",
			matchedT1
		);
	}

	matchedName = baselineProcessedFolder + "flair_corrected_matched.nii.gz";
	ProbabilityImage matchedFLAIR = brainio->ReadProbabilityImage(matchedName);
	if ((matchedFLAIR == (ProbabilityImage)NULL) && (correctedFLAIRBaseline != (ProbabilityImage)NULL) && (correctedFLAIRFollowup != (ProbabilityImage)NULL)) {
		matchedFLAIR = BrainPreprocessing<ProbabilityImageType>::MatchHistogram(correctedFLAIRBaseline, correctedFLAIRFollowup);	
		brainio->WriteProbabilityImage(
			baselineProcessedFolder + "flair_corrected_matched.nii.gz",
			matchedFLAIR
		);
	}

	delete(brainio);

    return EXIT_SUCCESS;

}
