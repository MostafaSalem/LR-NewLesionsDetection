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
    std::cout << "|           Tissuetool          |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baseFolderName, imagesFolderName, transformsFolderName, atlasFolderName, segmentationFolderName;
	switch (argc) {
		case 2:
			baseFolderName = argv[1];
            imagesFolderName = baseFolderName + "preprocessed/";
            transformsFolderName = baseFolderName + "transforms/";
            atlasFolderName = baseFolderName + "atlas/";
            segmentationFolderName = baseFolderName + "segmentation/";
			break;
		case 3:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
            transformsFolderName = baseFolderName + "transforms/";
            atlasFolderName = baseFolderName + "atlas/";
            segmentationFolderName = baseFolderName + "segmentation/";
			break;
		case 4:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
			transformsFolderName = baseFolderName + argv[3];
            atlasFolderName = baseFolderName + "atlas/";
            segmentationFolderName = baseFolderName + "segmentation/";
		case 5:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
			transformsFolderName = baseFolderName + argv[3];
			atlasFolderName = baseFolderName + argv[4];
            segmentationFolderName = baseFolderName + "segmentation/";
		case 6:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
			transformsFolderName = baseFolderName + argv[3];
			atlasFolderName = baseFolderName + argv[4];
			segmentationFolderName = baseFolderName + argv[5];
		default:
            std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: TissueTool base_folder [images_folder] [transforms_folder] [atlas_folder] [segmentation_folder]" << std::endl;
			return EXIT_FAILURE;
			break;
	}

    #ifdef WIN32
    CreateDirectory(segmentationFolderName.c_str(),NULL);
    #else
    mkdir(segmentationFolderName.c_str(), 0777);
    #endif
	
	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// Get the names from Baseline and read the files (mask + image for baseline image and mask + image for follow-up)
	// First we search for the files in this order: brainmask, pd, t1, t2 and flair
	// Get the names from input and read the files (mask + image for base image and mask + image for reference + transformation for base image)
	std::vector<FileNameType> toFindVector;

	//------> PD
    std::cout << "\t/- PD" << std::endl;
	FileNameType searchName = imagesFolderName + "*.nii.gz";
	toFindVector.clear();
	toFindVector.push_back("pd_corrected");
	FileNameType pdName = BrainIO::SearchForFile(searchName, toFindVector);
	if (!pdName.empty())
		pdName = imagesFolderName + pdName;
	else
	{
		std::cerr << "ERROR: There is no corrected PD image with any of these names in folder " << imagesFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
	ProbabilityImage pd = brainio->ReadProbabilityImage(pdName);

	//------> T2
    std::cout << "\t/- T2" << std::endl;
	toFindVector.clear();
	toFindVector.push_back("t2_corrected");
	FileNameType t2Name = BrainIO::SearchForFile(searchName, toFindVector);
	if (!t2Name.empty())
		t2Name = imagesFolderName + t2Name;
	else
	{
		std::cerr << "ERROR: There is no corrected T2 image with any of these names in folder " << imagesFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
	ProbabilityImage t2 = brainio->ReadProbabilityImage(t2Name);

	//------> FLAIR
    std::cout << "\t/- FLAIR" << std::endl;
	FileNameType flairTransformName = transformsFolderName + "affineFLAIRtoSpace.tfm";
	toFindVector.clear();
	toFindVector.push_back("flair_corrected");
	FileNameType flairName = BrainIO::SearchForFile(searchName, toFindVector);
	if (!flairName.empty())
		flairName = imagesFolderName + flairName;
	else
	{
		std::cerr << "ERROR: There is no corrected FLAIR image with any of these names in folder " << imagesFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
	ProbabilityImage flair = brainio->ReadProbabilityImage(flairName);
	if (pd != (ProbabilityImage)NULL) {
		AffineTransform flairAffine = brainio->ReadAffineTransform(flairTransformName);
		flair = BrainRegistration::ResampleImage<AffineTransformType>(pd, flair, flairAffine);
	}

	//------> T1
    std::cout << "\t/- T1" << std::endl;
	FileNameType t1TransformName = transformsFolderName + "affineT1toSpace.tfm";
	toFindVector.clear();
	toFindVector.push_back("t1_corrected");
	FileNameType t1Name = BrainIO::SearchForFile(searchName, toFindVector);
	if (!t1Name.empty())
		t1Name = imagesFolderName + t1Name;
	else
	{
		std::cerr << "ERROR: There is no corrected T1 image with any of these names in folder " << imagesFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
	ProbabilityImage t1 = brainio->ReadProbabilityImage(t1Name);
	AffineTransform t1Affine = brainio->ReadAffineTransform(t1TransformName);
	if (pd != (ProbabilityImage)NULL)
		t1 = BrainRegistration::ResampleImage<AffineTransformType>(pd, t1, t1Affine);
	else
		t1 = BrainRegistration::ResampleImage<AffineTransformType>(flair, t1, t1Affine);

	//------> Brainmask
    std::cout << "\t/- Brainmask" << std::endl;
	searchName = baseFolderName + "*.nii";
	toFindVector.clear();
	toFindVector.push_back("brain_mask");
	toFindVector.push_back("brainmask");
	toFindVector.push_back("BrainMask");
	FileNameType maskName = BrainIO::SearchForFile(searchName, toFindVector);
	if (!maskName.empty())
		maskName = baseFolderName + maskName;
	else
	{
		std::cerr << "ERROR: There is no brain mask with any of these names in folder " << imagesFolderName << ":" << std::endl;
		for (std::vector<FileNameType>::iterator it = toFindVector.begin(); it != toFindVector.end(); ++it)
			std::cerr << "\t <*" << *it << "*>" << std::endl;
	}
	MaskImage mask = brainio->ReadMaskImage(maskName);

    // Reading the atlases
	std::vector<ProbabilityImage> atlases;
	ProbabilityImage similarity = brainio->ReadProbabilityImage(atlasFolderName + "atlas_similarity.nii.gz");
	for (unsigned int atlas = 0; atlas < 3; atlas++) {
		std::stringstream out;
		out << atlas;
		atlases.push_back(brainio->ReadProbabilityImage(atlasFolderName + "atlas" + out.str() + ".nii.gz"));
	}

	/*--- Segmentation ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Segmentation          " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	ResultsImage brain = brainio->ReadResultsImage(segmentationFolderName + "brain.nii.gz");

	if (brain == (ResultsImage)NULL) {

		std::vector<int*> pv;
		int csfgm[2] = {0,1};
		int csfwm[2] = {0,2};
		int gmwm[2] = {1,2};
		pv.push_back(csfgm);

		// Tissue segmentation
		std::vector<ProbabilityImage> images;
		if (t1 != (ProbabilityImage)NULL)
			images.push_back(t1);
		if (t2 != (ProbabilityImage)NULL)
			images.push_back(t2);
		if (pd != (ProbabilityImage)NULL)
			images.push_back(pd);
		if ((pd == (ProbabilityImage)NULL) && (flair != (ProbabilityImage)NULL))
			images.push_back(flair);
		if (images.size()>0) {
			std::vector<ProbabilityImage> prTissue = BrainSegmentation::Generic_WEM_Atlas(atlases, images, similarity, mask, pv, 0.75, 10);
			for (unsigned int pr = 0; pr < prTissue.size(); pr++) {
				std::stringstream out;
				out << pr;
				brainio->WriteProbabilityImage(
					segmentationFolderName + "pr" + out.str() + ".nii.gz",
					prTissue[pr]
				);
			}
			ResultsImage tissue = BrainSegmentation::SolutionFromPr(prTissue, mask);

			// Lesion segmentation
			if (flair != (ProbabilityImage)NULL) {
				ProbabilityImage flairMasked = BrainPreprocessing<ProbabilityImageType>::MaskImage(flair, mask);
				MaskImage wml = BrainSegmentation::Lesion_New(tissue, flairMasked, 2, 0.8, 0.6, 10);
				brainio->WriteMaskImage(
					segmentationFolderName + "lesoin.nii.gz",
					wml
					);
				brain = BrainSegmentation::Relabel(tissue, wml);
			
				brainio->WriteResultsImage(
					segmentationFolderName + "brain.nii.gz",
					brain
				);
			}
			else
				brain = tissue;
		}
	}

	// WM mask for subtraction
	if (brain != (ResultsImage)NULL) {
		MaskImage wmMask = BrainSegmentation::GetWMMask(brain);
		brainio->WriteMaskImage(
			segmentationFolderName + "wm_mask.nii.gz",
			wmMask
		);
	}

	delete(brainio);

    return EXIT_SUCCESS;

}
