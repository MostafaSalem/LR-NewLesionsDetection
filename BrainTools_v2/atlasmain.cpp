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
    std::cout << "|           Atlastool           |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters
	FileNameType baseFolderName, imagesFolderName, searchName, transformsFolderName, atlasFolderName;
	switch (argc) {
		case 2:
			baseFolderName = argv[1];
            imagesFolderName = baseFolderName + "preprocessed/";
            transformsFolderName = baseFolderName + "transforms/";
            atlasFolderName = baseFolderName + "atlas/";
			break;
		case 3:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
            transformsFolderName = baseFolderName + "transforms/";
            atlasFolderName = baseFolderName + "atlas/";
			break;
		case 4:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
			transformsFolderName = baseFolderName + argv[3];
            atlasFolderName = baseFolderName + "atlas/";
		case 5:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName +  argv[2];
			transformsFolderName = baseFolderName + argv[3];
			atlasFolderName = baseFolderName + argv[4];
		default:
			std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: AtlasTool base_folder [images_folder] [transforms_folder] [atlas_folder]" << std::endl; 
			return EXIT_FAILURE;
			break;
	}
	// Get the names from input and read the files (mask + image for base image and mask + image for reference + transformation for base image)
    #ifdef WIN32
    CreateDirectory(atlasFolderName.c_str(),NULL);
    #else
    mkdir(atlasFolderName.c_str(), 0777);
    #endif
	std::vector<FileNameType> toFindVector;

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	//------>  T1
	searchName = imagesFolderName + "*.nii.gz";
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
		return EXIT_FAILURE;
	}
	IntensityImage t1 = brainio->ReadIntensityImage(t1Name);

	//------>  Image for T1 resampling
	toFindVector.clear();
	toFindVector.push_back("pd_corrected");
	FileNameType referenceName = BrainIO::SearchForFile(searchName, toFindVector);
	if (!referenceName.empty()) {
		referenceName = imagesFolderName + referenceName;
		IntensityImage pd = brainio->ReadIntensityImage(referenceName);
		AffineTransform affine = brainio->ReadAffineTransform(t1TransformName);
		t1 = BrainRegistration::ResampleImage<AffineTransformType>(pd, t1, affine);
	}
	else {
		toFindVector.clear();
		toFindVector.push_back("flair_corrected");
		FileNameType referenceName = BrainIO::SearchForFile(searchName, toFindVector);
		referenceName = imagesFolderName + referenceName;
		IntensityImage flair = brainio->ReadIntensityImage(referenceName);
		AffineTransform affine = brainio->ReadAffineTransform(t1TransformName);
		t1 = BrainRegistration::ResampleImage<AffineTransformType>(flair, t1, affine);
	}

	//------>  Mask to help during registration
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
		return EXIT_FAILURE;
	}

	/*--- Atlas registration ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t       Atlas registration      " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
	// Check if registration was already applied
	if (brainio->ReadProbabilityImage(atlasFolderName + "atlas_similarity.nii.gz") == (ProbabilityImage)NULL) {
		// Preparing T1
		BrainRegistration *brainreg = new BrainRegistration();
		MaskImage mask = brainio->ReadMaskImage(maskName);
		// Preparing the registration
		brainreg->SetFixed(BrainPreprocessing<IntensityImageType>::MaskImage(t1, mask));
		brainreg->SetMoving(brainio->ReadIntensityImage("./ICBM452/icbm452_atlas_template.nii"));
		brainreg->AddAtlas(brainio->ReadProbabilityImage("./ICBM452/icbm452_atlas_probability_csf.nii"));
		brainreg->AddAtlas(brainio->ReadProbabilityImage("./ICBM452/icbm452_atlas_probability_gray.nii"));
		brainreg->AddAtlas(brainio->ReadProbabilityImage("./ICBM452/icbm452_atlas_probability_white.nii"));
		brainreg->RigidRegistration();
		brainreg->AffineRegistration();
        //brainreg->BSplineMultiRegistration(3, 8, 100);
		ProbabilityImage similarity = brainreg->GetSimilarity();
		std::vector<ProbabilityImage> atlases = brainreg->GetAtlases();
		
		brainio->WriteProbabilityImage(
			atlasFolderName + "atlas_similarity.nii.gz",
			similarity
		);
		for (unsigned int atlas = 0; atlas < atlases.size(); atlas++) {
			std::stringstream out;
			out << atlas;
			brainio->WriteProbabilityImage(
				atlasFolderName + "atlas" + out.str() + ".nii.gz",
				atlases[atlas]
			);
		}

		delete(brainreg);
	}

	delete(brainio);

    return EXIT_SUCCESS;

}
