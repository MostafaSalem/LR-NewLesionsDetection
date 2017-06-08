#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainregistration.h>


using namespace std;

int main(int argc, char **argv)
{
	std::cout << "/-------------------------------\\" << std::endl;
    std::cout << "|           Coregtool           |" << std::endl;
    std::cout << "\\-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	FileNameType baseFolderName, imagesFolderName, transformsFolderName;
	switch (argc) {
		case 2:
			baseFolderName = argv[1];
            imagesFolderName = baseFolderName + "preprocessed/";
            transformsFolderName = baseFolderName + "transforms/";
			break;
		case 3:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName + argv[2];
            transformsFolderName = baseFolderName + "transforms/";
			break;
		case 4:
			baseFolderName = argv[1];
			imagesFolderName = baseFolderName + argv[2];
			transformsFolderName = baseFolderName + argv[3];
			break;
		default:
			std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: CoregTool base_folder [images_folder] [transformations_folder]" << std::endl; 
			return EXIT_FAILURE;
			break;
	}

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	FileNameType pdName = imagesFolderName + "pd_corrected.nii.gz";
    ProbabilityImage pd = brainio->ReadProbabilityImage(pdName);
	FileNameType t2Name = imagesFolderName + "t2_corrected.nii.gz";
    ProbabilityImage t2 = brainio->ReadProbabilityImage(t2Name);
	FileNameType t1Name = imagesFolderName + "t1_corrected.nii.gz";
    ProbabilityImage t1 = brainio->ReadProbabilityImage(t1Name);
	FileNameType flairName = imagesFolderName + "flair_corrected.nii.gz";
    ProbabilityImage flair = brainio->ReadProbabilityImage(flairName);
	FileNameType affineName;

    #ifdef WIN32
	CreateDirectory(transformsFolderName.c_str(),NULL);
    #else
    mkdir(transformsFolderName.c_str(), 0777);
    #endif
	//Affine Registration Shinking&Smoothing Parameters

	/*--- Coregister ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t           Coregister          " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	affineName = transformsFolderName + "affineT1toSpace.tfm";
	AffineTransform transformT1 = brainio->ReadAffineTransform(affineName);
	if ((t1 != (ProbabilityImage)NULL) && (pd != (ProbabilityImage)NULL) && (transformT1 == (AffineTransform)NULL)) {
		transformT1 = BrainRegistration::MultiAffineRegistration(pd, t1);
		brainio->WriteAffineTransform(
			transformsFolderName + "affineT1toSpace.tfm",
			transformT1
		);
		ProbabilityImage t1Moved = BrainRegistration::ResampleImage<AffineTransformType>(pd, t1, transformT1);
		brainio->WriteProbabilityImage(
			imagesFolderName + "t1_registered.nii.gz",
			t1Moved
		);
	}
	else if ((t1 != (ProbabilityImage)NULL) && (flair != (ProbabilityImage)NULL) && (transformT1 == (AffineTransform)NULL)) {
		transformT1 = BrainRegistration::MultiAffineRegistration(flair, t1);
		brainio->WriteAffineTransform(
			transformsFolderName + "affineT1toSpace.tfm",
			transformT1
		);
		ProbabilityImage t1Moved = BrainRegistration::ResampleImage<AffineTransformType>(flair, t1, transformT1);
		brainio->WriteProbabilityImage(
			imagesFolderName + "t1_registered.nii.gz",
			t1Moved
		);
	}

	affineName = transformsFolderName + "affineFLAIRtoSpace.tfm";
	AffineTransform transformFLAIR = brainio->ReadAffineTransform(affineName);
	if ((flair != (ProbabilityImage)NULL) && (pd != (ProbabilityImage)NULL) && (transformFLAIR == (AffineTransform)NULL)) {
		AffineTransform transformFLAIR = BrainRegistration::MultiAffineRegistration(pd, flair);
		brainio->WriteAffineTransform(
			transformsFolderName + "affineFLAIRtoSpace.tfm",
			transformFLAIR
		);
		ProbabilityImage flairMoved = BrainRegistration::ResampleImage<AffineTransformType>(pd, flair, transformFLAIR);
		brainio->WriteProbabilityImage(
			imagesFolderName + "flair_registered.nii.gz",
			flairMoved
		);
	}

	delete(brainio);

    return EXIT_SUCCESS;

}
