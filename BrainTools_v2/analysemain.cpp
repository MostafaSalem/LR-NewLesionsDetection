#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
    std::cout << "/-------------------------------/" << std::endl;
    std::cout << "|          Analysetool          |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;

    BrainIO *brainio = new BrainIO();

	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;

	if (argc <2) {
		std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: Analysetool folder" << std::endl; 
		return EXIT_FAILURE;
	}

	FileNameType folderName = argv[1];
	//------> PD
	FileNameType pdName = folderName + "pd_demons_positive_activity.nii.gz";
	MaskImage pd = brainio->ReadMaskImage(pdName);

	//------> T2
	FileNameType t2Name = folderName + "t2_demons_positive_activity.nii.gz";
	MaskImage t2 = brainio->ReadMaskImage(t2Name);

	//------> FLAIR
	FileNameType flairName = folderName + "flair_demons_positive_activity.nii.gz";
	MaskImage flair = brainio->ReadMaskImage(flairName);

	//------> JOINT
	FileNameType jointName = folderName + "joint_demons_positive_activity.nii.gz";
	MaskImage joint = brainio->ReadMaskImage(jointName);

	/*--- Analysis ---*/
	ConnectedPixelType pdLesions, t2Lesions, flairLesions, jointLesions;
	double pdVolume, t2Volume, flairVolume, jointVolume;
	pdLesions = BrainSegmentation::CountRegions(pd);
	t2Lesions = BrainSegmentation::CountRegions(t2);
	flairLesions = BrainSegmentation::CountRegions(flair);
	jointLesions = BrainSegmentation::CountRegions(joint);

	pdVolume = BrainSegmentation::ComputeVolume(pd);
	t2Volume = BrainSegmentation::ComputeVolume(t2);
	flairVolume = BrainSegmentation::ComputeVolume(flair);
	jointVolume = BrainSegmentation::ComputeVolume(joint);

	std::cout << "PD;" << pdLesions << ";" << pdVolume << ";" << std::endl;
	std::cout << "T2;" << t2Lesions << ";" << t2Volume << ";" << std::endl;
	std::cout << "FLAIR;" << flairLesions << ";" << flairVolume << ";" << std::endl;
	std::cout << "FINAL;" << jointLesions << ";" << jointVolume << ";" << std::endl;

	return EXIT_SUCCESS;
}
