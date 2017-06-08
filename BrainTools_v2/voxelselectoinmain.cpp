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
    std::cout << "|           Creating Images for Voxel Selectoins         |" << std::endl;
    std::cout << "/-------------------------------/" << std::endl;
    BrainIO *brainio = new BrainIO();


	/*--- Init ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t             Init              " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	// We check the input parameters	
	FileNameType followupFolderName;
	FileNameType subtractionFolderName,voxelSelectionFolderName;	

			followupFolderName = argv[1];            
            subtractionFolderName = followupFolderName + "subtraction/";
            voxelSelectionFolderName = followupFolderName + "voxelselection1.5/";
	
    #ifdef WIN32
    CreateDirectory(voxelSelectionFolderName.c_str(),NULL);
    #else
    mkdir(voxelSelectionFolderName.c_str(), 0777);
    #endif

	/*--- Image reading ---*/
	std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
	
	
	// WM Masked Subtraction 
	//ProbabilityImage gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	ProbabilityImage gaussWMMaskedSubtractionPD = brainio->ReadProbabilityImage(subtractionFolderName + "pd_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionPD == (ProbabilityImage)NULL)) 
	{
		std::cerr << " gaussWMMaskedSubtractionPD image is not found" << std::endl; 
		return EXIT_FAILURE;	
	}

	//ProbabilityImage gaussWMMaskedSubtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_wmmasked_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	ProbabilityImage gaussWMMaskedSubtractionT1 = brainio->ReadProbabilityImage(subtractionFolderName + "t1_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionT1 == (ProbabilityImage)NULL)) 
	{
		std::cerr << " gaussWMMaskedSubtractionT1 image is not found" << std::endl; 
		return EXIT_FAILURE;	
	}

	//ProbabilityImage gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	ProbabilityImage gaussWMMaskedSubtractionT2 = brainio->ReadProbabilityImage(subtractionFolderName + "t2_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionT2 == (ProbabilityImage)NULL)) 
	{
		std::cerr << " gaussWMMaskedSubtractionT2 image is not found" << std::endl; 
		return EXIT_FAILURE;	
	}

	//ProbabilityImage gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_0.5sigma_subtraction.nii.gz");
	//ProbabilityImage gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_1sigma_subtraction.nii.gz");
	ProbabilityImage gaussWMMaskedSubtractionFLAIR = brainio->ReadProbabilityImage(subtractionFolderName + "flair_wmmasked_smoothed_1.5sigma_subtraction.nii.gz");
	if ((gaussWMMaskedSubtractionFLAIR == (ProbabilityImage)NULL)) 
	{
		std::cerr << " gaussWMMaskedSubtractionFLAIR image is not found" << std::endl; 
		return EXIT_FAILURE;	
	}
	
	    MaskImage positiveActivityPD, positiveActivityT2, positiveActivityFLAIR;
		//WM Masking Voxel Selectoin
		positiveActivityPD = BrainSegmentation::GetVoxelSelection(gaussWMMaskedSubtractionPD,0);
		brainio->WriteMaskImage(voxelSelectionFolderName + "pd_wmmasked_voxel_selection_0sigma.nii.gz",positiveActivityPD);		

		positiveActivityT2 = BrainSegmentation::GetVoxelSelection(gaussWMMaskedSubtractionT2,0);
		brainio->WriteMaskImage(voxelSelectionFolderName + "t2_wmmasked_voxel_selection_0sigma.nii.gz",positiveActivityT2);		

		positiveActivityFLAIR = BrainSegmentation::GetVoxelSelection(gaussWMMaskedSubtractionFLAIR,0);
		brainio->WriteMaskImage(voxelSelectionFolderName + "flair_wmmasked_voxel_selection_0sigma.nii.gz",positiveActivityFLAIR);
		

	delete(brainio);

    return EXIT_SUCCESS;

}
