#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainregistration.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
	BrainIO *brainio = new BrainIO();
	FileNameType followupFolderName;
				 followupFolderName = argv[1];            
            

	MaskImage lesionMask = brainio->ReadMaskImage(followupFolderName + "lesionMask.nii");

	ConnectComponentFilterType::Pointer connectComp = ConnectComponentFilterType::New();
    connectComp->SetInput(lesionMask);
    connectComp->Update();

    ConnectedImage cImage=connectComp->GetOutput();
	brainio->WriteConnectedImage(followupFolderName + "lesionMaskComponent.nii",cImage);

	RelabelComponentFilterType::Pointer relabelComp = RelabelComponentFilterType::New();
	relabelComp->SetInput(connectComp->GetOutput());
	relabelComp->Update();

	ConnectedImage rImage=relabelComp->GetOutput();
	brainio->WriteConnectedImage(followupFolderName + "lesionMaskLabel.nii",rImage);

	delete(brainio);

    return EXIT_SUCCESS;
}