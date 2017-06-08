#include <cstdlib>
#include <cstdio>

#include <brainio.h>
#include <brainpreprocessing.h>
#include <brainregistration.h>
#include <brainsegmentation.h>

using namespace std;

int main(int argc, char **argv)
{
		typedef itk::RelabelComponentImageFilter< ConnectedImageType, ConnectedImageType > RelabelComponentFilterType;

	BrainIO *brainio = new BrainIO();
	FileNameType followupFolderName;
				 followupFolderName = argv[1];            
	FileNameType inputImageName;
				 inputImageName = argv[2];
    FileNameType outputImageName;
				 outputImageName = argv[3];

	MaskImage lesionMask = brainio->ReadMaskImage(followupFolderName + inputImageName);

	ConnectComponentFilterType::Pointer connectComp = ConnectComponentFilterType::New();
    connectComp->SetInput(lesionMask);    
	RelabelComponentFilterType::Pointer relabelComp = RelabelComponentFilterType::New();
	relabelComp->SetInput(connectComp->GetOutput());
	relabelComp->Update();

	ConnectedImage rImage=relabelComp->GetOutput();
	brainio->WriteConnectedImage(followupFolderName + outputImageName,rImage);

	delete(brainio);

    return EXIT_SUCCESS;
}