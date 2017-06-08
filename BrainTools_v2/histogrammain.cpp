#include <cstdlib>
#include <cstdio>
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
        std::cerr << "Incorrect number of parameters." << std::endl << "Correct usage: HistogramTool reference_image input_image" << std::endl;
        return EXIT_FAILURE;
    }

    /*--- Image reading ---*/
    std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t         Image reading         " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    // Get the names from Baseline and read the files (mask + image for baseline image and mask + image for follow-up)
    // First we search for the files in this order: brainmask, pd, t1, t2 and flair
    //-- Baseline
    FileNameType referenceImageName = argv[1];
    ProbabilityImage referenceImage = brainio->ReadProbabilityImage(referenceImageName);

    FileNameType inputImageName = argv[2];
    ProbabilityImage inputImage = brainio->ReadProbabilityImage(inputImageName);

    /*--- Histogram matching ---*/
    std::cout << "\t-------------------------------" << std::endl;
    std::cout << "\t       Histogram matching      " << std::endl;
    std::cout << "\t-------------------------------" << std::endl;
    // Histogram matching of all the images
    // If the images exist we do nothing
    FileNameType matchedImageName = inputImageName.substr(0,inputImageName.find(".")) + "_matched.nii.gz";
    ProbabilityImage matchedImage = BrainPreprocessing<ProbabilityImageType>::MatchHistogram(inputImage, referenceImage);
    brainio->WriteProbabilityImage(
        matchedImageName,
        matchedImage
    );

}
