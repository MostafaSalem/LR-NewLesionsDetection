#ifndef BRAINSEGMENTATION_H
#define BRAINSEGMENTATION_H


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cfloat>

#include <vnl/vnl_math.h>
#include <imagedefinitions.h>
#include <gaussianestimator1d.h>
#include <gaussianestimator2d.h>
#include <gaussianestimator3d.h>
#include <gaussianestimatornd.h>
#include <brainio.h>

#include <itkMinimumMaximumImageCalculator.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAndImageFilter.h>
#include <itkOrImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkNeighborhoodAllocator.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMultiplyByConstantImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkNotImageFilter.h>
#include <itkSubtractImageFilter.h>

class BrainSegmentation
{

    typedef itk::MinimumMaximumImageCalculator< ProbabilityImageType > MinimumMaximumImageCalculatorType;
    typedef itk::MinimumMaximumImageCalculator< ResultsImageType > MinimumMaximumLabelCalculatorType;
    typedef itk::MinimumMaximumImageCalculator< ConnectedImageType > MinimumMaximumLabelComponentType;

	typedef itk::RelabelComponentImageFilter< ConnectedImageType, ConnectedImageType > RelabelComponentFilterType;
    typedef itk::Statistics::ScalarImageToHistogramGenerator< IntensityImageType > HistogramGeneratorType;
    typedef HistogramGeneratorType::HistogramType  HistogramType;
    typedef itk::Statistics::ScalarImageToHistogramGenerator< ProbabilityImageType > PHistogramGeneratorType;
    typedef PHistogramGeneratorType::HistogramType  PHistogramType;
    typedef itk::Statistics::ScalarImageToHistogramGenerator< ConnectedImageType > CHistogramGeneratorType;
    typedef HistogramGeneratorType::HistogramType  CHistogramType;

    typedef itk::BinaryThresholdImageFilter< ProbabilityImageType, MaskImageType > ThresholdingFilterType;
    typedef itk::AndImageFilter< MaskImageType, MaskImageType, MaskImageType > AndFilterType;
	typedef itk::OrImageFilter< MaskImageType, MaskImageType, MaskImageType > OrFilterType;
    typedef itk::MultiplyImageFilter< ProbabilityImageType, ProbabilityImageType > MultiplyFilterType;
	typedef itk::DivideImageFilter< ProbabilityImageType, ProbabilityImageType, ProbabilityImageType > DivideFilterType;
    typedef itk::AddImageFilter< ProbabilityImageType, ProbabilityImageType > AddFilterType;
	typedef itk::NotImageFilter< MaskImageType, MaskImageType > NotFilterType;
    typedef itk::MultiplyByConstantImageFilter< ProbabilityImageType, ProbabilityPixelType, ProbabilityImageType > MultiplyConstantFilterType;
	typedef itk::SubtractImageFilter< ProbabilityImageType, ProbabilityImageType, ProbabilityImageType > SubtractionFilterType;

    typedef itk::NeighborhoodAllocator< ProbabilityPixelType > NeighborhoodAllocatorType;
    typedef itk::BinaryBallStructuringElement< ProbabilityPixelType, 3, NeighborhoodAllocatorType >  BinaryBallStructuringElementType;
    typedef itk::GrayscaleDilateImageFilter< ProbabilityImageType, ProbabilityImageType, BinaryBallStructuringElementType > GrayscaleDilateFilterType;
    typedef itk::GrayscaleErodeImageFilter< ProbabilityImageType, ProbabilityImageType, BinaryBallStructuringElementType > GrayscaleErodeFilterType;


public:

    /* Constructors/Destructors */
    BrainSegmentation() {};
    ~BrainSegmentation() {};

    /* MS segmentation methods */
    static ResultsImage Relabel(ResultsImage tissue, MaskImage wml);
    static MaskImage Lesion_New(ResultsImage tissue, ProbabilityImage flair,double alpha=3, double omegaT = 0.8, double omegaN=0.3, int minSize=5);
    static MaskImage Lesion_New(ResultsImage tissue, ProbabilityImage flair,std::vector<int*> pv,double alpha=3, double omegaT = 0.8, double omegaN=0.3, int minSize=5);

    /* Generic algorithms */
    // Segmentation methods
    static std::vector<ProbabilityImage> Generic_EM_Atlas(std::vector<ProbabilityImage> atlas, std::vector<ProbabilityImage> input, MaskImage bm, std::vector<int*> pv, ProbabilityPixelType th = 0.5, int maxIter=25);
    static std::vector<ProbabilityImage> Generic_WEM_Atlas(std::vector<ProbabilityImage> atlas, std::vector<ProbabilityImage> input, ProbabilityImage weights, MaskImage bm, std::vector<int*> pv, ProbabilityPixelType th = 0.5, int maxIter=25);

    // Image generation methods
    static MaskImage GetWMMask(ResultsImage brain, ResultsPixelType wmlab = 3, ResultsPixelType wmllab = 5);
    static ResultsImage SolutionFromPr(std::vector<ProbabilityImage> pr, MaskImage bm);
    static MaskImage RefineMask(MaskImage init, std::vector<ProbabilityImage> atlas);
    static ProbabilityImage IntegralImage(ProbabilityImage source);
    static ProbabilityImage IntegralSquaredImage(ProbabilityImage source);
	static MaskImage Intersection(MaskImage maskA, MaskImage maskB);
	static MaskImage Intersection(std::vector<MaskImage> masks);
	static MaskImage Union(MaskImage maskA, MaskImage maskB);
	static ProbabilityImage Ratio(ProbabilityImage numerator, ProbabilityImage denominator);
	static ProbabilityImage Subtraction(ProbabilityImage minuend, ProbabilityImage subtrahend);
	static MaskImage ThresholdImage(ProbabilityImage input, ProbabilityPixelType threshold); 
	static MaskImage GetPositiveActivity(ProbabilityImage subtraction, double alpha = 5);
	static MaskImage GetVoxelSelection(ProbabilityImage subtraction, double alpha = 1);  
	static MaskImage OnurPostProcessing(MaskImage activity, ProbabilityImage basal, ProbabilityImage following, int minSize = 4, double bnr_t = 0.95, double fnr_t = 1);
	static MaskImage DeformationPostProcessing(MaskImage activity, ProbabilityImage jacobian, ProbabilityImage divergence, int minSize = 4, double jacobian_t = 1, double divergence_t = 0);
	static ProbabilityImage GetCorrelation(std::vector<ProbabilityImage> images, unsigned long radiusSize=1);

	// Results analysis
	static ConnectedPixelType CountRegions(MaskImage mask);
	static double ComputeVolume(MaskImage mask);

private:
    /* Generic functions */
    // Contrast enhancement
    static ProbabilityImage ContrastEnhancement(ProbabilityImage input, unsigned long radius=1);
    static ProbabilityImage NeighborhoodAverage(ProbabilityImage input, unsigned long radiusSize=1);

};



#endif // BRAINSEGMENTATION_H
