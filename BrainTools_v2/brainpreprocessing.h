#ifndef BRAINPREPROCESSING_H
#define BRAINPREPROCESSING_H

#include <imagedefinitions.h>
#include <brainio.h>

#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkListSample.h>
#include <itkArray.h>
#include <itkVariableSizeMatrix.h>

#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

#include <itkHistogramMatchingImageFilter.h>
#include <itkMRIBiasFieldCorrectionFilter.h>
#include <itkN4BiasFieldCorrectionImageFilter.h>

#include <itkBinaryThresholdImageFilter.h>
#include <itkVotingBinaryHoleFillingImageFilter.h>

#include <itkMaskImageFilter.h>

template <class TImage>
class BrainPreprocessing
{
	// Denoising
	typedef itk::GradientAnisotropicDiffusionImageFilter< TImage, TImage > GradientAnisotropicDiffusionFilterType;
    typedef itk::CurvatureAnisotropicDiffusionImageFilter< TImage, TImage > CurvatureAnisotropicDiffusionFilterType;
	typedef itk::DiscreteGaussianImageFilter< TImage, TImage> GaussianFilterType;

	// Bias correction
    typedef itk::N4BiasFieldCorrectionImageFilter<TImage, MaskImageType, TImage> N4CorrectorType;
    typedef itk::N4BiasFieldCorrectionImageFilter<TImage, TImage> N4CorrectorType1;

	// Histogram matching
    typedef itk::HistogramMatchingImageFilter< TImage, TImage >   MatchingFilterType;

	// Masking
    typedef itk::MaskImageFilter< TImage, MaskImageType, TImage > MaskImageFilterType;

public:
	static typename TImage::Pointer GaussianFiltering(typename TImage::Pointer image, double sigma = 0.5) {
        typename GaussianFilterType::Pointer gaussianfilter = GaussianFilterType::New();
		gaussianfilter->SetInput(image);
		gaussianfilter->SetVariance(sigma*sigma);
		std::cout << "\t---------------------------" << std::endl;
		std::cout << "\tDenoising image (gaussian)" << std::endl;
		std::cout << "\t\\-- Sigma = " << gaussianfilter->GetVariance() << std::endl;
		gaussianfilter->Update();
		std::cout << "\t---------------------------" << std::endl;
		return(gaussianfilter->GetOutput());
	}
	static typename TImage::Pointer GradientAnisotropicDiffusion(typename TImage::Pointer image, int nIter = 100) {
        typename GradientAnisotropicDiffusionFilterType::Pointer gradientfilter = GradientAnisotropicDiffusionFilterType::New();
		gradientfilter->SetInput(image);
		gradientfilter->SetNumberOfIterations(nIter);
		std::cout << "\t---------------------------" << std::endl;
		std::cout << "\tDenoising image (gradient)" << std::endl;
		std::cout << "\t\\-- " << gradientfilter->GetNumberOfIterations() << " iterations" << std::endl;
		std::cout << "\t\\-- Time step " << gradientfilter->GetTimeStep() << "" << std::endl;
		gradientfilter->Update();
		std::cout << "\t---------------------------" << std::endl;
		return(gradientfilter->GetOutput());
	}
	static typename TImage::Pointer CurvatureAnisotropicDiffusion(typename TImage::Pointer image, int nIter = 100) {
        typename CurvatureAnisotropicDiffusionFilterType::Pointer curvaturefilter = CurvatureAnisotropicDiffusionFilterType::New();
		curvaturefilter->SetInput(image);
		curvaturefilter->SetNumberOfIterations(nIter);
		std::cout << "\t---------------------------" << std::endl;
		std::cout << "\tDenoising image (curvature)" << std::endl;
		std::cout << "\t\\-- " << curvaturefilter->GetNumberOfIterations() << " iterations" << std::endl;
		std::cout << "\t\\-- Time step " << curvaturefilter->GetTimeStep() << "" << std::endl;
		curvaturefilter->Update();
		std::cout << "\t---------------------------" << std::endl;
		return(curvaturefilter->GetOutput());
	}
	static typename TImage::Pointer N4BiasCorrection(typename TImage::Pointer image, MaskImage mask, int nIter = 400, int nLevels = 1) {
        typename N4CorrectorType::Pointer corrector = N4CorrectorType::New();
		corrector->SetInput(image);
		std::cout << "\t---------------------------" << std::endl;
		corrector->SetNumberOfFittingLevels(nLevels);
        typename N4CorrectorType::VariableSizeArrayType iterations = corrector->GetMaximumNumberOfIterations();
		iterations.SetSize(nLevels);
		iterations.Fill(nIter);
		corrector->SetMaximumNumberOfIterations(iterations);
		std::cout << "\tEstimating bias field (N4)" << std::endl;
		std::cout << "\t\\-- " << corrector->GetMaximumNumberOfIterations()[0] <<  " iterations" << std::endl;
		corrector->SetMaskImage(mask);
		try {
			corrector->Update();
		}
		catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
			return image;
        }
		std::cout << "\t---------------------------" << std::endl;
		return(corrector->GetOutput());
	}
		static typename TImage::Pointer N4BiasCorrection(typename TImage::Pointer image, int nIter = 400, int nLevels = 1) {
        typename N4CorrectorType1::Pointer corrector = N4CorrectorType1::New();
		corrector->SetInput(image);
		std::cout << "\t---------------------------" << std::endl;
		corrector->SetNumberOfFittingLevels(nLevels);
        typename N4CorrectorType::VariableSizeArrayType iterations = corrector->GetMaximumNumberOfIterations();
		iterations.SetSize(nLevels);
		iterations.Fill(nIter);
		corrector->SetMaximumNumberOfIterations(iterations);
		std::cout << "\tEstimating bias field (N4)" << std::endl;
		std::cout << "\t\\-- " << corrector->GetMaximumNumberOfIterations()[0] <<  " iterations" << std::endl;
		try {
			corrector->Update();
		}
		catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
			return image;
        }
		std::cout << "\t---------------------------" << std::endl;
		return(corrector->GetOutput());
	}
	static typename TImage::Pointer MatchHistogram(typename TImage::Pointer image, typename TImage::Pointer reference, int histogramLevels = 1024, int matchPoints = 7, bool meanIntensityOn = true) {
        typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();
		matcher->SetInput( image );
		matcher->SetReferenceImage( reference );
		matcher->SetNumberOfHistogramLevels( histogramLevels );
		matcher->SetNumberOfMatchPoints( matchPoints );
		if (meanIntensityOn)
			matcher->ThresholdAtMeanIntensityOn();
		matcher->Update();

		return(matcher->GetOutput());
	}
	static typename TImage::Pointer MaskImage(typename TImage::Pointer image, MaskImage mask) {
		std::cout << "\tBrainPreprocessing::MaskImage" << std::endl;
        typename MaskImageFilterType::Pointer maskFilter = MaskImageFilterType::New();
		maskFilter->SetInput(image);
		maskFilter->SetInput2(mask);
		maskFilter->Update();

		return(maskFilter->GetOutput());
	}

};

#endif // BRAINPREPROCESSING_H
