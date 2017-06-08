#ifndef IMAGEDEFINITIONS_H
#define IMAGEDEFINITIONS_H

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkConnectedComponentImageFilter.h>

#include <itkVector.h>
#include <itkMatrix.h>
#include <itkVariableSizeMatrix.h>
#include <itkCompose2DVectorImageFilter.h>
#include <itkCompose3DVectorImageFilter.h>
#include <itkImageToVectorImageFilter.h>
#include <itkTranslationTransform.h>
#include <itkCompositeTransform.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkMemoryProbesCollectorBase.h>
#include <itkExtractImageFilter.h>
#include <itkCenteredTransformInitializer.h>
#include <itkVersorRigid3DTransform.h>
#include <itkAffineTransform.h>
#include <itkBSplineTransform.h>

#ifdef ITK_USE_REVIEW
#define itkProbesCreate()  \
  itk::TimeProbesCollectorBase chronometer; \
  itk::MemoryProbesCollectorBase memorymeter
#define itkProbesStart( text ) memorymeter.Start( text ); chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text ); memorymeter.Stop( text  )
#define itkProbesReport( stream )  chronometer.Report( stream ); memorymeter.Report( stream  )
#else
#define itkProbesCreate()  \
  itk::TimeProbesCollectorBase chronometer
#define itkProbesStart( text ) chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text )
#define itkProbesReport( stream )  chronometer.Report( stream )
#endif


// Transformation definitions
typedef itk::TranslationTransform< double, 3 > TranslateTransformType;
typedef TranslateTransformType::Pointer TranslateTransform;
typedef itk::VersorRigid3DTransform< double > RigidTransformType;
typedef RigidTransformType::Pointer RigidTransform;
typedef itk::AffineTransform< double, 3 > AffineTransformType;
typedef AffineTransformType::Pointer AffineTransform;
typedef itk::CompositeTransform< double,3 >  CompositeTransformType;
typedef CompositeTransformType::Pointer CompositeTransform;
typedef itk::BSplineTransform<double,3,3>    DeformableTransformType;
typedef DeformableTransformType::Pointer DeformableTransform;

// Floating images
typedef  float ProbabilityPixelType;
typedef itk::Image< ProbabilityPixelType, 3 > ProbabilityImageType;
typedef itk::Image< ProbabilityPixelType, 2 > ProbabilitySliceType;
typedef ProbabilityImageType::Pointer ProbabilityImage;
typedef ProbabilitySliceType::Pointer ProbabilitySlice;
typedef ProbabilityImageType::RegionType ProbabilityRegionType;
typedef ProbabilityImageType::SizeType ProbabilitySizeType;
typedef ProbabilityImageType::IndexType ProbabilityIndexType;
typedef ProbabilityImageType::DirectionType ProbabilityDirectionType;
typedef itk::ImageRegionIterator<ProbabilityImageType> ProbabilityIterator;
typedef itk::NeighborhoodIterator<ProbabilityImageType> ProbabilityNeighborhoodIterator;
typedef itk::ImageRegionIterator<ProbabilitySliceType> ProbabilitySliceIterator;
typedef itk::NeighborhoodIterator<ProbabilitySliceType> ProbabilitySliceNeighborhoodIterator;

typedef float SignedIntensityPixelType;
typedef itk::Image< SignedIntensityPixelType, 3 > SignedImageType;
typedef SignedImageType::Pointer SignedImage;
typedef SignedImageType::RegionType SignedRegionType;
typedef SignedImageType::SizeType SignedSizeType;
typedef SignedImageType::DirectionType SignedDirectionType;
typedef SignedImageType::SizeType::SizeValueType SignedSizeValueType;

// Integer images
typedef unsigned short IntensityPixelType;
typedef itk::Image< IntensityPixelType, 3 > IntensityImageType;
typedef itk::Image< IntensityPixelType, 2 > IntensitySliceType;
typedef IntensityImageType::Pointer IntensityImage;
typedef IntensitySliceType::Pointer IntensitySlice;
typedef IntensityImageType::RegionType IntensityRegionType;
typedef IntensityImageType::SizeType IntensitySizeType;
typedef IntensityImageType::SpacingType IntensitySpacingType;
typedef IntensityImageType::PointType IntensityOriginType;
typedef IntensityImageType::DirectionType IntensityDirectionType;
typedef itk::ImageRegionIterator<IntensityImageType> IntensityIterator;
typedef itk::NeighborhoodIterator<IntensityImageType> IntensityNeighborhoodIterator;
typedef itk::ImageRegionIterator<IntensitySliceType> IntensitySliceIterator;
typedef itk::NeighborhoodIterator<IntensitySliceType> IntensitySliceNeighborhoodIterator;

typedef unsigned char ResultsPixelType;
typedef itk::Image< ResultsPixelType, 3 > ResultsImageType;
typedef ResultsImageType::Pointer ResultsImage;
typedef ResultsImageType::RegionType ResultsRegionType;
typedef ResultsImageType::SizeType ResultsSizeType;
typedef ResultsImageType::DirectionType ResultsDirectionType;
typedef itk::ImageRegionIterator<ResultsImageType> ResultsIterator;
typedef itk::NeighborhoodIterator<ResultsImageType> ResultsNeighborhoodIterator;

typedef bool MaskPixelType;
typedef itk::Image< MaskPixelType, 3 > MaskImageType;
typedef itk::Image< MaskPixelType, 2 > MaskSliceType;
typedef MaskImageType::Pointer MaskImage;
typedef MaskImageType::RegionType MaskRegionType;
typedef MaskImageType::SizeType MaskSizeType;
typedef MaskImageType::DirectionType MaskDirectionType;
typedef MaskImageType::SizeType::SizeValueType MaskSizeValueType;
typedef MaskImageType::IndexType MaskIndexType;
typedef itk::ImageRegionIterator<MaskImageType> MaskIterator;
typedef itk::NeighborhoodIterator<MaskImageType> MaskNeighborhoodIterator;

// Image definitions (connected iterators)
typedef unsigned int ConnectedPixelType;
typedef itk::Image< ConnectedPixelType, 3 > ConnectedImageType;
typedef ConnectedImageType::Pointer ConnectedImage;
typedef itk::ImageRegionIterator<ConnectedImageType> ConnectedIterator;
typedef itk::NeighborhoodIterator<ConnectedImageType> ConnectedNeighborhoodIterator;
typedef itk::ConnectedComponentImageFilter< MaskImageType, ConnectedImageType, MaskImageType > ConnectComponentFilterType;


// Vector images
typedef itk::Vector< float, 3 > VectorPixelType; 
typedef itk::Image< VectorPixelType, 3 > DeformationFieldType; 
typedef DeformationFieldType::Pointer DeformationField; 
typedef itk::ImageRegionIterator<DeformationFieldType> DeformationIterator;

//Rigid Transform Initializer
typedef itk::CenteredTransformInitializer< RigidTransformType, IntensityImageType, IntensityImageType > TransformInitializerType;
typedef itk::CenteredTransformInitializer< RigidTransformType, ProbabilityImageType, ProbabilityImageType > PTransformInitializerType;
typedef PTransformInitializerType::Pointer PTransformInitializer;

// Multispectral images
typedef itk::Image< itk::Vector< ProbabilityPixelType, 2 > , 3 > MultiSpectral2ImageType;
typedef MultiSpectral2ImageType::Pointer MultiSpectral2Image;
typedef itk::ImageRegionIterator<MultiSpectral2ImageType> Multi2Iterator;
typedef itk::Image< itk::Vector< ProbabilityPixelType, 3 > , 3 > MultiSpectral3ImageType;
typedef MultiSpectral3ImageType::Pointer MultiSpectral3Image;
typedef itk::ImageRegionIterator<MultiSpectral3ImageType> Multi3Iterator;
typedef itk::Image< itk::VariableSizeMatrix< ProbabilityPixelType > , 3 > MultiSpectralImageType;
typedef MultiSpectralImageType::Pointer MultiSpectralImage;
typedef itk::ImageRegionIterator<MultiSpectralImageType> MultiIterator;

typedef itk::Vector< ProbabilityPixelType, 2 > Vector2Type;
typedef itk::Vector< ProbabilityPixelType, 3 > Vector3Type;
typedef itk::VariableLengthVector< ProbabilityPixelType > VectorType;
typedef itk::Matrix< ProbabilityPixelType, 2, 2 > Matrix2x2Type;
typedef itk::Matrix< ProbabilityPixelType, 3, 3 > Matrix3x3Type;
typedef itk::VariableSizeMatrix< ProbabilityPixelType > MatrixType;
typedef itk::Compose2DVectorImageFilter<ProbabilityImageType,MultiSpectral2ImageType> ScalarToVector2Type;
typedef itk::Compose3DVectorImageFilter<ProbabilityImageType,MultiSpectral3ImageType> ScalarToVector3Type;


// Slicers
typedef itk::ExtractImageFilter< MaskImageType, MaskSliceType > MaskSlicerType;
typedef itk::ExtractImageFilter< ProbabilityImageType, ProbabilitySliceType > ProbabilitySlicerType;
typedef itk::ExtractImageFilter< IntensityImageType, IntensitySliceType > IntensitySlicerType;

// Others
typedef std::string FileNameType;

//enum { MI_MATTES, MI, MI_HIST, NMI, CORR, NORM_CORR, GRADIENT_DIFF, KULLBACK_LEIBLER, KAPPA, MEAN_SQUARES, MEAN_SQUARES_HIST, MEAN_SQUARES_RECIPROCAL };


// Min/max definitions
template <class T> const T& min ( const T& a, const T& b ) {
  return (a<b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T, class Tt> const T& min ( const T& a, const Tt& b ) {
  return (a<b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T> const T& max ( const T& a, const T& b ) {
  return (a>b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}
template <class T, class Tt> const T& max ( const T& a, const Tt& b ) {
  return (a>b)?a:b;     // or: return comp(a,b)?a:b; for the comp version
}


#endif // IMAGEDEFINITIONS_H
