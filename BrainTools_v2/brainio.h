#ifndef BRAINIO_H
#define BRAINIO_H

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>

#include <imagedefinitions.h>
#include <iostream>
#include <ctime>
#include <string>
#ifdef WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

class BrainIO
{

	// Intensity Images (T1 basically)
	typedef itk::ImageFileReader< IntensityImageType > IntensityReaderType;
	typedef itk::ImageFileWriter< IntensityImageType > IntensityWriterType;

	// Connected Images
	typedef itk::ImageFileReader< ConnectedImageType > ConnectedReaderType;
	typedef itk::ImageFileWriter< ConnectedImageType > ConnectedWriterType;

	// Probability atlases
	typedef itk::ImageFileReader< ProbabilityImageType > ProbabilityReaderType;
	typedef itk::ImageFileWriter< ProbabilityImageType > ProbabilityWriterType;

	// Segmentation results
	typedef itk::ImageFileReader< ResultsImageType > ResultsReaderType;
	typedef itk::ImageFileWriter< ResultsImageType > ResultsWriterType;

	// Deformation fields
	typedef itk::ImageFileReader< DeformationFieldType >  DeformationReaderType;
	typedef itk::ImageFileWriter< DeformationFieldType >  DeformationWriterType;

	// Transformation
	typedef itk::TransformFileReader TransformReaderType;
	typedef itk::TransformFileWriter TransformWriterType;

	// Masks
	typedef itk::ImageFileReader< MaskImageType > MaskReaderType;
	typedef itk::ImageFileWriter< MaskImageType > MaskWriterType;

	// Image caster
	typedef itk::CastImageFilter< MaskImageType, ResultsImageType> CastFilterType;

	// Connected components
	typedef itk::BinaryBallStructuringElement< MaskPixelType, 2> BallStructuringElementType;
	typedef itk::BinaryDilateImageFilter< MaskSliceType, MaskSliceType, BallStructuringElementType > DilateFilterType;
	typedef itk::ConnectedComponentImageFilter< MaskSliceType, IntensitySliceType, MaskSliceType > ConnectComponentFilterType;
	typedef itk::MinimumMaximumImageCalculator< IntensitySliceType > MinimumMaximumLabelComponentType;


public:
	BrainIO();

	// Read
	ProbabilityImage ReadProbabilityImage(FileNameType fileName);
	IntensityImage ReadIntensityImage(FileNameType fileName);
	ConnectedImage ReadConnectedImage(FileNameType fileName);
	ResultsImage ReadResultsImage(FileNameType fileName);
	MaskImage ReadMaskImage(FileNameType fileName);
	DeformationField ReadDeformationField(FileNameType fileName);
	AffineTransform ReadAffineTransform(FileNameType fileName);
	CompositeTransform ReadCompositeTransform(FileNameType fileName);

	// Write
	void WriteProbabilityImage(FileNameType fileName, ProbabilityImage fimage);
	void WriteIntensityImage(FileNameType fileName, IntensityImage intimage);
	void WriteConnectedImage(FileNameType fileName, ConnectedImage conimage);
	void WriteResultsImage(FileNameType fileName, ResultsImage resimage);
	void WriteMaskImage(FileNameType fileName, MaskImage bimage);
	void WriteDeformationField(FileNameType fileName, DeformationField dfimage);
	void WriteAffineTransform(FileNameType fileName, AffineTransform affine);
	void WriteRigidTransform(FileNameType fileName, RigidTransform rigidT);
	void WriteCompositeTransform(FileNameType fileName, CompositeTransform compT);

	// Filename manipulation
	static FileNameType SearchForFile(FileNameType folder, FileNameType toFind, bool verbose = false);
	static FileNameType SearchForFile(FileNameType folder, std::vector<FileNameType> toFindVector, bool verbose = false);

	template <class TInputImage, class TOutputImage>
	static typename TOutputImage::Pointer Initialise(typename TInputImage::Pointer input) {
        typename TOutputImage::Pointer output = TOutputImage::New();
		output->SetRegions( input->GetLargestPossibleRegion() );
		output->SetSpacing( input->GetSpacing() );
		output->SetOrigin( input->GetOrigin() );
		output->SetDirection( input->GetDirection() );
		output->Allocate();
		output->Update();

		return(output);
	}
};

#endif // BRAINIO_H
