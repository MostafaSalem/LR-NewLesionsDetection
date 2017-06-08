#ifndef REGISTRATIONLEVELUPDATE_H
#define REGISTRATIONLEVELUPDATE_H

#include <itkCommand.h>
#include <itkRegularStepGradientDescentBaseOptimizer.h>
#include <itkBSplineDeformableTransform.h>
#include <itkScaleTransform.h>
#include <itkBSplineResampleImageFunction.h>
#include <itkBSplineDecompositionImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkImageToImageMetric.h>
#include <itkDemonsRegistrationFilter.h>

#include <imagedefinitions.h>

template <typename TRegistration>
class RegistrationLevelUpdate: public itk::Command {

protected:
    RegistrationLevelUpdate() { };

public:
    typedef RegistrationLevelUpdate Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro( Self );

    typedef TRegistration RegistrationType;
    typedef RegistrationType * RegistrationPointer;
    typedef itk::RegularStepGradientDescentBaseOptimizer OptimizerType;
    typedef OptimizerType::ScalesType OptimizerScalesType;
    typedef itk::BSplineDeformableTransform< double, 3, 3 >  BSplineTransformType;
    typedef itk::IdentityTransform< double, 3 > IdentityTransformType;
    typedef itk::ScaleTransform< double, 3 > ScaleTransformType;
    typedef BSplineTransformType::RegionType RegionType;
    typedef BSplineTransformType::OriginType OriginType;
    typedef BSplineTransformType::SpacingType SpacingType;
    typedef BSplineTransformType::ParametersType ParametersType;
    typedef BSplineTransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
    typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
    typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType> DecompositionType;
    typedef itk::ImageToImageMetric<SignedImageType,SignedImageType> MetricType;
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
    typedef OptimizerType * OptimizerPointer;
    typedef BSplineTransformType * BSplineTransformPointer;
    typedef const SignedImageType * SignedImagePointer;
	typedef itk::DemonsRegistrationFilter< ProbabilityImageType, ProbabilityImageType, DeformationFieldType >   RegistrationFilterType;

    void Execute(itk::Object * object, const itk::EventObject & event)
    {
        int level;/*, gridnodesxdir;
        unsigned int dim;*/
        OptimizerPointer optimizer;

        if( !(itk::IterationEvent().CheckEvent( &event )) ) {
            return;
        }

        if (!std::string("MultiResolutionImageRegistrationMethod").compare(object->GetNameOfClass())) {
            RegistrationPointer registration;
            // We get the multi-resolution registration pointer and check the pyramid level
            registration = dynamic_cast<RegistrationPointer>( object );
            //level = registration->GetCurrentLevel();
            level = ((RegistrationPointer)object)->GetCurrentLevel();

            //std::cout << "Registration level " << level << " using " << registration->GetTransform()->GetNameOfClass() << std::endl;
            std::cout << "\t\\- Registration level " << level << " using " << ((RegistrationPointer)object)->GetTransform()->GetNameOfClass() << std::endl;

            // We also get the optimizer to change those parameters affected by the number of grid nodes
            optimizer = dynamic_cast< OptimizerPointer >(((RegistrationPointer)object)->GetOptimizer() );

            // We only need to interact with the BSpline transform, so we need to check the actual transform
            //if (!std::string("BSplineDeformableTransform").compare(registration->GetTransform()->GetNameOfClass())) {
            if (!std::string("BSplineDeformableTransform").compare(((RegistrationPointer)object)->GetTransform()->GetNameOfClass())) {

                // We get the fixed image and the transform
                SignedImagePointer fixed = dynamic_cast<SignedImagePointer>(registration->GetFixedImagePyramid()->GetOutput(level));
                BSplineTransformPointer bspline = dynamic_cast<BSplineTransformPointer>(registration->GetTransform());
                unsigned int num_parameters,gridnodesxdir,dim;

                // We double the nodes per grid at each level as the resolution increases and change the B-Splines parameters
                RegionType gridreg;
                RegionType::SizeType gridsize = bspline->GetGridRegion().GetSize();
                if (level>0) {
                    gridnodesxdir = gridsize[0]-bspline->SplineOrder;
                    gridsize.Fill(static_cast<long unsigned int>(gridnodesxdir*2+bspline->SplineOrder));
                    gridreg.SetSize(gridsize);
                }


                // We also update the origin and the spacing of the new grid
                SpacingType spacing = fixed->GetSpacing();
                OriginType  origin  = fixed->GetOrigin();
                SignedDirectionType griddir = fixed->GetDirection();
                SignedSizeType fixedsize = fixed->GetLargestPossibleRegion().GetSize();

                for(dim=0; dim<3; dim++) {
                    spacing[dim] *= static_cast<double>(fixedsize[dim] - 1) / static_cast<double>(gridsize[dim] - 1);
                }
                origin -= griddir * spacing;

                // We finally set the number of parameters
                num_parameters = gridsize[0] * gridsize[1] * gridsize[2] * 3;

                // We upsample the already computed grid parameters
                if (level>0) {
                    // We prepare the new parameters
                    ParametersType new_parameters( num_parameters );
                    new_parameters.Fill( 0.0 );

                    unsigned int counter = 0;

                    // We create the resampling objects
                    ResamplerType::Pointer upsampler = ResamplerType::New();
                    FunctionType::Pointer function = FunctionType::New();
                    DecompositionType::Pointer decomposition = DecompositionType::New();

                    ScaleTransformType::Pointer scale = ScaleTransformType::New();
                    ScaleTransformType::ScaleType scale_vector;
                    scale_vector.Fill(2.0);
                    scale->SetScale(scale_vector);

                    for ( dim = 0; dim < 3; dim++ ) {

                        // We prepare the grid upsampling
                        upsampler->SetInput( bspline->GetCoefficientImage()[dim] );
                        upsampler->SetInterpolator( function );
                        upsampler->SetTransform( scale );
                        upsampler->SetSize( gridsize );
                        upsampler->SetOutputSpacing( spacing );
                        upsampler->SetOutputOrigin( origin );

                        // We prepare the new grid
                        decomposition->SetSplineOrder( bspline->SplineOrder );
                        decomposition->SetInput( upsampler->GetOutput() );
                        decomposition->Update();

                        // We copy the new coefficients
                        ParametersImageType::Pointer new_coefficients = decomposition->GetOutput();
                        Iterator it( new_coefficients, bspline->GetGridRegion() );
                        while ( !it.IsAtEnd() ) {
                            new_parameters[ counter++ ] = it.Get()*2;
                            ++it;
                        }
                    }

                    // We set the final grid parameters
                    bspline->SetGridRegion( gridreg );
                    bspline->SetGridSpacing( spacing );
                    bspline->SetGridOrigin( origin );
					
                    bspline->SetParameters( new_parameters );
                    ((RegistrationPointer)object)->SetInitialTransformParametersOfNextLevel( new_parameters );
                    ((RegistrationPointer)object)->GetMetric()->SetTransformParameters( new_parameters );
					((RegistrationPointer)object)->SetTransform(bspline);
                }
				
                // Finally we prepare the optimizer and the metric
                OptimizerScalesType optimizer_scales = OptimizerScalesType( num_parameters );
                optimizer_scales.Fill( 1.0 );
                optimizer->SetScales( optimizer_scales );
            }
        }
        else if (std::string(object->GetNameOfClass()).find("Optimizer")) {

            // We get the optimizer to show registration progress
            optimizer = dynamic_cast< OptimizerPointer >(object);

            //std::cout << optimizer->GetNameOfClass() << " iteration: ";
            std::cout << "\t\\--- " << ((OptimizerPointer)object)->GetNameOfClass() << " iteration: ";
            //std::cout << optimizer->GetCurrentIteration();
            std::cout << ((OptimizerPointer)object)->GetCurrentIteration();
            //std::cout << " (" << optimizer->GetCostFunction()->GetNameOfClass() << " = ";
            std::cout << " (" << ((OptimizerPointer)object)->GetCostFunction()->GetNameOfClass() << " = ";
            //std::cout << -optimizer->GetValue() << ")";
            std::cout << -((OptimizerPointer)object)->GetValue() << ")";
            std::cout << std::endl;
        }
        else if (std::string(object->GetNameOfClass()).find("Demons")) {
            RegistrationFilterType * filter = dynamic_cast< RegistrationFilterType * >( object );
			std::cout << "\t\\--- " << filter->GetNameOfClass() << "iteration: ";
            std::cout << filter->GetElapsedIterations();
            std::cout << "(" << filter->GetMetric() << ")" << std::endl;
        }
    }
    void Execute(const itk::Object * , const itk::EventObject & ) { return; }

};

#endif // REGISTRATIONLEVELUPDATE_H
