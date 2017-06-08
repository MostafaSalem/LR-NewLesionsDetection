#include <brainregistration.h>
#include <itkMultiResolutionImageRegistrationMethod.h>


int  BrainRegistration::numberOfLevels=3;
int  BrainRegistration::steps=300;
float BrainRegistration::shinrkingFactors[3]={3,2,1};
float BrainRegistration::smoothingSigmas[3]={2,1,0};


BrainRegistration::BrainRegistration()
{
int num_Levels=3; int num_Steps=300;
float shrinking[3]={3,2,1};float smmothing[3]={2,1,0};
initialize(num_Levels,num_Steps,num_Levels,num_Steps,shrinking,smmothing,shrinking,smmothing); 
}
BrainRegistration::BrainRegistration(int rigidLevels,int rigidSteps,int affineLevels,int affineSteps,float* rigidShrinking,float*rigidSmoothing,float* affineShrinking,float*affineSmoothing) 
{
    initialize(rigidLevels, rigidSteps, affineLevels, affineSteps, rigidShrinking, rigidSmoothing, affineShrinking, affineSmoothing);
}
   
void BrainRegistration::initialize(int rigidLevels,int rigidSteps,int affineLevels,int affineSteps,float* rigidShrinking,float*rigidSmoothing,float* affineShrinking,float*affineSmoothing) {
    std::cout<<"Argument Constructor Constructor"<<std::endl;
    rigid = false;
    affine = false;
    bspline = false;

    this->rigidlevels = rigidLevels;
    this->affinelevels = affineLevels;
    this->rigidSteps = rigidSteps;
    this->affineSteps = affineSteps;
    bsplinenodes = 0;

    //Default values for Shrinking&Smoothing Factors
    //rigidShrinkingFactors = new float[rigidLevels]; 
    rigidShrinkFactorsPerLevel.SetSize( rigidLevels );    
    rigidSmoothingSigmasPerLevel.SetSize( rigidLevels );
    affineShrinkFactorsPerLevel.SetSize( affineLevels );    
    affineSmoothingSigmasPerLevel.SetSize( affineLevels );
    //rigidSmoothingFactors= new float[rigidLevels];
    for(int i=0;i<rigidLevels;i++)
    {
        rigidShrinkFactorsPerLevel[i]=rigidShrinking[i];
        rigidSmoothingSigmasPerLevel[i]=rigidSmoothing[i];
    }    
    for(int i=0;i<affineLevels;i++)
    {
        affineShrinkFactorsPerLevel[i]=affineShrinking[i];
        affineSmoothingSigmasPerLevel[i]=affineSmoothing[i];
    }


    // Transforms init
    idtransf = IdentityTransformType::New();

    rigidTransform = RigidTransformType::New();
    axis[0] = 0.0;    axis[1] = 0.0;    axis[2] = 1.0;    rotation.Set( axis, 0.0 );    
    rigidTransform->SetRotation( rotation );

    affineTransform = AffineTransformType::New();

    bsplinetransf = DeformableTransformType::New();
    //Rigid Intializer
    initializer = PTransformInitializerType::New();
    initializer->SetTransform(rigidTransform);
    initializer->MomentsOn();  
    //Observers
    optimizer_observer = CommandIterationUpdate::New();        
    registrationRigid_observer = UnTemplateRegistrationObserverType::New();
    registrationAffine_observer = UnTemplateRegistrationObserverType::New();
    
    // Metric init    
    rigid_mattesmi_metric = PMattesMIMetricType::New();
    rigid_mattesmi_metric->SetNumberOfHistogramBins(50);
    affine_mattesmi_metric = PMattesMIMetricType::New();
    affine_mattesmi_metric->SetNumberOfHistogramBins(50);
    
    scalesEstimatorRigid = ScalesEstimatorType::New();
    scalesEstimatorRigid->SetMetric( rigid_mattesmi_metric );
    scalesEstimatorRigid->SetTransformForward( true );
    scalesEstimatorAffine = ScalesEstimatorType::New();
    scalesEstimatorAffine->SetMetric( affine_mattesmi_metric );
    scalesEstimatorAffine->SetTransformForward( true );
    // Optimizer init
    rigidOptimizer = OptimizerType::New();
    rigidOptimizer->AddObserver( itk::IterationEvent(), optimizer_observer );
    rigidOptimizer->SetScalesEstimator( scalesEstimatorRigid );
    rigidOptimizer->SetNumberOfIterations( rigidSteps );
    rigidOptimizer->SetLearningRate( 0.2 );
    rigidOptimizer->SetMinimumStepLength( 0.0001 );
    rigidOptimizer->SetReturnBestParametersAndValue(true);   


    affineOptimizer = OptimizerType::New();
    affineOptimizer->AddObserver( itk::IterationEvent(), optimizer_observer );
    affineOptimizer->SetScalesEstimator( scalesEstimatorAffine );
    affineOptimizer->SetNumberOfIterations( affineSteps );
    affineOptimizer->SetLearningRate( 0.2 );
    affineOptimizer->SetMinimumStepLength( 0.0001 );
    affineOptimizer->SetReturnBestParametersAndValue(true);    

    // Interpolator init
    interpolator = InterpolatorType::New();
    binterpolator = BSplineInterpolatorType::New();
    bintinterpolator = BSplineImageInterpolatorType::New();
    bprinterpolator = BSplineAtlasInterpolatorType::New();
    ninterpolator = NNInterpolatorType::New();
    nanainterpolator = NNAnatomicInterpolatorType::New();

    // Registration framework 
    rigidRegistration = UnTemplateRegistrationType::New();
    rigidRegistration->SetMetric(rigid_mattesmi_metric);   
    rigidRegistration->SetOptimizer(rigidOptimizer); 
    rigidRegistration->SetInitialTransform( rigidTransform );
    rigidRegistration->InPlaceOn();    
    rigidRegistration->AddObserver( itk::MultiResolutionIterationEvent(), registrationRigid_observer );    
    rigidRegistration->SetObjectName("RigidRegistration");
    rigidRegistration->SetNumberOfLevels( rigidlevels );
    rigidRegistration->SetSmoothingSigmasPerLevel( rigidSmoothingSigmasPerLevel );
    rigidRegistration->SetShrinkFactorsPerLevel( rigidShrinkFactorsPerLevel );

    affineRegistration = UnTemplateRegistrationType::New();
    affineRegistration->SetMetric(affine_mattesmi_metric);   
    affineRegistration->SetOptimizer(affineOptimizer); 
    affineRegistration->SetInitialTransform( affineTransform );
    affineRegistration->InPlaceOn();
    affineRegistration->AddObserver( itk::MultiResolutionIterationEvent(), registrationAffine_observer );    
    affineRegistration->SetObjectName("AffineRegistration");
    affineRegistration->SetNumberOfLevels( affinelevels );
    affineRegistration->SetSmoothingSigmasPerLevel( affineSmoothingSigmasPerLevel );
    affineRegistration->SetShrinkFactorsPerLevel( affineShrinkFactorsPerLevel );

    // Resampling init
    imresample = ResampleImageFilterType::New();
    imresample->SetTransform(this->idtransf);
    imresample->ReleaseDataFlagOn();
    atlasresample = ResampleAtlasFilterType::New();
    atlasresample->SetTransform(this->idtransf);
    atlasresample->ReleaseDataFlagOn();
    maskresample = ResampleMaskFilterType::New();
    maskresample->SetTransform(this->idtransf);
    maskresample->SetInterpolator(ninterpolator);
    maskresample->ReleaseDataFlagOn();
    anatomicresample = ResampleResultsFilterType::New();
    anatomicresample->SetTransform(this->idtransf);
    anatomicresample->SetInterpolator(nanainterpolator);
    anatomicresample->ReleaseDataFlagOn();   

    // Observer init
  //  observer = ObserverType::New();
  //  optimizer->AddObserver( itk::IterationEvent(), observer );
  //  versor_optimizer->AddObserver( itk::IterationEvent(), observer );
  //  multires_registration->AddObserver( itk::IterationEvent(), observer );

    // Casting init
    fixint2signed = IntensityToSignedCastType::New();
    movint2signed = IntensityToSignedCastType::New();
    transfint2signed = IntensityToSignedCastType::New();

    // Preprocessing init
    matcher = MatchingFilterType::New();

    //Others
    minmax = MinimumMaximumImageCalculatorType::New();
}
//#####################################################################################################################################
RigidTransform BrainRegistration::GetRigidTransform(){return this->rigidTransform;}
AffineTransform BrainRegistration::GetAffineTransform(){return this->affineTransform;}
//#####################################################################################################################################
void BrainRegistration::SetFixed(IntensityImage fix) {
    // Setting fixed image   
    this->fiximg=fix;

    // Casting to a signed type for registration    
    fixint2signed->SetInput(this->fiximg);
    fixint2signed->Update();
    
    // Registration parameters
    fixreg = this->fiximg->GetBufferedRegion();
    fixsize = this->fixreg.GetSize();
    origin = this->fiximg->GetOrigin();
    spacing = this->fiximg->GetSpacing();
    griddir = this->fiximg->GetDirection();

    //RigidInitializer    
    initializer->SetFixedImage(fixint2signed->GetOutput());
   
    //Registration Methods
    rigidRegistration->SetFixedImage(fixint2signed->GetOutput());
    affineRegistration->SetFixedImage(fixint2signed->GetOutput());    

    // Intensity thresholds
    minmax->SetImage(fiximg);
    minmax->Compute();
    //mattesmi_metric->SetFixedImageSamplesIntensityThreshold(minmax->GetMaximum()*0.05);


    // No registration done yet
    affine = false;
    rigid = false;
}
//####################################################################################################################
void BrainRegistration::SetMoving(IntensityImage mov) {
    std::cout<<"Setting Moving Image...."<<std::endl;
     // Setting moving image
    this->movimg=mov;

    // Retrieving image information
    movreg = this->fiximg->GetBufferedRegion();
    movsize = this->fixreg.GetSize();
    movorigin = this->fiximg->GetOrigin();
    movspacing = this->fiximg->GetSpacing();

    // Casting to a signed type for registration
    movint2signed->SetInput(this->movimg);
    movint2signed->Update();

    //RigidInitializer
    initializer->SetMovingImage(movint2signed->GetOutput());
    //Registration Methods    
    rigidRegistration->SetMovingImage(movint2signed->GetOutput());
    affineRegistration->SetMovingImage(movint2signed->GetOutput());
   
    // No registration done yet
    affine = false;
    rigid = false;
   
}
//#####################################################################################################################################
void BrainRegistration::SetMask(MaskImage mask) {
    this->maskimg=mask;
    this->maskresample->SetInput(this->maskimg);
}
//#####################################################################################################################################
void BrainRegistration::SetAnatomic(ResultsImage anatomic) {
    this->anaimg=anatomic;
    this->anatomicresample->SetInput(this->anaimg);
}
//#####################################################################################################################################
void BrainRegistration::AddAtlas(ProbabilityImage atlas) {
    this->movatlas.push_back(atlas);
}
//#####################################################################################################################################
void BrainRegistration::SetAtlases(std::vector<ProbabilityImage> atlases) {
    this->movatlas.clear();
    for (unsigned int i=0;i<atlases.size();i++) {
        this->movatlas.push_back(atlases[i]);
    }
}
//#####################################################################################################################################
void BrainRegistration::ClearAtlases() {
    this->movatlas.clear();
}
//#####################################################################################################################################
//#####################################################################################################################################
void BrainRegistration::RigidRegistration() {
    std::cout << "\tBrainRegistration::RigidRegistration" << std::endl;
    if (!rigid) 
    {
        this->bspline=false;
        this->affine=false;

    this->initializer->InitializeTransform();
    try {
        this->rigidRegistration->Update();
        this->rigid = true;
    	}
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return;        
    }    
    std::cout << "\t\\-- Rigid transform completed" << std::endl;
    std::cout << std::endl;
    }
    imresample->SetTransform(this->rigidTransform);
    atlasresample->SetTransform(this->rigidTransform);
    maskresample->SetTransform(this->rigidTransform);
    anatomicresample->SetTransform(this->rigidTransform);
    }
//#####################################################################################################################################
    void BrainRegistration::AffineRegistration() {
    std::cout << "\tBrainRegistration::AffineRegistration" << std::endl;
    if (!affine )
     {
        this->bspline=false;

         if (!rigid) 
         {
            this->RigidRegistration();
         }

        affineTransform->SetCenter(rigidTransform->GetCenter());
        affineTransform->SetTranslation(rigidTransform->GetTranslation());
        affineTransform->SetMatrix(rigidTransform->GetMatrix());

    try {
        this->affineRegistration->Update();
        affine = true;
        }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return;        
    }    
    std::cout << "\t\\-- Affine transform completed" << std::endl;
    std::cout << std::endl;
    }
    imresample->SetTransform(this->affineTransform);
    atlasresample->SetTransform(this->affineTransform);
    maskresample->SetTransform(this->affineTransform);
    anatomicresample->SetTransform(this->affineTransform);
    }

//#####################################################################################################################################
//#####################################################################################################################################
IntensityImage BrainRegistration::GetMoving() {

    // We resample the image to the last transform
    imresample->SetInput(this->movimg);
    imresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    imresample->SetOutputOrigin(this->fiximg->GetOrigin());
    imresample->SetOutputSpacing(this->fiximg->GetSpacing());
    imresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        imresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return imresample->GetOutput();

}
//#####################################################################################################################################
MaskImage BrainRegistration::GetMask() {

    // We resample the image to the last transform
    maskresample->SetInput(this->maskimg);
    maskresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    maskresample->SetOutputOrigin(this->fiximg->GetOrigin());
    maskresample->SetOutputSpacing(this->fiximg->GetSpacing());
    maskresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        maskresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return maskresample->GetOutput();
}
//#####################################################################################################################################
ResultsImage BrainRegistration::GetAnatomic() {

    // We resample the image to the last transform
    anatomicresample->SetInput(this->anaimg);
    anatomicresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    anatomicresample->SetOutputOrigin(this->fiximg->GetOrigin());
    anatomicresample->SetOutputSpacing(this->fiximg->GetSpacing());
    anatomicresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        anatomicresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    return anatomicresample->GetOutput();

}
//#####################################################################################################################################
IntensityImage BrainRegistration::GetFixed() {
    return(this->fiximg);
}
//#####################################################################################################################################
std::vector<ProbabilityImage> BrainRegistration::GetAtlases() {

    std::vector<ProbabilityImage> atlases;

    // We resample all probabilistic atlases to the last transform
    for (unsigned int i=0;i<this->movatlas.size();i++) {
        atlasresample = ResampleAtlasFilterType::New();
        if (bspline) atlasresample->SetTransform(this->bsplinetransf);
        else if (affine) atlasresample->SetTransform(this->affineTransform);
        else if (rigid) atlasresample->SetTransform(this->rigidTransform);
        else atlasresample->SetTransform(this->idtransf);
        atlasresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
        atlasresample->SetOutputOrigin(this->fiximg->GetOrigin());
        atlasresample->SetOutputSpacing(this->fiximg->GetSpacing());
        atlasresample->SetOutputDirection(this->fiximg->GetDirection());
        atlasresample->SetInput(this->movatlas[i]);
        try {
            atlasresample->Update();
            atlases.push_back(atlasresample->GetOutput());
        }
        catch( itk::ExceptionObject & err ) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return atlases;
        }
    }

    return atlases;
}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::GetSimilarity(unsigned long radiusSize) {
    unsigned int i;
    ProbabilityImage similarity = ProbabilityImageType::New();
    similarity->SetRegions( this->fiximg->GetLargestPossibleRegion() );
    similarity->SetSpacing( fiximg->GetSpacing() );
    similarity->SetOrigin( fiximg->GetOrigin() );
    similarity->SetDirection( fiximg->GetDirection() );
    similarity->Allocate();
    similarity->FillBuffer(0);
    similarity->Update();

    // We resample the image to the last transform
    imresample->SetInput(this->movimg);
    imresample->SetSize(this->fiximg->GetLargestPossibleRegion().GetSize());
    imresample->SetOutputOrigin(this->fiximg->GetOrigin());
    imresample->SetOutputSpacing(this->fiximg->GetSpacing());
    imresample->SetOutputDirection(this->fiximg->GetDirection());
    try {
        imresample->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return NULL;
    }

    // We prepare the nighborhood iterators for both images
    IntensityImage movedimg = imresample->GetOutput();
    IntensityNeighborhoodIterator::RadiusType radius;
    for (i = 0; i < IntensityImageType::ImageDimension; ++i)
        radius[i] = radiusSize;
    IntensityNeighborhoodIterator movIt = IntensityNeighborhoodIterator(radius,movedimg,movedimg->GetRequestedRegion());
    IntensityNeighborhoodIterator fixIt = IntensityNeighborhoodIterator(radius,fiximg,fiximg->GetRequestedRegion());
    ProbabilityIterator simIt = ProbabilityIterator(similarity,similarity->GetRequestedRegion());
    ProbabilityPixelType movMean, fixMean, ncc, fixStd, movStd;

    for (movIt.GoToBegin();!movIt.IsAtEnd();++movIt,++fixIt,++simIt) {
        if (movIt.GetCenterPixel()>0 && fixIt.GetCenterPixel()>0) {
            // Mean computation
            movMean = 0;
            fixMean = 0;
            movStd = 0;
            fixStd = 0;
            ncc = 0;
            for (i=0; i<movIt.Size(); i++) {
                movMean += movIt.GetPixel(i);
                fixMean += fixIt.GetPixel(i);
            }
            movMean = movMean/movIt.Size();
            fixMean = fixMean/fixIt.Size();

            // Normalised cross correlation computation
            for (i=0; i<movIt.Size(); i++) {
                fixStd += (fixIt.GetPixel(i) - fixMean)*(fixIt.GetPixel(i) - fixMean);
                movStd += (movIt.GetPixel(i) - movMean)*(movIt.GetPixel(i) - movMean);
                ncc += (movIt.GetPixel(i) - movMean)*(fixIt.GetPixel(i) - fixMean);
            }

            simIt.Set(fabs(ncc/(sqrt(fixStd)*sqrt(movStd))));
        }
        else
            simIt.Set(0);
    }

    return(similarity);
}
//#####################################################################################################################################
//#####################################################################################################################################
RigidTransform BrainRegistration::MultiRigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int numSteps) {
    return BrainRegistration::MultiRigidRegistration(fiximg, movimg, numberOfLevels, shinrkingFactors, smoothingSigmas , numSteps);
}
//#####################################################################################################################################    
RigidTransform BrainRegistration::MultiRigidRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int numberOfLevels, float* shinrkingFactors, float* smoothingSigmas ,int steps) {
std::cout << "\tBrainRegistration::RigidRegistration" << std::endl;
     /* Init */
     // Observer
    OptimizerIterationObserver optimizer_observer = CommandIterationUpdate::New();    
    UnTemplateRegistrationObserver registration_observer = UnTemplateRegistrationObserverType::New();    
    // Metric parameters
    PMattesMIMetricType::Pointer rigid_mattesmi_metric = PMattesMIMetricType::New();
    rigid_mattesmi_metric->SetNumberOfHistogramBins(50);    
    // Transformations
    RigidTransform rigidTransform  = RigidTransformType::New();
    VersorType rotation;
    VectorType axis;    
    // Optimizers
    Optimizer rigidOptimizer = OptimizerType::New();
    rigidOptimizer->AddObserver( itk::IterationEvent(), optimizer_observer );
    // Registration framework         
    UnTemplateRegistration rigidRegistration = UnTemplateRegistrationType::New();
    rigidRegistration->SetMovingImage(movimg);
    rigidRegistration->SetFixedImage(fiximg);
    rigidRegistration->SetMetric(rigid_mattesmi_metric);   
    rigidRegistration->SetOptimizer(rigidOptimizer);         
    rigidRegistration->AddObserver( itk::MultiResolutionIterationEvent(), registration_observer );    
    rigidRegistration->SetObjectName("RigidRegistration");
    
    /* Rigid */    
    std::cout << std::endl;
    std::cout << "\t\\-- Starting Rigid transform" << std::endl;
    std::cout << std::endl;
    // Transform parameters    
    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;
    rotation.Set( axis, 0.0 );
    rigidTransform->SetRotation( rotation );
    PTransformInitializer initializer = PTransformInitializerType::New();
    initializer->SetFixedImage(fiximg);
    initializer->SetMovingImage(movimg);
    initializer->SetTransform(rigidTransform);
    initializer->MomentsOn();    
    initializer->InitializeTransform();
    
    rigidRegistration->SetInitialTransform( rigidTransform );
    rigidRegistration->InPlaceOn();  

    // Optimizer parameters
    ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( rigid_mattesmi_metric );
    scalesEstimator->SetTransformForward( true );

    rigidOptimizer->SetScalesEstimator( scalesEstimator );
    rigidOptimizer->SetNumberOfIterations( steps );
    rigidOptimizer->SetLearningRate( 0.2 );
    rigidOptimizer->SetRelaxationFactor(0.5);    
    rigidOptimizer->SetMinimumStepLength( 0.0001 );    
    rigidOptimizer->SetReturnBestParametersAndValue(true);    
    
    //Multiresolution Rigid    
    UnTemplateRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize( numberOfLevels );
    UnTemplateRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( numberOfLevels );

    for(int i=0;i<numberOfLevels;i++)
    {
        shrinkFactorsPerLevel[i] = shinrkingFactors[i];
        smoothingSigmasPerLevel[i] = smoothingSigmas[i];
    }

    rigidRegistration->SetNumberOfLevels( numberOfLevels );
    rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
    rigidRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

    // Registration    
    try {
        rigidRegistration->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(rigidTransform);
    } 
    std::cout << "\t\\-- Rigid transform completed" << std::endl;
    std::cout << std::endl;

    return(rigidTransform);
}
//#####################################################################################################################################
AffineTransform BrainRegistration::MultiAffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg,int numSteps) {
    return BrainRegistration::MultiAffineRegistration(fiximg, movimg, numberOfLevels, shinrkingFactors, smoothingSigmas , numSteps);
}
//#####################################################################################################################################    
AffineTransform BrainRegistration::MultiAffineRegistration(ProbabilityImage fiximg, ProbabilityImage movimg,  int numberOfLevels, float* shinrkingFactors, float* smoothingSigmas ,int steps) {
    std::cout << "\tBrainRegistration::AffineRegistration" << std::endl;
     /* Init */
     // Observer
    OptimizerIterationObserver optimizer_observer = CommandIterationUpdate::New();    
    UnTemplateRegistrationObserver registration_observer = UnTemplateRegistrationObserverType::New();        
    // Metric parameters    
    PMattesMIMetricType::Pointer affine_mattesmi_metric = PMattesMIMetricType::New();
    affine_mattesmi_metric->SetNumberOfHistogramBins(50);    
    // Transformations    
    AffineTransform affineTransform  = AffineTransformType::New();    
    RigidTransform rigidTransform  = BrainRegistration::MultiRigidRegistration(fiximg, movimg);
    // Optimizers    
    Optimizer affineOptimizer = OptimizerType::New();
    affineOptimizer->AddObserver( itk::IterationEvent(), optimizer_observer );
    // Registration framework         
    UnTemplateRegistration affineRegistration = UnTemplateRegistrationType::New();
    affineRegistration->SetMovingImage(movimg);
    affineRegistration->SetFixedImage(fiximg);
    affineRegistration->SetMetric(affine_mattesmi_metric);   
    affineRegistration->SetOptimizer(affineOptimizer);         
    affineRegistration->AddObserver( itk::MultiResolutionIterationEvent(), registration_observer );    
    affineRegistration->SetObjectName("AffineRegistration");

    /* Affine */    
    std::cout << std::endl;
    std::cout << "\t\\-- Starting affine transform" << std::endl;
    std::cout << std::endl;
    // Transform parameters
    affineTransform->SetCenter(rigidTransform->GetCenter());
    affineTransform->SetTranslation(rigidTransform->GetTranslation());
    affineTransform->SetMatrix(rigidTransform->GetMatrix());
    
    affineRegistration->SetInitialTransform( affineTransform );
    affineRegistration->InPlaceOn();  

    // Optimizer parameters
    ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( affine_mattesmi_metric );
    scalesEstimator->SetTransformForward( true );

    affineOptimizer->SetScalesEstimator( scalesEstimator );
    affineOptimizer->SetNumberOfIterations( steps );
    affineOptimizer->SetLearningRate( 0.2 );
    affineOptimizer->SetRelaxationFactor(0.5);    
    affineOptimizer->SetMinimumStepLength( 0.0001 );    
    affineOptimizer->SetReturnBestParametersAndValue(true);    
    
    //Multiresolution Affine    
    UnTemplateRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize( numberOfLevels );
    UnTemplateRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( numberOfLevels );

    for(int i=0;i<numberOfLevels;i++)
    {
        shrinkFactorsPerLevel[i] = shinrkingFactors[i];
        smoothingSigmasPerLevel[i] = smoothingSigmas[i];
    }

    affineRegistration->SetNumberOfLevels( numberOfLevels );
    affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
    affineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

    // Registration    
    try {
        affineRegistration->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(affineTransform);
    }    
    std::cout << "\t\\-- Affine transform completed" << std::endl;
    std::cout << std::endl;

    return(affineTransform);
}
//#####################################################################################################################################
void BrainRegistration::BSplineRegistration(unsigned int ngridnodesdim,int steps) {
   
}
//#####################################################################################################################################
void BrainRegistration::BSplineMultiRegistration(int levels, unsigned int ngridnodesdim,int steps) {
	
}

//#####################################################################################################################################
void BrainRegistration::ResampleBSpline() {
   

}
//#####################################################################################################################################
DeformableTransform BrainRegistration::BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, AffineTransform affinetransf, int levels, unsigned int ngridnodesdim, int steps) {
	std::cout << "\tBrainRegistration::BSplineMultiRegistration" << std::endl;
	
    return(NULL);
}
//#####################################################################################################################################
DeformableTransform BrainRegistration::BSplineMultiRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, unsigned int numberOfGridNodesInOneDimension, int steps) {
	/* Init */
    const unsigned int ImageDimension = 3;
    const unsigned int SpaceDimension = ImageDimension;
    const unsigned int SplineOrder = 3;
     // Observer
    OptimizerIterationObserver optimizer1_observer = CommandIterationUpdate::New();    
    UnTemplateRegistrationObserver registration1_observer = UnTemplateRegistrationObserverType::New();    
    // Metric parameters
    PMSMetricType::Pointer ms_metric = PMSMetricType::New();    
    // Transformations
    DeformableTransform outputBSplineTransform = DeformableTransformType::New();
    BSplineInitializer transformInitializer = BSplineInitializerType::New();
    // Optimizers
    LBFGSBOptimizer optimizer = LBFGSBOptimizerType::New();    
    // Registration framework         
    UnTemplateRegistration registration = UnTemplateRegistrationType::New();
    registration->SetMovingImage(movimg);
    registration->SetFixedImage(fiximg);
    registration->SetMetric(ms_metric);   
    registration->SetOptimizer(optimizer); 
    registration->SetInitialTransform( outputBSplineTransform );
    registration->InPlaceOn();
    registration->AddObserver( itk::MultiResolutionIterationEvent(), registration1_observer );
    registration->SetObjectName("BSplineRegistration");
    //Initalize the BSpline transform
    DeformableTransformType::MeshSizeType     meshSize;
    meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder ); 
    
    transformInitializer->SetTransform( outputBSplineTransform );
    transformInitializer->SetImage( fiximg );
    transformInitializer->SetTransformDomainMeshSize( meshSize );
    transformInitializer->InitializeTransform();
    // Set transform to identity    
    const unsigned int numberOfParameters =
               outputBSplineTransform->GetNumberOfParameters();
    DeformableTransformType::ParametersType parameters( numberOfParameters );
    parameters.Fill( 0.0 );
    outputBSplineTransform->SetParameters( parameters );
    //Multiresolution BSpline
    const unsigned int numberOfLevels = 1;
    UnTemplateRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize( numberOfLevels );
    shrinkFactorsPerLevel[0] = 1;

    UnTemplateRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( numberOfLevels );
    smoothingSigmasPerLevel[0] = 0;
    
    registration->SetNumberOfLevels( numberOfLevels );
    registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
    registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
    //Using Adaptor
    UnTemplateRegistrationType::TransformParametersAdaptorsContainerType adaptors;
    // First, get fixed image physical dimensions
    DeformableTransformType::PhysicalDimensionsType      fixedPhysicalDimensions;
    for( unsigned int i=0; i< ImageDimension; i++ )
    {
    fixedPhysicalDimensions[i] = fiximg->GetSpacing()[i] *
    static_cast<double>(
      fiximg->GetLargestPossibleRegion().GetSize()[i] - 1 );
    }

    // Create the transform adaptors specific to B-splines
    for( unsigned int level = 0; level < numberOfLevels; level++ )
    {    
    PShrinkFilter shrinkFilter = PShrinkFilterType::New();
    shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
    shrinkFilter->SetInput( fiximg );
    shrinkFilter->Update();

    // A good heuristic is to double the b-spline mesh resolution at each level
    //
    DeformableTransformType::MeshSizeType requiredMeshSize;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      requiredMeshSize[d] = meshSize[d] << level;
      }

    BSplineAdaptor bsplineAdaptor = BSplineAdaptorType::New();
    bsplineAdaptor->SetTransform( outputBSplineTransform );
    bsplineAdaptor->SetRequiredTransformDomainMeshSize( requiredMeshSize );
    bsplineAdaptor->SetRequiredTransformDomainOrigin(
      shrinkFilter->GetOutput()->GetOrigin() );
    bsplineAdaptor->SetRequiredTransformDomainDirection(
      shrinkFilter->GetOutput()->GetDirection() );
    bsplineAdaptor->SetRequiredTransformDomainPhysicalDimensions(
      fixedPhysicalDimensions );

    adaptors.push_back( bsplineAdaptor.GetPointer() );
    }

    registration->SetTransformParametersAdaptorsPerLevel( adaptors );
    // Configure the LBFGSB Optimizer
    const unsigned int numParameters =
    outputBSplineTransform->GetNumberOfParameters();
    LBFGSBOptimizerType::BoundSelectionType boundSelect( numParameters );
    LBFGSBOptimizerType::BoundValueType upperBound( numParameters );
    LBFGSBOptimizerType::BoundValueType lowerBound( numParameters );

    boundSelect.Fill( LBFGSBOptimizerType::UNBOUNDED );
    upperBound.Fill( 0.0 );
    lowerBound.Fill( 0.0 );

    optimizer->SetBoundSelection( boundSelect );
    optimizer->SetUpperBound( upperBound );
    optimizer->SetLowerBound( lowerBound );

    optimizer->SetCostFunctionConvergenceFactor( 1e+12 );
    optimizer->SetGradientConvergenceTolerance( 1.0e-35 );
    optimizer->SetNumberOfIterations( steps );
    optimizer->SetMaximumNumberOfFunctionEvaluations( 500 );
    optimizer->SetMaximumNumberOfCorrections( 5 );      

    /* BSpline */    
    std::cout << std::endl;
    std::cout << "\t\\-- Starting MultiBSpline Registration" << std::endl;
    std::cout << std::endl;
// Registration    
    try {
        registration->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
      //  return(NULL);
    }
    std::cout << "\t\\-- BSpline transform completed" << std::endl;
    std::cout << std::endl;

    //return(outputBSplineTransform);
}

//#####################################################################################################################################
DeformationField BrainRegistration::DemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int steps) {
   std::cout << "\tBrainRegistration::DemonsRegistration" << std::endl;
    /* Init */
    DemonsRegistrationFilterType::Pointer registration = DemonsRegistrationFilterType::New();    
    DeformationField deformation = DeformationFieldType::New();

    // Parameters
    registration->SetFixedImage( fiximg );
    registration->SetMovingImage( movimg );
    registration->SetNumberOfIterations( steps );
    registration->SetStandardDeviations( 1.0 );
    //registration->AddObserver( itk::IterationEvent(), observer );

    // Registration
    try {
        registration->Update();
    }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(deformation);
    }
    std::cout << "\t\\-- Demons Registration completed" << std::endl;
    std::cout << std::endl;
    deformation = registration->GetOutput();

    return(deformation);

}
//#####################################################################################################################################
DeformationField BrainRegistration::MultiDemonsRegistration(ProbabilityImage fiximg, ProbabilityImage movimg, int levels, int steps) {
	std::cout << "\tBrainRegistration::MultiDemonsRegistration" << std::endl;
    /* Init */
    unsigned int *nIterations = new unsigned int[levels];
    for(int i=0; i<levels; i++) nIterations[i] = steps;        
  
    MultiDemonsRegistrationFilterType::Pointer registration = MultiDemonsRegistrationFilterType::New();
    DemonsRegistrationFilterType::Pointer demons = DemonsRegistrationFilterType::New();    
    DeformationField deformation = DeformationFieldType::New();
    // Parameters
    demons->SetStandardDeviations( 1.0 );

    registration->SetRegistrationFilter( demons );
    registration->SetNumberOfLevels( levels  );    
    registration->SetFixedImage( fiximg );
    registration->SetMovingImage( movimg );
    registration->SetNumberOfIterations( nIterations );
    //Observers  
    DemonsRegistrationFilterObserver demonsLevelObserver = DemonsRegistrationFilterObserverType::New();
  demons->AddObserver( itk::IterationEvent(), demonsLevelObserver );

  MultiDemonsRegistrationFilterObserver multiDemonsLevelObserver = MultiDemonsRegistrationFilterObserverType::New();
  registration->AddObserver( itk::IterationEvent(), multiDemonsLevelObserver );
    
    // Registration
    try {        
        registration->Update();        
        }
    catch( itk::ExceptionObject & err ) {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return(deformation);
    }
    std::cout << "\t\\-- MultiDemons Registration completed" << std::endl;
    std::cout << std::endl;

    deformation = registration->GetOutput();

    //delete [] nIterations;

    return(deformation);
}
///#####################################################################################################################################
//#####################################################################################################################################
ProbabilityImage BrainRegistration::Subtraction(ProbabilityImage minuend, ProbabilityImage subtrahend, RigidTransform transformation) {
	std::cout << "\tBrainRegistration::Subtraction" << std::endl;
	// Init
	RigidTransformType::OutputVectorType half_translation;
	RigidTransformType::OutputVectorType inverse_half_translation;
	VersorType half_rotation;
	VersorType inverse_half_rotation;
	RigidTransformType::OutputVectorType::ValueType tx;
	RigidTransformType::OutputVectorType::ValueType ty;
	RigidTransformType::OutputVectorType::ValueType tz;
	VersorType::ValueType angle_radian;
	VersorType::VectorType versor_axis;
	double half_angle;
	double inverse_half_angle;


	// Transformation data
	std::cout << "\t\\- Getting the versor" << std::endl;
	VersorType versor =  transformation->GetVersor();
	AffineTransformType::OutputVectorType translation = transformation->GetTranslation();

	// Center
	std::cout << "\t\tCenter" << std::endl;
	BaseImageCalculatorType::Pointer m_BaseCalculator = BaseImageCalculatorType::New();
	m_BaseCalculator->SetImage(  minuend );		    
	m_BaseCalculator->Compute();
	BaseImageCalculatorType::VectorType baseCenter  = m_BaseCalculator->GetCenterOfGravity();

	AffineTransformType::CenterType centerFollowingImage;
	centerFollowingImage[0] = baseCenter[0];
	centerFollowingImage[1] = baseCenter[1];
	centerFollowingImage[2] = baseCenter[2];

	// Translation
	std::cout << "\t\tTranslation" << std::endl;
	tx = translation[0];
	ty = translation[1];
	tz = translation[2];

	half_translation[0] = tx/2;
	half_translation[1] = ty/2;
	half_translation[2] = tz/2;

	inverse_half_translation[0] = -half_translation[0];
	inverse_half_translation[1] = -half_translation[1];
	inverse_half_translation[2] = -half_translation[2];

	// Rotation
	std::cout << "\t\tRotation" << std::endl;
	angle_radian = versor.GetAngle();
	versor_axis = versor.GetAxis();
	VersorType::VectorType versor_axis_t;
	versor_axis_t[0] = versor_axis[0];
	versor_axis_t[1] = versor_axis[1];
	versor_axis_t[2] = versor_axis[2];
	std::cout << "\t\\--- Axis: " << versor_axis << std::endl;
	std::cout << "\t\\--- Angle: " << angle_radian << std::endl;
	std::cout << "\t\\--- Translation : " << half_translation << std::endl;
	half_angle = angle_radian/2;
	inverse_half_angle = -half_angle;
	

	half_rotation.Set(  versor_axis_t, half_angle  );
	inverse_half_rotation.Set( versor_axis_t,inverse_half_angle );		

	// Transforms
	std::cout << "\t\\- Transformations" << std::endl;
	RigidTransform halfRigid = RigidTransformType::New();
	halfRigid->SetCenter(transformation->GetCenter());
	halfRigid->SetRotation( half_rotation );
	halfRigid->SetTranslation( half_translation );
	
	RigidTransform inverseHalfRigid = RigidTransformType::New();
	inverseHalfRigid->SetCenter(transformation->GetCenter());
	inverseHalfRigid->SetRotation( inverse_half_rotation );
	inverseHalfRigid->SetTranslation(inverse_half_translation);

	// Halfway resample and subtraction
	std::cout << "\t\\- Halfway resample and subtraction" << std::endl;
	BSplineAtlasInterpolatorType::Pointer interpolator = BSplineAtlasInterpolatorType::New();

	ResampleAtlasFilterType::Pointer basalResample = ResampleAtlasFilterType::New();
	basalResample = ResampleAtlasFilterType::New();
	basalResample->SetInput(subtrahend);
	basalResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	basalResample->SetOutputOrigin(minuend->GetOrigin());
	basalResample->SetOutputSpacing(minuend->GetSpacing());
	basalResample->SetOutputDirection(minuend->GetDirection());
	basalResample->SetInterpolator(interpolator);
	basalResample->SetDefaultPixelValue( 0 );
	basalResample->SetTransform(halfRigid);
	try {
		basalResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	ResampleAtlasFilterType::Pointer followResample = ResampleAtlasFilterType::New();
	followResample = ResampleAtlasFilterType::New();
	followResample->SetInput(minuend);
	followResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	followResample->SetOutputOrigin(minuend->GetOrigin());
	followResample->SetOutputSpacing(minuend->GetSpacing());
	followResample->SetOutputDirection(minuend->GetDirection());
	followResample->SetInterpolator(interpolator);
	followResample->SetDefaultPixelValue( 0 );
	followResample->SetTransform(inverseHalfRigid);
	try {
		followResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	// We resample the subtraction
	std::cout << "\t\\- Subtraction resampling" << std::endl;
	SubtractionFilterType::Pointer subtractor = SubtractionFilterType::New();
	subtractor->SetInput1( followResample->GetOutput() );
	subtractor->SetInput2( basalResample->GetOutput() );
	subtractor->Update();

	ResampleAtlasFilterType::Pointer subResample = ResampleAtlasFilterType::New();
	subResample = ResampleAtlasFilterType::New();
	subResample->SetInput(subtractor->GetOutput());
	subResample->SetSize(minuend->GetLargestPossibleRegion().GetSize());
	subResample->SetOutputOrigin(minuend->GetOrigin());
	subResample->SetOutputSpacing(minuend->GetSpacing());
	subResample->SetOutputDirection(minuend->GetDirection());
	subResample->SetInterpolator(interpolator);
	subResample->SetDefaultPixelValue( 0 );
	subResample->SetTransform(halfRigid);
	try {
		subResample->Update();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	return(subResample->GetOutput());
}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::Warp(ProbabilityImage fiximg, ProbabilityImage movimg, DeformationField deformation) {

	WarperType::Pointer warper = WarperType::New(); 
	PInterpolatorType::Pointer interpolator = PInterpolatorType::New(); 
	
	warper->SetInput( movimg ); 
	warper->SetInterpolator( interpolator ); 
	warper->SetOutputSpacing( fiximg->GetSpacing() ); 
	warper->SetOutputOrigin( fiximg->GetOrigin() );
	warper->SetOutputDirection( movimg->GetDirection() );
	warper->SetDeformationField( deformation );
	warper->Update();

	return( warper->GetOutput() );
}
//#####################################################################################################################################
MaskImage BrainRegistration::Warp(MaskImage fiximg, MaskImage movimg, DeformationField deformation) {

	WarperMaskType::Pointer warper = WarperMaskType::New(); 
	NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New(); 
	
	warper->SetInput( fiximg ); 
	warper->SetInterpolator( interpolator ); 
	warper->SetOutputSpacing( movimg->GetSpacing() ); 
	warper->SetOutputOrigin( movimg->GetOrigin() );
	warper->SetOutputDirection( movimg->GetDirection() );
	warper->SetDeformationField( deformation );
	warper->Update();

	return( warper->GetOutput() );
}
//#####################################################################################################################################
DeformationField BrainRegistration::Sum(std::vector<DeformationField> deformations) {
	std::cout << "\tBrainRegistration::Sum" << std::endl;
	VectorPixelType zero;
	zero.Fill(0);
	DeformationField average = DeformationFieldType::New();
    average->SetRegions( deformations[0]->GetLargestPossibleRegion() );
    average->SetSpacing( deformations[0]->GetSpacing() );
    average->SetOrigin( deformations[0]->GetOrigin() );
    average->SetDirection( deformations[0]->GetDirection() );
    average->Allocate();
    average->FillBuffer(zero);
    average->Update();
	std::vector<DeformationField>::iterator it;
	MultiplyConstantFilterType::Pointer constFilter = MultiplyConstantFilterType::New();
	
	it = deformations.begin();
	std::cout << "\t\\-- Sum of the images" << std::endl;
	for (it = deformations.begin(); it != deformations.end(); ++it) {
		AddFilterType:: Pointer add = AddFilterType::New();
		add->SetInput1(average);
        add->SetInput2(*it);
        add->Update();
		average = add->GetOutput();
	}
	
	return(average);
}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::Jacobian(DeformationField deformation) {

	JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
	jacobianFilter->SetInput( deformation );
	jacobianFilter->Update();

	return( jacobianFilter->GetOutput() );
}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::Norm(DeformationField deformation) {

   /* Init */
    ProbabilityImage norm = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
    norm->FillBuffer(0);
    ProbabilityIterator itNorm = ProbabilityIterator(norm, norm->GetLargestPossibleRegion());

    DeformationIterator itDef = DeformationIterator(deformation, deformation->GetLargestPossibleRegion());
    for (itDef.GoToBegin(); !itDef.IsAtEnd(); ++itDef, ++itNorm) 
        itNorm.Set(  sqrt( itDef.Get()[0]*itDef.Get()[0] + itDef.Get()[1]*itDef.Get()[1] + itDef.Get()[2]*itDef.Get()[2])  );
        
    return( norm );

}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::Divergence(DeformationField deformation) {

	/* Init */
	ProbabilityImage divergence = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	divergence->FillBuffer(0);
	ProbabilityIterator itDiv = ProbabilityIterator(divergence, divergence->GetLargestPossibleRegion());
	

	/* Vector decomposition */
	ProbabilityImage vx = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vx->FillBuffer(0);
	ProbabilityImage vy = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vy->FillBuffer(0);
	ProbabilityImage vz = BrainIO::Initialise<DeformationFieldType, ProbabilityImageType>(deformation);
	vz->FillBuffer(0);

	DeformationIterator itDef = DeformationIterator(deformation, deformation->GetLargestPossibleRegion());
	ProbabilityIterator itVx = ProbabilityIterator(vx, vx->GetLargestPossibleRegion());
	ProbabilityIterator itVy = ProbabilityIterator(vy, vy->GetLargestPossibleRegion());
	ProbabilityIterator itVz = ProbabilityIterator(vz, vz->GetLargestPossibleRegion());

	for (itDef.GoToBegin(); !itDef.IsAtEnd(); ++itDef, ++itVx, ++itVy, ++itVz) {
		itVx.Set(itDef.Get()[0]);
		itVy.Set(itDef.Get()[1]);
		itVz.Set(itDef.Get()[2]);
	}

	/* Gradients */
	GradientFilterType::Pointer GradientVx = GradientFilterType::New();
	GradientVx->SetInput(vx);
	GradientVx->Update();
	GradientImage gx = GradientVx->GetOutput();
	GradientIterator itGx = GradientIterator(gx, gx->GetLargestPossibleRegion());
	GradientFilterType::Pointer GradientVy = GradientFilterType::New();
	GradientVy->SetInput(vy);
	GradientVy->Update();
	GradientImage gy = GradientVy->GetOutput();
	GradientIterator itGy = GradientIterator(gy, gy->GetLargestPossibleRegion());
	GradientFilterType::Pointer GradientVz = GradientFilterType::New();
	GradientVz->SetInput(vz);
	GradientVz->Update();
	GradientImage gz = GradientVz->GetOutput();
	GradientIterator itGz = GradientIterator(gz, gz->GetLargestPossibleRegion());

	for (itDiv.GoToBegin(); !itDiv.IsAtEnd(); ++itDiv, ++itGx, ++itGy, ++itGz)
		itDiv.Set(itGx.Get()[0]+itGy.Get()[1]+itGz.Get()[2]);

	return( divergence );
}
//#####################################################################################################################################
ProbabilityImage BrainRegistration::NormMulDivergence(ProbabilityImage norm, ProbabilityImage divergence) {

    MultiplyImageFilterType::Pointer multiplier = MultiplyImageFilterType::New();

    multiplier->SetInput1( norm );
    multiplier->SetInput2( divergence );
    multiplier->Update();

    return(multiplier->GetOutput());
}
//#####################################################################################################################################