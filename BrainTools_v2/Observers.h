#include "itkImageRegistrationMethodv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkCommand.h"

template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand() {};

public:
  typedef   TRegistration      RegistrationType;
  typedef   RegistrationType * RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
  typedef   OptimizerType * OptimizerPointer;

  void Execute( itk::Object * object,
                const itk::EventObject & event) ITK_OVERRIDE
    {
    if( !(itk::MultiResolutionIterationEvent().CheckEvent( &event ) ) )
      {
      return;
      }
    RegistrationPointer registration =
      static_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer =  static_cast< OptimizerPointer >(
        registration->GetModifiableOptimizer() );

    unsigned int currentLevel = registration->GetCurrentLevel();
    typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
      registration->GetShrinkFactorsPerDimension( currentLevel );
    typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
      registration->GetSmoothingSigmasPerLevel();

    std::cout << "-------------------------------------" << std::endl;
    std::cout << " registration Transform = " <<  registration->GetObjectName() << std::endl;
    std::cout << " Current level = " << currentLevel << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    std::cout << std::endl;

/*    if ( registration->GetCurrentLevel() == 0 )
      {
      optimizer->SetLearningRate( 16.00 );
      optimizer->SetMinimumStepLength( 2.5 );
      }
    else
      {
      optimizer->SetLearningRate( optimizer->GetCurrentStepLength() );
      optimizer->SetMinimumStepLength(
        optimizer->GetMinimumStepLength() * 0.2 );
      }*/

    }

  void Execute(const itk::Object * , const itk::EventObject & ) ITK_OVERRIDE
    {
    return;
    }
};
//#############################################################################################################
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate(): m_CumulativeIterationIndex(0) {};

public:
  typedef   itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
  typedef   const OptimizerType *                               OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
  {
  Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
  {
  OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
  if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
    return;
    }
  std::cout << optimizer->GetCurrentIteration() << "   ";
  std::cout << optimizer->GetValue() << "   ";
  std::cout << optimizer->GetCurrentPosition() << "   ";
  std::cout << m_CumulativeIterationIndex++ << std::endl;
  }
private:
  unsigned int m_CumulativeIterationIndex;
};
//#############################################################################################################
template <typename TRegistration>
class DemonsCommandIterationUpdate : public itk::Command
{
public:
  typedef  DemonsCommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  DemonsCommandIterationUpdate() {};
  
public:

  void Execute(const itk::Object *, const itk::EventObject & ) ITK_OVERRIDE
    {
    std::cout << "Warning: The const Execute method shouldn't be called" << std::endl;
    }

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
       TRegistration * filter = static_cast<  TRegistration * >( caller );

       if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      if(filter)
        {
        //filter->SetMaximumRMSError(MaxRmsE[RmsCounter]);
        std::cout << filter->GetMetric() <<  "  RMS Change: " << filter->GetRMSChange() << std::endl;

         std::cout << "Level Tolerance=  "<<filter->GetMaximumRMSError ()<<std::endl;
    }

}
};
//#############################################################################################################
class CommandResolutionLevelUpdate : public itk::Command
{
public:
  typedef  CommandResolutionLevelUpdate   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  CommandResolutionLevelUpdate() {};

public:
  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
    Execute( (const itk::Object *)caller, event);
    }
  void Execute(const itk::Object *, const itk::EventObject & ) ITK_OVERRIDE
    {
    std::cout << "----------------------------------" << std::endl;
    //RmsCounter = RmsCounter + 1;
    std::cout << "----------------------------------" << std::endl;
    }
};

