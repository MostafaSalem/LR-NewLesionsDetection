#ifndef GAUSSIANESTIMATOR_H
#define GAUSSIANESTIMATOR_H

#include <itkMath.h>
#include <vnl/vnl_inverse.h>

#include <imagedefinitions.h>
using namespace itk::Math;

class GaussianEstimator
{

public:
    GaussianEstimator() {}
    virtual void SetData(std::vector<ProbabilityImage> data) = 0;
    //virtual void UpdateParameters(ProbabilityImage pr) = 0;
    virtual void UpdateParameters(ProbabilityImage pr,ProbabilityPixelType threshold = 0.0) = 0;
    virtual ProbabilityImage GetCurrentLikelihood() = 0;
    virtual ProbabilityImage GetCurrentLogLikelihood() = 0;
    virtual ProbabilityImage GetMahalanobisMap() = 0;
    virtual ProbabilityImage GetDistanceMap() {return ProbabilityImageType::New();};
    virtual void SetParameters(ProbabilityPixelType mu_t, ProbabilityPixelType sigma_t) {}
    virtual void SetParameters(Vector2Type mu_t, Matrix2x2Type sigma_t) {}
    virtual void SetParameters(Vector3Type mu_t, Matrix3x3Type sigma_t) {}
    virtual void SetParameters(vnl_vector<ProbabilityPixelType> mu_t, vnl_matrix<ProbabilityPixelType> sigma_t) {}
    virtual vnl_vector<ProbabilityPixelType> GetMu() {return vnl_vector<ProbabilityPixelType>();}
    virtual vnl_matrix<ProbabilityPixelType> GetSigma() {return vnl_matrix<ProbabilityPixelType>();}
    virtual void PrintParameters() = 0;
};

#endif // GAUSSIANESTIMATOR_H
