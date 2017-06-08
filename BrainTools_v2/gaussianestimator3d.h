#ifndef GAUSSIANESTIMATOR3D_H
#define GAUSSIANESTIMATOR3D_H

#include <gaussianestimator.h>

class GaussianEstimator3D : public GaussianEstimator
{

public:
    GaussianEstimator3D();
    GaussianEstimator3D(Vector3Type mu_t, Matrix3x3Type sigma_t);
    GaussianEstimator3D(std::vector<ProbabilityImage> data);
    GaussianEstimator3D(std::vector<ProbabilityImage> data, ProbabilityImage pr);
    GaussianEstimator3D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold);
    GaussianEstimator3D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c);
    void SetData(std::vector<ProbabilityImage> data);
    void SetParameters(Vector3Type mu_t, Matrix3x3Type sigma_t);
    virtual void SetParameters(vnl_vector<ProbabilityPixelType> mu_t, vnl_matrix<ProbabilityPixelType> sigma_t);
    vnl_vector<ProbabilityPixelType> GetMu();
    vnl_matrix<ProbabilityPixelType> GetSigma();
    void UpdateParameters(ProbabilityImage pr,ProbabilityPixelType threshold = 0.0);
    ProbabilityImage GetCurrentLikelihood();
    ProbabilityImage GetCurrentLogLikelihood();
    ProbabilityImage GetMahalanobisMap();
    void PrintParameters();

private:

    Vector3Type mu;
    Matrix3x3Type sigma;
    MultiSpectral3Image data;
    ScalarToVector3Type::Pointer sc3vect3;

};

#endif // GAUSSIANESTIMATOR3D_H
