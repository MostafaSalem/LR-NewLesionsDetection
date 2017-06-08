#ifndef GAUSSIANESTIMATOR2D_H
#define GAUSSIANESTIMATOR2D_H

#include <gaussianestimator.h>

class GaussianEstimator2D : public GaussianEstimator
{

public:
    GaussianEstimator2D();
    GaussianEstimator2D(Vector2Type mu_t, Matrix2x2Type sigma_t);
    GaussianEstimator2D(std::vector<ProbabilityImage> data);
    GaussianEstimator2D(std::vector<ProbabilityImage> data, ProbabilityImage pr);
    GaussianEstimator2D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold);
    GaussianEstimator2D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c);
    void SetData(std::vector<ProbabilityImage> data);
    void SetParameters(Vector2Type mu_t, Matrix2x2Type sigma_t);
    virtual void SetParameters(vnl_vector<ProbabilityPixelType> mu_t, vnl_matrix<ProbabilityPixelType> sigma_t);
    vnl_vector<ProbabilityPixelType> GetMu();
    vnl_matrix<ProbabilityPixelType> GetSigma();
    void UpdateParameters(ProbabilityImage pr, ProbabilityPixelType threshold = 0.0);
    ProbabilityImage GetCurrentLikelihood();
    ProbabilityImage GetCurrentLogLikelihood();
    ProbabilityImage GetMahalanobisMap();
    void PrintParameters();

private:

    Vector2Type mu;
    Matrix2x2Type sigma;
    MultiSpectral2Image data;
    ScalarToVector2Type::Pointer sc2vect2;

};

#endif // GAUSSIANESTIMATOR2D_H
