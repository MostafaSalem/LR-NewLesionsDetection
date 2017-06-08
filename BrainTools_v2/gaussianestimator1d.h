#ifndef GAUSSIANESTIMATOR1D_H
#define GAUSSIANESTIMATOR1D_H

#include <gaussianestimator.h>

class GaussianEstimator1D : public GaussianEstimator
{
public:
    GaussianEstimator1D();
    GaussianEstimator1D(ProbabilityPixelType mu_t, ProbabilityPixelType sigma_t);
    GaussianEstimator1D(std::vector<ProbabilityImage> data);
    GaussianEstimator1D(std::vector<ProbabilityImage> data, ProbabilityImage pr);
    GaussianEstimator1D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold);
    GaussianEstimator1D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c);
    void SetData(ProbabilityImage data);
    void SetData(std::vector<ProbabilityImage> data);
    void SetParameters(ProbabilityPixelType mu_t, ProbabilityPixelType sigma_t);
    virtual void SetParameters(vnl_vector<ProbabilityPixelType> mu_t, vnl_matrix<ProbabilityPixelType> sigma_t);
    vnl_vector<ProbabilityPixelType> GetMu();
    vnl_matrix<ProbabilityPixelType> GetSigma();
    void UpdateParameters(ProbabilityImage pr,ProbabilityPixelType threshold = 0.0);
    ProbabilityImage GetCurrentLikelihood();
    ProbabilityImage GetCurrentLogLikelihood();
    ProbabilityImage GetMahalanobisMap();
    ProbabilityImage GetPositiveOutlierMap();
    ProbabilityImage GetPositiveOutlierMap(ProbabilityImage prior,bool t1=false);
    void PrintParameters();

private:

    double mu, sigma;
    ProbabilityImage data;

};

#endif // GAUSSIANESTIMATOR1D_H
