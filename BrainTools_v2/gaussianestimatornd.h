#ifndef GAUSSIANESTIMATORND_H
#define GAUSSIANESTIMATORND_H

#include <gaussianestimator.h>

class GaussianEstimatorND : public GaussianEstimator
{
public:
    GaussianEstimatorND();
    GaussianEstimatorND(MatrixType mu_t, MatrixType sigma_t);
    GaussianEstimatorND(std::vector<ProbabilityImage> data);
    GaussianEstimatorND(std::vector<ProbabilityImage> data, ProbabilityImage pr);
    GaussianEstimatorND(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold);
    GaussianEstimatorND(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c);
    void SetData(std::vector<ProbabilityImage> data);
    void SetParameters(MatrixType mu_t, MatrixType sigma_t);
    virtual void SetParameters(vnl_vector<ProbabilityPixelType> mu_t, vnl_matrix<ProbabilityPixelType> sigma_t);
    void UpdateParameters(ProbabilityImage pr,ProbabilityPixelType threshold = 0.0);
    vnl_vector<ProbabilityPixelType> GetMu();
    vnl_matrix<ProbabilityPixelType> GetSigma();
    ProbabilityImage GetCurrentLikelihood();
    ProbabilityImage GetCurrentLogLikelihood();
    ProbabilityImage GetMahalanobisMap();
    ProbabilityImage GetDistanceMap();
    ProbabilityImage GetMeanDistanceMap();
    ProbabilityImage GetOutlierMap();
    ProbabilityImage GetPositiveOutlierMap();
    ProbabilityImage GetPositiveOutlierMap(ProbabilityImage prior,bool t1=false);
    void PrintParameters();

private:

    MatrixType mu;
    MatrixType sigma;
    MultiSpectralImage image;
    unsigned int D;
};

#endif // GAUSSIANESTIMATORND_H
