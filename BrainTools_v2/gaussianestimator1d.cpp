#include "gaussianestimator1d.h"

GaussianEstimator1D::GaussianEstimator1D() : GaussianEstimator()
{
    
}

GaussianEstimator1D::GaussianEstimator1D(ProbabilityPixelType mu_t, ProbabilityPixelType sigma_t) : GaussianEstimator()
{
    mu = mu_t;
    sigma = sigma_t;
}

GaussianEstimator1D::GaussianEstimator1D(std::vector<ProbabilityImage> data) : GaussianEstimator()
{
    this->data = data[0];
}

GaussianEstimator1D::GaussianEstimator1D(std::vector<ProbabilityImage> data, ProbabilityImage pr) : GaussianEstimator()
{
    this->data = data[0];
    this->UpdateParameters(pr);
}

GaussianEstimator1D::GaussianEstimator1D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold) : GaussianEstimator()
{
    this->data = data[0];
    this->UpdateParameters(pr,threshold);
}

GaussianEstimator1D::GaussianEstimator1D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c) : GaussianEstimator()
{
    this->data = data[0];
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( labels->GetLargestPossibleRegion() );
    pr->SetSpacing( labels->GetSpacing() );
    pr->SetOrigin( labels->GetOrigin() );
    pr->SetDirection( labels->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    ResultsIterator itLab( labels, labels->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    while (!itPr.IsAtEnd()) {

        itPr.Set((ProbabilityPixelType)(itLab.Get()==c));

        ++itPr;
        ++itLab;
    }

    this->UpdateParameters(pr);
}

void GaussianEstimator1D::SetData(ProbabilityImage data) {
    this->data = data;
}

void GaussianEstimator1D::SetData(std::vector<ProbabilityImage> data) {
    this->data = data[0];
}

void GaussianEstimator1D::SetParameters(ProbabilityPixelType mu_t, ProbabilityPixelType sigma_t) {
    mu = mu_t;
    sigma = sigma_t;
    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = " << mu << std::endl;
    //std::cout << "\t\tStandard deviation = " << sigma << std::endl;
    //std::cout << "\t\t------------------------------" << std::endl;
}

void GaussianEstimator1D::SetParameters(vnl_vector<ProbabilityPixelType>mu_t, vnl_matrix<ProbabilityPixelType>sigma_t) {
    mu = mu_t.get(0);
    sigma = sigma_t.get(0,0);
    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = " << mu << std::endl;
    //std::cout << "\t\tStandard deviation = " << sigma << std::endl;
    //std::cout << "\t\t------------------------------" << std::endl;
}

vnl_vector<ProbabilityPixelType> GaussianEstimator1D::GetMu() {
    vnl_vector<ProbabilityPixelType> mu_t(1);
    mu_t.put(0,mu);
    return(mu_t);
}

vnl_matrix<ProbabilityPixelType> GaussianEstimator1D::GetSigma() {
    vnl_matrix<ProbabilityPixelType> sigma_t(1,1);
    sigma_t.put(0,0,sigma);
    return(sigma_t);
}

void GaussianEstimator1D::UpdateParameters(ProbabilityImage pr, ProbabilityPixelType threshold) {
    /* Init */
    ProbabilityPixelType pixVal;
    ProbabilityPixelType sumPr = 0;
    mu = 0;
    sigma = 0;

    /* Iterators for tissue and images */
    ProbabilityIterator it( this->data, this->data->GetRequestedRegion() );
    ProbabilityIterator prIt( pr, pr->GetRequestedRegion() );

    /* Mean estimation */
    while ( !it.IsAtEnd() ) {
        if (prIt.Get()>=threshold) {
            mu += it.Get()*prIt.Get();
            sumPr += prIt.Get();
        }

        ++it;
        ++prIt;
    }
    mu /= sumPr;

    /* Standard deviation */
    it.GoToBegin();
    prIt.GoToBegin();
    while ( !it.IsAtEnd() ) {
        if (prIt.Get()>=threshold) {
            pixVal = it.Get();
            sigma += (pixVal-mu)*(pixVal-mu)*prIt.Get();
        }

        ++it;
        ++prIt;
    }
    sigma = vcl_sqrt(sigma / sumPr);

    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = " << mu << std::endl;
    //std::cout << "\t\tStandard deviation = " << sigma << std::endl;
    //std::cout << "\t\t------------------------------" << std::endl;


}

ProbabilityImage GaussianEstimator1D::GetCurrentLikelihood() {
    /* Init */
    ProbabilityPixelType pixVal;
    ProbabilityPixelType right;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    ProbabilityIterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Unidimensional Gaussian */
    while ( !it.IsAtEnd() ) {
        pixVal = it.Get();
        right = -0.5*(pixVal-mu)*(pixVal-mu)/(2*sigma*sigma);
        itPr.Set( (one_over_sqrt2pi / sigma) * vcl_pow((double)e,(double)right) );
        ++it;
        ++itPr;
    }

    return(pr);

}

ProbabilityImage GaussianEstimator1D::GetCurrentLogLikelihood() {
    /* Init */
    ProbabilityPixelType pixVal;
    ProbabilityPixelType left;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    ProbabilityIterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Unidimensional Gaussian */
    left = 0.5*vcl_log(2.0*pi*sigma*sigma);
    while ( !it.IsAtEnd() ) {
        pixVal = it.Get();
        pixVal = left + (pixVal-mu)*(pixVal-mu) / (2*sigma*sigma);
        if (pixVal<0)
            pixVal = (pixVal-mu)*(pixVal-mu) / (2*sigma*sigma);
        itPr.Set( pixVal );
        ++it;
        ++itPr;
    }

    return(pr);
}

ProbabilityImage GaussianEstimator1D::GetMahalanobisMap() {
    /* Init */
    ProbabilityPixelType maxPr=0,currentPr;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    ProbabilityIterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    while ( !it.IsAtEnd() ) {
        currentPr=(it.Get()-mu)*(it.Get()-mu)/(sigma*sigma);
        if (currentPr>maxPr)
            maxPr=currentPr;
        itPr.Set( currentPr );
        ++it;
        ++itPr;
    }
    itPr.GoToBegin();
    while ( !itPr.IsAtEnd() ) {
        itPr.Set( itPr.Get()/maxPr );
        ++itPr;
    }

    return(pr);
}

ProbabilityImage GaussianEstimator1D::GetPositiveOutlierMap() {
    /* Init */
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    ProbabilityIterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr ) {
        if (mu<it.Get()) {
             itPr.Set(vcl_sqrt((mu-it.Get())*(mu-it.Get())/sigma));
         }
    }

    return(pr);
}
ProbabilityImage GaussianEstimator1D::GetPositiveOutlierMap(ProbabilityImage prior, bool t1) {
    /* Init */
    int i;
    double right;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    ProbabilityIterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );
    ProbabilityIterator itPi( prior, prior->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr,++itPi ) {
        i = 0;
       if (((t1)&(mu>=it.Get()))||((!t1)&(mu+2*sigma<it.Get()))) {
            right = vcl_sqrt((mu-it.Get())*(mu-it.Get())/sigma);
            itPr.Set( right*itPi.Get() );
        }
    }

    return(pr);
}

void GaussianEstimator1D::PrintParameters() {
    std::cout << "\t\tDistribution parameters: " << std::endl;
    std::cout << "\t\tMean = " << mu << std::endl;
    std::cout << "\t\tStandard deviation = " << sigma << std::endl;
    std::cout << "\t\t------------------------------" << std::endl;
}
