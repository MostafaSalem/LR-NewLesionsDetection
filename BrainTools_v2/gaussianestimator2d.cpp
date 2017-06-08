#include "gaussianestimator2d.h"

GaussianEstimator2D::GaussianEstimator2D() : GaussianEstimator()
{
    sc2vect2 = ScalarToVector2Type::New();
}

GaussianEstimator2D::GaussianEstimator2D(Vector2Type mu_t, Matrix2x2Type sigma_t) : GaussianEstimator()
{
    mu = mu_t;
    sigma = sigma_t;
}

GaussianEstimator2D::GaussianEstimator2D(std::vector<ProbabilityImage> data) : GaussianEstimator()
{
    sc2vect2 = ScalarToVector2Type::New();
    sc2vect2->SetInput1(data[0]);
    sc2vect2->SetInput2(data[1]);
    sc2vect2->Update();
    this->data = sc2vect2->GetOutput();
}

GaussianEstimator2D::GaussianEstimator2D(std::vector<ProbabilityImage> data, ProbabilityImage pr) : GaussianEstimator()
{
    sc2vect2 = ScalarToVector2Type::New();
    sc2vect2->SetInput1(data[0]);
    sc2vect2->SetInput2(data[1]);
    sc2vect2->Update();
    this->data = sc2vect2->GetOutput();
    this->UpdateParameters(pr);
}

GaussianEstimator2D::GaussianEstimator2D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold) : GaussianEstimator()
{
    sc2vect2 = ScalarToVector2Type::New();
    sc2vect2->SetInput1(data[0]);
    sc2vect2->SetInput2(data[1]);
    sc2vect2->Update();
    this->data = sc2vect2->GetOutput();
    this->UpdateParameters(pr, threshold);
}

GaussianEstimator2D::GaussianEstimator2D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c) : GaussianEstimator()
{
    sc2vect2 = ScalarToVector2Type::New();
    sc2vect2->SetInput1(data[0]);
    sc2vect2->SetInput2(data[1]);
    sc2vect2->Update();
    this->data = sc2vect2->GetOutput();

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

void GaussianEstimator2D::SetData(std::vector<ProbabilityImage> data) {
    sc2vect2 = ScalarToVector2Type::New();
    sc2vect2->SetInput1(data[0]);
    sc2vect2->SetInput2(data[1]);
    sc2vect2->Update();
    this->data = sc2vect2->GetOutput();
}

void GaussianEstimator2D::SetParameters(Vector2Type mu_t, Matrix2x2Type sigma_t) {
    mu = mu_t;
    sigma = sigma_t;
    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = (" << mu[0] << "," << mu[1] << ")" << std::endl;
    //std::cout << "\t\tCovariance matrix" << std::endl;
    //std::cout << "\t\t|" << sigma(0,0) << "\t" << sigma(0,1) << "|" << std::endl;
    //std::cout << "\t\t|" << sigma(1,0) << "\t" << sigma(1,1) << "|" << std::endl;
}

void GaussianEstimator2D::SetParameters(vnl_vector<ProbabilityPixelType>mu_t, vnl_matrix<ProbabilityPixelType>sigma_t) {
    mu.SetVnlVector(mu_t);
    sigma = sigma_t;
    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = (" << mu[0] << "," << mu[1] << ")" << std::endl;
    //std::cout << "\t\tCovariance matrix" << std::endl;
    //std::cout << "\t\t|" << sigma(0,0) << "\t" << sigma(0,1) << "|" << std::endl;
    //std::cout << "\t\t|" << sigma(1,0) << "\t" << sigma(1,1) << "|" << std::endl;
}

vnl_vector<ProbabilityPixelType> GaussianEstimator2D::GetMu() {
    return(mu.GetVnlVector());
}

vnl_matrix<ProbabilityPixelType> GaussianEstimator2D::GetSigma() {
    return(sigma.GetVnlMatrix());
}

void GaussianEstimator2D::UpdateParameters(ProbabilityImage pr, ProbabilityPixelType threshold) {
    /* Init */
    itk::Vector<ProbabilityPixelType,2> value;
    ProbabilityPixelType sumPr=0;
    mu.Fill(0);
    sigma.Fill(0);

    /* Iterators for tissue and images */
    Multi2Iterator it( this->data, this->data->GetRequestedRegion() );
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
            value = it.Get()-mu;
            sigma(0,0) += value[0]*value[0]*prIt.Get();
            sigma(0,1) += value[0]*value[1]*prIt.Get();
            sigma(1,0) += value[1]*value[0]*prIt.Get();
            sigma(1,1) += value[1]*value[1]*prIt.Get();
        }

        ++it;
        ++prIt;
    }
    sigma = sigma / sumPr;

    //std::cout << "\t\tDistribution parameters: " << std::endl;
    //std::cout << "\t\tMean = (" << mu[0] << "," << mu[1] << ")" << std::endl;
    //std::cout << "\t\tCovariance matrix" << std::endl;
    //std::cout << "\t\t|" << sigma(0,0) << "\t" << sigma(0,1) << "|" << std::endl;
    //std::cout << "\t\t|" << sigma(1,0) << "\t" << sigma(1,1) << "|" << std::endl;


}

ProbabilityImage GaussianEstimator2D::GetCurrentLikelihood() {
    /* Init */
    Vector2Type mu_data;
    Matrix2x2Type sigmaInv(sigma.GetInverse());
    ProbabilityPixelType right, left;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi2Iterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Unidimensional Gaussian */
    left = one_over_pi/(4 * vnl_determinant(sigma.GetVnlMatrix()));
    while ( !it.IsAtEnd() ) {
        mu_data = (mu-it.Get());
        right = -0.5*(mu_data[0]*(mu_data[0]*sigmaInv(0,0)+mu_data[1]*sigmaInv(1,0))+mu_data[1]*(mu_data[0]*sigmaInv(0,1)+mu_data[1]*sigmaInv(1,1)));
        itPr.Set( left * vcl_pow((double)e,(double)right) );
        ++it;
        ++itPr;
    }

    return(pr);

}

ProbabilityImage GaussianEstimator2D::GetCurrentLogLikelihood() {
    /* Init */
    Vector2Type mu_data;
    Matrix2x2Type sigmaInv(sigma.GetInverse());
    ProbabilityPixelType right, left, pixVal;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi2Iterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    left = vcl_log(4*pi*pi*vnl_determinant(sigma.GetVnlMatrix()));
    while ( !it.IsAtEnd() ) {
        mu_data = (mu-it.Get());
        right = mu_data[0]*(mu_data[0]*sigmaInv(0,0)+mu_data[1]*sigmaInv(1,0))+mu_data[1]*(mu_data[0]*sigmaInv(0,1)+mu_data[1]*sigmaInv(1,1));
        pixVal = 0.5*(left + right);
        if (pixVal<0)
            pixVal = right;
        itPr.Set( pixVal );
        ++it;
        ++itPr;
    }

    return(pr);
}

ProbabilityImage GaussianEstimator2D::GetMahalanobisMap() {
    /* Init */
    ProbabilityPixelType maxPr=0,currentPr;
    vnl_vector<ProbabilityPixelType> mu_data;
    vnl_matrix<ProbabilityPixelType> sigma_inv(sigma.GetInverse());
    vnl_matrix<ProbabilityPixelType> t(1,2);
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi2Iterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    while ( !it.IsAtEnd() ) {
        mu_data = (it.Get().GetVnlVector()-mu.GetVnlVector());
        t.set_row(0,sigma_inv*mu_data);
        currentPr=(t*mu_data)[0];
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

void GaussianEstimator2D::PrintParameters() {
    std::cout << "\t\tDistribution parameters: " << std::endl;
    std::cout << "\t\tMean = (" << mu[0] << "," << mu[1] << ")" << std::endl;
    std::cout << "\t\tCovariance matrix" << std::endl;
    std::cout << "\t\t|" << sigma(0,0) << "\t" << sigma(0,1) << "|" << std::endl;
    std::cout << "\t\t|" << sigma(1,0) << "\t" << sigma(1,1) << "|" << std::endl;
}
