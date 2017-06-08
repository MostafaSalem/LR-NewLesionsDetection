#include "gaussianestimator3d.h"

GaussianEstimator3D::GaussianEstimator3D() : GaussianEstimator()
{
    sc3vect3 = ScalarToVector3Type::New();
}

GaussianEstimator3D::GaussianEstimator3D(Vector3Type mu_t, Matrix3x3Type sigma_t) : GaussianEstimator()
{
    mu = mu_t;
    sigma = sigma_t;
}

GaussianEstimator3D::GaussianEstimator3D(std::vector<ProbabilityImage> data) : GaussianEstimator()
{
    sc3vect3 = ScalarToVector3Type::New();
    sc3vect3->SetInput1(data[0]);
    sc3vect3->SetInput2(data[1]);
    sc3vect3->SetInput3(data[2]);
    sc3vect3->Update();
    this->data = sc3vect3->GetOutput();
}

GaussianEstimator3D::GaussianEstimator3D(std::vector<ProbabilityImage> data, ProbabilityImage pr) : GaussianEstimator()
{
    sc3vect3 = ScalarToVector3Type::New();
    sc3vect3->SetInput1(data[0]);
    sc3vect3->SetInput2(data[1]);
    sc3vect3->SetInput3(data[2]);
    sc3vect3->Update();
    this->data = sc3vect3->GetOutput();
    this->UpdateParameters(pr);
}

GaussianEstimator3D::GaussianEstimator3D(std::vector<ProbabilityImage> data, ProbabilityImage pr, ProbabilityPixelType threshold) : GaussianEstimator()
{
    sc3vect3 = ScalarToVector3Type::New();
    sc3vect3->SetInput1(data[0]);
    sc3vect3->SetInput2(data[1]);
    sc3vect3->SetInput3(data[2]);
    sc3vect3->Update();
    this->data = sc3vect3->GetOutput();
    this->UpdateParameters(pr, threshold);
}

GaussianEstimator3D::GaussianEstimator3D(std::vector<ProbabilityImage> data, ResultsImage labels, ResultsPixelType c) : GaussianEstimator()
{
    sc3vect3 = ScalarToVector3Type::New();
    sc3vect3->SetInput1(data[0]);
    sc3vect3->SetInput2(data[1]);
    sc3vect3->SetInput3(data[2]);
    sc3vect3->Update();
    this->data = sc3vect3->GetOutput();
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

void GaussianEstimator3D::SetData(std::vector<ProbabilityImage> data) {
    sc3vect3 = ScalarToVector3Type::New();
    sc3vect3->SetInput1(data[0]);
    sc3vect3->SetInput2(data[1]);
    sc3vect3->SetInput3(data[2]);
    sc3vect3->Update();
    this->data = sc3vect3->GetOutput();
}

void GaussianEstimator3D::SetParameters(Vector3Type mu_t, Matrix3x3Type sigma_t) {
    mu = mu_t;
    sigma = sigma_t;
    //std::cout << "\tDistribution parameters: " << std::endl;
    //std::cout << "\tMean = (" << mu[0] << "," << mu[1] << "," << mu[2] << ")" << std::endl;
    //std::cout << "\tCovariance matrix" << std::endl;
    //std::cout << "\t|" << sigma(0,0) << "\t" << sigma(0,1) << "\t" << sigma(0,2) << "|" << std::endl;
    //std::cout << "\t|" << sigma(1,0) << "\t" << sigma(1,1) << "\t" << sigma(1,2) << "|" <<  std::endl;
    //std::cout << "\t|" << sigma(2,0) << "\t" << sigma(2,1) << "\t" << sigma(2,2) << "|" <<  std::endl;
}

void GaussianEstimator3D::SetParameters(vnl_vector<ProbabilityPixelType>mu_t, vnl_matrix<ProbabilityPixelType>sigma_t) {
    mu.SetVnlVector(mu_t);
    sigma = sigma_t;
    //std::cout << "\tDistribution parameters: " << std::endl;
    //std::cout << "\tMean = (" << mu[0] << "," << mu[1] << "," << mu[2] << ")" << std::endl;
    //std::cout << "\tCovariance matrix" << std::endl;
    //std::cout << "\t|" << sigma(0,0) << "\t" << sigma(0,1) << "\t" << sigma(0,2) << "|" << std::endl;
    //std::cout << "\t|" << sigma(1,0) << "\t" << sigma(1,1) << "\t" << sigma(1,2) << "|" <<  std::endl;
    //std::cout << "\t|" << sigma(2,0) << "\t" << sigma(2,1) << "\t" << sigma(2,2) << "|" <<  std::endl;
}

vnl_vector<ProbabilityPixelType> GaussianEstimator3D::GetMu() {
    return(mu.GetVnlVector());
}

vnl_matrix<ProbabilityPixelType> GaussianEstimator3D::GetSigma() {
    return(sigma.GetVnlMatrix());
}

void GaussianEstimator3D::UpdateParameters(ProbabilityImage pr, ProbabilityPixelType threshold) {
    /* Init */
    itk::Vector<ProbabilityPixelType,3> value;
    ProbabilityPixelType sumPr=0;
    mu.Fill(0);
    sigma.Fill(0);;

    /* Iterators for tissue and images */
    Multi3Iterator it( this->data, this->data->GetRequestedRegion() );
    ProbabilityIterator prIt( pr, pr->GetRequestedRegion() );

    /* Mean estimation */
    while ( !it.IsAtEnd() ) {
        if (prIt.Get()>threshold) {
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
        if (prIt.Get()>threshold) {
            value = it.Get()-mu;
            sigma(0,0) += value[0]*value[0]*prIt.Get();
            sigma(0,1) += value[0]*value[1]*prIt.Get();
            sigma(0,2) += value[0]*value[2]*prIt.Get();
            sigma(1,0) += value[1]*value[0]*prIt.Get();
            sigma(1,1) += value[1]*value[1]*prIt.Get();
            sigma(1,2) += value[1]*value[2]*prIt.Get();
            sigma(2,0) += value[2]*value[0]*prIt.Get();
            sigma(2,1) += value[2]*value[1]*prIt.Get();
            sigma(2,2) += value[2]*value[2]*prIt.Get();
        }

        ++it;
        ++prIt;
    }
    sigma = sigma / sumPr;

    //std::cout << "\tDistribution parameters: " << std::endl;
    //std::cout << "\tMean = (" << mu[0] << "," << mu[1] << "," << mu[2] << ")" << std::endl;
    //std::cout << "\tCovariance matrix" << std::endl;
    //std::cout << "\t|" << sigma(0,0) << "\t" << sigma(0,1) << "\t" << sigma(0,2) << "|" << std::endl;
    //std::cout << "\t|" << sigma(1,0) << "\t" << sigma(1,1) << "\t" << sigma(1,2) << "|" <<  std::endl;
    //std::cout << "\t|" << sigma(2,0) << "\t" << sigma(2,1) << "\t" << sigma(2,2) << "|" <<  std::endl;


}

ProbabilityImage GaussianEstimator3D::GetCurrentLikelihood() {
    /* Init */
    Vector3Type mu_data;
    Matrix3x3Type sigmaInv(sigma.GetInverse());
    ProbabilityPixelType right, right1, right2, right3, left;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi3Iterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Unidimensional Gaussian */
    left = 1/(8*pi*pi*pi* vcl_sqrt(vnl_determinant(sigma.GetVnlMatrix())));
    while ( !it.IsAtEnd() ) {
        mu_data = (mu-it.Get());
        right1 = mu_data[0]*(mu_data[0]*sigmaInv(0,0)+mu_data[1]*sigmaInv(1,0)+mu_data[2]*sigmaInv(2,0));
        right2 = mu_data[1]*(mu_data[0]*sigmaInv(0,1)+mu_data[1]*sigmaInv(1,1)+mu_data[2]*sigmaInv(2,1));
        right3 = mu_data[2]*(mu_data[0]*sigmaInv(0,2)+mu_data[1]*sigmaInv(1,2)+mu_data[2]*sigmaInv(2,2));
        right = -0.5*(right1+right2+right3);
        itPr.Set( left * vcl_pow((double)e,(double)right) );
        ++it;
        ++itPr;
    }

    return(pr);

}

ProbabilityImage GaussianEstimator3D::GetCurrentLogLikelihood() {
    /* Init */
    Vector3Type mu_data;
    Matrix3x3Type sigmaInv(sigma.GetInverse());
    ProbabilityPixelType right, right1, right2, right3, left, pixVal;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi3Iterator it( data, data->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    left = vcl_log(8*pi*pi*pi*vnl_determinant(sigma.GetVnlMatrix()));
    while ( !it.IsAtEnd() ) {
        mu_data = (mu-it.Get());
        right1 = mu_data[0]*(mu_data[0]*sigmaInv(0,0)+mu_data[1]*sigmaInv(1,0)+mu_data[2]*sigmaInv(2,0));
        right2 = mu_data[1]*(mu_data[0]*sigmaInv(0,1)+mu_data[1]*sigmaInv(1,1)+mu_data[2]*sigmaInv(2,1));
        right3 = mu_data[2]*(mu_data[0]*sigmaInv(0,2)+mu_data[1]*sigmaInv(1,2)+mu_data[2]*sigmaInv(2,2));
        right = -0.5*(right1+right2+right3);
        pixVal = 0.5*(left + right);
        if (pixVal<0)
            pixVal = right;
         itPr.Set( pixVal );
        ++it;
        ++itPr;
    }

    return(pr);
}

ProbabilityImage GaussianEstimator3D::GetMahalanobisMap() {
    /* Init */
    ProbabilityPixelType maxPr=0,currentPr;
    vnl_vector<ProbabilityPixelType> mu_data;
    vnl_matrix<ProbabilityPixelType> sigma_inv(sigma.GetInverse());
    vnl_matrix<ProbabilityPixelType> t(1,3);
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( data->GetLargestPossibleRegion() );
    pr->SetSpacing( data->GetSpacing() );
    pr->SetOrigin( data->GetOrigin() );
    pr->SetDirection( data->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    Multi3Iterator it( data, data->GetRequestedRegion() );
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
void GaussianEstimator3D::PrintParameters() {
    std::cout << "\tDistribution parameters: " << std::endl;
    std::cout << "\tMean = (" << mu[0] << "," << mu[1] << "," << mu[2] << ")" << std::endl;
    std::cout << "\tCovariance matrix" << std::endl;
    std::cout << "\t|" << sigma(0,0) << "\t" << sigma(0,1) << "\t" << sigma(0,2) << "|" << std::endl;
    std::cout << "\t|" << sigma(1,0) << "\t" << sigma(1,1) << "\t" << sigma(1,2) << "|" <<  std::endl;
    std::cout << "\t|" << sigma(2,0) << "\t" << sigma(2,1) << "\t" << sigma(2,2) << "|" <<  std::endl;
}
