#include <gaussianestimatornd.h>

GaussianEstimatorND::GaussianEstimatorND() : GaussianEstimator()
{
}

GaussianEstimatorND::GaussianEstimatorND(MatrixType mu_t, MatrixType sigma_t) : GaussianEstimator()
{
    mu = mu_t;
    sigma = sigma_t;
}

GaussianEstimatorND::GaussianEstimatorND(std::vector<ProbabilityImage>data) : GaussianEstimator()
{
    SetData(data);

}

GaussianEstimatorND::GaussianEstimatorND(std::vector<ProbabilityImage>data, ProbabilityImage pr) : GaussianEstimator()
{
    SetData(data);
    UpdateParameters(pr);

}

GaussianEstimatorND::GaussianEstimatorND(std::vector<ProbabilityImage>data, ProbabilityImage pr, ProbabilityPixelType threshold) : GaussianEstimator()
{
    SetData(data);
    UpdateParameters(pr,threshold);

}

GaussianEstimatorND::GaussianEstimatorND(std::vector<ProbabilityImage>data, ResultsImage labels, ResultsPixelType c) : GaussianEstimator()
{
    SetData(data);

    /* From labels to "probabilities" */
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
    for (itPr.GoToBegin();!itPr.IsAtEnd();++itPr,++itLab)
        itPr.Set((ProbabilityPixelType)(itLab.Get()==c));

    UpdateParameters(pr);

}


void GaussianEstimatorND::SetParameters(MatrixType mu_t, MatrixType sigma_t) {
    mu = mu_t;
    sigma = sigma_t;
}

void GaussianEstimatorND::SetParameters(vnl_vector<ProbabilityPixelType>mu_t, vnl_matrix<ProbabilityPixelType>sigma_t) {
    vnl_matrix<ProbabilityPixelType> mu_m;
    mu_m.set_size(1,mu_t.size());
    mu_m.set_row(0,mu_t);
    mu = mu_m;
    sigma = sigma_t;
}

void GaussianEstimatorND::SetData(std::vector<ProbabilityImage> data) {
    /* Init */
    D = data.size();
    unsigned int image_n;
    // Intermediate pixel containers
    vnl_vector<ProbabilityPixelType> pixel_v;
    pixel_v.set_size(data.size());
    vnl_matrix<ProbabilityPixelType> pixel_m;
    pixel_m.set_size(1,D);
    MatrixType pixel;
    pixel.SetSize(1,D);
    // Image init
    image = MultiSpectralImageType::New();
    image->SetRegions( data[0]->GetLargestPossibleRegion() );
    image->SetSpacing( data[0]->GetSpacing() );
    image->SetOrigin( data[0]->GetOrigin() );
    image->SetDirection( data[0]->GetDirection() );
    image->Allocate();
    image->Update();
    // Iterator definition
    MultiIterator imIt = MultiIterator( image, image->GetRequestedRegion() );
    std::vector<ProbabilityIterator> imIts;
    for (image_n=0;image_n<D;image_n++)
        imIts.push_back(ProbabilityIterator(data[image_n], data[image_n]->GetRequestedRegion()));
    /* Image computation */
    for (imIt.GoToBegin();!imIt.IsAtEnd();++imIt) {
        for (image_n=0;image_n<D;image_n++) {
            pixel_v[image_n] = imIts[image_n].Get();
            ++imIts[image_n];
        }
        pixel_m.set_row(0,pixel_v);
        pixel = pixel_m;
        imIt.Set(pixel);
    }
}

void GaussianEstimatorND::UpdateParameters(ProbabilityImage pr,ProbabilityPixelType threshold) {
    /* Init */
    MatrixType mu_data, sigma_t;
    ProbabilityPixelType sumPr = 0;
    mu.SetSize(1,D);
    mu.Fill(0.0);
    sigma_t.SetSize(D,D);
    sigma.SetSize(D,D);
    sigma.Fill(0.0);

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Mean estimation */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        if (itPr.Get()>threshold) {
            mu += it.Get()*itPr.Get();
            sumPr += itPr.Get();
        }
    }
    mu /= sumPr;

    /* Standard deviation */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        if (itPr.Get()>threshold) {
            mu_data = (mu - it.Get());
            sigma_t = mu_data.GetTranspose();
            sigma+=sigma_t*mu_data*itPr.Get();
        }
    }
    sigma = sigma / sumPr;

}

ProbabilityImage GaussianEstimatorND::GetCurrentLikelihood() {
    /* Init */
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    ProbabilityPixelType right, left;
    // Image creation
    ProbabilityImage likelihood = ProbabilityImageType::New();
    likelihood->SetRegions( image->GetLargestPossibleRegion() );
    likelihood->SetSpacing( image->GetSpacing() );
    likelihood->SetOrigin( image->GetOrigin() );
    likelihood->SetDirection( image->GetDirection() );
    likelihood->Allocate();
    likelihood->FillBuffer(0.0);
    likelihood->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( likelihood, likelihood->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    left = 1/(vcl_sqrt(vcl_pow(2.0*pi,D)*vnl_determinant(sigma.GetVnlMatrix())));
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        mu_data = (mu-it.Get());
        mu_data *= sigmaInv*mu_data.GetTranspose();
        right = -0.5*mu_data(0,0);
        itPr.Set( left * vcl_pow((double)e,(double)right) );
    }
    return(likelihood);
}

ProbabilityImage GaussianEstimatorND::GetCurrentLogLikelihood() {
    /* Init */
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    ProbabilityPixelType right, left;
    // Image creation
    ProbabilityImage likelihood = ProbabilityImageType::New();
    likelihood->SetRegions( image->GetLargestPossibleRegion() );
    likelihood->SetSpacing( image->GetSpacing() );
    likelihood->SetOrigin( image->GetOrigin() );
    likelihood->SetDirection( image->GetDirection() );
    likelihood->Allocate();
    likelihood->FillBuffer(0.0);
    likelihood->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( likelihood, likelihood->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    left = vcl_log(vcl_pow(2.0*pi,D)* vnl_determinant(sigma.GetVnlMatrix()));
    if (vnl_math_isinf(left))
        left = vcl_pow(2.0,20.0);

    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        mu_data = (mu-it.Get());
        mu_data *= sigmaInv*mu_data.GetTranspose();
        right = 0.5*mu_data(0,0);
        itPr.Set( left + right );
    }
    return(likelihood);
}

ProbabilityImage GaussianEstimatorND::GetMahalanobisMap() {
    /* Init */
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    // Image creation
    ProbabilityImage mahalanobis = ProbabilityImageType::New();
    mahalanobis->SetRegions( image->GetLargestPossibleRegion() );
    mahalanobis->SetSpacing( image->GetSpacing() );
    mahalanobis->SetOrigin( image->GetOrigin() );
    mahalanobis->SetDirection( image->GetDirection() );
    mahalanobis->Allocate();
    mahalanobis->FillBuffer(0.0);
    mahalanobis->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( mahalanobis, mahalanobis->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        mu_data = (mu-it.Get());
        mu_data *= sigmaInv*mu_data.GetTranspose();
        itPr.Set( mu_data(0,0) );
    }
    return(mahalanobis);
}

ProbabilityImage GaussianEstimatorND::GetDistanceMap() {
    /* Init */
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    // Image creation
    ProbabilityImage distance = ProbabilityImageType::New();
    distance->SetRegions( image->GetLargestPossibleRegion() );
    distance->SetSpacing( image->GetSpacing() );
    distance->SetOrigin( image->GetOrigin() );
    distance->SetDirection( image->GetDirection() );
    distance->Allocate();
    distance->FillBuffer(0.0);
    distance->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( distance, distance->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr) {
        mu_data = (mu-it.Get())*sigmaInv;
        itPr.Set( mu_data.GetVnlMatrix().get_row(0).inf_norm() );
    }
    return(distance);
}

ProbabilityImage GaussianEstimatorND::GetMeanDistanceMap() {
    /* Init */
    ProbabilityPixelType sumDist;
    ProbabilityPixelType sumPr = 0;
    MatrixType dVector;
    // Image creation
    ProbabilityImage distance = ProbabilityImageType::New();
    distance->SetRegions( image->GetLargestPossibleRegion() );
    distance->SetSpacing( image->GetSpacing() );
    distance->SetOrigin( image->GetOrigin() );
    distance->SetDirection( image->GetDirection() );
    distance->Allocate();
    distance->FillBuffer(0.0);
    distance->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    MultiIterator itData( image, image->GetRequestedRegion() );
    ProbabilityIterator itDist( distance, distance->GetRequestedRegion() );
    /* Current likelihood probabilities */
    ProbabilityImage pr = GetCurrentLikelihood();
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    // Probability sum
    for (itPr.GoToBegin();!itPr.IsAtEnd();++itPr) {
        sumPr += itPr.Get();
    }

    /* Mean distance to all points */
    for (it.GoToBegin(),itDist.GoToBegin();!it.IsAtEnd(),!itDist.IsAtEnd();++it,++itDist) {
        sumDist = 0;
        for (itData.GoToBegin(),itPr.GoToBegin();!itData.IsAtEnd(),!itPr.IsAtEnd();++itData,++itPr) {
            dVector = (it.Get()-itData.Get());
            sumDist += itPr.Get()*dVector.GetVnlMatrix().get_row(0).two_norm();
        }
        itDist.Set(sumDist/sumPr);
    }
    return(distance);
}

ProbabilityImage GaussianEstimatorND::GetOutlierMap() {
    /* Init */
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    double right;
    double K = vcl_log(vcl_pow(2*pi,D)* vnl_determinant(sigma.GetVnlMatrix()));
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( image->GetLargestPossibleRegion() );
    pr->SetSpacing( image->GetSpacing() );
    pr->SetOrigin( image->GetOrigin() );
    pr->SetDirection( image->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr ) {
        mu_data = (mu-it.Get());
        mu_data *= sigmaInv*mu_data.GetTranspose();
        right = -0.5*mu_data(0,0);
        itPr.Set(K - K * vcl_pow(e,right) );
    }

    return(pr);
}

ProbabilityImage GaussianEstimatorND::GetPositiveOutlierMap() {
    /* Init */
    MatrixType pix;
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    int i;
    bool positive = true;
    double right;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( image->GetLargestPossibleRegion() );
    pr->SetSpacing( image->GetSpacing() );
    pr->SetOrigin( image->GetOrigin() );
    pr->SetDirection( image->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr ) {
        i = 0;
        do {
                positive = mu(0,i)<=it.Get()(0,i);
            i++;
        } while (positive&&i<D);
        if (positive) {
            mu_data = (mu-it.Get());
            mu_data *= sigmaInv*mu_data.GetTranspose();
            right = vcl_sqrt(mu_data(0,0));
            itPr.Set( right );
        }
    }

    return(pr);
}

ProbabilityImage GaussianEstimatorND::GetPositiveOutlierMap(ProbabilityImage prior, bool t1) {
    /* Init */
    MatrixType pix;
    MatrixType mu_data;
    MatrixType sigmaInv;
    sigmaInv = sigma.GetInverse();
    int i;
    bool positive = true;
    double right;
    ProbabilityImage pr = ProbabilityImageType::New();
    pr->SetRegions( image->GetLargestPossibleRegion() );
    pr->SetSpacing( image->GetSpacing() );
    pr->SetOrigin( image->GetOrigin() );
    pr->SetDirection( image->GetDirection() );
    pr->Allocate();
    pr->FillBuffer(0.0);
    pr->Update();

    /* Iterators for tissue and images */
    MultiIterator it( image, image->GetRequestedRegion() );
    ProbabilityIterator itPr( pr, pr->GetRequestedRegion() );
    ProbabilityIterator itPi( prior, prior->GetRequestedRegion() );

    /* Multidimensional Gaussian */
    for (it.GoToBegin(),itPr.GoToBegin();!it.IsAtEnd(),!itPr.IsAtEnd();++it,++itPr,++itPi ) {
        i = 0;
        do {
            if ((!t1) || (i<D-1))
                positive = mu(0,i)<=it.Get()(0,i);
            else
                positive = mu(0,i)+2*vcl_sqrt(sigma(i,i))>=it.Get()(0,i);
            i++;
        } while ((positive)&&(i<D));
       if (positive) {
            mu_data = (mu-it.Get());
            mu_data *= sigmaInv*mu_data.GetTranspose();
            right = vcl_sqrt(mu_data(0,0));
            itPr.Set( right*itPi.Get() );
        }
    }

    return(pr);
}

void GaussianEstimatorND::PrintParameters() {
    unsigned int d;
    std::cout << "\t\tMean = (" << mu.GetVnlMatrix().get(0,0);
    for (d = 1; d < D; d++)
        std::cout << "," << mu.GetVnlMatrix().get(0,d);
    std::cout << ")" << std::endl;

    std::cout << "\t\tCovariance matrix" << std::endl;
    for (d=0; d < D; d++)
        std::cout << "\t\t|" << sigma.GetVnlMatrix().get_row(d) << "|" << std::endl;
}

vnl_vector<ProbabilityPixelType> GaussianEstimatorND::GetMu() {
    return(mu.GetVnlMatrix().get_row(0));
}

vnl_matrix<ProbabilityPixelType> GaussianEstimatorND::GetSigma() {
    return(sigma.GetVnlMatrix());
}
