#include <brainsegmentation.h>
#include <sstream>

/*************************************************************************************************************************************************/
//
// Function to segment lesions based on de Boer et al. and Souplet et al. approaches. This lesion uses a previous tissue segmentation and
// FLAIR images.
//
// Inputs:
//  tissue = Segmentation image.
//   flair = FLAIR MRI image to be segmented.
//   alpha = Number of deviations defining the lesions. The thresholded is computed using the following formula: T = mean_gm + alpha*std_gm,
//           where mean_gm and std_gm are the mean and standard deviation of the GM mask in the FLAIR image.
//    beta = This parameter is used to define the minimum percentage of WM voxels surrounding a lesion. Lesions with values below this
//           thresholded are discarded.
//
/*************************************************************************************************************************************************/
MaskImage BrainSegmentation::Lesion_New(ResultsImage tissue, ProbabilityImage flair,double alpha, double omegaT, double omegaN, int minSize) {
	std::cout << "\tBrainSegmentation::Lesion_New" << std::endl;
    /* Init */
	ResultsImageType::SizeType imgSize = tissue->GetLargestPossibleRegion().GetSize();
    ResultsImageType::IndexType index;
    std::cout << "\t\\- Init" << std::endl;
    ProbabilityImage flair_e = ContrastEnhancement(flair,1);
    BrainIO* brainio = new BrainIO();
    brainio->WriteProbabilityImage("contrasted_flair.nii",flair_e);
    unsigned int i;

    /* Threshold estimation */
    std::cout << "\t\\- Threshold estimation" << std::endl;
    ResultsIterator itLabel = ResultsIterator(tissue,tissue->GetRequestedRegion());
    ProbabilityIterator itFLAIR = ProbabilityIterator(flair_e,flair_e->GetRequestedRegion());
    for (itLabel.GoToBegin();!itLabel.IsAtEnd();++itLabel,++itFLAIR) {
        if (itLabel.Get()!=2)
            itFLAIR.Set(-1);
    }

    MinimumMaximumImageCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumImageCalculatorType::New();
    //maxIntensityCalculator->SetImage(flair_e);
    maxIntensityCalculator->SetImage(flair);
    maxIntensityCalculator->Compute();
    ProbabilityPixelType maxIntensity = maxIntensityCalculator->GetMaximum();

    PHistogramGeneratorType::Pointer histogramGenerator = PHistogramGeneratorType::New();
    //histogramGenerator->SetInput(flair_e);
    histogramGenerator->SetInput(flair);
    histogramGenerator->SetNumberOfBins( 256 );
    histogramGenerator->SetHistogramMin( 0 );
    histogramGenerator->SetHistogramMax( maxIntensity );
    histogramGenerator->Compute();

    const PHistogramType * histogram = histogramGenerator->GetOutput();
    const unsigned int histogramSize = histogram->Size();
    unsigned int currentValue, maxBin=1, maxBinValue=0;
    ProbabilityPixelType mu = 0, sigma;
    // Mean as peak
    for(currentValue=1; currentValue < histogramSize; currentValue++ )
    {

        if (histogram->GetFrequency( currentValue, 0 )>maxBinValue) {
            mu = histogram->GetMeasurement( currentValue, 0 );
            maxBinValue = histogram->GetFrequency( currentValue, 0 );
            maxBin = currentValue;
        }
    }


    // Sigma via full half width maximum estimation
    // Lower band
    unsigned int x1=maxBin-1,x2=maxBin+1;
    while (double(histogram->GetFrequency( x1, 0))/maxBinValue > 0.5 )
        x1--;
    double hx1,hx1m1, mx1, mx1m1, loband;
    hx1=double(histogram->GetFrequency( x1, 0))/maxBinValue;
    hx1m1=double(histogram->GetFrequency( x1+1, 0))/maxBinValue;
    mx1=double(histogram->GetMeasurement( x1, 0));
    mx1m1=double(histogram->GetMeasurement( x1+1, 0));
    loband = mx1+(0.5-hx1)*(mx1m1-mx1)/(hx1m1-hx1);
    // Higher band
    while (double(histogram->GetFrequency( x2, 0))/maxBinValue > 0.5 )
        x2++;
    double hx2,hx2m1,mx2,mx2m1, hiband;
    hx2=double(histogram->GetFrequency( x2, 0))/maxBinValue;
    hx2m1=double(histogram->GetFrequency( x2-1, 0))/maxBinValue;
    mx2=double(histogram->GetMeasurement( x2, 0));
    mx2m1=double(histogram->GetMeasurement( x2-1, 0));
    hiband = mx2m1 -(0.5-hx2)*(mx2m1-mx2)/(hx2m1-hx2);
    sigma=(hiband-loband)/(2*vcl_sqrt(2*vcl_log(2.0)));

    ProbabilityPixelType t = mu + alpha*sigma;
    std::cout << "\tThreshold: " << t << " (" << mu << " + " << alpha << " * " << sigma << ")" << std::endl;
    std::cout << "\t\\- Lower band: " << loband << " (" << mx1 << ")" << std::endl;
    std::cout << "\t\\- Higher band: " << hiband << " (" << mx2 << ")" << std::endl;

    /* Thresholding and refinement */
    std::cout << "\t\\- Thresholding and refinement" << std::endl;

    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    thresholder->SetLowerThreshold(t);
    //thresholder->SetInput(flair_e);
    thresholder->SetInput(flair);
    thresholder->Update();

    // Thresholding
    MaskImage wml = thresholder->GetOutput();


    // Lesion neighborhood refinement (WM label = 3)
    MaskIterator itWML = MaskIterator(wml,wml->GetRequestedRegion());

    ConnectComponentFilterType::Pointer connectComp = ConnectComponentFilterType::New();
    connectComp->SetInput(wml);
    connectComp->Update();

    ConnectedImage components = connectComp->GetOutput();
    MinimumMaximumLabelComponentType::Pointer maxComponentCalculator = MinimumMaximumLabelComponentType::New();
    maxComponentCalculator->SetImage(components);
    maxComponentCalculator->Compute();
    ConnectedPixelType maxComponent = maxComponentCalculator->GetMaximum();

    if (maxComponent > 0) {

        double *yesNeighbors = new double[maxComponent];
        double *noNeighbors = new double[maxComponent];
        double *yesTissue = new double[maxComponent];
        double *noTissue = new double[maxComponent];
        int *nVoxels = new int[maxComponent];
        bool *checked = new bool[maxComponent];
        bool *centered = new bool[maxComponent];
        std::fill_n(yesNeighbors,maxComponent,0);
        std::fill_n(noNeighbors,maxComponent,0);
        std::fill_n(yesTissue,maxComponent,0);
        std::fill_n(noTissue,maxComponent,0);
        std::fill_n(nVoxels,maxComponent,0);
        std::fill_n(checked,maxComponent,false);
        std::fill_n(centered,maxComponent,false);

        ConnectedIterator itComp = ConnectedIterator(components,components->GetRequestedRegion());
        ConnectedNeighborhoodIterator::RadiusType radius;
        for (i = 0; i < IntensityImageType::ImageDimension; ++i)
            radius[i] = 1;
        ConnectedNeighborhoodIterator itNeighbors = ConnectedNeighborhoodIterator(radius,components,components->GetRequestedRegion());
        ConnectedPixelType lesionId;

        // The idea is to check at each voxel if it is a lesion neighbour. If that is the case, we check if the voxel is a permitted tissue voxel (WM) or not.
        // This will allow us to compute the WM fraction value.
        for (itLabel.GoToBegin();!itLabel.IsAtEnd();++itLabel,++itComp,++itNeighbors) {
            if (itComp.Get()==0) {
                for (i = 0; i < itNeighbors.Size(); ++i) {
                    lesionId = itNeighbors.GetPixel(i);
                    if ( lesionId > 0 && !checked[lesionId-1]) {
                        checked[lesionId-1] = true;
                        switch (itLabel.Get()) {
                        // WM
                        case 3:
                            yesNeighbors[lesionId-1]++;
                            break;
                        // Other tissues
                        default:
                            noNeighbors[lesionId-1]++;
                            break;
                        }
                    }
                }
            }
            else {
                lesionId = itComp.Get();
                switch (itLabel.Get()) {
                // WM + GM-based
                case 2:
                    yesTissue[lesionId-1]++;
                    break;
                case 3:
                    yesTissue[lesionId-1]++;
                    break;
                case 4:
                    yesTissue[lesionId-1]++;
                    break;
                // Other tissues
                default:
                    noTissue[lesionId-1]++;
                    break;
                }
                nVoxels[lesionId-1]++;
                index = itComp.GetIndex();
                centered[lesionId-1] = centered[lesionId-1] |\
                                       (((index[0] > imgSize[0]/2-3) & (index[0] < imgSize[0]/2+4)) &\
                                       ((index[1] > imgSize[1]/2-3) & (index[1] < imgSize[1]/2+4)));
            }
            std::fill_n(checked,maxComponent,false);
        }

        // Finally we relabel those lesions with a WM fraction value lower than beta
        bool neighborsCondition, tissueCondition, sizeCondition, centeredCondition;
        for (itComp.GoToBegin();!itComp.IsAtEnd();++itComp,++itWML) {
            lesionId = itComp.Get();
            if (lesionId>0) {
                // Neighbouring condition
                neighborsCondition = yesNeighbors[lesionId-1] > (yesNeighbors[lesionId-1] + noNeighbors[lesionId-1])*omegaN;
                // Tissue condition
                tissueCondition = yesTissue[lesionId-1] > (yesTissue[lesionId-1] + noTissue[lesionId-1])*omegaT;
                // Size condition
                sizeCondition = nVoxels[lesionId-1]>=minSize;
                // Centered condition
                centeredCondition = !centered[lesionId-1];

                itWML.Set(  neighborsCondition & tissueCondition & sizeCondition & centeredCondition );
                //itWML.Set( true );

            }
        }

    }

	delete(brainio);

    return(wml);
}

MaskImage BrainSegmentation::Lesion_New(ResultsImage tissue, ProbabilityImage flair,std::vector<int*> pv,double alpha, double omegaT, double omegaN, int minSize) {
	std::cout << "\tBrainSegmentation::Lesion_New" << std::endl;
    /* Init */
    ResultsImageType::SizeType imgSize = tissue->GetLargestPossibleRegion().GetSize();
    ResultsImageType::IndexType index;
    std::cout << "\t\\- Init" << std::endl;
    ProbabilityImage flair_e = ContrastEnhancement(flair,1);
    BrainIO* brainio = new BrainIO();
    brainio->WriteProbabilityImage("contrasted_flair.nii",flair_e);
    unsigned int i;

    /* Threshold estimation */
    std::cout << "\t\\- Threshold estimation" << std::endl;
    ResultsIterator itLabel = ResultsIterator(tissue,tissue->GetRequestedRegion());
    ProbabilityIterator itFLAIR = ProbabilityIterator(flair_e,flair_e->GetRequestedRegion());
    for (itLabel.GoToBegin();!itLabel.IsAtEnd();++itLabel,++itFLAIR) {
        if (itLabel.Get()!=2)
            itFLAIR.Set(-1);
    }

    MinimumMaximumImageCalculatorType::Pointer maxIntensityCalculator = MinimumMaximumImageCalculatorType::New();
    //maxIntensityCalculator->SetImage(flair_e);
    maxIntensityCalculator->SetImage(flair);
    maxIntensityCalculator->Compute();
    ProbabilityPixelType maxIntensity = maxIntensityCalculator->GetMaximum();

    PHistogramGeneratorType::Pointer histogramGenerator = PHistogramGeneratorType::New();
    //histogramGenerator->SetInput(flair_e);
    histogramGenerator->SetInput(flair);
    histogramGenerator->SetNumberOfBins( 256 );
    histogramGenerator->SetHistogramMin( 0 );
    histogramGenerator->SetHistogramMax( maxIntensity );
    histogramGenerator->Compute();

    const PHistogramType * histogram = histogramGenerator->GetOutput();
    const unsigned int histogramSize = histogram->Size();
    unsigned int currentValue, maxBin=1, maxBinValue=0;
    ProbabilityPixelType mu = 0, sigma;
    // Mean as peak
    for(currentValue=1; currentValue < histogramSize; currentValue++ )
    {

        if (histogram->GetFrequency( currentValue, 0 )>maxBinValue) {
            mu = histogram->GetMeasurement( currentValue, 0 );
            maxBinValue = histogram->GetFrequency( currentValue, 0 );
            maxBin = currentValue;
        }
    }


    // Sigma via full half width maximum estimation
    // Lower band
    unsigned int x1=maxBin-1,x2=maxBin+1;
    while (double(histogram->GetFrequency( x1, 0))/maxBinValue > 0.5 )
        x1--;
    double hx1,hx1m1, mx1, mx1m1, loband;
    hx1=double(histogram->GetFrequency( x1, 0))/maxBinValue;
    hx1m1=double(histogram->GetFrequency( x1+1, 0))/maxBinValue;
    mx1=double(histogram->GetMeasurement( x1, 0));
    mx1m1=double(histogram->GetMeasurement( x1+1, 0));
    loband = mx1+(0.5-hx1)*(mx1m1-mx1)/(hx1m1-hx1);
    // Higher band
    while (double(histogram->GetFrequency( x2, 0))/maxBinValue > 0.5 )
        x2++;
    double hx2,hx2m1,mx2,mx2m1, hiband;
    hx2=double(histogram->GetFrequency( x2, 0))/maxBinValue;
    hx2m1=double(histogram->GetFrequency( x2-1, 0))/maxBinValue;
    mx2=double(histogram->GetMeasurement( x2, 0));
    mx2m1=double(histogram->GetMeasurement( x2-1, 0));
    hiband = mx2m1 -(0.5-hx2)*(mx2m1-mx2)/(hx2m1-hx2);
    sigma=(hiband-loband)/(2*vcl_sqrt(2*vcl_log(2.0)));

    ProbabilityPixelType t = mu + alpha*sigma;
    std::cout << "\tThreshold: " << t << " (" << mu << " + " << alpha << " * " << sigma << ")" << std::endl;
    std::cout << "\t\\- Lower band: " << loband << " (" << mx1 << ")" << std::endl;
    std::cout << "\t\\- Higher band: " << hiband << " (" << mx2 << ")" << std::endl;

    /* Thresholding and refinement */
    std::cout << "\t\\- Thresholding and refinement" << std::endl;

    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
    thresholder->SetLowerThreshold(t);
    //thresholder->SetInput(flair_e);
    thresholder->SetInput(flair);
    thresholder->Update();

    // Thresholding
    MaskImage wml = thresholder->GetOutput();


    // Lesion neighborhood refinement (WM label = 3)
    MaskIterator itWML = MaskIterator(wml,wml->GetRequestedRegion());

    ConnectComponentFilterType::Pointer connectComp = ConnectComponentFilterType::New();
    connectComp->SetInput(wml);
    connectComp->Update();

    ConnectedImage components = connectComp->GetOutput();
    MinimumMaximumLabelComponentType::Pointer maxComponentCalculator = MinimumMaximumLabelComponentType::New();
    maxComponentCalculator->SetImage(components);
    maxComponentCalculator->Compute();
    ConnectedPixelType maxComponent = maxComponentCalculator->GetMaximum();

    if (maxComponent > 0) {

        double *yesNeighbors = new double[maxComponent];
        double *noNeighbors = new double[maxComponent];
        double *yesTissue = new double[maxComponent];
        double *noTissue = new double[maxComponent];
        int *nVoxels = new int[maxComponent];
        bool *checked = new bool[maxComponent];
        bool *centered = new bool[maxComponent];
        std::fill_n(yesNeighbors,maxComponent,0);
        std::fill_n(noNeighbors,maxComponent,0);
        std::fill_n(yesTissue,maxComponent,0);
        std::fill_n(noTissue,maxComponent,0);
        std::fill_n(nVoxels,maxComponent,0);
        std::fill_n(checked,maxComponent,false);
        std::fill_n(centered,maxComponent,false);

        ConnectedIterator itComp = ConnectedIterator(components,components->GetRequestedRegion());
        ConnectedNeighborhoodIterator::RadiusType radius;
        for (i = 0; i < IntensityImageType::ImageDimension; ++i)
            radius[i] = 1;
        ConnectedNeighborhoodIterator itNeighbors = ConnectedNeighborhoodIterator(radius,components,components->GetRequestedRegion());
        ConnectedPixelType lesionId;

        // The idea is to check at each voxel if it is a lesion neighbour. If that is the case, we check if the voxel is a permitted tissue voxel (WM) or not.
        // This will allow us to compute the WM fraction value.
        for (itLabel.GoToBegin();!itLabel.IsAtEnd();++itLabel,++itComp,++itNeighbors) {
            if (itComp.Get()==0) {
                for (i = 0; i < itNeighbors.Size(); ++i) {
                    lesionId = itNeighbors.GetPixel(i);
                    if ( lesionId > 0 && !checked[lesionId-1]) {
                        checked[lesionId-1] = true;
                        switch (itLabel.Get()) {
                        // WM
                        case 3:
                            yesNeighbors[lesionId-1]++;
                            break;
                        // Other tissues
                        default:
                            if ((itLabel.Get()-4>=0) && (pv.at(itLabel.Get()-4)[0]==3||pv.at(itLabel.Get()-4)[1]==3))
                                yesTissue[lesionId-1]++;
                            else
                                noNeighbors[lesionId-1]++;
                            break;
                        }
                    }
                }
            }
            else {
                lesionId = itComp.Get();
                switch (itLabel.Get()) {
                // WM + GM-based
                case 2:
                case 3:
                    yesTissue[lesionId-1]++;
                    break;
                // Other tissues
                default:
                    if ((itLabel.Get()-4>=0) && (pv.at(itLabel.Get()-4)[0]==2|pv.at(itLabel.Get()-4)[0]==3|pv.at(itLabel.Get()-4)[1]==2|pv.at(itLabel.Get()-4)[1]==3))
                        yesTissue[lesionId-1]++;
                    else
                        noTissue[lesionId-1]++;
                    break;
                }
                nVoxels[lesionId-1]++;
                index = itComp.GetIndex();
                centered[lesionId-1] = centered[lesionId-1] |\
                                       (((index[0] > imgSize[0]/2-3) & (index[0] < imgSize[0]/2+4)) &\
                                       ((index[1] > imgSize[1]/2-3) & (index[1] < imgSize[1]/2+4)));
            }
            std::fill_n(checked,maxComponent,false);
        }

        // Finally we relabel those lesions with a WM fraction value lower than beta
        bool neighborsCondition, tissueCondition, sizeCondition, centeredCondition;
        for (itComp.GoToBegin();!itComp.IsAtEnd();++itComp,++itWML) {
            lesionId = itComp.Get();
            if (lesionId>0) {
                // Neighbouring condition
                neighborsCondition = yesNeighbors[lesionId-1] > (yesNeighbors[lesionId-1] + noNeighbors[lesionId-1])*omegaN;
                // Tissue condition
                tissueCondition = yesTissue[lesionId-1] > (yesTissue[lesionId-1] + noTissue[lesionId-1])*omegaT;
                // Size condition
                sizeCondition = nVoxels[lesionId-1]>=minSize;
                // Centered condition
                centeredCondition = !centered[lesionId-1];

                itWML.Set(  neighborsCondition & tissueCondition & sizeCondition & centeredCondition );
                //itWML.Set( true );

            }
        }

    }

	delete(brainio);

    return(wml);
}


/*************************************************************************************************************************************************/
//
// Function to segment lesions based on de Boer et al. and Souplet et al. approaches. This method uses a previous tissue segmentation and
// lesion mask.
//
// Inputs:
//  tissue = Tissue segmentation image.
//     wml = Lesion segmentation image.
//
/*************************************************************************************************************************************************/
ResultsImage BrainSegmentation::Relabel(ResultsImage tissue, MaskImage wml) {
	std::cout << "\tBrainSegmentation::Relabel" << std::endl;
    // Init
    ResultsImage solution = ResultsImageType::New();
    solution->SetRegions( tissue->GetLargestPossibleRegion() );
    solution->SetSpacing( tissue->GetSpacing() );
    solution->SetOrigin( tissue->GetOrigin() );
    solution->SetDirection( tissue->GetDirection() );
    solution->Allocate();
    solution->FillBuffer(0);
    solution->Update();
    ResultsIterator itLabel = ResultsIterator(tissue,tissue->GetRequestedRegion());
    MaskIterator itWML = MaskIterator(wml,wml->GetRequestedRegion());

    /* Image relabeling */
    std::cout << "\t\\- Relabeling to include lesions" << std::endl;
    MinimumMaximumLabelCalculatorType::Pointer maxLabelCalculator = MinimumMaximumLabelCalculatorType::New();
    maxLabelCalculator->SetImage(tissue);
    maxLabelCalculator->Compute();
    ResultsPixelType maxLabel = maxLabelCalculator->GetMaximum();

    ResultsIterator itSolution = ResultsIterator(solution,solution->GetRequestedRegion());
    for (itWML.GoToBegin(),itLabel.GoToBegin();!itWML.IsAtEnd();++itWML,++itSolution,++itLabel) {
        if (itWML.Get())
            itSolution.Set(maxLabel+1);
        else
            itSolution.Set(itLabel.Get());
    }

    std::cout << "\tRelabeling done" << std::endl;
    return(solution);
}

/*************************************************************************************************************************************************/
//
// Function used to segment Brain MRI images into pure tissues and partial volumes using an Expectation-Maximisation and Gaussian mixture models.
// This algorithm is initialised with atlases. Those atlases are also used during the probability computation as priors.
//
// Inputs:
//   atlas = Vector containing all atlas probability maps for pure tissues (GM, WM and CSF).
//   input = Vector containing all images used during the tissue segmentation step (T1, T2 and PD commonly).
//      bm = Brain mask containing the region of interest to be segmented.
//      pv = Vector cointaining pairs of pure tissue. Example: If we have an atlas vector with the following atlases:
//           [CSF GM  WM] in this order, a partial volume between CSF and GM would be defined by the pair array [1 2].
//      th = Threshold used during the initialisation of the pure tissues classes. Probabilities lower than th won't
//           be used to estimate the Gaussian parameters (mean and standard deviation/covariance matrix).
// maxIter = Number of maximum iteration for the Expectation-Maximisation algorithm.
//
/*************************************************************************************************************************************************/
std::vector<ProbabilityImage> BrainSegmentation::Generic_EM_Atlas(std::vector<ProbabilityImage> atlas, std::vector<ProbabilityImage> input, MaskImage bm, std::vector<int*> pv, ProbabilityPixelType th, int maxIter) {
	std::cout << "\tBrainSegmentation::Generic_EM_Atlas" << std::endl;
    /* Image initialisation */
     std::cout << "\t\\- Init" << std::endl;
    int nPV = pv.size();
    int nPure = atlas.size();
    int nClasses = nPure+nPV;
    ProbabilityImage* ipr = new ProbabilityImage[nClasses]; // Gaussian likelihood at step i
    double *pr = new double[nClasses];
    double prSum, sumlog =-100/10,sumlog_ant=-100;
    int i=0,c,c1,c2;
    int nImages = input.size();

    /* Brain mask iterator */
    MaskIterator itBM(bm, bm->GetRequestedRegion());

    /* Atlas iterators */
    ProbabilityIterator *prior = new ProbabilityIterator[nPure];

    /* Gaussian memberships */
    GaussianEstimator** gauss = new GaussianEstimator*[nClasses];
    for (c=0;c<nPure;c++) {
        switch(nImages) {
            case 1:
                    gauss[c] = new GaussianEstimator1D(input,atlas[c],th);
                    break;
            case 2:
                    gauss[c] = new GaussianEstimator2D(input,atlas[c],th);
                    break;
            case 3:
                    gauss[c] = new GaussianEstimator3D(input,atlas[c],th);
                    break;
            default:
                    gauss[c] = new GaussianEstimatorND(input,atlas[c],th);
                    break;
        }
        //gauss[c]->PrintParameters();
    }
    for (c=0;c<nPV;c++) {
        c1 = pv.at(c)[0];
        c2 = pv.at(c)[1];
        switch(nImages) {
            case 1:
                    gauss[c+nPure] = new GaussianEstimator1D();
                    break;
            case 2:
                    gauss[c+nPure] = new GaussianEstimator2D();
                    break;
            case 3:
                    gauss[c+nPure] = new GaussianEstimator3D();
                    break;
            default:
                    gauss[c+nPure] = new GaussianEstimatorND();
                    break;
        }
        gauss[c+nPure]->SetData(input);
        if (nImages==1)
            gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                  gauss[c1]->GetSigma()*(1/2.0)+gauss[c2]->GetSigma()*(1/2.0));
        else
            gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                      gauss[c1]->GetSigma()*(1/4.0)+gauss[c2]->GetSigma()*(1/4.0));
    }
    ProbabilityIterator *conditional = new ProbabilityIterator[nClasses];

     std::cout << "\t\\- Segmentation" << std::endl;

    while ((i<maxIter)&&(sumlog!=sumlog_ant)) {
        /* Expectation (E-step) */
        for (c=0;c<nClasses;c++) {
            if (c<nPure)
                prior[c] = ProbabilityIterator(atlas[c],atlas[c]->GetRequestedRegion());
            ipr[c] = gauss[c]->GetCurrentLikelihood();
            conditional[c] = ProbabilityIterator(ipr[c],ipr[c]->GetRequestedRegion());
        }

        itBM.GoToBegin();
        sumlog_ant=sumlog;
        sumlog = 0;
        while(!itBM.IsAtEnd()) {
            // Inside the mask
            if (itBM.Get()) {
                prSum = 0;
                // Posteriori probabilities (Conditional*Priors)
                for (c=0;c<nClasses;c++) {
                    if (c<nPure)
                        pr[c] = conditional[c].Get()*prior[c].Get();
                    else
                        pr[c] = conditional[c].Get()*0.5*(prior[pv.at(c-nPure)[0]].Get()+prior[pv.at(c-nPure)[1]].Get());
                    prSum += pr[c];
                }

                /* Probability normalisation */
                if (prSum<=0) {
                    // Same probability for each voxel
                    for (c=0;c<nClasses;c++) {
                        conditional[c].Set(1.0/nClasses);
                    }
                }
                else {
                    // Normalized probability
                    for (c=0;c<nClasses;c++) {
                        conditional[c].Set(pr[c]/prSum);
                    }
                    sumlog += vcl_log(prSum);
                }
            }
            // Not brain
            else {
                for (c=0;c<nClasses;c++) {
                    conditional[c].Set(0);
                }
            }

            ++itBM;
            for (c=0;c<nClasses;c++) {
                ++conditional[c];
                if (c<nPure)
                    ++prior[c];
            }

        }

        /* Maximisation (M-step) */
        for (c=0;c<nPure;c++) {
            gauss[c]->UpdateParameters(ipr[c],th);
        }
        for (c=0;c<nPV;c++) {
            c1 = pv.at(c)[0];
            c2 = pv.at(c)[1];
            if (nImages==1)
                gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                      gauss[c1]->GetSigma()*(1/2.0)+gauss[c2]->GetSigma()*(1/2.0));
            else
                gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                          gauss[c1]->GetSigma()*(1/4.0)+gauss[c2]->GetSigma()*(1/4.0));
        }


         std::cout << "\t\\--- Step " << i << " log-likelihood = " <<  sumlog << std::endl;
        i++;
    }

    std::vector<ProbabilityImage> prFinal;
    for (c=0;c<nClasses;c++) {
        //gauss[c]->PrintParameters();
        prFinal.push_back(ipr[c]);
    }

     std::cout << "\t\\--- Segmentation finished" << std::endl;

    /* Memory cleansing */
    for (c=0;c<nClasses;c++)
        delete gauss[c];
    delete [] gauss;
    delete [] pr;
    delete [] ipr;
    delete [] prior;
    delete [] conditional;
    //delete [] atlasPV;

    return (prFinal);
}

/*************************************************************************************************************************************************/
//
// Function used to segment Brain MRI images into pure tissues and partial volumes using an Expectation-Maximisation and Gaussian mixture models.
// This algorithm is initialised with atlases. Those atlases are also used during the probability computation as priors. A weighting image is used
// with the atlases. This image determines how good is the atlas value at each point
//
// Inputs:
//    atlas = Vector containing all atlas probability maps for pure tissues (GM, WM and CSF).
//    input = Vector containing all images used during the tissue segmentation step (T1, T2 and PD commonly).
//   weight = Weighting image for the atlases.
//       bm = Brain mask containing the region of interest to be segmented.
//       pv = Vector cointaining pairs of pure tissue. Example: If we have an atlas vector with the following atlases:
//            [CSF GM  WM] in this order, a partial volume between CSF and GM would be defined by the pair array [1 2].
//       th = Threshold used during the initialisation of the pure tissues classes. Probabilities lower than th won't
//            be used to estimate the Gaussian parameters (mean and standard deviation/covariance matrix).
//  maxIter = Number of maximum iteration for the Expectation-Maximisation algorithm.
//
/*************************************************************************************************************************************************/
std::vector<ProbabilityImage> BrainSegmentation::Generic_WEM_Atlas(std::vector<ProbabilityImage> atlas, std::vector<ProbabilityImage> input, ProbabilityImage weights, MaskImage bm, std::vector<int*> pv, ProbabilityPixelType th, int maxIter) {
	std::cout << "\tBrainSegmentation::Generic_WEM_Atlas" << std::endl;
    /* Image initialisation */
    std::cout << "\t\\- Init" << std::endl;
    int nPV = pv.size();
    int nPure = atlas.size();
    int nClasses = nPure+nPV;
    double *pr = new double[nClasses];
    double prSum, sumlog =-100/10,sumlog_ant=-100;
    int i=0,c,c1,c2;
    int nImages = input.size();

    ProbabilityImage* ipr = new ProbabilityImage[nClasses]; // Gaussian likelihood at step i
    ProbabilityImage* npr = new ProbabilityImage[nClasses]; // Mean of the neighboring probabilities
    ProbabilityImage* ppr = new ProbabilityImage[nClasses]; // Posterior probability
    for (c = 0; c < nClasses; c++) {
        ppr[c] = ProbabilityImageType::New();
        ppr[c]->SetRegions( bm->GetLargestPossibleRegion() );
        ppr[c]->SetSpacing( bm->GetSpacing() );
        ppr[c]->SetOrigin( bm->GetOrigin() );
        ppr[c]->SetDirection( bm->GetDirection() );
        ppr[c]->Allocate();
        ppr[c]->FillBuffer(0);
        ppr[c]->Update();
    }


    /* Brain mask iterator */
    MaskIterator itBM(bm, bm->GetLargestPossibleRegion());

    /* Prior iterator */
    ProbabilityIterator w = ProbabilityIterator(weights, weights->GetRequestedRegion());
    ProbabilityIterator *nprior = new ProbabilityIterator[nClasses];
    ProbabilityIterator *aprior = new ProbabilityIterator[nPure];
    ProbabilityPixelType prior;

    /* Gaussian memberships */
    GaussianEstimator** gauss = new GaussianEstimator*[nClasses];
    for (c=0;c<nPure;c++) {
        switch(nImages) {
            case 1:
                    gauss[c] = new GaussianEstimator1D(input,atlas[c],th);
                    break;
            case 2:
                    gauss[c] = new GaussianEstimator2D(input,atlas[c],th);
                    break;
            case 3:
                    gauss[c] = new GaussianEstimator3D(input,atlas[c],th);
                    break;
            default:
                    gauss[c] = new GaussianEstimatorND(input,atlas[c],th);
                    break;
        }
        //gauss[c]->PrintParameters();
        npr[c] = NeighborhoodAverage(atlas[c],2);
        /*sprintf(name,"cprinic%d.nii",c);
        brainio->WriteProbabilityImage(name,gauss[c]->GetCurrentLikelihood());
        sprintf(name,"nprinic%d.nii",c);
        brainio->WriteProbabilityImage(name,npr[c]);*/
    }
    for (c=0;c<nPV;c++) {
        AddFilterType::Pointer sum = AddFilterType::New();
        MultiplyConstantFilterType::Pointer const1 = MultiplyConstantFilterType::New();
        MultiplyConstantFilterType::Pointer const2 = MultiplyConstantFilterType::New();

        c1 = pv.at(c)[0];
        c2 = pv.at(c)[1];
        switch(nImages) {
            case 1:
                    gauss[c+nPure] = new GaussianEstimator1D();
                    break;
            case 2:
                    gauss[c+nPure] = new GaussianEstimator2D();
                    break;
            case 3:
                    gauss[c+nPure] = new GaussianEstimator3D();
                    break;
            default:
                    gauss[c+nPure] = new GaussianEstimatorND();
                    break;
        }
        gauss[c+nPure]->SetData(input);
        if (nImages==1)
            gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                  gauss[c1]->GetSigma()*(1/2.0)+gauss[c2]->GetSigma()*(1/2.0));
        else
            gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                      gauss[c1]->GetSigma()*(1/4.0)+gauss[c2]->GetSigma()*(1/4.0));
        const1->SetInput(atlas[c1]);
        const1->SetConstant(0.5);
        const1->Update();
        const2->SetInput(atlas[c2]);
        const2->SetConstant(0.5);
        const2->Update();
        sum->SetInput1(const1->GetOutput());
        sum->SetInput2(const2->GetOutput());
        sum->Update();

        npr[c+nPure] = NeighborhoodAverage(sum->GetOutput(),2);
    }
    ProbabilityIterator *conditional = new ProbabilityIterator[nClasses];
    ProbabilityIterator *posterior = new ProbabilityIterator[nClasses];

    std::cout << "\t\\- Segmentation" << std::endl;

    while ((i<maxIter)&&(sumlog!=sumlog_ant)) {
        /* Expectation (E-step) */
        for (c=0;c<nClasses;c++) {
            if (c<nPure)
                aprior[c] = ProbabilityIterator(atlas[c],atlas[c]->GetRequestedRegion());
            nprior[c] = ProbabilityIterator(npr[c],npr[c]->GetRequestedRegion());
            ipr[c] = gauss[c]->GetCurrentLikelihood();
            conditional[c] = ProbabilityIterator(ipr[c],ipr[c]->GetRequestedRegion());
            posterior[c] = ProbabilityIterator(ppr[c],ppr[c]->GetRequestedRegion());
        }

        sumlog_ant=sumlog;
        sumlog = 0;
        for(itBM.GoToBegin(),w.GoToBegin();!itBM.IsAtEnd();++itBM,++w) {
            // Inside the mask
            if (itBM.Get()) {
                prSum = 0;
                // Posteriori probabilities (Conditional*Priors)
                for (c=0;c<nClasses;c++) {
                    if (c<nPure)
                        prior = (aprior[c].Get()*w.Get() + (1-w.Get())*nprior[c].Get());
                    else
                        prior = (w.Get()*(0.5*(aprior[pv.at(c-nPure)[0]].Get()+aprior[pv.at(c-nPure)[1]].Get())) +\
                                (1-w.Get())*nprior[c].Get());
                    pr[c] = conditional[c].Get()*prior;
                    prSum += pr[c];
                }

                /* Probability normalisation */
                if (prSum<=0) {
                    // Same probability for each voxel
                    for (c=0;c<nClasses;c++) {
                        posterior[c].Set(1.0/nClasses);
                    }
                }
                else {
                    // Normalized probability
                    for (c=0;c<nClasses;c++) {
                        posterior[c].Set(pr[c]/prSum);
                    }
                    sumlog += vcl_log(prSum);
                }
            }
//            // Not brain
//            else {
//                for (c=0;c<nClasses;c++) {
//                    posterior[c].Set(0);
//                }
//            }

            for (c=0;c<nClasses;c++) {
                ++conditional[c];
                ++posterior[c];
                ++nprior[c];
                if (c<nPure)
                    ++aprior[c];
            }

        }

        /* Maximisation (M-step) */
        for (c = 0; c < nClasses; c++) {
            npr[c] = NeighborhoodAverage(ppr[c],2);
            //npr[c] = NeighborhoodAverage(ipr[c],2);
            if (c < nPure)
                gauss[c]->UpdateParameters(ppr[c],th);
            //gauss[c]->PrintParameters();
        }
        for (c=0;c<nPV;c++) {
            c1 = pv.at(c)[0];
            c2 = pv.at(c)[1];
            if (nImages==1)
                gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                      gauss[c1]->GetSigma()*(1/2.0)+gauss[c2]->GetSigma()*(1/2.0));
            else
                gauss[c+nPure]->SetParameters(gauss[c1]->GetMu()*(1/2.0)+gauss[c2]->GetMu()*(1/2.0),
                                          gauss[c1]->GetSigma()*(1/4.0)+gauss[c2]->GetSigma()*(1/4.0));
        }

        std::cout << "\t\\--- Step " << i << " log-likelihood = " <<  sumlog << std::endl;
        i++;
    }

    std::vector<ProbabilityImage> prFinal;
    for (c=0;c<nClasses;c++) {
        gauss[c]->PrintParameters();
        prFinal.push_back(ppr[c]);
    }

    std::cout << "\t\\--- Segmentation finished" << std::endl;

    /* Memory cleansing */
    for (c=0;c<nClasses;c++)
        delete gauss[c];
    delete [] gauss;
    delete [] pr;
	delete [] ppr;
	delete [] npr;
    delete [] ipr;
    delete [] aprior;
    delete [] nprior;
    delete [] conditional;
	delete [] posterior;

    return (prFinal);
}

/*************************************************************************************************************************************************/
//
// Simple function to obtain a hard segmentation from a set of soft segmentations (probabilities).
//
// Inputs:
// pr = Soft segmentations (probability images). It is assumed that for any voxel, the sum of all probabilities equals to 1.
// bm = Brain mask.
//
/*************************************************************************************************************************************************/
ResultsImage BrainSegmentation::SolutionFromPr(std::vector<ProbabilityImage> pr, MaskImage bm) {
    /* Labeled image preparation */
    double prMax, prActual;
    int i, labelActual, labelMax, classes=pr.size();
    ResultsImage solution = ResultsImageType::New();
    solution->SetRegions( bm->GetLargestPossibleRegion() );
    solution->SetSpacing( bm->GetSpacing() );
    solution->SetOrigin( bm->GetOrigin() );
    solution->SetDirection( bm->GetDirection() );
    solution->Allocate();
    solution->FillBuffer(0);
    solution->Update();
    /* Brain mask */
    MaskIterator itBM(bm, bm->GetRequestedRegion());

    /* Probability maps */
    ProbabilityIterator* its = new ProbabilityIterator[classes];
    for (i=0;i<classes;i++) {
        its[i] = ProbabilityIterator(pr[i],pr[i]->GetRequestedRegion());
    }

    ResultsIterator it( solution, solution->GetRequestedRegion() );
    while (!it.IsAtEnd()) {
        prMax = 0;
        labelMax = 0;
        for (labelActual=1;labelActual<=classes;labelActual++) {
            if (itBM.Get()) {
                prActual = its[labelActual-1].Get();
                if (prActual>prMax) {
                    labelMax = labelActual;
                    prMax = prActual;
                }
            }
            ++its[labelActual-1];
        }
        it.Set(labelMax);
        ++it;
        ++itBM;
    }

	delete [] its;
    return(solution);
}

/*************************************************************************************************************************************************/
//
// Simple function to create a WM mask from a tissue segmentation and a lesion mask.
//
// Inputs:
//   brain = Tissue segmentation.
//   wmlab = Label for the white matter tissue.
//   wmlab = Label for the white matter lesion label.
//
/*************************************************************************************************************************************************/
MaskImage BrainSegmentation::GetWMMask(ResultsImage brain, ResultsPixelType wmlab, ResultsPixelType wmllab) {
	std::cout << "\tBrainSegmentation::GetWMMask" << std::endl;
    ResultsPixelType maxLabel, currentLabel;
    int maxPixelSize=0;
    MinimumMaximumLabelCalculatorType::Pointer maxLabelCalculator = MinimumMaximumLabelCalculatorType::New();
    maxLabelCalculator->SetImage(brain);
    maxLabelCalculator->Compute();
    maxLabel = maxLabelCalculator->GetMaximum();

    /* We compute the CSF/BCK mask */
    MaskImage notWM = MaskImageType::New();
    notWM->SetRegions( brain->GetLargestPossibleRegion() );
    notWM->SetSpacing( brain->GetSpacing() );
    notWM->SetOrigin( brain->GetOrigin() );
    notWM->SetDirection( brain->GetDirection() );
    notWM->Allocate();
    notWM->FillBuffer(false);
    notWM->Update();

    MaskIterator notWMIt = MaskIterator(notWM, notWM->GetLargestPossibleRegion());
    ResultsIterator brainIt = ResultsIterator(brain, brain->GetLargestPossibleRegion());
	for (brainIt.GoToBegin(); !brainIt.IsAtEnd(); ++notWMIt, ++brainIt) {
        currentLabel = brainIt.Get();
        notWMIt.Set(currentLabel!=wmlab);
    }

    /* We get the components of the image in order to find the biggest one */
    ConnectComponentFilterType::Pointer connectComp1 = ConnectComponentFilterType::New();
    connectComp1->SetInput(notWM);
    connectComp1->Update();

	RelabelComponentFilterType::Pointer relabelComp1 = RelabelComponentFilterType::New();
	relabelComp1->SetInput(connectComp1->GetOutput());
	relabelComp1->Update();

    /* We now compute the final mask negating the not WM mask */
    MaskImage WM = MaskImageType::New();
    WM->SetRegions( brain->GetLargestPossibleRegion() );
    WM->SetSpacing( brain->GetSpacing() );
    WM->SetOrigin( brain->GetOrigin() );
    WM->SetDirection( brain->GetDirection() );
    WM->Allocate();
    WM->FillBuffer(false);
    WM->Update();
    MaskIterator WMIt = MaskIterator(WM, WM->GetLargestPossibleRegion());
    ConnectedIterator componentIt = ConnectedIterator(relabelComp1->GetOutput(), relabelComp1->GetOutput()->GetLargestPossibleRegion());
	for (componentIt.GoToBegin(), brainIt.GoToBegin(); !componentIt.IsAtEnd(); ++componentIt, ++WMIt, ++brainIt) {
		WMIt.Set((componentIt.Get() != 1) | (brainIt.Get() == wmllab));
    }

	/* Finally we get the biggest component for the WM mask (to reduce the effect of oversegmentation) */
	ConnectComponentFilterType::Pointer connectComp2 = ConnectComponentFilterType::New();
    connectComp2->SetInput(WM);
    connectComp2->Update();

	RelabelComponentFilterType::Pointer relabelComp2 = RelabelComponentFilterType::New();
	relabelComp2->SetInput(connectComp2->GetOutput());
	relabelComp2->Update();
	
	WMIt = MaskIterator(WM, WM->GetLargestPossibleRegion());
	componentIt = ConnectedIterator(relabelComp2->GetOutput(), relabelComp2->GetOutput()->GetLargestPossibleRegion());
	for (componentIt.GoToBegin(), brainIt.GoToBegin(); !componentIt.IsAtEnd(); ++componentIt, ++WMIt, ++brainIt) {
		WMIt.Set(componentIt.Get() == 1);
    }

    return(WM);
}

/*************************************************************************************************************************************************/
//
// Simple function to enhance the contrast of an image using morphological operations. This algorithm is based on the approach of Souplet et al.
//
// Inputs:
//  input = Image to be enhanced. Usually a FLAIR image is used.
// radius = Radius used for the morphological operations. Higher radiuses cause more blur.
//
/*************************************************************************************************************************************************/
ProbabilityImage BrainSegmentation::ContrastEnhancement(ProbabilityImage input, unsigned long radius) {
    ProbabilityImage solution = ProbabilityImageType::New();
    solution->SetRegions( input->GetLargestPossibleRegion() );
    solution->SetSpacing( input->GetSpacing() );
    solution->SetOrigin( input->GetOrigin() );
    solution->SetDirection( input->GetDirection() );
    solution->Allocate();
    solution->FillBuffer(0.0);
    solution->Update();

    BinaryBallStructuringElementType se;
    se.SetRadius(radius);
    se.CreateStructuringElement();

    GrayscaleDilateFilterType::Pointer dilate = GrayscaleDilateFilterType::New();
    dilate->SetKernel(se);
    dilate->SetInput(input);
    dilate->Update();
    GrayscaleErodeFilterType::Pointer erode = GrayscaleErodeFilterType::New();
    erode->SetKernel(se);
    erode->SetInput(input);
    erode->Update();

    ProbabilityImage Dima = dilate->GetOutput();
    ProbabilityImage Eima = erode->GetOutput();

    ProbabilityIterator itIma = ProbabilityIterator(input, input->GetRequestedRegion());
    ProbabilityIterator itEIma = ProbabilityIterator(Eima, Eima->GetRequestedRegion());
    ProbabilityIterator itDIma = ProbabilityIterator(Dima, Dima->GetRequestedRegion());
    ProbabilityIterator itCIma = ProbabilityIterator(solution, solution->GetRequestedRegion());

    for (itIma.GoToBegin();!itIma.IsAtEnd();++itIma,++itEIma,++itDIma,++itCIma) {
        if (itDIma.Get() - itIma.Get() <= itIma.Get() - itEIma.Get()) {
            itCIma.Set(itDIma.Get());
        }
        else if (itIma.Get() - itEIma.Get() <= itDIma.Get() - itIma.Get()) {
            itCIma.Set(itEIma.Get());
        }
        else {
            itCIma.Set(itIma.Get());
        }
    }

    return(solution);
}


/*************************************************************************************************************************************************/
//
// Simple function to compute the average probability of a neighborhood.
//
// Inputs:
//  input = Probability values for a tissue class
// radius = Radius used for the morphological operations. Higher radiuses cause more blur.
//
/*************************************************************************************************************************************************/
ProbabilityImage BrainSegmentation::NeighborhoodAverage(ProbabilityImage input, unsigned long radiusSize) {
    unsigned int i;
    ProbabilityImage solution = ProbabilityImageType::New();
    solution->SetRegions( input->GetLargestPossibleRegion() );
    solution->SetSpacing( input->GetSpacing() );
    solution->SetOrigin( input->GetOrigin() );
    solution->SetDirection( input->GetDirection() );
    solution->Allocate();
    solution->FillBuffer(0.0);
    solution->Update();
    
    ProbabilityNeighborhoodIterator::RadiusType radius;
    for (i = 0; i < ProbabilityImageType::ImageDimension; ++i)
        radius[i] = radiusSize;
    ProbabilityNeighborhoodIterator neighbors = ProbabilityNeighborhoodIterator(radius,input,input->GetRequestedRegion());
    ProbabilityIterator it = ProbabilityIterator(solution,solution->GetRequestedRegion());
    ProbabilityPixelType sum;

    for (neighbors.GoToBegin();!neighbors.IsAtEnd();++neighbors,++it) {
        sum = 0;
        for (i = 0; i < neighbors.Size(); i++) {
            if (i != ceil(neighbors.Size()/2.0))
                sum += neighbors.GetPixel(i);
        }
        it.Set(sum/(neighbors.Size()-1));
    }

    return(solution);
}

/*************************************************************************************************************************************************/
//
// Simple function to compute the intersection between the masks of the fixed (patient image) and moving (atlas) images
//
// Inputs:
//  init = Mask image for the patient.
// atlas = Atlas probabilities.
//
/*************************************************************************************************************************************************/
MaskImage BrainSegmentation::RefineMask(MaskImage init, std::vector<ProbabilityImage> atlas) {
    /* Init */
    unsigned int c;
    ProbabilityIterator* atlas_pr = new ProbabilityIterator[atlas.size()];
    MaskIterator bm = MaskIterator(init, init->GetRequestedRegion());
    ProbabilityPixelType pixelSum;

    for ( c = 0; c < atlas.size(); c++ )
        atlas_pr[c] = ProbabilityIterator( atlas[c], atlas[c]->GetRequestedRegion() );

    for ( bm.GoToBegin(); !bm.IsAtEnd(); ++bm) {
        pixelSum = 0;
        for ( c = 0; c < atlas.size(); c++ ) {
            pixelSum += atlas_pr[c].Get();
            ++atlas_pr[c];
        }
        bm.Set(bm.Get() & (pixelSum > 0.5));
    }

	delete [] atlas_pr;
    return( init );
}

/*************************************************************************************************************************************************/
//
// Function to compute the integral image of a source image. This image is used to compute fast means of regions on a volume.
//
// Inputs:
//  source = Source image. It can be any type of image (probabilistic, map, intensity...).
//
/*************************************************************************************************************************************************/
ProbabilityImage BrainSegmentation::IntegralImage(ProbabilityImage source) {
    ProbabilityImage integral = ProbabilityImageType::New();
    integral->SetRegions( source->GetLargestPossibleRegion() );
    integral->SetSpacing( source->GetSpacing() );
    integral->SetOrigin( source->GetOrigin() );
    integral->SetDirection( source->GetDirection() );
    integral->Allocate();
    integral->FillBuffer(1.0);
    integral->Update();

    ProbabilityIterator itSource = ProbabilityIterator(source, source->GetRequestedRegion());
    ProbabilityIterator it = ProbabilityIterator(integral, integral->GetRequestedRegion());
    ProbabilityImageType::IndexType pixel;
    ProbabilityPixelType x000,x001,x010,x011,x100,x101,x110;

    int i, j, k;
    for (it.GoToBegin();!it.IsAtEnd();++it,++itSource) {
        i = it.GetIndex()[0];
        j = it.GetIndex()[1];
        k = it.GetIndex()[2];
        x000 = 0;
        x001 = 0;
        x010 = 0;
        x011 = 0;
        x100 = 0;
        x101 = 0;
        x110 = 0;
        pixel[0] = i-1;
        pixel[1] = j-1;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x000 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j-1;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x001 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x010 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x011 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j-1;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x100 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j-1;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x101 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x110 = integral->GetPixel(pixel);

        it.Set(itSource.Get()+x110+x101+x011+x000-x001-x010-x100);
    }


    return(integral);
}


/*************************************************************************************************************************************************/
//
// Function to compute the integral image of a source image. This image is used to compute fast standard deviations of regions on a volume.
//
// Inputs:
//  source = Source image. It can be any type of image (probabilistic, map, intensity...).
//
/*************************************************************************************************************************************************/
ProbabilityImage BrainSegmentation::IntegralSquaredImage(ProbabilityImage source) {
    ProbabilityImage integral = ProbabilityImageType::New();
    integral->SetRegions( source->GetLargestPossibleRegion() );
    integral->SetSpacing( source->GetSpacing() );
    integral->SetOrigin( source->GetOrigin() );
    integral->SetDirection( source->GetDirection() );
    integral->Allocate();
    integral->FillBuffer(1.0);
    integral->Update();

    ProbabilityIterator itSource = ProbabilityIterator(source, source->GetRequestedRegion());
    ProbabilityIterator it = ProbabilityIterator(integral, integral->GetRequestedRegion());
    ProbabilityImageType::IndexType pixel;
    ProbabilityPixelType x000,x001,x010,x011,x100,x101,x110;

    int i, j, k;
    for (it.GoToBegin();!it.IsAtEnd();++it,++itSource) {
        i = it.GetIndex()[0];
        j = it.GetIndex()[1];
        k = it.GetIndex()[2];
        x000 = 0;
        x001 = 0;
        x010 = 0;
        x011 = 0;
        x100 = 0;
        x101 = 0;
        x110 = 0;
        pixel[0] = i-1;
        pixel[1] = j-1;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x000 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j-1;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x001 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x010 = integral->GetPixel(pixel);
        pixel[0] = i-1;
        pixel[1] = j;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x011 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j-1;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x100 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j-1;
        pixel[2] = k;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x101 = integral->GetPixel(pixel);
        pixel[0] = i;
        pixel[1] = j;
        pixel[2] = k-1;
        if (pixel[0]>=0 && pixel[1]>=0 && pixel[2]>=0)
            x110 = integral->GetPixel(pixel);

        it.Set(itSource.Get()*itSource.Get()+x110+x101+x011+x000-x001-x010-x100);
    }


    return(integral);
}

MaskImage BrainSegmentation::Intersection(MaskImage maskA, MaskImage maskB) {
	std::cout << "\tBrainSegmentation::Intersection" << std::endl;
	AndFilterType::Pointer andFilter = AndFilterType::New();

	andFilter->SetInput1(maskA);
	andFilter->SetInput2(maskB);

	andFilter->Update();

	return(andFilter->GetOutput());
}

MaskImage BrainSegmentation::Intersection(std::vector<MaskImage> masks) {
	std::cout << "\tBrainSegmentation::Intersection" << std::endl;
	std::vector<MaskImage>::iterator it;
    MaskImage intersection = BrainIO::Initialise<MaskImageType, MaskImageType>(masks[0]);
	intersection->FillBuffer(true);
	intersection->Update();

	for (it = masks.begin(); it != masks.end(); ++it) {
		AndFilterType::Pointer andFilter = AndFilterType::New();

		andFilter->SetInput1(intersection);
		andFilter->SetInput2(*it);

		andFilter->Update();

		intersection = andFilter->GetOutput();
	}

	return(intersection);
}

MaskImage BrainSegmentation::Union(MaskImage maskA, MaskImage maskB) {
	OrFilterType::Pointer orFilter = OrFilterType::New();

	orFilter->SetInput1(maskA);
	orFilter->SetInput2(maskB);

	orFilter->Update();

	return(orFilter->GetOutput());
}

ProbabilityImage BrainSegmentation::Subtraction(ProbabilityImage minuend, ProbabilityImage subtrahend) {
	SubtractionFilterType::Pointer subtractor = SubtractionFilterType::New();

	subtractor->SetInput1( minuend );
	subtractor->SetInput2( subtrahend );

	subtractor->Update();

	return(subtractor->GetOutput());

}

ProbabilityImage BrainSegmentation::Ratio(ProbabilityImage numerator, ProbabilityImage denominator) {
	DivideFilterType::Pointer divider = DivideFilterType::New();

	divider->SetInput1( numerator );
	divider->SetInput2( denominator );

	divider->Update();

	return(divider->GetOutput());

}

MaskImage BrainSegmentation::ThresholdImage(ProbabilityImage input, ProbabilityPixelType threshold) {
    ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
	thresholder->SetLowerThreshold(threshold);
    thresholder->SetInput(input);
    thresholder->Update();

	return(thresholder->GetOutput());
}

MaskImage BrainSegmentation::GetPositiveActivity(ProbabilityImage subtraction, double alpha) {
    std::cout << "\tBrainSegmentation::GetPositiveActivity" << std::endl;
	std::vector<ProbabilityPixelType> pixels;
	ProbabilityPixelType mean = 0, std = 0;
	ProbabilityIterator subIt = ProbabilityIterator(subtraction, subtraction->GetLargestPossibleRegion());

	// We compute the mean and standard deviation of the positive activity.
	// This values define the threshold (by default 5 sigmas over the mean).
	for (subIt.GoToBegin(); !subIt.IsAtEnd(); ++subIt) {
		if (subIt.Get()>0) {
			mean += subIt.Get();
			pixels.push_back(subIt.Get());
		}
	}

	if (pixels.size()>0)
		mean /= pixels.size();
	
	for (std::vector<ProbabilityPixelType>::iterator pixIt = pixels.begin() ; pixIt != pixels.end(); ++pixIt)
		std += (*pixIt - mean)*(*pixIt - mean);

	if (pixels.size()>0) {
		std /= pixels.size();
		std = sqrt(std);
	}

	double threshold = mean + alpha * std;

	return(ThresholdImage(subtraction, threshold));
}
MaskImage BrainSegmentation::GetVoxelSelection(ProbabilityImage subtraction, double alpha) {
    std::cout << "\tBrainSegmentation::GetPositiveActivity" << std::endl;
    std::vector<ProbabilityPixelType> pixels;
    ProbabilityPixelType mean = 0, std = 0;
    ProbabilityIterator subIt = ProbabilityIterator(subtraction, subtraction->GetLargestPossibleRegion());

    // We compute the mean and standard deviation of the positive activity.
    // This values define the threshold (by default 5 sigmas over the mean).
    for (subIt.GoToBegin(); !subIt.IsAtEnd(); ++subIt) {
        if (subIt.Get()>0) {
            mean += subIt.Get();
            pixels.push_back(subIt.Get());
        }
    }

    if (pixels.size()>0)
        mean /= pixels.size();
    
    for (std::vector<ProbabilityPixelType>::iterator pixIt = pixels.begin() ; pixIt != pixels.end(); ++pixIt)
        std += (*pixIt - mean)*(*pixIt - mean);

    if (pixels.size()>0) {
        std /= pixels.size();
        std = sqrt(std);
    }

    double threshold = mean + alpha * std;

    return(ThresholdImage(subtraction, threshold));
}

MaskImage BrainSegmentation::OnurPostProcessing(MaskImage activity, ProbabilityImage basal, ProbabilityImage following, int minSize, double bnr_t, double fnr_t) {
    typedef std::vector<ProbabilityPixelType> ProbabilityVector;
    std::cout << "\tBrainSegmentation::OnurPostProcessing" << std::endl;
    MaskImage final = BrainIO::Initialise<MaskImageType,MaskImageType>(activity);
    final->FillBuffer(false);
    final->Update();

    /* -- Lesion labeling -- */
    ConnectComponentFilterType::Pointer connectedFilter = ConnectComponentFilterType::New();
    connectedFilter->SetInput(activity);
    connectedFilter->Update();

    ConnectedImage lesions = connectedFilter->GetOutput();

    MinimumMaximumLabelComponentType::Pointer maxComponentCalculator = MinimumMaximumLabelComponentType::New();
    maxComponentCalculator->SetImage(lesions);
    maxComponentCalculator->Compute();
    ConnectedPixelType maxLabel = maxComponentCalculator->GetMaximum();

    /* -- Intensity compiling and lesion size computation -- */
    // Mean, standard deviation and size per lesion
    double global_mean = 0, global_lesion_size = 0, global_std = 0;
    int *lesion_size = new int[maxLabel];
    double *mean_basal = new double[maxLabel];
    double *mean_neigh_basal = new double[maxLabel];
    double *mean_following = new double[maxLabel];
    double *mean_neigh_following = new double[maxLabel];
    double *bnr = new double[maxLabel];
    double *fnr = new double[maxLabel];
    bool *checked = new bool[maxLabel];
    std::fill_n(lesion_size, maxLabel, 0);
    std::fill_n(mean_basal, maxLabel, 0);
    std::fill_n(mean_neigh_basal, maxLabel, 0);
    std::fill_n(mean_following, maxLabel, 0);
    std::fill_n(mean_neigh_following, maxLabel, 0);
    std::fill_n(bnr, maxLabel, 0);
    std::fill_n(fnr, maxLabel, 0);
    std::fill_n(checked, maxLabel, false);

    // Intensity vector
    ProbabilityVector intensities;

    // Iterators
    ProbabilityIterator itBasal = ProbabilityIterator(basal, basal->GetLargestPossibleRegion());
    ProbabilityIterator itFollowing = ProbabilityIterator(following, following->GetLargestPossibleRegion());
    ConnectedIterator itComp = ConnectedIterator(lesions, lesions->GetLargestPossibleRegion());
    ConnectedNeighborhoodIterator::RadiusType radius;
    for (int i = 0; i < IntensityImageType::ImageDimension; ++i)
        radius[i] = 1;
    ConnectedNeighborhoodIterator itNeighbors = ConnectedNeighborhoodIterator(radius, lesions, lesions->GetLargestPossibleRegion());
    
    ConnectedPixelType label, neigh_label;

    std::cout << "\t\\- Mean computation + intensity storage" << std::endl;
    // Mean computation + intensity storage
    for (itComp.GoToBegin(); !itComp.IsAtEnd(); ++itComp, ++itBasal, ++itFollowing, ++itNeighbors) {
        label = itComp.Get();
        // We check the lesion intensities
        if (label > 0) {
            intensities.push_back(itBasal.Get());

            mean_basal[label-1] += itBasal.Get();
            mean_following[label-1] += itFollowing.Get();

            global_lesion_size++;
            global_mean += itBasal.Get();

            lesion_size[label-1]++;
        }
        // We check the neighbours intensity
        else {
            for (int i = 0; i < itNeighbors.Size(); ++i) {
                neigh_label = itNeighbors.GetPixel(i);
                if ( neigh_label > 0 && !checked[neigh_label-1]) {
                    checked[neigh_label-1] = true;
                    mean_neigh_basal[neigh_label-1] += itBasal.Get();
                    mean_neigh_following[neigh_label-1] += itFollowing.Get();
                }
            }
        }

        std::fill_n(checked, maxLabel, false);
    }

    global_mean /= global_lesion_size;

    for (label = 0; label < maxLabel; label++) {
        mean_basal[label] /= lesion_size[label];
        mean_neigh_basal[label] /= lesion_size[label];
        bnr[label] = mean_basal[label] / mean_neigh_basal[label];
        mean_following[label] /= lesion_size[label];
        mean_neigh_following[label] /= lesion_size[label];
        fnr[label] = mean_following[label] / mean_neigh_following[label];
    }

    for (ProbabilityVector::iterator voxel = intensities.begin(); voxel != intensities.end(); voxel++)
        global_std += (*voxel - global_mean) * (*voxel - global_mean);

    global_std /= global_lesion_size;
    global_std = sqrt(global_std);
    

    /*-- Thresholding (size, BNR, FNR and intensity) -- */
    std::cout << "\t\\- Thresholding" << std::endl;
    MaskIterator itFinal = MaskIterator(final, final->GetLargestPossibleRegion());
    bool bnr_condition, fnr_condition, size_condition, intensity_condition;
    for (itComp.GoToBegin(); !itComp.IsAtEnd(); ++itComp, ++itFinal) {
        label = itComp.Get();
        if (label > 0) {
            //bnr_condition = true;
            bnr_condition = bnr[label-1] > bnr_t;
            //fnr_condition = true;
            fnr_condition = fnr[label-1] > fnr_t;
            //size_condition = true;
            size_condition = lesion_size[label-1] > minSize;
            //intensity_condition = true;
            intensity_condition = mean_basal[label-1] > (global_mean - 2*global_std);
            itFinal.Set(bnr_condition & fnr_condition & size_condition & intensity_condition);
        }
    }
    std::cout << "\t\\- Done" << std::endl;

    delete[] lesion_size;
    delete[] mean_basal;
    delete[] mean_neigh_basal;
    delete[] mean_following;
    delete[] mean_neigh_following;
    delete[] bnr;
    delete[] fnr;
    delete[] checked;

    return(final);
}

MaskImage BrainSegmentation::DeformationPostProcessing(MaskImage activity, ProbabilityImage jacobian, ProbabilityImage divergence, int minSize, double jacobian_t, double divergence_t) {
	typedef std::vector<ProbabilityPixelType> ProbabilityVector;
    std::cout << "\tBrainSegmentation::OnurPostProcessing" << std::endl;
	MaskImage final = BrainIO::Initialise<MaskImageType,MaskImageType>(activity);
    final->FillBuffer(false);
    final->Update();

	/* -- Lesion labeling -- */
	ConnectComponentFilterType::Pointer connectedFilter = ConnectComponentFilterType::New();
	connectedFilter->SetInput(activity);
    connectedFilter->Update();

    ConnectedImage lesions = connectedFilter->GetOutput();

    MinimumMaximumLabelComponentType::Pointer maxComponentCalculator = MinimumMaximumLabelComponentType::New();
    maxComponentCalculator->SetImage(lesions);
    maxComponentCalculator->Compute();
    ConnectedPixelType maxLabel = maxComponentCalculator->GetMaximum();

	/* -- Intensity compiling and lesion size computation -- */
	// Mean, standard deviation and size per lesion
	int *lesion_size = new int[maxLabel];
	double *mean_jacobian = new double[maxLabel];
	double *mean_divergence = new double[maxLabel];
	std::fill_n(lesion_size, maxLabel, 0);
	std::fill_n(mean_jacobian, maxLabel, 0);
	std::fill_n(mean_divergence, maxLabel, 0);

	// Intensity vector
	ProbabilityVector intensities;

	// Iterators
	ProbabilityIterator itJacobian = ProbabilityIterator(jacobian, jacobian->GetLargestPossibleRegion());
	ProbabilityIterator itDivergence = ProbabilityIterator(divergence, divergence->GetLargestPossibleRegion());
	ConnectedIterator itComp = ConnectedIterator(lesions, lesions->GetLargestPossibleRegion());
	ConnectedPixelType label;

	std::cout << "\t\\- Mean computation + intensity storage" << std::endl;
	// Mean computation + intensity storage
	for (itComp.GoToBegin(); !itComp.IsAtEnd(); ++itComp, ++itJacobian, ++itDivergence) {
		label = itComp.Get();
		// We check the lesion intensities
		if (label > 0) {
			mean_jacobian[label-1] += itJacobian.Get();
			mean_divergence[label-1] += itDivergence.Get();
			lesion_size[label-1]++;
		}

	}

	for (label = 0; label < maxLabel; label++) {
		mean_jacobian[label] /= lesion_size[label];
		mean_divergence[label] /= lesion_size[label];
		std::cout << "\t\tLesion " << label+1 << ": (Jacobian = " << mean_jacobian[label] << ", Divergence = " << mean_divergence[label] << ")" << std::endl;
	}


	/*-- Thresholding (size, BNR, FNR and intensity) -- */
	std::cout << "\t\\- Thresholding" << std::endl;
	MaskIterator itFinal = MaskIterator(final, final->GetLargestPossibleRegion());
	bool jacobian_condition, divergence_condition, size_condition;
	for (itComp.GoToBegin(); !itComp.IsAtEnd(); ++itComp, ++itFinal) {
		label = itComp.Get();
		if (label > 0) {
			jacobian_condition = mean_jacobian[label-1] < jacobian_t;
			divergence_condition = mean_divergence[label-1] < divergence_t;
			size_condition = lesion_size[label-1] > minSize;
			itFinal.Set(jacobian_condition & divergence_condition & size_condition);
		}
	}

	delete[] lesion_size;
	delete[] mean_jacobian;
	delete[] mean_divergence;

	return(final);
}

ConnectedPixelType BrainSegmentation::CountRegions(MaskImage mask) {
	/* -- Lesion labeling -- */
	ConnectComponentFilterType::Pointer connectedFilter = ConnectComponentFilterType::New();
	connectedFilter->SetInput(mask);
    connectedFilter->Update();

    ConnectedImage lesions = connectedFilter->GetOutput();

    MinimumMaximumLabelComponentType::Pointer maxComponentCalculator = MinimumMaximumLabelComponentType::New();
    maxComponentCalculator->SetImage(lesions);
    maxComponentCalculator->Compute();
    ConnectedPixelType nRegions = maxComponentCalculator->GetMaximum();

	return(nRegions);
}

double BrainSegmentation::ComputeVolume(MaskImage mask) {
	/* -- Init -- */
	int nPixels = 0;
	MaskImageType::SpacingType spacing = mask->GetSpacing();
	double pixVolume = spacing[0]*spacing[1]*spacing[2];

	/* -- Pixel counting -- */
	MaskIterator maskIt = MaskIterator(mask, mask->GetLargestPossibleRegion());
	for (maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt) {
		if (maskIt.Get())
			nPixels++;
	}

	return(nPixels*pixVolume);
}
