#include <brainio.h>

BrainIO::BrainIO()
{
}

IntensityImage BrainIO::ReadIntensityImage(FileNameType fileName) {
    // Init
    IntensityReaderType::Pointer intreader = IntensityReaderType::New();;
    IntensityImageType::Pointer intimage = NULL;

    if (fileName.compare("")!=0) {
        intreader->SetFileName( fileName );
        // Reading the image.
        try {
            intreader->Update();
            intimage = intreader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(intimage);
}

ProbabilityImage BrainIO::ReadProbabilityImage(FileNameType fileName) {
    // Init
    ProbabilityReaderType::Pointer freader = ProbabilityReaderType::New();
    ProbabilityImageType::Pointer fimage = NULL;

    if (fileName.compare("")!=0) {
        freader->SetFileName( fileName );
        // Reading the image.
        try {
            freader->Update();
            fimage = freader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(fimage);
}

ConnectedImage BrainIO::ReadConnectedImage(FileNameType fileName) {
    // Init
    ConnectedReaderType::Pointer conreader = ConnectedReaderType::New();
    ConnectedImageType::Pointer conimage = NULL;

    if (fileName.compare("")!=0) {
        conreader->SetFileName( fileName );
        // Reading the image.
        try {
            conreader->Update();
            conimage = conreader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(conimage);
}

ResultsImage BrainIO::ReadResultsImage(FileNameType fileName) {
    // Init
    ResultsReaderType::Pointer rreader = ResultsReaderType::New();
    ResultsImageType::Pointer rimage = NULL;

    if (fileName.compare("")!=0) {
        rreader->SetFileName( fileName );
        // Reading the image.
        try {
            rreader->Update();
            rimage = rreader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(rimage);
}

MaskImage BrainIO::ReadMaskImage(FileNameType fileName) {
    // Init
    MaskReaderType::Pointer breader = MaskReaderType::New();
    MaskImageType::Pointer bimage = NULL;

    if (fileName.compare("")!=0) {
        breader->SetFileName( fileName );
        // Reading the image.
        try {
            breader->Update();
            bimage = breader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(bimage);
}

DeformationField BrainIO::ReadDeformationField(FileNameType fileName) {
    // Init
	DeformationReaderType::Pointer dfreader = DeformationReaderType::New();
	DeformationFieldType::Pointer dfimage = NULL;

    if (fileName.compare("")!=0) {
        dfreader->SetFileName( fileName );
        // Reading the image.
        try {
            dfreader->Update();
            dfimage = dfreader->GetOutput();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(dfimage);
}

AffineTransform BrainIO::ReadAffineTransform(FileNameType fileName) {
    // Init
	TransformReaderType::Pointer affreader = TransformReaderType::New();
	AffineTransform affine = NULL;

    if (fileName.compare("")!=0) {
        affreader->SetFileName( fileName );
        // Reading the image.
        try {
            affreader->Update();
			affine = static_cast< AffineTransformType* > ((*affreader->GetTransformList()->begin()).GetPointer());
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(affine);
}
CompositeTransform BrainIO::ReadCompositeTransform(FileNameType fileName) {
    // Init
    TransformReaderType::Pointer compositereader = TransformReaderType::New();
    CompositeTransform compositeTransform = NULL;

    if (fileName.compare("")!=0) {
        compositereader->SetFileName( fileName );
        // Reading the image.
        try {
            compositereader->Update();
            compositeTransform = static_cast< CompositeTransformType* > ((*compositereader->GetTransformList()->begin()).GetPointer());
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file reader " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
    return(compositeTransform);
}

void BrainIO::WriteIntensityImage(FileNameType fileName,IntensityImage intimage) {
    // Init
    IntensityWriterType::Pointer intwriter = IntensityWriterType::New();

    if (fileName.compare("")!=0) {
        intwriter->SetFileName( fileName );
        intwriter->SetInput( intimage );
        // Writing the image.
        try {
            intwriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteConnectedImage(FileNameType fileName,ConnectedImage resimage) {
    // Init
    ConnectedWriterType::Pointer conwriter = ConnectedWriterType::New();

    if (fileName.compare("")!=0) {
        conwriter->SetFileName( fileName );
        conwriter->SetInput( resimage );
        // Writing the image.
        try {
            conwriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteResultsImage(FileNameType fileName,ResultsImage resimage) {
    // Init
    ResultsWriterType::Pointer reswriter = ResultsWriterType::New();

    if (fileName.compare("")!=0) {
        reswriter->SetFileName( fileName );
        reswriter->SetInput( resimage );
        // Writing the image.
        try {
            reswriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteProbabilityImage(FileNameType fileName,ProbabilityImage fimage) {
    // Init
    ProbabilityWriterType::Pointer fwriter = ProbabilityWriterType::New();

    if (fileName.compare("")!=0) {
        fwriter->SetFileName( fileName );
        fwriter->SetInput( fimage );
        // Writing the image.
        try {
            fwriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteMaskImage(FileNameType fileName,MaskImage bimage) {
    // Init
    ResultsWriterType::Pointer mwriter = ResultsWriterType::New();
    CastFilterType::Pointer caster = CastFilterType::New();

    caster->SetInput(bimage);

    if (fileName.compare("")!=0) {
        mwriter->SetFileName( fileName );
        mwriter->SetInput( caster->GetOutput() );
        // Writing the image.
        try {
            mwriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}

void BrainIO::WriteDeformationField(FileNameType fileName,DeformationField dfimage) {
    // Init
	DeformationWriterType::Pointer dfwriter = DeformationWriterType::New();

    if (fileName.compare("")!=0) {
        dfwriter->SetFileName( fileName );
        dfwriter->SetInput( dfimage );
        // Writing the image.
        try {
            dfwriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}
void BrainIO::WriteRigidTransform(FileNameType fileName, RigidTransform rigidT) {
    // Init
    TransformWriterType::Pointer rigidWriter = TransformWriterType::New();

    if (fileName.compare("")!=0) {
        rigidWriter->SetFileName( fileName );
        rigidWriter->SetInput( rigidT );
        // Writing the image.
        try {
            rigidWriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}
void BrainIO::WriteAffineTransform(FileNameType fileName, AffineTransform affine) {
    // Init
	TransformWriterType::Pointer affWriter = TransformWriterType::New();

    if (fileName.compare("")!=0) {
        affWriter->SetFileName( fileName );
        affWriter->SetInput( affine );
        // Writing the image.
        try {
            affWriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}
void BrainIO::WriteCompositeTransform(FileNameType fileName, CompositeTransform compT) {
    // Init
    TransformWriterType::Pointer compTWriter = TransformWriterType::New();

    if (fileName.compare("")!=0) {
        compTWriter->SetFileName( fileName );
        compTWriter->SetInput( compT );
        // Writing the image.
        try {
            compTWriter->Update();
        }
        catch (itk::ExceptionObject & e) {
            std::cerr << "exception in file writer " << std::endl;
            std::cerr << e.GetDescription() << std::endl;
            std::cerr << e.GetLocation() << std::endl;
        }
    }
}
#ifdef WIN32
FileNameType BrainIO::SearchForFile(FileNameType folder, FileNameType toFind, bool verbose) {
	/*--- Init ---*/
	FileNameType currentFile, file;
	bool found = false;
	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;

	/*--- Search ---*/
	if (verbose)
		std::cout << "Search: " << folder << std::endl;

	if((hFind = FindFirstFile(folder.c_str(), &FindFileData)) != INVALID_HANDLE_VALUE){
		do {
			currentFile = FindFileData.cFileName;
			if (currentFile.find(toFind) != FileNameType::npos) {
				found = true;
				file = currentFile;
				
				if (verbose)
					std::cout << "\tFound file: " << file << std::endl;
			}
		} while(FindNextFile(hFind, &FindFileData) & !found);
		
		FindClose(hFind);
	}
	else
		std::cerr << "ERROR: Invalid handle value for folder " << folder << std::endl;

	return(file);
}
#else
FileNameType BrainIO::SearchForFile(FileNameType folder, FileNameType toFind, bool verbose) {
	/*--- Init ---*/
	FileNameType currentFile, file;
	bool found = false;

	DIR *dir;
	struct dirent *ent;
    if ((dir = opendir (folder.substr(0,folder.find("*")).c_str())) != NULL) {
	  /* print all the files and directories within directory */
      while (((ent = readdir (dir)) != NULL) && !found) {
          currentFile = ent->d_name;
          if (currentFile.find(toFind) != FileNameType::npos) {
              found = true;
              file = currentFile;

              if (verbose)
                  std::cout << "\tFound file: " << file << std::endl;
          }
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  perror ("");
	}
	return(file);
}
#endif

#ifdef WIN32
FileNameType BrainIO::SearchForFile(FileNameType folder, std::vector<FileNameType> toFindVector, bool verbose) {
	/*--- Init ---*/
	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;
	FileNameType currentFile, file;
	bool found = false;
	std::vector<FileNameType>::iterator it;

	/*--- Search ---*/
	if (verbose)
		std::cout << "Search: " << folder << std::endl;

	if((hFind = FindFirstFile(folder.c_str(), &FindFileData)) != INVALID_HANDLE_VALUE){
		do {
			currentFile = FindFileData.cFileName;
			it = toFindVector.begin();
			do {
				if (currentFile.find(*it) != FileNameType::npos) {
					found = true;
					file = currentFile;
				
					if (verbose)
						std::cout << "\tFound file: " << file << std::endl;
				}
				++it;
			} while(!found & it != toFindVector.end());
		} while(FindNextFile(hFind, &FindFileData) & !found);
		
		FindClose(hFind);
	}
	else
		std::cerr << "ERROR: Invalid handle value for folder " << folder << std::endl;

	return(file);
}
#else
FileNameType BrainIO::SearchForFile(FileNameType folder, std::vector<FileNameType> toFindVector, bool verbose) {
	/*--- Init ---*/
	FileNameType currentFile, file;
	bool found = false;
	std::vector<FileNameType>::iterator it;
	DIR *dir;
	struct dirent *ent;
    if ((dir = opendir (folder.substr(0,folder.find("*")).c_str())) != NULL) {
	  /* print all the files and directories within directory */
      while (((ent = readdir (dir)) != NULL) && !found) {
          currentFile = ent->d_name;
          it = toFindVector.begin();
          do {
              if (currentFile.find(*it) != FileNameType::npos) {
                  found = true;
                  file = currentFile;

                  if (verbose)
                      std::cout << "\tFound file: " << file << std::endl;
              }
              ++it;
          } while(!found & it != toFindVector.end());
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  perror ("");
	}
	return(file);

}
#endif


