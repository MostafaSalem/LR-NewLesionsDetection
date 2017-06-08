#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#ifdef WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif


using namespace std;

typedef std::string FileNameType;

#ifdef WIN32
FileNameType SearchForFile(FileNameType folder, FileNameType toFind, bool verbose = false) {
	/*--- Init ---*/
	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;
	FileNameType currentFile, file;
	bool found = false;

	/*--- Search ---*/
	if (verbose)
		printf("Search: %s\n", folder.c_str());
	if((hFind = FindFirstFile(folder.c_str(), &FindFileData)) != INVALID_HANDLE_VALUE){
		do {
			currentFile = FindFileData.cFileName;
			if (currentFile.find(toFind) != string::npos) {
				found = true;
				file = currentFile;
				
				if (verbose)
					printf("\t%s\n", currentFile.c_str());
			}
		} while(FindNextFile(hFind, &FindFileData) & !found);
		
		FindClose(hFind);
	}
	else
		printf("ERROR: Invalid handle value for folder %s",folder.c_str());

	return(file);
};
#else
FileNameType SearchForFile(FileNameType folder, FileNameType toFind, bool verbose = false) {
    /*--- Init ---*/
    FileNameType currentFile, file;
    bool found = false;

    return(file);
};
#endif

#ifdef WIN32
FileNameType SearchForFile(FileNameType folder, std::vector<FileNameType> toFindVector, bool verbose) {
	/*--- Init ---*/
	HANDLE hFind;
	WIN32_FIND_DATA FindFileData;
	FileNameType currentFile, file;
	bool found = false;
	std::vector<FileNameType>::iterator it;

	/*--- Search ---*/
	if (verbose)
		printf("Search: %s\n", folder.c_str());

	if((hFind = FindFirstFile(folder.c_str(), &FindFileData)) != INVALID_HANDLE_VALUE){
		do {
			currentFile = FindFileData.cFileName;
			it = toFindVector.begin();
			do {
				if (currentFile.find(*it) != FileNameType::npos) {
					found = true;
					file = currentFile;
				
					if (verbose)
						printf("\t%s\n", currentFile.c_str());
				}
				++it;
			} while(!found & it != toFindVector.end());
		} while(FindNextFile(hFind, &FindFileData) & !found);
		
		FindClose(hFind);
	}
	else
		printf("ERROR: Invalid handle value for folder %s",folder.c_str());

	return(file);
}
#else
FileNameType SearchForFile(FileNameType folder, std::vector<FileNameType> toFindVector, bool verbose) {
    /*--- Init ---*/
    FileNameType currentFile, file;
    bool found = false;
    std::vector<FileNameType>::iterator it;
    return(file);
}
#endif

int main(int argc, char **argv)
{

	/*--- Init ---*/
	printf("Init: %d parameters\n", argc);
	for (int i = 0; i < argc; i++)
		printf("\tParameter %d/%d: %s\n", i, argc, argv[i]);
	FileNameType baselineFolderName = argv[1];
	FileNameType followingFolderName = argv[2];
	FileNameType baselineSearchName = baselineFolderName + "*.nii";
	FileNameType followingSearchName = followingFolderName + "*.nii";

	FileNameType toFind = "dp";
	std::vector<FileNameType> toFindVector;
	toFindVector.push_back("DP");
	toFindVector.push_back("PD");
	toFindVector.push_back("pd");
	toFindVector.push_back("dp");

	printf("Simple search test\n");
	FileNameType testSimple = SearchForFile(baselineSearchName, toFind, true);
	if (!testSimple.empty())
        printf("Simple search results: %s", testSimple.c_str());
	printf("Vector search test\n");
	FileNameType testVector = SearchForFile(baselineSearchName, toFindVector, true);
	if (!testVector.empty())
        printf("Vector search results: %s", testVector.c_str());

	FileNameType testFolder = baselineFolderName + "processed\\";

    #ifdef WIN32
	CreateDirectory(testFolder.c_str(),NULL);
    #else
    mkdir(testFolder.c_str(), 0777);
    #endif


    return EXIT_SUCCESS;
}

