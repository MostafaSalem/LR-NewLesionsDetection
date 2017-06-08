PROJECT(BrainToolsCARM)

cmake_minimum_required(VERSION 2.8.9)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()



FIND_PACKAGE(ITK REQUIRED)
include(${ITK_USE_FILE})




SET (SEGMENTATION_SOURCES
    brainsegmentation.cpp
    gaussianestimator1d.cpp
    gaussianestimator2d.cpp
    gaussianestimator3d.cpp
    gaussianestimatornd.cpp
)
SET (SEGMENTATION_HEADERS
    brainsegmentation.h
    gaussianestimator.h
    gaussianestimator1d.h
    gaussianestimator2d.h
    gaussianestimator3d.h
    gaussianestimatornd.h
)
SET (SEGMENTATION_FILES
	${SEGMENTATION_SOURCES}
	${SEGMENTATION_HEADERS}
)

SET (REGISTRATION_SOURCES
    brainregistration.cpp
)
SET (REGISTRATION_HEADERS
    brainregistration.h
)
SET (REGISTRATION_FILES
	${REGISTRATION_SOURCES}
	${REGISTRATION_HEADERS}
)

SET (PREPROCESSING_SOURCES
    brainpreprocessing.cpp
)
SET (PREPROCESSING_HEADERS
    brainpreprocessing.h
)
SET (PREPROCESSING_FILES
	${PREPROCESSING_SOURCES}
	${PREPROCESSING_HEADERS}
)

SET (OTHER_SOURCES
    brainio.cpp
)
SET (OTHER_HEADERS
    imagedefinitions.h
    brainio.h
)
SET (OTHER_FILES
	${OTHER_SOURCES}
	${OTHER_HEADERS}
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)



SET(PRE_FILES ${PREPROCESSING_FILES} ${OTHER_FILES})
SET(COREG_FILES ${REGISTRATION_FILES} ${OTHER_FILES})
SET(ATLAS_FILES ${REGISTRATION_FILES} ${OTHER_FILES})
SET(TISSUE_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(SUB_FILES ${REGISTRATION_FILES} ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(DEFO_FILES ${PREPROCESSING_FILES} ${REGISTRATION_FILES} ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(POST_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})
SET(ANALYSE_FILES ${SEGMENTATION_FILES} ${OTHER_FILES})

ADD_EXECUTABLE( PreTool  premain.cpp ${PRE_FILES})
ADD_EXECUTABLE( CoregTool  coregmain.cpp ${COREG_FILES})
ADD_EXECUTABLE( AtlasTool  atlasmain.cpp ${ATLAS_FILES})
ADD_EXECUTABLE( TissueTool  tissuemain.cpp ${TISSUE_FILES})
ADD_EXECUTABLE( RelabelTool  labelmain.cpp ${TISSUE_FILES})
ADD_EXECUTABLE( SubtractionTool  submain.cpp ${SUB_FILES})
ADD_EXECUTABLE( VoxelSelectionTool  voxelselectoinmain.cpp ${SUB_FILES})
ADD_EXECUTABLE( DeformableTool defomain.cpp ${DEFO_FILES})
ADD_EXECUTABLE( PostTool postmain.cpp ${POST_FILES})
ADD_EXECUTABLE( AnalyseTool analysemain.cpp ${ANALYSE_FILES})
ADD_EXECUTABLE( HistogramTool histogrammain.cpp ${PRE_FILES})
ADD_EXECUTABLE( Test testmain.cpp )


TARGET_LINK_LIBRARIES( PreTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( CoregTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( AtlasTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( TissueTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( RelabelTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( SubtractionTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( VoxelSelectionTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( DeformableTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( PostTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( AnalyseTool ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES( HistogramTool ${ITK_LIBRARIES})



