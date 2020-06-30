PROJECT (VisIVOServer_New072012)
# Well technicall you also need:
# http://cmake.org/Bug/bug.php?op=show&bugid=3582
cmake_minimum_required(VERSION 2.8.5)
#-----------------------------------------------------------------------------
OPTION (VSBIGENDIAN "You should select on if your sistem is big endian" OFF)
OPTION (VSAPILIB "You must select ON to create Library" OFF)
OPTION (VSMAC "You must select ON to compile on a MAC system" ON)
OPTION (LIGHT "You must select OFF to compile on a full features" ON)
OPTION (GLITE "You must select ON to compile on gLite grid nodes" OFF)
OPTION (AHF "You must select ON to compile with AHF features" OFF)
OPTION (SPLVISIVO "You must select ON to compile with Splotch" ON)
OPTION (SPLCUDA "You must select ON to compile with CUDA support" OFF)
OPTION (VSMPI "You must select ON to compile with MPI support" OFF)


SET(CFITSIO_INCLUDE_DIR "/opt/local/" CACHE PATH "")
SET(CFITSIO_LIB_DIR "/opt/local/" CACHE PATH "")
SET(CURL_DIR "/opt/local/" CACHE PATH "")
SET(HDF5_INCLUDE_DIR "/opt/local/" CACHE PATH "")
SET(HDF5_LIB_DIR "/opt/local/" CACHE PATH "")
SET(VTK_DIR "/opt/vtk6.20/lib/cmake/vtk-6.2/")

#SET(VTK_DIR "/opt/vtk5.10/lib/vtk-5.10/")
#SET(CFITSIO_DIR "SET_CFITSIO_DIR" CACHE PATH "")
#SET(CURL_DIR "SET_CURL_DIR" CACHE PATH "")
#SET(HDF5_DIR "SET_HDF5_DIR" CACHE PATH "")

if(VSMPI)
	ADD_DEFINITIONS(-DVSMPI)
    SET(CMAKE_CXX_COMPILER "Users/ubesoft/mpich2-install/bin/mpicxx" CACHE PATH "")
endif(VSMPI)

if(SPLCUDA)
	ADD_DEFINITIONS(-DCUDA)
	ADD_DEFINITIONS(-DNO_WIN_THREAD)
	find_package(CUDA)
endif(SPLCUDA)

#ADD_DEFINITIONS(-DHDF5)


FIND_PACKAGE(VTK REQUIRED)


#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
#CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#Searching CUDA
FIND_PACKAGE(CUDA)

#Adding src directory to the project
#ADD_SUBDIRECTORY(src)

INCLUDE(${VTK_USE_FILE})

#-----------------------------------------------------------------------------
SET (EXECUTABLE_OUTPUT_PATH ./bin CACHE PATH "Single output directory for building all executables.")
SET (LIBRARY_OUTPUT_PATH ${VisIVOApiTest_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)



#---CONFIGURE_FILE (${PROJECT_SOURCE_DIR}/VisIVOApiTestConfigure.h.in ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Configure.h @ONLYIMMEDIATE)


SET(SRCS_VisIVOFilterMain
    mainFilter.cpp
)

SET(SRCS_VisIVOImporterMain
    mainImporter.cpp
)
SET(SRCS_VisIVOViewerMain
    mainViewer.cpp
)

SET(SRCS_VisIVOUtilsMain
    mainUtils.cpp
)



SET(SRCS_VisIVOApi

    visivoset.cpp
    visivoview.cpp
    visivoimp.cpp
    visivofil.cpp
    )

if(GLITE)
 ADD_DEFINITIONS(-DGLITE) 
endif(GLITE)



if(AHF)
SET(SRCS_ahf
VisIVOFilters/Filters/vsvbt2ahf.cpp
)
 ADD_DEFINITIONS(-DAHF) 
endif(AHF)

if(LIGHT)
 SET(DIRS  VisIVOImporter/Importers  Utils VisIVOViewer/pipe  VisIVOFilters/Filters VisIVOUtils )
 ADD_DEFINITIONS(-DLIGHT) 
else(LIGHT)   

SET(DIRS VisIVOImporter/External VisIVOImporter/Importers  Utils VisIVOViewer/pipe  VisIVOFilters/Filters VisIVOUtils )

endif(LIGHT)

if(SPLVISIVO)
 if(SPLCUDA)
 	SET(SPLDIRS VisIVOViewer/Splotch/kernel VisIVOViewer/Splotch/cxxsupport VisIVOViewer/Splotch/cuda VisIVOViewer/Splotch/reader VisIVOViewer/Splotch/splotch VisIVOViewer/Splotch/booster /VisIVOViewer/Splotch/c_utils)
 else(SPLCUDA)
	 SET(SPLDIRS VisIVOViewer/Splotch/kernel VisIVOViewer/Splotch/cxxsupport VisIVOViewer/Splotch/reader VisIVOViewer/Splotch/splotch VisIVOViewer/Splotch/booster /VisIVOViewer/Splotch/c_utils)
endif(SPLCUDA)
FOREACH(DIR_VAL ${SPLDIRS})
        INCLUDE(./${DIR_VAL}/CMakeLists.txt)
        FOREACH(VAL ${PRJ_SRC})
                SET(SOURCESSpl ${SOURCESSpl} ./${DIR_VAL}/${VAL})
                INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${DIR_VAL})
        ENDFOREACH(VAL)
ENDFOREACH(DIR_VAL)
 ADD_DEFINITIONS(-DSPLVISIVO) 

endif(SPLVISIVO)

FOREACH(DIR_VAL ${DIRS})
        INCLUDE(./${DIR_VAL}/CMakeLists.txt)
        FOREACH(VAL ${PRJ_SRC})
                SET(SOURCES ${SOURCES} ./${DIR_VAL}/${VAL})
                INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${DIR_VAL})
                #FILE(APPEND tmp.txt ./${DIR_VAL}/${VAL})
        ENDFOREACH(VAL)
ENDFOREACH(DIR_VAL)

if(AHF)
SET(AHFDIRS VisIVOFilters/Filters/AHFstep  VisIVOFilters/Filters/AHFstep  VisIVOFilters/Filters/AHFstep/libahf  VisIVOFilters/Filters/AHFstep/libgravity  VisIVOFilters/Filters/AHFstep/libgrids  VisIVOFilters/Filters/AHFstep/libio  VisIVOFilters/Filters/AHFstep/libio_serial  VisIVOFilters/Filters/AHFstep/libmhd  VisIVOFilters/Filters/AHFstep/libparticles  VisIVOFilters/Filters/AHFstep/libsfc  VisIVOFilters/Filters/AHFstep/libutility)
ADD_SUBDIRECTORY( VisIVOFilters/Filters/AHFstep)
FOREACH(DIR_VAL ${AHFDIRS})
        INCLUDE(./${DIR_VAL}/CMakeLists.txt)
        FOREACH(VAL ${PRJ_SRC})
                SET(SOURCESAhf ${SOURCESAhf} ./${DIR_VAL}/${VAL})
                INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/${DIR_VAL})
        ENDFOREACH(VAL)
ENDFOREACH(DIR_VAL)
endif(AHF)



SET(SRCS_VisIVOApiTest
    ${SRCS_VisIVOApi}
${SOURCES} ${SOURCESAhf} ${SRCS_ahf}  ${SOURCESSpl}

)

SET(SRCS_VisIVOFilters
    ${SRCS_VisIVOFilterMain} 
${SOURCES} ${SOURCESAhf} ${SRCS_ahf} ${SOURCESSpl})

SET(SRCS_VisIVOImporter
    ${SRCS_VisIVOImporterMain} 
${SOURCES} ${SOURCESSpl})

SET(SRCS_VisIVOViewer
    ${SRCS_VisIVOViewerMain} 
${SOURCES}  ${SOURCESSpl} )

SET(SRCS_VisIVOUtils
    ${SRCS_VisIVOUtilsMain} 
${SOURCES} ${SOURCESSpl})

if(VSMAC)
ADD_DEFINITIONS(-DMAC)
endif(VSMAC)

if(AHF)
 SET_SOURCE_FILES_PROPERTIES(${SRCS_ahf} PROPERTIES COMPILE_FLAGS "-DHAVE_STDINT_H -DHAVE_STDBOOL_H -DNONC99=0 -DNEWSTARTRUN -DSTANDARD -DAHFstep -DAHF")
 SET_SOURCE_FILES_PROPERTIES(${SOURCESAhf} PROPERTIES COMPILE_FLAGS "-std=c99  -O2  -DHAVE_STDINT_H -DHAVE_STDBOOL_H -DNONC99=0 -DNEWSTARTRUN -DWITH_AMIGA_input -DSTANDARD -DAHFstep -DAHF")
 INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VisIVOFilters/Filters/AHFstep/replace)
endif(AHF)

SET_SOURCE_FILES_PROPERTIES(${SOURCESSpl} PROPERTIES COMPILE_FLAGS "-c -Wno-long-long -Wfatal-errors -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast -O2 -g")
#SET_SOURCE_FILES_PROPERTIES(${SOURCESSpl} PROPERTIES COMPILE_FLAGS "-c -std=c++98 -Wno-long-long -Wfatal-errors -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast -O2 -g")
#SET_SOURCE_FILES_PROPERTIES(${SOURCESSpl} PROPERTIES COMPILE_FLAGS "-c -std=c++98 -pedantic -Wno-long-long -Wfatal-errors -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast -O2 -g")


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VisIVOViewer/Splotch)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VisIVOViewer/Splotch/cxxsupport)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VisIVOViewer/Splotch/kernel)


if(LIGHT)
#INCLUDE_DIRECTORIES(/opt/local/var/macports/software/cfitsio/3.280_1/opt/local/include/)
#INCLUDE_DIRECTORIES(/opt/local/var/macports/software/curl/7.21.7_0+ssl/opt/local/include)
#INCLUDE_DIRECTORIES(/opt/local/var/macports/software/hdf5-18/1.8.7_1/opt/local/include)

#LINK_DIRECTORIES(/opt/local/var/macports/software/cfitsio/3.280_1/opt/local/lib)
#LINK_DIRECTORIES(/opt/local/var/macports/software/curl/7.21.7_0+ssl/opt/local/lib/)
#LINK_DIRECTORIES(/opt/local/var/macports/software/hdf5-18/1.8.7_1/opt/local/lib)


INCLUDE_DIRECTORIES(${CFITSIO_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CURL_DIR}/include)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

LINK_DIRECTORIES(${CFITSIO_LIB_DIR})
LINK_DIRECTORIES(${CURL_DIR}/lib)
LINK_DIRECTORIES(${HDF5_LIB_DIR})

else(LIGHT)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/VisIVOImporter/CodeSynthesisXSD2.3/include)
INCLUDE_DIRECTORIES(/opt/xerces2_8_0/include)
INCLUDE_DIRECTORIES(/opt/cfitsio3_1_0/include)
#INCLUDE_DIRECTORIES(/opt/curl-7.19.4/include)
INCLUDE_DIRECTORIES(/Developer/SDKs/MacOSX10.6.sdk/usr/include/curl/)

INCLUDE_DIRECTORIES(/opt/hdf5-1.8.5-patch1/include)

LINK_DIRECTORIES(/opt/xerces2_8_0/lib)
LINK_DIRECTORIES(/opt/cfitsio3_1_0/lib)
#LINK_DIRECTORIES(/opt/curl-7.19.4/lib)
LINK_DIRECTORIES(/Developer/SDKs/MacOSX10.6.sdk/usr/lib/)
LINK_DIRECTORIES(/opt/libssh2_1/lib)
LINK_DIRECTORIES(/opt/hdf5-1.8.5-patch1/lib)
LINK_LIBRARIES(xerces-c)
LINK_LIBRARIES(xerces-depdom)

LINK_LIBRARIES(ssh2)
#LINK_LIBRARIES(ldap)
LINK_LIBRARIES(ssl)
LINK_LIBRARIES(crypto)
LINK_LIBRARIES(idn)
#LINK_LIBRARIES(rt)
LINK_LIBRARIES(ssl)
#LINK_LIBRARIES(curl)

endif(LIGHT)

INCLUDE_DIRECTORIES(/usr/include)
LINK_DIRECTORIES(/usr/lib64)

if(GLITE)
 INCLUDE_DIRECTORIES(/opt/lcg/include)
 LINK_DIRECTORIES(/opt/lcg/lib64)
 LINK_LIBRARIES(lcg_util)
endif(GLITE)

LINK_LIBRARIES(cfitsio)
LINK_LIBRARIES(curl)
LINK_LIBRARIES(hdf5)


LINK_LIBRARIES(vtkRendering)
LINK_LIBRARIES(vtkHybrid)
LINK_LIBRARIES(vtkVolumeRendering)


SET(GUI_EXECUTABLE )

### NOTE:
# The installation paths (INSTALL_TARGETS, INSTALL_FILES) are relative to CMAKE_INSTALL_PREFIX that is /usr/local by defautl
#----------- Alternative------
if(VSAPILIB)
SET(HEAD_VisIVOApi
 ${PROJECT_SOURCE_DIR}/visivo.h
 ${PROJECT_SOURCE_DIR}/visivodef.h
 ${PROJECT_SOURCE_DIR}/visivoserver.h
 )
ADD_LIBRARY(VisIVOApiTest   ${SRCS_VisIVOApiTest})
INSTALL_TARGETS(/lib VisIVOApiTest)
INSTALL_FILES(/include .h ${HEAD_VisIVOApi})
TARGET_LINK_LIBRARIES( VisIVOApiTest vtkRendering)
TARGET_LINK_LIBRARIES( VisIVOApiTest vtkHybrid)
TARGET_LINK_LIBRARIES( VisIVOApiTest vtkVolumeRendering)

else(VSAPILIB)
#----------------------------------------------------------


ADD_EXECUTABLE(VisIVOImporter ${GUI_EXECUTABLE} ${SRCS_VisIVOImporter})
#INSTALL_TARGETS(/bin VisIVOImporter)

#ADD_EXECUTABLE(VisIVOFilters ${GUI_EXECUTABLE} ${SRCS_VisIVOFilters})
#INSTALL_TARGETS(/bin VisIVOFilters)

include_directories(/Developer/GPU\ Computing/CUDALibraries/common/inc/ ) 
if(SPLCUDA)
CUDA_ADD_EXECUTABLE(VisIVOViewer ${GUI_EXECUTABLE} ${SRCS_VisIVOViewer})
else(SPLCUDA)
#ADD_EXECUTABLE(VisIVOViewer ${GUI_EXECUTABLE} ${SRCS_VisIVOViewer})
endif(SPLCUDA)
#INSTALL_TARGETS(/bin VisIVOViewer)

#ADD_EXECUTABLE(VisIVOUtils ${GUI_EXECUTABLE} ${SRCS_VisIVOUtils})
#INSTALL_TARGETS(/bin VisIVOUtils)

endif(VSAPILIB)
#-----------------------------------------------------
