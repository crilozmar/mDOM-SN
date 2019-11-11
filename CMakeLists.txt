# Find the manual here:
# http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/InstallationGuide/html/ch03s02.html

# check for cmake version, initialize project -> sets name and paths
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mdom)

# choose option
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

# locate geant4 and set all required variables
if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

# try to use ROOT libraries
set (CMAKE_MODULE_PATH /Applic.ZIV/CentOS7/software/p0kapp/software/geant4/10.4.02.p01/lib64/Geant4-10.2.1/Modules/)
#find_package(ROOT)

# sets search paths for include
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIR})

# create list of all the files to compile
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# manually add argtable
LINK_DIRECTORIES(/Applic.ZIV/CentOS7/software/p0kapp/software/lib)

# details for the finished binary
add_executable(mdom mdom.cc ${sources} ${headers})

target_link_libraries(mdom ${Geant4_LIBRARIES} argtable2 ${ROOT_LIBRARIES})


# which auxiliary files out of the source directory should be copied to the directory of the binary
set(mdom_aux
	aux/vis.mac
	aux/vis.ogl
	refraction.cfg
	mdom.gps
	QuantumEfficiency.cfg
)
foreach(_script ${mdom_aux})
	configure_file(
		${PROJECT_SOURCE_DIR}/${_script}
		${PROJECT_BINARY_DIR}/${_script}
	COPYONLY
	)
endforeach()

#auxiliary files for SN and solar fluxes
set(SN_aux
	Flux_Nu_Heavy_ls220.cfg
	Flux_Nubar_Heavy_ls220.cfg
	Flux_Nu_light_ls220.cfg
	Flux_Nubar_light_ls220.cfg
	solarneutrinospectrum.cfg
	nu_DDT.cfg
	nubar_DDT.cfg
	nu_GCD.cfg
	nubar_GCD.cfg
)
foreach(_script ${SN_aux})
	configure_file(
		${PROJECT_SOURCE_DIR}/SNandSolarStuff/${_script}
		${PROJECT_BINARY_DIR}/${_script}
	COPYONLY
	)
endforeach()

