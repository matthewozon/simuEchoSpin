cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(mriSimulation)#Spline)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

#Inclusion de Boost
include_directories(. lib/include)
link_directories(. lib/lib)

#Configuration de l'exécutable
file(
	GLOB
	source_files
	*.h *.cpp
)
add_executable(
        mriSimulation #Spline
	${source_files}
)

#Configuration de l'édition de liens
target_link_libraries(
	mriSimulation
        gif # for gif images
        tiff #for tiff images
        png #read png images
        X11 #needed by png (different on windows)
        pthread #needed by png (different on windows)
	)

