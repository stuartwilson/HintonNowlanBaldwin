#
# This is an example CMakeLists.txt file to compile a program against
# libmorphologica (which must already have been built and installed).
#

# cmake version 3.1 provides the set(CMAKE_CXX_STANDARD 17) feature
cmake_minimum_required(VERSION 3.1)

# Give your project a name
project(HintonNolan)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or __WIN__ are required
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(EXTRA_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(EXTRA_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
  set(EXTRA_HOST_DEFINITION "-D__OSX__")
elseif(CMAKE_SYSTEM MATCHES Win.*)
  set(EXTRA_HOST_DEFINITION "-D__WIN__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

# morphologica uses c++-17 language features
set(CMAKE_CXX_STANDARD 17)

# Common cflags
set(CMAKE_CXX_FLAGS "${EXTRA_HOST_DEFINITION} -Wall -Wfatal-errors -Wno-sign-compare -g")

# OS dependent cflags
if (APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
  # This assumes a gcc compiler (or a gcc mimic like Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -Wno-unknown-pragmas -march=native -O3")
endif()

# Tell clang to be quiet about brace initialisers:
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# Correct way to determine if OpenMP is available
find_package(OpenMP)
if(OpenMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Additional GL compiler flags.
#
# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")
# Add GL defines to the CXX_FLAGS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

# Tell the program where the morph fonts are, to compile them into the binary
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMORPH_FONTS_DIR=\"\\\"${PROJECT_SOURCE_DIR}/morphologica/fonts\\\"\"")

# Find the libraries which will be needed
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(glfw3 REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Freetype REQUIRED)

# Define collections of includes that have to be made for morphologica
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
set(MORPH_INC_CORE ${JSON_INC_PATH} ${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
set(MORPH_INC_GL ${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS})
include_directories(${MORPH_INC_CORE} ${MORPH_INC_GL})

# Assuming that you installed morphologica in-tree (i.e. 'next to' XXX.cpp).
set(MORPH_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/morphologica" CACHE PATH "The path to morphologica")
include_directories(BEFORE ${MORPH_INCLUDE_PATH}/include)
include_directories(BEFORE ${MORPH_INCLUDE_PATH})

#
# All the library finding is now done. It's time to compile the program.
#



# Define which code makes up the recurrent net program
#add_executable(hinton hinton.cpp)
#add_executable(hintonBatch hintonBatch.cpp)
#add_executable(hintonBatchCritical hintonBatchCritical.cpp)
#add_executable(hintonBatchForced hintonBatchForced.cpp)
#add_executable(hintonBatchVolatile hintonBatchVolatile.cpp)
#add_executable(hintonBatchOpposing hintonBatchOpposing.cpp)


# Define how model will be compiled
#target_compile_definitions(hinton PUBLIC FLT=float COMPILE_PLOTTING)
#target_compile_definitions(hintonBatch PUBLIC FLT=float COMPILE_PLOTTING)
#target_compile_definitions(hintonBatchForced PUBLIC FLT=float COMPILE_PLOTTING)
#target_compile_definitions(hintonBatchCritical PUBLIC FLT=float COMPILE_PLOTTING)
#target_compile_definitions(hintonBatchVolatile PUBLIC FLT=float COMPILE_PLOTTING)
#target_compile_definitions(hintonBatchOpposing PUBLIC FLT=float COMPILE_PLOTTING)


# Morphologica code requires a number of libraries, collected into 'CORE' and 'GL'.
set(MORPH_LIBS_CORE ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${HDF5_C_LIBRARIES} ${LAPACK_LIBRARIES} jsoncpp_lib)
set(MORPH_LIBS_GL ${OpenCV_LIBS} OpenGL::GL Freetype::Freetype glfw)

#foreach(MODEL hinton hintonBatch hintonBatchForced hintonBatchCritical hintonBatchVolatile hintonBatchOpposing hintonBatchVolatileConfigurable)
foreach(MODEL hinton hintonSeparate)
	set(MODELNAME ${MODEL})
	add_executable(${MODELNAME} ${MODELNAME}.cpp)
	target_compile_definitions(${MODELNAME} PUBLIC FLT=float COMPILE_PLOTTING)
	target_link_libraries(${MODELNAME} ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL})
	#target_link_libraries(${MODELNAME} OpenMP::OpenMP_CXX)
endforeach()

#if(OpenMP_FOUND)
#  target_link_libraries(hinton OpenMP::OpenMP_CXX)
#  target_link_libraries(hintonBatch OpenMP::OpenMP_CXX)
#  target_link_libraries(hintonBatchForced OpenMP::OpenMP_CXX)
#  target_link_libraries(hintonBatchCritical OpenMP::OpenMP_CXX)
  #target_link_libraries(hintonBatchVolatile OpenMP::OpenMP_CXX)
#  target_link_libraries(hintonBatchOpposing OpenMP::OpenMP_CXX)
#endif()

# For debugging of variables:
option(DEBUG_VARIABLES OFF)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
