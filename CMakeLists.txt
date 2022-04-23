cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(VirtualEndoscopy)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/VirtualEndoscopy")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Yi Gao (Shenzhen University), Haofan Huang (Shenzhen University)")
set(EXTENSION_DESCRIPTION "This extension performs the virtual endoscopy, including: centerline extraction, lumen segmentation, and virtual fly-through.")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/VirtualEndoscopy.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/VirtualEndoscopy/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})



#--------------------------------------------------------------------------------
# OMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


#--------------------------------------------------------------------------------
# find GTH818N path
set(GTH818N_PATH Libs/gth818n/)

include_directories(${GTH818N_PATH})

#--------------------------------------------------------------------------------
# my libs
add_library(GTH818N_LIBRARIES
  ${GTH818N_PATH}/fibheap.cpp
  )


#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(VirtualEndoscopy)
add_subdirectory(ComputeVesselness)
add_subdirectory(ComputeAxisFromVesselness)
add_subdirectory(SegmentLumenFromAxis)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
