cmake_minimum_required(VERSION 3.10)

project(MPCNS LANGUAGES CXX)

option(MPCNS_ENABLE_MPI "Build with MPI support" OFF)
option(MPCNS_ENABLE_METIS "Build with METIS graph partitioning support" OFF)
option(MPCNS_ENABLE_SAFE_COUT "Use a mutex for parallel cout" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(MAIN
    main.cpp
    src/MPCNS_Pre_App.cpp
    src/MPCNS_Pre_Data.cpp
    src/MPCNS_Pre_Datatrans.cpp
    src/MPCNS_Pre_Group.cpp
    src/MPCNS_Pre_IOinfo.cpp
    src/MPCNS_Pre_Parameter.cpp
    src/MPCNS_Pre_Split.cpp
    src/Pre_ACANSOutput.cpp
    src/Pre_MPCNSOutput.cpp
)

target_include_directories(MAIN PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

if(MPCNS_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_link_libraries(MAIN PRIVATE MPI::MPI_CXX)
    target_compile_definitions(MAIN PRIVATE MPCNS_IF_MPI=1)
else()
    target_compile_definitions(MAIN PRIVATE MPCNS_IF_MPI=0)
endif()

if(MPCNS_ENABLE_METIS)
    find_path(METIS_INCLUDE_DIR metis.h
        PATHS ${PROJECT_SOURCE_DIR}/metis
        NO_DEFAULT_PATH
    )
    find_library(METIS_LIBRARY metis
        PATHS ${PROJECT_SOURCE_DIR}/metis
        NO_DEFAULT_PATH
    )

    if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
        find_path(METIS_INCLUDE_DIR metis.h)
        find_library(METIS_LIBRARY metis)
    endif()

    if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
        message(FATAL_ERROR "MPCNS_ENABLE_METIS=ON, but METIS was not found. Set METIS_INCLUDE_DIR and METIS_LIBRARY, or turn it OFF.")
    endif()

    target_include_directories(MAIN PRIVATE ${METIS_INCLUDE_DIR})
    target_link_libraries(MAIN PRIVATE ${METIS_LIBRARY})
    target_compile_definitions(MAIN PRIVATE MPCNS_IF_METIS=1)
else()
    target_compile_definitions(MAIN PRIVATE MPCNS_IF_METIS=0)
endif()

if(MPCNS_ENABLE_SAFE_COUT)
    target_compile_definitions(MAIN PRIVATE MPCNS_Para_COUT=1)
else()
    target_compile_definitions(MAIN PRIVATE MPCNS_Para_COUT=0)
endif()
