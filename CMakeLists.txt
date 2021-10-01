cmake_minimum_required (VERSION 3.5)
cmake_policy(VERSION 3.5...3.20)
project(IsosurfaceCpp)

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c src/*.h src/*.hpp)

find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_library(IsosurfaceCpp SHARED ${SOURCES})
else()
    add_library(IsosurfaceCpp ${SOURCES})
endif()

target_link_libraries(IsosurfaceCpp PRIVATE glm::glm)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
