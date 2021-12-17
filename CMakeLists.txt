cmake_minimum_required (VERSION 3.7)
cmake_policy(VERSION 3.7...3.20)
project(IsosurfaceCpp)

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c src/*.h src/*.hpp)

find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_library(IsosurfaceCpp SHARED ${SOURCES})
else()
    add_library(IsosurfaceCpp ${SOURCES})
endif()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11")
    target_link_libraries(IsosurfaceCpp PRIVATE glm::glm)
else()
    target_include_directories(IsosurfaceCpp PRIVATE ${GLM_INCLUDE_DIR})
endif()

# According to https://devblogs.microsoft.com/cppblog/improved-openmp-support-for-cpp-in-visual-studio/,
# support for LLVM OpenMP was added with Visual Studio 2019 version 16.9. According to
# https://docs.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-170,
# both Visual Studio 2019 version 16.8 and 16.9 have _MSC_VER / MSVC_VERSION 1928, so we need to check
# for a minimum version 1929.
find_package(OpenMP REQUIRED)
if(MSVC)
    option(USE_MSVC_OPENMP_LLVM "Use LLVM OpenMP with MSVC (available since Visual Studio 2019 version 16.9)." ON)
    option(USE_MSVC_OPENMP_EXPERIMENTAL "Use experimental OpenMP with MSVC (available since Visual Studio 2019 version 16.9)." ON)

    if(${MSVC_VERSION} GREATER_EQUAL 1929 AND (USE_MSVC_OPENMP_LLVM OR USE_MSVC_OPENMP_EXPERIMENTAL))
        if(USE_MSVC_OPENMP_LLVM)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp:llvm")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp:llvm")
        endif()
        if(USE_MSVC_OPENMP_EXPERIMENTAL)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp:experimental")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp:experimental")
        endif()
        if(USE_MSVC_OPENMP_LLVM AND USE_MSVC_OPENMP_EXPERIMENTAL)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /w44005")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w44005")
        endif()
    endif()
else()
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()
