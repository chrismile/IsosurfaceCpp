cmake_minimum_required(VERSION 3.10...4.0)
project(IsosurfaceCpp)
set(CMAKE_CXX_STANDARD 17)

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c src/*.h src/*.hpp)

if (MSVC)
    set(DEFAULT_USE_TBB ON)
else()
    set(DEFAULT_USE_TBB OFF)
endif()

option(USE_TBB "Build with TBB threading support instead of using OpenMP." ${DEFAULT_USE_TBB})

find_package(glm REQUIRED)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    add_library(IsosurfaceCpp SHARED ${SOURCES})
else()
    add_library(IsosurfaceCpp ${SOURCES})
endif()

if (TARGET glm::glm)
    target_link_libraries(IsosurfaceCpp PRIVATE glm::glm)

    # Homebrew sets GLM_INCLUDE_DIRS to the wrong directory on Linux for some reason.
    if (GLM_INCLUDE_DIRS MATCHES "lib$")
        target_include_directories(IsosurfaceCpp PRIVATE "${GLM_INCLUDE_DIRS}/../include")
    endif()
elseif (TARGET glm)
    target_link_libraries(IsosurfaceCpp PRIVATE glm)
else()
    target_include_directories(IsosurfaceCpp PRIVATE ${GLM_INCLUDE_DIR})
endif()

if (${USE_TBB})
    if (NOT TARGET TBB::tbb)
        # Add guard, as parent project may use OpenVDB, which uses a custom FindTBB script.
        find_package(TBB REQUIRED)
    endif()
    target_link_libraries(IsosurfaceCpp PRIVATE TBB::tbb)
    target_compile_definitions(IsosurfaceCpp PRIVATE USE_TBB)
else()
    # According to https://devblogs.microsoft.com/cppblog/improved-openmp-support-for-cpp-in-visual-studio/,
    # support for LLVM OpenMP was added with Visual Studio 2019 version 16.9. According to
    # https://docs.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-170,
    # both Visual Studio 2019 version 16.8 and 16.9 have _MSC_VER / MSVC_VERSION 1928, so we need to check
    # for a minimum version 1929.
    if (APPLE)
        # There are some problems with Homebrew OpenMP as of 2022-10-20.
        find_package(OpenMP QUIET)
        if (OpenMP_FOUND)
            MESSAGE(STATUS "Found OpenMP version ${OpenMP_CXX_VERSION}. Enabling OpenMP threading support.")
        else()
            MESSAGE(STATUS "OpenMP not found. Disabling OpenMP threading support.")
        endif()
    else()
        find_package(OpenMP REQUIRED)
    endif()
    if(OpenMP_FOUND)
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
        elseif(NOT TARGET OpenMP::OpenMP_CXX)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        else()
            target_link_libraries(IsosurfaceCpp PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()
endif()


if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/IsosurfaceCppConfig.cmake.in" CONTENTS)
    file(WRITE "${CMAKE_BINARY_DIR}/IsosurfaceCppConfig.cmake.tmp" "${CONTENTS}")

    # For creating IsosurfaceCppConfig.cmake.in
    include(GenerateExportHeader)
    generate_export_header(IsosurfaceCpp)
    install(
            TARGETS IsosurfaceCpp EXPORT IsosurfaceCppTargets
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin
            INCLUDES DESTINATION include/IsosurfaceCpp
    )

    install(
            DIRECTORY ${CMAKE_SOURCE_DIR}/src/
            DESTINATION include/IsosurfaceCpp
            FILES_MATCHING PATTERN "*.h*"
    )

    include(CMakePackageConfigHelpers)

    export(
            EXPORT IsosurfaceCppTargets
            FILE "${CMAKE_CURRENT_BINARY_DIR}/IsosurfaceCpp/IsosurfaceCppTargets.cmake"
    )
    configure_file(
            "${CMAKE_BINARY_DIR}/IsosurfaceCppConfig.cmake.tmp"
            "${CMAKE_CURRENT_BINARY_DIR}/IsosurfaceCpp/IsosurfaceCppConfig.cmake"
            COPYONLY
    )

    set(ConfigPackageLocation lib/cmake/IsosurfaceCpp)
    install(
            EXPORT IsosurfaceCppTargets
            FILE IsosurfaceCppTargets.cmake
            DESTINATION ${ConfigPackageLocation}
    )
    install(
            FILES "${CMAKE_CURRENT_BINARY_DIR}/IsosurfaceCpp/IsosurfaceCppConfig.cmake"
            DESTINATION ${ConfigPackageLocation}
    )
    if (WIN32 AND NOT ${BUILD_STATIC_LIBRARY} AND ((CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
            AND NOT (MSYS OR MINGW OR (${CMAKE_GENERATOR} STREQUAL "MinGW Makefiles") OR (${CMAKE_GENERATOR} STREQUAL "MSYS Makefiles")))
        install(
                FILES $<TARGET_PDB_FILE:IsosurfaceCpp>
                DESTINATION bin OPTIONAL
        )
    endif()
endif()
