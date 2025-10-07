cmake_minimum_required(VERSION 3.15)
project(QISS VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------
# 1. Collect your project sources
# -----------------------------
file(GLOB_RECURSE QISS_SOURCES src/*.cpp)

add_executable(QISS ${QISS_SOURCES})

# -----------------------------
# 2. Include directories
# -----------------------------
# Your headers
target_include_directories(QISS PRIVATE include)

# FFTW headers
set(FFTW_ROOT "${PROJECT_SOURCE_DIR}/extern/FFTW")
set(FFTW_INCLUDE_DIR "${FFTW_ROOT}/include")
target_include_directories(QISS PRIVATE ${FFTW_INCLUDE_DIR})

# Matplot++ headers
target_include_directories(QISS PRIVATE extern/Matplotplusplus/source)

# ImGui + ImPlot headers
target_include_directories(QISS PRIVATE extern/ImGui extern/ImPlot)

# -----------------------------
# 3. Add external sources manually
# -----------------------------
# ImGui sources
file(GLOB IMGUI_SOURCES
    extern/ImGui/*.cpp
)

# ImPlot sources
file(GLOB IMPLOT_SOURCES
    extern/ImPlot/*.cpp
)

# Add them to your executable
target_sources(QISS PRIVATE ${IMGUI_SOURCES} ${IMPLOT_SOURCES})

# -----------------------------
# 4. Add Matplot++ subdirectory
# -----------------------------
add_subdirectory(extern/Matplotplusplus)
target_link_libraries(QISS PRIVATE matplot)

# -----------------------------
# 5. FFTW linking
# -----------------------------
set(FFTW_LIB_DIR "${FFTW_ROOT}/lib")
set(FFTW_LIB "${FFTW_LIB_DIR}/fftw3.lib")

if(NOT EXISTS ${FFTW_LIB})
    message(FATAL_ERROR "FFTW import library not found: ${FFTW_LIB}. You must run lib.exe on the .def file.")
endif()

link_directories(${FFTW_LIB_DIR})
target_link_libraries(QISS PRIVATE fftw3)

# Link OpenGL (required by ImGui + OpenGL backend)
find_package(OpenGL REQUIRED)
target_link_libraries(QISS PRIVATE OpenGL::GL)

# Copy runtime DLLs
add_custom_command(TARGET QISS POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${FFTW_LIB_DIR}/libfftw3-3.dll
        $<TARGET_FILE_DIR:QISS>
)

# -----------------------------
# 6. GLFW linking
# -----------------------------
set(GLFW_ROOT "${PROJECT_SOURCE_DIR}/extern/GLFW")
add_subdirectory(${GLFW_ROOT})
target_link_libraries(QISS PRIVATE glfw)

# Add ImGui backends (only the .cpp files)
set(IMGUI_BACKENDS
    ${CMAKE_SOURCE_DIR}/extern/ImGui/backends/imgui_impl_glfw.cpp
    ${CMAKE_SOURCE_DIR}/extern/ImGui/backends/imgui_impl_opengl3.cpp
)

target_sources(QISS PRIVATE ${IMGUI_BACKENDS})

# Include directories for the backend headers
target_include_directories(QISS PRIVATE
    ${CMAKE_SOURCE_DIR}/extern/ImGui/backends
)

# -----------------------------
# 7. Optional compile definitions
# -----------------------------
target_compile_definitions(QISS PRIVATE _CRT_SECURE_NO_WARNINGS)
