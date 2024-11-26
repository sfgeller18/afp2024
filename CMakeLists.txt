cmake_minimum_required(VERSION 3.18)
project(AFP)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the HDF5 package (both C++ and C API)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Find Python 3.12 (adjust the version number if needed)
find_package(Python3 3.12 REQUIRED COMPONENTS Interpreter Development NumPy)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Add the executable
add_executable(afp ${SOURCES})

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_compile_options(afp PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(afp PRIVATE OpenMP::OpenMP_CXX)
endif()

# Link the libraries (HDF5, Python, and math library)
target_link_libraries(afp 
    PRIVATE
    HDF5::HDF5  # This will automatically link both C and C++ parts of the HDF5 library
    Python3::Python  # Links Python3.12
    m  # Math library
)

# Set include directories (assuming the "include" directory is in your project)
target_include_directories(afp
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Python3_INCLUDE_DIRS}  # Include directories for Python headers
    ${Python3_NumPy_INCLUDE_DIRS}  # Include directories for NumPy headers
)

# Optionally, specify the Python library directory if necessary
target_link_directories(afp
    PRIVATE
    ${Python3_LIBRARY_DIRS}  # Python library directory
)