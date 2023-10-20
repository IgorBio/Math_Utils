# Math Utils

`MathUtils` is a C library that provides basic mathematical operations on floating-point numbers. It implements functions from the standard `math.h` library using Taylor series expansions.

## Features

- Basic mathematical functions (`abs`, `fabs`, `ceil`, `floor`, `trunc`, `sqrt`)
- Trigonometric functions (`sin`, `cos`, `tan`, `asin`, `acos`, `atan`)
- Exponential functions (`exp`, `pow`)
- Logarithmic functions (`log`)
- Special constants (`π`, `e`, square roots, golden ratio, Catalan's constant, Cahen's constant)

## Usage

To use `MathUtils` in your C project, follow these steps:

1. Clone the repository:
To get started with `MathUtils`, you can clone the repository to your local machine using the following command:

```bash
git clone https://github.com/IgorBio/Math_Utils.git
```
2. Building the Library:
After cloning the repository, you can build the library using the provided `Makefile`. Navigate to the repository directory in your terminal and run:

```bash
make build
```
This command will compile the source files and generate the static library `libMathUtils.a` in the lib folder.

3. Integrating `MathUtils` in CMake Projects
   To include `MathUtils` in your CMake project, follow these steps:

- Copy the `libMathUtils.a` library file to your project's library folder, for example, libs.

- Add the following lines to your project's `CMakeLists.txt` file:

```cmake
# Specify the path to the MathUtils library
set(MATH_UTILS_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/libs/libMathUtils.a)

# Add the MathUtils library to your project
add_library(MathUtils STATIC IMPORTED)
set_target_properties(MathUtils PROPERTIES IMPORTED_LOCATION ${MATH_UTILS_LIBRARY_PATH})

# Link the MathUtils library to your target
target_link_libraries(<your_target_name> MathUtils)
```
Replace <your_target_name> with the actual name of your CMake target.

4. Including `MathUtils` Header

In your C source files where you want to use `MathUtils` functions, include the `math_utils.h` header file at the beginning of the file:

```c
#include "math_utils.h"
```

## Documentation

Check the library documentation for specific function details and usage examples.

```bash
make dvi
```