# cp2-feynman-ppiep

To compile the code your system needs the following applications:

- c compiler (gcc or clang recommended)
- cmake (3.13.4 or higher)
- make (4.2.1 or higher)
- fftw3 (explicitly set path in `source/cmake/Findfftw3.cmake`)

Clone the repository and navigate to the project root.
Build with `cmake -S source -B build && make`. The main executable will be `build/main[.exe]` and
the tests will be located in `build/tests/`.
