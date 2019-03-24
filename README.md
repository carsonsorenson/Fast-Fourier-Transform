## Fast Fourier Transform

This program implements a recursive version of the fast fourier transform algorithm in order to compute polynomial multiplication between very large numbers.

The program first runs a test to confirm that the FFT is getting the same results as the naive version of polynomial multiplication and the slightly better using divide and conquer to solve three-subproblems.

Next the program is set on an infinite loop, where n (size of the polynomial) starts at 32 and will keep doubling until killed. The data from the program is outputted to a .csv file that was used to create the graphs.

The output of the program was expected, the FFT algorithm started off slower than the naive algorithm, however as the polynomials grew the FFT algorithm outperformed the naive algorithm.

![GitHub Logo](/graphs/linear_fit.png)



This program was written in c++ and was compiled using g++, I used the following commands to compile and run the program:

### Compile
    * g++ -std=c++14 -O3 -c main.cpp
    * g++ -std=c++14 -O3 main.o

### Run
    * ./a.out