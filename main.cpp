#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

//function to generate a random double between -1, and 1
vector<long double> randomDouble(int n) {
    vector<long double> v;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);
    for (int i = 0; i < n; i++) {
        v.push_back(dis(gen));
    }
    return v;
}

//kept in to compare times, and to make sure the fft algorithm is correct
//double for loop implements the "high school" polynomial multiplication. returns the product array pq
vector<long double> polynomialMult(vector<long double>& p, vector<long double>& q) {
    size_t size = p.size() + q.size() - 1;
    vector<long double>pq(size);

    for (unsigned int i = 0; i < p.size(); i++) {
        for (unsigned int j = 0; j < q.size(); j++) {
            pq[i + j] = p[i] * q[j] + pq[i + j];
        }
    }
    return pq;
}

//the three sub problem algorithm
vector<long double> polynomialMult2(vector<long double>&p, vector<long double>&q, int n) {
    if (n == 1) {
        return vector<long double>{p[0] * q[0]};
    }
    int half = n / 2;
    vector<long double> pFirstHalf(half);
    vector<long double> pSecondHalf(half);
    vector<long double> qFirstHalf(half);
    vector<long double> qSecondHalf(half);
    for (int i = 0; i < half; i++) {
        pFirstHalf[i] = p[i];
        pSecondHalf[i] = p[i + half];
        qFirstHalf[i] = q[i];
        qSecondHalf[i] = q[i + half];
    }

    vector<long double> pTermsAdded(half);
    vector<long double> qTermsAdded(half);

    for (int i = 0; i < half; i++) {
        pTermsAdded[i] = pFirstHalf[i] + pSecondHalf[i];
        qTermsAdded[i] = qFirstHalf[i] + qSecondHalf[i];
    }

    //make the 3 recursive calls
    vector<long double> pqFirstHalf = polynomialMult2(pFirstHalf, qFirstHalf, half);
    vector<long double> middlePart = polynomialMult2(pTermsAdded, qTermsAdded, half);
    vector<long double> pqSecondHalf = polynomialMult2(pSecondHalf, qSecondHalf, half);

    //store the final solution in pq
    vector<long double> pq(2 * n - 1);

    //pq=A+(B+C)x^n/2+Dx^n
    //=> p0q0 + ((p0q0)-(p1q1)-p0q0-p1q1) + p1q1
    for (int i = 0; i < n - 1; i++) {
        pq[i] += pqFirstHalf[i];
        pq[i + half] += middlePart[i] - pqFirstHalf[i] - pqSecondHalf[i];
        pq[i + (2 * half)] += pqSecondHalf[i];
    }
    return pq;
}

//function to compute omega so we don't need to include it in the timing of our algorithm
vector<complex<long double>> computeOmega(int n) {
    vector<complex<long double>> V;
    for (int i = 0; i < 2 * n; i++) {
        complex<long double> w_i(cos((2 * M_PI * i) / (2 * n)), sin((2 * M_PI * i) / (2 * n)));
        V.push_back(w_i);
    }
    return V;
}

//recursive section of the FFT algorithm
vector<complex<long double>> FFT(vector<complex<long double>>& P, vector<complex<long double>>& V, int n) {
    if (n == 1) {
        return P;
    }

    //calculate even and odd sections
    vector<complex<long double>> pEven(n / 2);
    vector<complex<long double>> pOdd(n / 2);
    for (int i = 0; i < (n / 2); i++) {
        pEven[i] = P[2 * i];
        pOdd[i] = P[2 * i + 1];
    }

    vector<complex<long double>> VSquared(n / 2);
    for (int i = 0; i < (n / 2); i++) {
        VSquared[i] = V[i] * V[i];
    }

    //recursive calls
    vector<complex<long double>> solEven = FFT(pEven, VSquared, n / 2);
    vector<complex<long double>> solOdd = FFT(pOdd, VSquared, n / 2);

    vector<complex<long double>> sol(n);
    for (int i = 0; i < (n / 2); i++) {
        sol[i] = solEven[i] + V[i] * solOdd[i];
        sol[i + n / 2] = solEven[i] - V[i] * solOdd[i];
    }

    return sol;
}

//Main section of the fft algorithm
vector<long double> polyMultFFT(vector<complex<long double>>&P, vector<complex<long double>>&Q, vector<complex<long double>>&V, vector<complex<long double>>&VInverse, int n) {
    //pad the two vectors with 0's
    for (int i = 0; i < n; i++) {
        P.push_back(0);
        Q.push_back(0);
    }

    //compute the two solutions
    vector<complex<long double>> solP = FFT(P, V, 2 * n);
    vector<complex<long double>> solQ = FFT(Q, V, 2 * n);

    vector<complex<long double>> solPQ(2 * n);
    for (int i = 0; i < 2 * n; i++) {
        solPQ[i] = solP[i] * solQ[i];
    }

    vector<complex<long double>> PQ = FFT(solPQ, VInverse, 2 * n);
    vector<long double> finalSolution;

    //PQ is the inverse so we take the real part of the vector and divide it by (2 * n) this returns the correct value
    for (int i = 0; i < 2 * n; i++) {
        finalSolution.push_back(real(PQ[i]) / double(2 * n));
    }

    //this is the solution with only the real values, gets rid of the imaginary values
    return finalSolution;

}


/* EMPIRACLE STUDIES SECTION */
void compileDataForTheThreeAlgorithms() {
    int n = 32;
    double totalOne = 0;
    double totalTwo = 0;
    double totalFFT = 0;
    ofstream outStream;
    outStream.open("data.csv");
    outStream << "n,Naive(nanoseconds),3 sub problem(nanoseconds),FFT(nanoseconds)\n";
    outStream.close();

    while (true) {
        outStream.open("data.csv", ios_base::app);
        cout << "2^" << log2(n) << endl;
        for (int i = 0; i < 10; i++) {
            vector<long double> p = randomDouble(n); //create 2 random vectors
            vector<long double> q = randomDouble(n);

            //fft algorithm
            vector<complex<long double>> fft_p;
            vector<complex<long double>> fft_q;
            for (int j = 0; j < n; j++) { //cast the array of doubles to an array of complex doubles, the imaginary part is just 0 to start
                fft_p.push_back(complex<long double>(p[j], 0));
                fft_q.push_back(complex<long double>(q[j], 0));
            }
            vector<complex<long double>> omega = computeOmega(n);
            //compute the inverse of omega by using the conj() function
            vector<complex<long double>> omegaInverse(n * 2);
            for (int i = 0; i < 2 * n; i++) {
                omegaInverse[i] = conj(omega[i]);
            }
            auto startFFT = chrono::steady_clock::now(); //start the time after we computeOmega
            vector<long double> solFFT = polyMultFFT(fft_p, fft_q, omega, omegaInverse, n);
            auto endFFT = chrono::steady_clock::now();
            auto totalTimeFFT = chrono::duration_cast<chrono::nanoseconds>(endFFT - startFFT).count();
            totalFFT += totalTimeFFT;
            cout << totalFFT << " ";
            solFFT.clear();

            //n log 3 algorithm
            auto startTwo = chrono::steady_clock::now();
            vector<long double> solTwo = polynomialMult2(p, q, n);
            auto endTwo = chrono::steady_clock::now();
            auto totalTimeTwo = chrono::duration_cast<chrono::nanoseconds>(endTwo - startTwo).count();
            totalTwo += totalTimeTwo;
            cout << totalTwo << " ";
            solTwo.clear();

            //n^2 algorithm
            auto startOne = chrono::steady_clock::now(); //start the time right before we call the polynomialMultiplication function
            vector<long double> solOne = polynomialMult(p, q);
            auto endOne = chrono::steady_clock::now();
            auto totalTimeOne = chrono::duration_cast<chrono::nanoseconds>(endOne - startOne).count(); //calculate the time the fucntion took up
            totalOne += totalTimeOne; //keep track of the total time it takes to run the 10 test
            cout << totalOne << endl;
            solOne.clear();

            p.clear();
            q.clear();
        }
        outStream << n << "," << totalOne << "," << totalTwo << "," << totalFFT << "\n";
        n *= 2;
        totalOne = 0;
        totalTwo = 0;
        totalFFT = 0;
        outStream.close();
    }
}

int main()
{

    bool runStudies = true;

    //confirm the algorithm is running the same as the simple n^2 algorithm.

    vector<long double> p1 { 0, 1, 2, 3};
    vector<long double> q1 = { 10, 11, 12, 13};
    int n = p1.size();

    vector<long double> result1 = polynomialMult(p1, q1);
    vector<complex<long double>> P;
    vector<complex<long double>> Q;
    for (int i = 0; i < n; i++) { //cast the array of doubles to an array of complex doubles, the imaginary part is just 0 to start
        P.push_back(complex<long double>(p1[i], 0));
        Q.push_back(complex<long double>(q1[i], 0));
    }
    vector<complex<long double>> omega = computeOmega(n);
    //compute the inverse of omega by using the conj() function
    vector<complex<long double>> omegaInverse(n * 2);
    for (int i = 0; i < 2 * n; i++) {
        omegaInverse[i] = conj(omega[i]);
    }
    vector<long double> result2 = polyMultFFT(P, Q, omega, omegaInverse, n);

    cout << "Simple polynomial algorithm result:\n";
    for (unsigned int i = 0; i < result1.size(); i++) {
        cout << result1[i] << " ";
    }
    cout << "\n\nFFT result:\n";
    for (unsigned int i = 0; i < result2.size() -1; i++) {
        cout << result2[i] << " ";
    }
    cout << "\n\n";

    if (runStudies){
        compileDataForTheThreeAlgorithms(); //create the raw data files to test when the fft algorithm becomes faster than the other two
    }
    return 0;
}
