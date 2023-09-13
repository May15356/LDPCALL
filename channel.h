
#ifndef channel_hpp
#define channel_hpp

#include <iostream>
#include <random>
#include <string>
#include <math.h>
#include <chrono>
#include "Matrix.h"
#include "randomc.h"


using namespace std;

#define PI 3.14159
struct CN {
    int degree;             // the number of connected edges
    int* edgeIndex;
    int* vnIndex;
};

struct VN {
    int degree;         // the number of connected edges
    int* edgeIndex;
    int* cnIndex;
};

class codeSetting
{
public:
    int nMaxIter;
    int rep;
    const char* path;
    const char* path_code;
    int nCN;
    int nVN;
    int nMessage;
    struct CN* cnArray;
    struct VN* vnArray;
    long int nMostFrame;
    int nRep;
    double R;
    int nLeastErrorFrame;
    int nEbNo;
    double* EbN0;
    void initChannelSetting();
    //void ResetChannelSetting();
};

void codeSetting::initChannelSetting() {
    double EbNo[] = { 1.8,2,2.1,2.5,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.9,4.1,4.3,4.5,4.7,4.8,4.9,5,5.1 };
    nMaxIter = 150;
    rep = 1;
    nLeastErrorFrame = 50;
    nEbNo = sizeof(EbNo) / sizeof(double);
    EbN0 = new double[nEbNo];
    memcpy(EbN0,EbNo, sizeof(EbNo));
    path = "codeFile//H576_144.code";
    //path_code = "codeFile//code.txt";

    FILE* fp;
    fp=fopen(path, "r");
    fscanf(fp, "%d", &nVN);
    fscanf(fp, "%d", &nCN);
    vnArray = new VN[nVN];
    cnArray = new CN[nCN];
    int nullint;
    fscanf(fp, "%d", &nullint);
    fscanf(fp, "%d", &nullint);
    for (int i = 0; i < nVN; i++) {
        fscanf(fp, "%d", &vnArray[i].degree);
        vnArray[i].cnIndex = new int[vnArray[i].degree];
    }
    for (int i = 0; i < nCN; i++) {
        fscanf(fp, "%d", &cnArray[i].degree);
        cnArray[i].vnIndex = new int[cnArray[i].degree];
    }
    int temp;
    for (int thisVN = 0; thisVN < nVN; thisVN++) {
        for (int j = 0; j < vnArray[thisVN].degree; j++)
        {
            fscanf(fp, "%d", &temp);
            temp--;
            vnArray[thisVN].cnIndex[j] = temp;
        }
    }
    for (int thisCN = 0; thisCN < nCN; thisCN++) {
        for (int j = 0; j < cnArray[thisCN].degree; j++)
        {
            fscanf(fp, "%d", &temp);
            temp--;
            cnArray[thisCN].vnIndex[j] = temp;
        }
    }
    fclose(fp);
    nMessage = nVN - nCN;
    nMostFrame = 1000000000;
    nRep = nVN * rep;
    R = double(nMessage) / double(nVN) / double(rep);
    printf("EbNo,   wer,       rawber,       ber,    hardber,nErrorFrame,duration,averItr,hardaverItr \n");
    return;
}

/*Genrate 0-1 sequence of length k*/
void BinarySource(int nRow, int nCol, struct boolMatrix* mat) {

    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count(); // ns
    default_random_engine e(seed);
    bernoulli_distribution distribution;
    // extern CRandomMersenne rng;
    mat->nCol = nCol;
    mat->nRow = nRow;

    for (int i = 0; i < nRow; i++) {
        for (int j = 0; j < nCol; j++) {
            //mat->val[i][j] = rng.Random()<0;
            mat->val[i][j] = distribution(e);
        }
    }
    return;
}

/* BPSK modulate */
void BPSK(bool* message, double* bpsk, int n) {

    for (int j = 0; j < n; j++) {
        if (message[j])
            bpsk[j] = -1;
        else
            bpsk[j] = +1;
    }

    return;
}

// generate gaussian random noise
double gaussG(double sigma) {
    double x;
    double x1;
    static double x2;
    static int x2Valid = 0;
    extern CRandomMersenne rng;

    if (x2Valid) {
        x2Valid = 0;
        x2 *= sigma;
        return x2;
    }

    do {
        x1 = 2.0 * rng.Random() - 1.0;
        x2 = 2.0 * rng.Random() - 1.0;
        x = x1 * x1 + x2 * x2;
    } while (x >= 1.0 || x == 0);

    x1 *= sqrt((-2.0) * log(x) / x);
    x2 *= sqrt((-2.0) * log(x) / x);

    x2Valid = 1;

    x1 *= sigma;

    return x1;
}


/* Cross the awgn channel 0 mean and var var */
void awgn(double stddev, struct doubleMatrix* input, struct doubleMatrix* output) {
    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count(); // ns
    default_random_engine e(seed);
    normal_distribution<double> distribution(0.0, stddev);

    // double noise;
    for (int i = 0; i < output->nRow; i++) {
        for (int j = 0; j < output->nCol; j++) {
            // noise = gaussG(stddev);
            // output->val[i][j] = input->val[i][j] + noise;
            output->val[i][j] = input->val[i][j] + distribution(e);
        }
    }
    return;
}







#endif /* channel_hpp */
