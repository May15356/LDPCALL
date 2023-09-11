#include<stdio.h>
#include<iostream>
using namespace std;

#include <math.h>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <thread>

#include "channel.h"
#include "Matrix.h"
#include "Decoder.h"
#include "randomc.h"
#include "BasicFunction.h"
CRandomMersenne rng(1);

void simulateSPA(codeSetting* channel, double iEbNo, boolMatrix* parityc_mat, int* exchange);

int main() {
    double EbNo[ ] = { 1.8,2,2.1,2.5,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.9,4.1,4.3,4.5,4.7,4.8,4.9,5,5.1 };
    int nEbNo = sizeof(EbNo) / sizeof(double);
    codeSetting* channel = new codeSetting;
    channel->initChannelSetting();
    boolMatrix* parityc_mat = new boolMatrix;
    boolMatrixInit(channel->nCN, channel->nVN, parityc_mat);
    int* exchange;     // row exchange recoding in gauss elimination for 
    int rank;          // used to store the rank of parity check matirx
    exchange = (int*)malloc(sizeof(int) * (channel->nVN));
    for (int i = 0; i < channel->nVN; i++)
        exchange[i] = i;
    initGauss(channel, parityc_mat, exchange);
    countGauss(channel, parityc_mat, exchange, rank);
    checkGauss(channel, parityc_mat, rank);
    for (int i = 0; i < nEbNo; i++)
    {
        double iEbNo = EbNo[i];
        thread th0(simulateSPA, channel, iEbNo, parityc_mat,exchange);
        thread th1(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th2(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th3(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th4(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th5(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th6(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th7(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th8(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        thread th9(simulateSPA, channel, iEbNo, parityc_mat, exchange);
        /**/
        th0.join();
        th1.join();
        th2.join();
        th3.join();
        th4.join();
        th5.join();
        th6.join();
        th7.join();
        th8.join();
        th9.join();
        /**/
        //printf("errorframe= %d", errorframe[i]);
    }
    
    delete channel;
    system ("pause");
    return 0;
}

void simulateSPA(codeSetting* channel, double EbNo, boolMatrix* parityc_mat, int* exchange)
{
    // Init decoder 
    GallagerSPA* decoder=new GallagerSPA;
    decoderSetting* varSet = new decoderSetting;
    varSet->init();
    double sigma;
    double var;
    sigma = 1 / sqrt(2 * channel->R) * pow(10, -EbNo/ 20);
    var = pow(sigma, 2);
    boolMatrix* message=new boolMatrix;
    boolMatrix* encoded = new boolMatrix;
    doubleMatrix* transimit = new doubleMatrix;
    doubleMatrix* received = new doubleMatrix;
    doubleMatrix* hardreceived = new doubleMatrix;
    boolMatrix* decoded = new boolMatrix;

    boolMatrixInit(1, channel->nVN, message);
    boolMatrixInit(1, channel->nVN, encoded);
    doubleMatrixInit(1, channel->nRep, transimit);
    doubleMatrixInit(1, channel->nRep, received);
    boolMatrixInit(1, channel->nVN, decoded);
    doubleMatrixInit(1, channel->nRep, hardreceived);


    initdecoder(decoder, channel->path, channel->nMaxIter);

    for (int iFrame = 0; iFrame < channel->nMostFrame; iFrame++)
    {
        bool flagerror = false;
        //encoder
        initMessage(message,channel->nVN);
        encoder(channel, parityc_mat, exchange, message->val[0], encoded->val[0]);
        if (checkCode(channel, encoded->val[0]) == false) {
            printf("Code Invalid!!!\n");
            break;
        }
        //channel
        BPSK(encoded->val[0], transimit->val[0], channel->nVN);
        awgn(sigma, transimit, received);
        //rawber
        for (int i = 0; i < channel->nVN; i++) {
            hardreceived->val[0][i] = received->val[0][i] > 0 ? 0 : 1;
            if ((char)hardreceived->val[0][i] != encoded->val[0][i])
            {
                varSet->rawbiterror++;
            }
        }
        for (int i = 0; i < channel->nVN; i++) {
            hardreceived->val[0][i] = received->val[0][i] > 0 ? 1 : -1;
        }
        //ber
        decoder->SPAdecode(received->val[0], var, decoded->val[0], channel->rep);
        //varSet->averItr= varSet->averItr+(double)decoder->MSdecodeshuffle(received->val[0], var, decoded->val[0], channel->rep);
        flagerror=!isCorrectFrame(encoded->val[0], decoded->val[0], channel->nVN, varSet->biterror);
        //hardber
        varSet->hardaverItr = varSet->hardaverItr + (double)decoder->MSdecodeshuffle(hardreceived->val[0], var, decoded->val[0], channel->rep);
        isCorrectFrame(encoded->val[0], decoded->val[0], channel->nVN, varSet->hardbiterror);
        if (flagerror) {
            varSet->nErrorFrame++;
        }
        if (varSet->nErrorFrame == channel->nLeastErrorFrame || iFrame == (channel->nMostFrame - 1)) {
            varSet->computation(channel,iFrame,EbNo);
            break;
        }
    }    
    delete message;
    delete encoded;
    delete received;
    delete hardreceived;
    delete decoded;
    delete transimit;
    delete decoder;
    delete varSet;
    return;
} 


