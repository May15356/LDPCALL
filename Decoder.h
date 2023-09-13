// GallagerSPA.hpp
#ifndef GallagerSPA_hpp
#define GallagerSPA_hpp
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Matrix.h"
#include "Encoder.h"
using namespace std;
void quantifyUniform(double* llr, int length, int q, int f);
double min(double* arr, double* edgeAbsArray, int* cn_edge, int z, double VNedge);

class GallagerSPA {
public:
    int nCN;
    int nVN;
    int maxCNDegree;
    int nEdge;
    double* edgeAbsArray;
    bool* edgeSgnArray;
    bool* edgeCheck;
    struct CN* cnArray;
    struct VN* vnArray;
    int nMaxIter;
    double maxVal = 20.0;
    double minVal = 0.001;
    // register
    double* llr;
    double* llrTotal;
    double tmpLogSum;
    bool tmpSgn;
    double tmpval;
    bool thisSydrome;
    bool isCodeword;
    double sgnllrTotal;

    // function
    void probeOverflow(double val);
    void init(boolMatrix* H, int number);
    void SPAdecode(double* recieved, double gaussVar, bool* decoded, int rep);
    int MSdecodeshuffle(double* recieved, double gaussVar, bool* decoded, int rep);
    void TBMPdecode(double* recieved, double gaussVar, bool* decoded, int rep);
};

void GallagerSPA::probeOverflow(double val) {
    if (abs(val) > 20) {
        printf("%.2f", val);
    }
}

void GallagerSPA::init(boolMatrix* H, int number)
{
    nCN = H->nRow;
    nVN = H->nCol;
    nMaxIter = number;
    // Init edges
    doubleMatrix edgeIdMatrix;

    edgeIdMatrix.val = new double* [nCN];
    if (edgeIdMatrix.val == NULL) {
        printf("edgeIdMatrix.val == NULL!");
    }
    nEdge = 0;
    for (int i = 0; i < nCN; i++) {
        edgeIdMatrix.val[i] = new double[nVN];
        for (int j = 0; j < nVN; j++) {
            if (H->val[i][j]) {
                edgeIdMatrix.val[i][j] = nEdge;
                nEdge++;
            }
        }
    }
    edgeAbsArray = new double[nEdge];
    for (int i = 0; i < nEdge; i++)
        edgeAbsArray[i] = 0.0;
    edgeSgnArray = new bool[nEdge];
    edgeCheck = new bool[nEdge];
    // Init CN
    cnArray = new CN[nCN];
    int thisEdge;
    for (int i = 0; i < nCN; i++) {
        cnArray[i].degree = 0;
        // caculate degree of i-th CN
        for (int j = 0; j < nVN; j++) {
            if (H->val[i][j]) {
                cnArray[i].degree++;
            }
        }
        // store the location of the edges 
        cnArray[i].edgeIndex = new int[cnArray[i].degree];
        thisEdge = 0;
        for (int j = 0; j < nVN; j++) {
            if (H->val[i][j]) {
                cnArray[i].edgeIndex[thisEdge] = (int)edgeIdMatrix.val[i][j];
                thisEdge++;
            }
        }
    }
    // Init VN
    vnArray = new VN[nVN];
    for (int j = 0; j < nVN; j++) {
        // caculate the degree of j-th vn
        vnArray[j].degree = 0;
        for (int i = 0; i < nCN; i++) {
            if (H->val[i][j]) {
                vnArray[j].degree++;
            }
        }
        // store the location of the edges
        vnArray[j].edgeIndex = new int[vnArray[j].degree];
        thisEdge = 0;
        for (int i = 0; i < nCN; i++) {
            if (H->val[i][j]) {
                vnArray[j].edgeIndex[thisEdge] = (int)edgeIdMatrix.val[i][j];
                thisEdge++;
            }
        }
    }
    // Init the registers
    llr = new double[nVN];
    llrTotal = new double[nVN];
    return;
}

// See Algorithm5.1 in the book channel codes:classic and modern by Ryan and lin
void GallagerSPA::SPAdecode(double* recieved, double gaussVar, bool* decoded, int rep)
{
    // Initialization
    for (int j = 0; j < nVN; j++) {
        llr[j] = 0;
        for (int iRep = 0; iRep < rep; iRep++) {
            llr[j] += 2 * recieved[j * rep + iRep] / gaussVar;
        }
        for (int e = 0; e < vnArray[j].degree; e++)
        {
            edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(llr[j]);
            edgeSgnArray[vnArray[j].edgeIndex[e]] = (llr[j] < 0);
        }
    }
    //quantify(edgeAbsArray, nEdge);
    //for (int j = 0; j < nEdge; j++) {
    //    printf("%0.4f ", edgeAbsArray[j]);
    //}
    for (int nIter = 0; nIter < nMaxIter; nIter++)
    {
        // CN update
        int cn = 0;
        //quantify(edgeAbsArray, nEdge);
        // VN update  LLR total  
        for (int j = 0; j < nVN; j++) {
            //CN update
            for (int e = 0; e < vnArray[j].degree; e++) {
                cn = vnArray[j].cnIndex[e];
                tmpLogSum = 0;
                tmpSgn = 0;
                for (int cnDeg = 0; cnDeg < cnArray[cn].degree; cnDeg++) {
                    tmpLogSum = tmpLogSum + log(tanh(0.5 * edgeAbsArray[cnArray[cn].edgeIndex[cnDeg]]));
                    tmpSgn = tmpSgn ^ edgeSgnArray[cnArray[cn].edgeIndex[cnDeg]];
                }
                edgeAbsArray[vnArray[j].edgeIndex[e]] = 2 * atanh(pow(2.718, tmpLogSum - log(tanh(0.5 * edgeAbsArray[vnArray[j].edgeIndex[e]]))));
                edgeSgnArray[vnArray[j].edgeIndex[e]] = tmpSgn ^ edgeSgnArray[vnArray[j].edgeIndex[e]];
            }
            //VN update
            llrTotal[j] = llr[j];
            for (int e = 0; e < vnArray[j].degree; e++) {
                llrTotal[j] = llrTotal[j] + pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
            }
            sgnllrTotal = (llrTotal[j]) < 0 ? -1 : 1;
            llrTotal[j] = (abs(llrTotal[j]) > maxVal) ? (sgnllrTotal * maxVal) : llrTotal[j];
            llrTotal[j] = (abs(llrTotal[j]) < minVal) ? (sgnllrTotal * minVal) : llrTotal[j];
            decoded[j] = (llrTotal[j] < 0);
            //update
            for (int e = 0; e < vnArray[j].degree; e++) {
                tmpval = llrTotal[j] - pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
                edgeSgnArray[vnArray[j].edgeIndex[e]] = (tmpval < 0);
                edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(tmpval);
                edgeCheck[vnArray[j].edgeIndex[e]] = (decoded[j]);
            }
        }
        //quantify(edgeAbsArray, nEdge);
         // Stopping criteria
        bool isCheck;
        isCheck = false;
        for (int i = 0; i < nCN; i++) {
            for (int e = 0; e < cnArray[i].degree; e++) {
                isCheck = isCheck ^ edgeCheck[cnArray[i].edgeIndex[e]];
            }
            if (isCheck)
                break;
        }
        if (!isCheck) {
            break;
        }

    } // End of one iteration
    return;
}
int GallagerSPA::MSdecodeshuffle(double* recieved, double gaussVar, bool* decoded, int rep)
{
    double* arr;
    arr = new double[maxCNDegree];
    double delta = 0.625;
    // Initialization
    for (int j = 0; j < nVN; j++) {
        llr[j] = 0;
        for (int iRep = 0; iRep < rep; iRep++) {
            llr[j] += 2 * recieved[j * rep + iRep] / gaussVar;
        }
        for (int e = 0; e < vnArray[j].degree; e++)
        {
            edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(llr[j]);
            edgeSgnArray[vnArray[j].edgeIndex[e]] = (llr[j] < 0);
        }
    }
    //quantify(edgeAbsArray, nEdge);
    //for (int j = 0; j < nEdge; j++) {
    //    printf("%0.4f ", edgeAbsArray[j]);
    //}
    int nIter;
    for (nIter = 0; nIter < nMaxIter; nIter++)
    {
        // CN update
        double tmp_x=0;
        int cn = 0;
        //quantify(edgeAbsArray, nEdge);
        // VN update  LLR total  
        for (int j = 0; j < nVN; j++) {
            //CN update
            for (int e = 0; e < vnArray[j].degree; e++) {
                cn = vnArray[j].cnIndex[e];
                tmpLogSum = 0;
                tmpSgn = 0;
                tmp_x= min(arr,edgeAbsArray, cnArray[cn].edgeIndex, cnArray[cn].degree, edgeAbsArray[vnArray[j].edgeIndex[e]]);
                edgeAbsArray[vnArray[j].edgeIndex[e]] = delta * tmp_x;
                for (int cnDeg = 0; cnDeg < cnArray[cn].degree; cnDeg++) {     
                    tmpSgn = tmpSgn ^ edgeSgnArray[cnArray[cn].edgeIndex[cnDeg]];
                }
                edgeSgnArray[vnArray[j].edgeIndex[e]] = tmpSgn ^ edgeSgnArray[vnArray[j].edgeIndex[e]];
            }
            //VN update
            llrTotal[j] = llr[j];
            for (int e = 0; e < vnArray[j].degree; e++) {
                llrTotal[j] = llrTotal[j] + pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
            }
            sgnllrTotal = (llrTotal[j]) < 0 ? -1 : 1;
            llrTotal[j] = (abs(llrTotal[j]) > maxVal) ? (sgnllrTotal * maxVal) : llrTotal[j];
            llrTotal[j] = (abs(llrTotal[j]) < minVal) ? (sgnllrTotal * minVal) : llrTotal[j];
            decoded[j] = (llrTotal[j] < 0);
            //update
            for (int e = 0; e < vnArray[j].degree; e++) {
                tmpval = llrTotal[j] - pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
                edgeSgnArray[vnArray[j].edgeIndex[e]] = (tmpval < 0);
                edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(tmpval);
                edgeCheck[vnArray[j].edgeIndex[e]] = (decoded[j]);
            }
        }
        //quantify(edgeAbsArray, nEdge);
         // Stopping criteria
        bool isCheck;
        isCheck = false;
        for (int i = 0; i < nCN; i++) {
            for (int e = 0; e < cnArray[i].degree; e++) {
                isCheck = isCheck ^ edgeCheck[cnArray[i].edgeIndex[e]];
            }
            if (isCheck)
                break;
        }
        if (!isCheck) {
            break;
        }

    } // End of one iteration
    delete[] arr;
    return nIter;
}

void GallagerSPA::TBMPdecode(double* recieved, double gaussVar, bool* decoded, int rep)
{
    double* arr;
    arr = new double[28];
    double delta = 0.625;
    // Initialization
    for (int j = 0; j < nVN; j++) {
        llr[j] = 0;
        for (int iRep = 0; iRep < rep; iRep++) {
            llr[j] += 2 * recieved[j * rep + iRep] / gaussVar;
        }
        for (int e = 0; e < vnArray[j].degree; e++)
        {
            edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(llr[j]);
            edgeSgnArray[vnArray[j].edgeIndex[e]] = (llr[j] < 0);
        }
    }
    //quantify(edgeAbsArray, nEdge);
    //for (int j = 0; j < nEdge; j++) {
    //    printf("%0.4f ", edgeAbsArray[j]);
    //}
    int nIter;
    for (nIter = 0; nIter < nMaxIter; nIter++)
    {
        // CN update
        double tmp_x = 0;
        int cn = 0;
        //quantify(edgeAbsArray, nEdge);
        // VN update  LLR total  
        for (int j = 0; j < nVN; j++) {
            //CN update
            for (int e = 0; e < vnArray[j].degree; e++) {
                cn = vnArray[j].cnIndex[e];
                tmpLogSum = 0;
                tmpSgn = 0;
                tmp_x = min(arr, edgeAbsArray, cnArray[cn].edgeIndex, cnArray[cn].degree, edgeAbsArray[vnArray[j].edgeIndex[e]]);
                edgeAbsArray[vnArray[j].edgeIndex[e]] = delta * tmp_x;
                for (int cnDeg = 0; cnDeg < cnArray[cn].degree; cnDeg++) {
                    tmpSgn = tmpSgn ^ edgeSgnArray[cnArray[cn].edgeIndex[cnDeg]];
                }
                edgeSgnArray[vnArray[j].edgeIndex[e]] = tmpSgn ^ edgeSgnArray[vnArray[j].edgeIndex[e]];
            }
            //VN update
            llrTotal[j] = llr[j];
            for (int e = 0; e < vnArray[j].degree; e++) {
                llrTotal[j] = llrTotal[j] + pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
            }
            sgnllrTotal = (llrTotal[j]) < 0 ? -1 : 1;
            llrTotal[j] = (abs(llrTotal[j]) > maxVal) ? (sgnllrTotal * maxVal) : llrTotal[j];
            llrTotal[j] = (abs(llrTotal[j]) < minVal) ? (sgnllrTotal * minVal) : llrTotal[j];
            decoded[j] = (llrTotal[j] < 0);
            //update
            for (int e = 0; e < vnArray[j].degree; e++) {
                tmpval = llrTotal[j] - pow(-1, edgeSgnArray[vnArray[j].edgeIndex[e]]) * edgeAbsArray[vnArray[j].edgeIndex[e]];
                edgeSgnArray[vnArray[j].edgeIndex[e]] = (tmpval < 0);
                edgeAbsArray[vnArray[j].edgeIndex[e]] = abs(tmpval);
                edgeCheck[vnArray[j].edgeIndex[e]] = (decoded[j]);
            }
        }
        //quantify(edgeAbsArray, nEdge);
         // Stopping criteria
        bool isCheck;
        isCheck = false;
        for (int i = 0; i < nCN; i++) {
            for (int e = 0; e < cnArray[i].degree; e++) {
                isCheck = isCheck ^ edgeCheck[cnArray[i].edgeIndex[e]];
            }
            if (isCheck)
                break;
        }
        if (!isCheck) {
            break;
        }

    } // End of one iteration
    delete[] arr;
}

void initdecoder(GallagerSPA* decoder, const char* path, int nMaxIter) {

    FILE* fp;

    fp=fopen(path, "r");
    fscanf(fp, "%d", &decoder->nVN);
    fscanf(fp, "%d", &decoder->nCN);

    decoder->vnArray = new VN[decoder->nVN];
    decoder->cnArray = new CN[decoder->nCN];

    int nullint;
    fscanf(fp, "%d", &nullint);
    fscanf(fp, "%d", &decoder->maxCNDegree);

    for (int i = 0; i < decoder->nVN; i++) {
        fscanf(fp, "%d", &decoder->vnArray[i].degree);
        decoder->vnArray[i].edgeIndex = new int[decoder->vnArray[i].degree];
        decoder->vnArray[i].cnIndex = new int[decoder->vnArray[i].degree];
    }
    for (int i = 0; i < decoder->nCN; i++) {
        fscanf(fp, "%d", &decoder->cnArray[i].degree);
        decoder->cnArray[i].edgeIndex = new int[decoder->cnArray[i].degree];
    }

    decoder->llr = new double[decoder->nVN];
    decoder->llrTotal = new double[decoder->nVN];

    int indexEdge = 0; // index of edge 
    int* cnthisedge;
    cnthisedge = new int[decoder->nCN];
    for (int i = 0; i < decoder->nCN; i++) {
        cnthisedge[i] = 0;
    }
    int thisCN;
    for (int thisVN = 0; thisVN < decoder->nVN; thisVN++) {
        for (int j = 0; j < decoder->vnArray[thisVN].degree; j++)
        {

            fscanf(fp, "%d", &thisCN);
            thisCN--;
            decoder->vnArray[thisVN].cnIndex[j] = thisCN;
            decoder->vnArray[thisVN].edgeIndex[j] = indexEdge;
            decoder->cnArray[thisCN].edgeIndex[cnthisedge[thisCN]] = indexEdge;
            cnthisedge[thisCN]++;
            indexEdge++;
        }
    }
    decoder->nEdge = indexEdge + 1;
    decoder->edgeAbsArray = new double[decoder->nEdge];
    for (int i = 0; i < decoder->nEdge; i++)
        decoder->edgeAbsArray[i] = 0.0;
    decoder->edgeSgnArray = new bool[decoder->nEdge];
    decoder->edgeCheck = new bool[decoder->nEdge];

    delete[]cnthisedge;

    decoder->nMaxIter = nMaxIter;
    fclose(fp);
    return;
}
double min(double* arr,double* edgeAbsArray, int* cn_edge, int z, double VNedge)
{
    for (int i = 0; i < z; i++)
    {
        arr[i] = edgeAbsArray[cn_edge[i]];
    }
    double temp = 10000;
    int m = z, n = z;
    for (int i = 0; i < z; i++)
   {    
        if (arr[i] == VNedge)
            continue;
        temp = (temp < arr[i]) ? temp : arr[i];
    }
    return temp;
}
void quantifyUniform(double* llr, int length, int q, int f) {
    double maxnumber = pow(2, q - 1) - 1;
    double temp = 0.0;
    int b = 0;
    int c = 0;
    double inv_delta;
    double delta;
    inv_delta = pow(2, f);
    delta = 1 / inv_delta;
    for (int i = 0; i < length; i++)
    {
        temp = llr[i] * inv_delta;
        if (temp > maxnumber)
            temp = maxnumber;
        if (llr[i] > 0)
            c = (int)(temp + 0.5);
        else
            c = (int)(temp - 0.5);
        temp = (double)c;
        llr[i] = temp * delta;
        if (abs(llr[i]) < 0.001)
            llr[i] = (llr[i] < 0) ? -0.001 : 0.001;
    }
    return;
}

#endif