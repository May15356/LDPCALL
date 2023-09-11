#ifndef BoolMatrix_hpp
#define BoolMatrix_hpp

#include<stdio.h>
#include<iostream>
#include <stdlib.h>
#include <random>
#include <string>
#include <math.h>
#include <chrono>
using namespace std;

struct boolMatrix
{
    int nRow;
    int nCol;
    bool** val;
};

struct doubleMatrix
{
    int nRow;
    int nCol;
    double** val;
};

void boolMatrixInit(int nRow, int nCol, struct boolMatrix* mat) {
    mat->nCol = nCol;
    mat->nRow = nRow;
    mat->val = new bool* [nRow];
    for (int i = 0; i < nRow; i++)
        mat->val[i] = new bool[nCol];

    return;
}

void doubleMatrixInit(int nRow, int nCol, struct doubleMatrix* mat) {
    mat->nCol = nCol;
    mat->nRow = nRow;
    mat->val = new double* [nRow];
    for (int i = 0; i < nRow; i++)
        mat->val[i] = new double[nCol];

    return;
}

void doubleMatrixDelete(struct doubleMatrix* mat) {
    for (int i = 0; i < mat->nRow; i++) {
        delete[]mat->val[i];
    }
    delete[]mat->val;
    return;
}

void boolMatrixDelete(struct boolMatrix* mat) {
    for (int i = 0; i < mat->nRow; i++) {
        delete[]mat->val[i];
    }
    delete[]mat->val;
    return;
}

void sortDescent(double* thisArray, int n, int* order) {

    double tmpval;
    double* tmpArray;
    tmpArray = new double[n];
    for (int i = 0; i < n; i++) {
        tmpArray[i] = thisArray[i];
    }

    for (int i = 0; i < n; i++) {

        for (int j = (i + 1); j < n; j++) {
            if (thisArray[i] < thisArray[j]) {
                tmpval = thisArray[i];
                thisArray[i] = thisArray[j];
                thisArray[j] = tmpval;

            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (thisArray[i] == tmpArray[j]) {
                tmpArray[j] = NULL;
                order[i] = j;
                break;
            }
        }
    }
    delete[]tmpArray;
    return;
}

/* mat1 X mat2 = mat3*/
void boolMatrixProduct(struct boolMatrix* mat1, struct boolMatrix* mat2, struct boolMatrix* mat3)
{
    bool isValid;
    isValid = (mat1->nCol == mat2->nRow);
    if (!isValid) {
        printf("Sizes of matrixes do not match!");
        return;
    }

    mat3->nCol = mat2->nCol;
    mat3->nRow = mat1->nRow;
    for (int i = 0; i < mat3->nRow; i++) {
        for (int j = 0; j < mat3->nCol; j++) {
            mat3->val[i][j] = false;
            for (int k = 0; k < mat1->nCol; k++)
                mat3->val[i][j] = mat3->val[i][j] ^ (mat1->val[i][k] & mat2->val[k][j]);
        }
    }

    return;

}

void boolMatrixInv(struct boolMatrix* mat, struct boolMatrix* matinv) {
    matinv->nCol = mat->nRow;
    matinv->nRow = mat->nCol;
    boolMatrixInit(matinv->nRow, matinv->nCol, matinv);
    for (int i = 0; i < matinv->nRow; i++) {
        for (int j = 0; j < matinv->nCol; j++) {
            matinv->val[i][j] = mat->val[j][i];
        }
    }

    return;
}

bool isFalseMatrix(struct boolMatrix* mat) {
    for (int i = 0; i < mat->nRow; i++) {
        for (int j = 0; j < mat->nCol; j++) {
            if (mat->val[i][j]) {
                cout << i << ' ' << j << endl;
                return false;
            }
        }
    }
    return true;
}


#endif