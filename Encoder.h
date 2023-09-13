

#ifndef Encoder_hpp
#define Encoder_hpp
#include "channel.h"
#include "Matrix.h"
using namespace std;


void initGauss(codeSetting* parityCheck, boolMatrix* gauss, int* exchange) {

	for (int i = 0; i < parityCheck->nCN; i++)
		for (int j = 0; j < parityCheck->nVN; j++)
			gauss->val[i][j] = 0;

	for (int i = 0; i < parityCheck->nCN; i++)
		for (int j = 0; j < parityCheck->cnArray[i].degree; j++)
			gauss->val[i][parityCheck->cnArray[i].vnIndex[j]] = 1;
}


void countGauss(codeSetting* parityCheck, boolMatrix* gauss, int* exchange, int& rank) {
	int  i, j, k, temp, col, changeCol, swap;
	bool unfullrank;

	int n = parityCheck->nVN;
	int m = parityCheck->nCN;
	for (col = n - 1; col >= n - m; col--) {
		temp = -1;
		unfullrank = false;

		for (i = 0; i <= m - n + col; i++) {
			if (gauss->val[i][col] == 1) {
				temp = i;
				break;
			}
		}
		if (temp == -1) {
			changeCol = -1;

			for (j = col - 1; j >= 0; j--) {
				for (k = 0; k <= m - n + col; k++) {
					if (gauss->val[k][j] == 1) {
						changeCol = j;
						break;
					}
				}
				if (changeCol != -1)
					break;
			}

			if (changeCol == -1)
				unfullrank = true;
			else {
				for (j = 0; j <= m - 1; j++) {
					swap = gauss->val[j][changeCol];
					gauss->val[j][changeCol] = gauss->val[j][col];
					gauss->val[j][col] = swap;
				}
				swap = exchange[changeCol];
				exchange[changeCol] = exchange[col];
				exchange[col] = swap;
				col++;
			}
		}
		else if (temp == m - n + col) {
			for (i = 0; i < m; i++) {
				if (i != temp) {
					if (gauss->val[i][col] == 1) {
						for (j = 0; j < n; j++)
							gauss->val[i][j] = (gauss->val[i][j] + gauss->val[temp][j]) % 2;
					}
				}
			}
		}
		else {
			for (i = 0; i < m; i++) {
				if (i != temp && i != m - n + col) {
					if (gauss->val[i][col] == 1) {
						for (j = 0; j < n; j++)
							gauss->val[i][j] = (gauss->val[i][j] + gauss->val[temp][j]) % 2;
					}
				}
			}

			if (gauss->val[m - n + col][col] == 1) {
				for (j = 0; j < n; j++)
					gauss->val[temp][j] = (gauss->val[m - n + col][j] + gauss->val[temp][j]) % 2;
			}
			else {
				for (j = 0; j < n; j++) {
					swap = gauss->val[temp][j];
					gauss->val[temp][j] = gauss->val[m - n + col][j];
					gauss->val[m - n + col][j] = swap;
				}
			}
		}

		if (unfullrank) {
			rank = n - col - 1;
			printf("%s", "Rank: ");
			printf("%d", rank);
			break;
		}
		else rank = parityCheck->nCN;
	}
	parityCheck->nMessage = parityCheck->nMessage + parityCheck->nCN - rank;
	//printf("\r%s%d", "Information Bits: ", parityCheck->nMessage);
}


// check whether the eliminated H matrix is right, the matrix should be like [p I] mode if it's rankfull or like modified [p I] mode with all-zero rows in the first several row if it's not rankfull.
void checkGauss(codeSetting* parityCheck, boolMatrix* gauss, int& rank) {

	int i, j;
	bool one, zero;
	printf("\n");

	if (rank == parityCheck->nCN) {
		one = true;
		for (j = parityCheck->nVN - parityCheck->nCN; j < parityCheck->nVN; j++) {
			i = j - parityCheck->nVN + parityCheck->nCN;
			if (gauss->val[i][j] == 0) one = false;
		}

		if (one) printf("matrix 1 right; ");
		else printf("matrix 1 false; ");

		printf("\n");
		zero = true;

		for (j = parityCheck->nVN - parityCheck->nCN + 1; j < parityCheck->nVN; j++)
		{
			for (i = 0; i < (j - parityCheck->nVN + parityCheck->nCN); i++)
				if (gauss->val[i][j] == 1) zero = false;
		}

		if (zero) printf("matrix 0 right;\n ");
		else printf("matrix 0 false; ");
		printf("\n");
	}
	else {
		one = true;
		for (j = parityCheck->nVN - rank; j < parityCheck->nVN; j++) {
			i = j - parityCheck->nVN + parityCheck->nCN;
			if (gauss->val[i][j] == 0) one = false;//cout<<i<<j;
		}
		if (one) printf("matrix 1 right; ");
		else printf("matrix 1 false; ");

		zero = true;
		for (j = 0; j < parityCheck->nVN; j++) {
			for (i = 0; i < (parityCheck->nCN - rank); i++) {
				if (gauss->val[i][j] == 1) zero = false;//cout<<i<<j;
			}
		}
		if (zero) printf("matrix 0 right;\n ");
		else printf("matrix 0 false; ");
	}
}

void initMessage(boolMatrix* message, int n)
{
	for(int i=0;i<n;i++)
		message->val[0][i] = rand() % 2;
}

void encoder(codeSetting* parityCheck, boolMatrix* parityc_mat, int* exchange, bool* message, bool* code) {

	char* codeSwap;
	int m, temp;

	codeSwap = (char*)malloc(sizeof(char) * (parityCheck->nVN));

	for (int i = 0; i < parityCheck->nVN; i++)
		codeSwap[i] = 0;

	//code equals to message in message bit
	for (int i = 0; i < parityCheck->nMessage; i++) {
		if (message[i] == 0)
			code[i] = 0;
		else
			code[i] = 1;
	}

	// encode using the relationship between H matrix and G matrix
	for (int j = parityCheck->nMessage; j < parityCheck->nVN; j++) {
		temp = 0;
		m = j - parityCheck->nVN + parityCheck->nCN;
		for (int col = 0; col < parityCheck->nMessage; col++)
			temp = (temp + parityc_mat->val[m][col] * code[col]) % 2;
		code[j] = temp;
	}


	//exchage the changed columns in gauss elimination
	for (int i = 0; i < parityCheck->nVN; i++)
		codeSwap[exchange[i]] = code[i];
	for (int i = 0; i < parityCheck->nVN; i++)
		code[i] = codeSwap[i];

	free(codeSwap);

}

// check if the encoded code satifies the parity check matrix
bool checkCode(codeSetting* parityCheck, bool* code) {

	bool codeValid = true;
	int temp;

	for (int i = 0; i < parityCheck->nCN; i++) {
		temp = 0;
		for (int j = 0; j < parityCheck->cnArray[i].degree; j++)
			temp = temp + code[parityCheck->cnArray[i].vnIndex[j]];
		temp = temp % 2;
		if (temp > 0)
		{
			codeValid = false;
			break;
		}
	}
	return codeValid;
}
#endif