#ifndef BasicFunction_hpp
#define BasicFunction_hpp
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Encoder.h"
#include "Matrix.h"
int nchoosek(int n, int k) {
	int ans = 1;
	for (int i = 0; i < k; i++) {
		ans = ans * (n - i);
	}
	for (int i = 1; i < (k + 1); i++) {
		ans = ans / i;
	}
	return ans;
}

void nchoosekInList(int n, int k, int nComb, int** combList) {
	/* same function as nchoosek(0:(n-1), k) in matlab */

		// The first combination
	int p = 0; // pointer 
	for (int i = 0; i < k; i++) {
		combList[p][i] = i;
	}

	//

	while (true) {
		p = p + 1;
		for (int i = k - 1; i > -1; i--) {
			if (combList[p - 1][i] < (n - k + i)) {
				combList[p][i] = combList[p - 1][i] + 1;
				for (int j = i + 1; j < k; j++) {
					combList[p][j] = combList[p][j - 1] + 1;
				}
				for (int j = 0; j < i; j++) {
					combList[p][j] = combList[p - 1][j];
				}
				break;
			}
		}
		if (combList[p][0] == n - k) {
			break;
		}
	}
	return;
}


bool isCorrectFrame(bool* encoded, bool* decoded, int n, int &biterror)
{
	bool flag = true;
	for (int i = 0; i < n; i++) {
		if (encoded[i] != decoded[i])
		{
			biterror++;
			flag= false;
		}
	}
	return flag;
}

void initEncode(const char* path, int n, boolMatrix* encoded) {
	FILE* fp2;
	fp2=fopen(path, "r");
	char c;
	for (int i = 0; i < 1; i++) {
		for (int j = 0; j < n; j++) {
			while (true) {
				fscanf(fp2, "%c", &c);
				if (c == '1' || c == '0') {
					break;
				}
			}
			if (c == '1')
				encoded->val[i][j] = 1;
			else
				encoded->val[i][j] = 0;
		}
	}
	fclose(fp2);
}

class decoderSetting {
public:
	int nErrorFrame;
	double wer;
	clock_t startTime, endTime;
	double duration;
	int biterror;
	double ber;
	int rawbiterror;
	double rawber;
	int hardbiterror;
	double hardber;
	int BP;
	double BPwer;
	double averItr;
	double hardaverItr;
	void init();
	void computation(codeSetting* channel, double iFrame, double EbNo);
};

void decoderSetting::init()
{
	nErrorFrame = 0;
	biterror = 0;
	ber = 0;
	rawbiterror = 0;
	rawber = 0;
	hardbiterror = 0;
	hardber = 0;
	BP = 0;
	BPwer = 0;
	averItr = 0;
	hardaverItr = 0;
	startTime = clock();
	return;
}

void decoderSetting::computation(codeSetting* channel, double iFrame, double EbNo)
{
	averItr = (double)averItr / (double)(iFrame + 1);
	hardaverItr = (double)hardaverItr / (double)(iFrame + 1);
	wer = (double)nErrorFrame / (double)(iFrame + 1);
	ber = (double)biterror / (double)((iFrame + 1) * channel->nVN);
	rawber = (double)rawbiterror / (double)((iFrame + 1) * channel->nVN);
	hardber = (double)hardbiterror / (double)((iFrame + 1) * channel->nVN);
	endTime = clock();
	duration = ((double)endTime - (double)startTime) / CLOCKS_PER_SEC;
	printf(" %.2f   %.3e  %.3e  %.3e %.3e  %d    %.2f     %.2f     %.2f \n", EbNo, wer, rawber, ber, hardber, nErrorFrame, duration, averItr, hardaverItr);
	return;
}
#endif