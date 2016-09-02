#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "unafold.h"

int main(void) {

	srand(time(NULL));

	int N = 500, range = 10000;

	int local_Qprime_ip[500][500] = {0};
	int local_Q[500][500] = {0};
	int local_Qprime[500][500] = {0};
	int local_QM[500][500] = {0};

	int i, j;
	for (i = 0 ; i < N; i++) {
		for (j = 0; j < N; j++) {
			local_Qprime_ip[i][j] = rand() % range;
			local_Q[i][j] = rand() % range;
			local_Qprime[i][j] = rand() % range;
			local_QM[i][j] = rand() % range;
		}
	}

	fillMatrices1_unafold(	N, 
							99999,
							local_Qprime_ip,
							local_Q,
							local_Qprime,
							local_QM);

	// fillMatrices1_unafold(int N, 
	// 						int MAXLOOP, 
	// 						int* _local_Qprime_ip, 
	// 						int* _local_Q, 
	// 						int* _local_Qprime, 
	// 						int* _local_QM);


	return 0;
}