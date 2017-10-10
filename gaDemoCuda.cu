#include "cuda_runtime.h"
#include "device_launch_parameters.h" 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

#define BOOL int
#define TRUE 1
#define FALSE 0
#define populationSize 100
#define chromosomeSize 10
#define maxGeneration 500
#define crossRate 0.8
#define mutationRate 0.01
#define eliteCount 0.05*populationSize


//typedef float float;
float LB[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; //lower bound
float UB[10] = {5, 4, 5, 4, 5, 5, 5, 5, 5, 4}; //upper bound
float *a;  //Tzaihe
float *aa;  //yingliK
float *aaa; //Tyingli
int aRow;
int aaaRow;
float Dysum[9];

__device__ float c_LB[10]; //lower bound
__device__ float c_UB[10]; //upper bound
__device__ float *c_a;  //Tzaihe
__device__ float *c_aa;  //yingliK
__device__ float *c_aaa; //Tyingli
__device__ int c_aRow;
__device__ int c_aaaRow;
__device__ float c_Dysum[9];
// float Dzsum[9];
// float populationArray[populationSize][chromosomeSize];  //种群数组
// float fitness[populationSize]; //每个种群的适应度
// float populationPro[populationSize]; //每个种群在select时被选中的概率
float bestFitnessOfGen; //每一代的最优适应度
int bestIndexOfGen; //每一代的最优适应度位置
float aveFitnessOfGen[maxGeneration]; //每一代的平均最优适应度
float X_10[chromosomeSize]; //最优适应度对应的x值
float fval; //最终最优适应度
int G; //取得最终最优适应度的迭代次数
// BOOL elitism = TRUE; //是否精英选择

float *createMatrix(int rows, int cols) {
    float *matrix = (float*)malloc(rows * cols * sizeof(float));
    return matrix;
}

//Get data from files
BOOL getData(const char *fileName, float *x, int rows, int cols) {
    // open file to read data
    FILE *fp;
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Open file %s error!!\n", fileName);
        return FALSE;
    }

    // read data
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fscanf(fp, "%f", (x + i * cols + j));
        }
    }
    return TRUE;
}

// //sig2ext in rainflow
__host__ __device__ float *sig2ext(float *sigy, float *dty, int rows, int *lenOfArray) {

    int n = 0, m = 0, k = 0, l = 0;

    //w=diff(sig);
    //w=logical([1;(w(1:end-1).*w(2:end))<=0;1]);
    float *w = (float*)malloc(rows * sizeof(float));
    for (int i = 1; i < rows; i++) {
        w[i - 1] = sigy[i] - sigy[i - 1];
    }

    for (int i = rows - 2; i > 0; i--) {
        float tmp = w[i] * w[i - 1];
        if (tmp <= 0) {
            w[i] = 1;
        } else {
            w[i] = 0;
            n++;
        }
    }
    w[0] = w[rows - 1] = 1;

    //ext=sigy(w);  exttime=dty(w);
    float *ext = (float*)malloc((rows - n) * sizeof(float));
    float *exttime = (float*)malloc((rows - n) * sizeof(float));
    for (int i = 0, j = 0; i < rows - n && j < rows;) {
        if (w[j] == 0) {
            j++;
        } else {
            ext[i] = sigy[j];
            exttime[i] = dty[j];
            i++;
            j++;
        }
    }

    //w=diff(ext);
    //w=~logical([0; w(1:end-1)==0 & w(2:end)==0; 0]);
    for (int i = 1; i < rows - n; i++) {
        w[i - 1] = ext[i] - ext[i - 1];
    }

    for (int i = rows - n - 2; i > 0; i--) {
        if (w[i - 1] == 0 && w[i] == 0) {
            w[i] = 0;
            m++;
        } else {
            w[i] = 1;
        }
    }
    w[0] = w[rows - n - 1] = 1;

    //ext=ext(w); exttime=exttime(w);
    for (int i = 0, j = 0; i < rows - n - m && j < rows - n;) {
        if (w[j] == 0) {
            j++;
        } else {
            ext[i] = ext[j];
            exttime[i] = exttime[j];
            i++;
            j++;
        }
    }

    //w=~logical([0; ext(1:end-1)==ext(2:end)]);
    for (int i = 1; i < rows - n - m; i++) {
        if (ext[i - 1] == ext[i]) {
            w[i] = 0;
            k++;
        } else {
            w[i] = 1;
        }
    }
    w[0] = 1;

    //ext=ext(w);
    for (int i = 0, j = 0; i < rows - n - m - k && j < rows - n - m;) {
        if (w[j] == 0) {
            j++;
        } else {
            ext[i] = ext[j];
            i++;
            j++;
        }
    }
    //w2=(exttime(2:end)-exttime(1:end-1))./2
    //exttime=[exttime(1:end-1)+w2.*~w(2:end); exttime(end)];
    //exttime=exttime(w);

    float *w2 = (float*)malloc((rows - n - m - 1) * sizeof(float));
    for (int i = 1; i < rows - n - m; i++) {
        w2[i - 1] = (exttime[i] - exttime[i - 1]) / 2.00;
    }

    for (int i = 0, j = 1; i < rows - n - m - 1 && j < rows - n - m;) {
        if (w[j] == 0) {
            exttime[i] = w2[i] * 1.00 + exttime[i];
            i++;
            j++;
        } else {
            exttime[i] = w2[i] * 0.00 + exttime[i];
            i++;
            j++;
        }
    }

    for (int i = 0, j = 0; i < rows - n - m - k && j < rows - n - m;) {
        if (w[j] == 0) {
            j++;
        } else {
            exttime[i] = exttime[j];
            i++;
            j++;
        }
    }

    //length(ext)>2,  w=diff(ext); w=logical([1; w(1:end-1).*w(2:end)<0; 1]);
    //ext4=ext(w); exttime=exttime(w);
    float *ext4 = NULL;
    *lenOfArray = 0;
    if (rows - n - m - k > 2) {
        for (int i = 1; i < rows - n - m - k; i++) {
            w[i - 1] = ext[i] - ext[i - 1];
        }

        for (int i = rows - n - m - k - 2; i > 0; i--) {
            if (w[i - 1] * w[i] < 0) {
                w[i] = 1;
            } else {
                w[i] = 0;
                l++;
            }
        }
        w[0] = 1;
        w[rows - n - m - k - 1] = 1;

        *lenOfArray = rows - n - m - k - l;
        ext4 = (float*)malloc(2 * (*lenOfArray) * sizeof(float));

        for (int i = 0, j = 0; i < rows - n - m - k - l && j < rows - n - m - k;) {
            if (w[j] == 0) {
                j++;
            } else {
                ext4[i] = ext[j];
                ext4[i + (*lenOfArray)] = exttime[j];
                i++;
                j++;
            }
        }
    }

    free(w);
    free(w2);
    free(ext);
    free(exttime);
    return ext4;
}

// //rainFlow in rainflow
__host__ __device__ float *rainFlow(float *ext, float *exttime, int lenOfSig2ext, int *lenOfRainflow) {
    float *rfy = NULL, *rfyResult = NULL;

    //function rfy5
    float a[100], t[100], ampl, mean, period, atime;
    int cNr = 1;
    int j = -1;

    //create 2D rfy(5 * (lenOfSig2ext -1))
    rfy = (float*)malloc(5 * (lenOfSig2ext - 1) * sizeof(float));

    int columnId = 0;
    int pointId = 0;

    for (int i = 0; i < lenOfSig2ext; i++) {
        a[++j] = *(ext + pointId);
        t[j] = *(exttime + pointId);
        while ((j >= 2) && (fabs(a[j - 1] - a[j - 2]) <= fabs(a[j] - a[j - 1]))) {
            ampl = fabs((a[j - 1] - a[j - 2]) / 2);
            switch (j) {
                case 0: {
                    break;
                }
                case 1: {
                    break;
                }
                case 2: {
                    mean = (a[0] + a[1]) / 2;
                    period = (t[1] - t[0]) * 2;
                    atime = t[0];
                    a[0] = a[1];
                    a[1] = a[2];
                    t[0] = t[1];
                    t[1] = t[2];
                    j = 1;
                    if (ampl > 0) {
                        *(rfy + columnId*5 + 0) = ampl;
                        *(rfy + columnId*5 + 1) = mean;
                        *(rfy + columnId*5 + 2) = 0.50;
                        *(rfy + columnId*5 + 3) = atime;
                        *(rfy + columnId*5 + 4) = period;
                        columnId++;
                    }
                    break;
                }
                default: {
                    mean = (a[j - 1] + a[j - 2]) / 2;
                    period = (t[j - 1] - t[j - 2]) * 2;
                    atime = t[j - 2];
                    a[j - 2] = a[j];
                    t[j - 2] = t[j];
                    j = j - 2;
                    if (ampl > 0) {
                        *(rfy + columnId*5 + 0) = ampl;
                        *(rfy + columnId*5 + 1) = mean;
                        *(rfy + columnId*5 + 2) = 1.00;
                        *(rfy + columnId*5 + 3) = atime;
                        *(rfy + columnId*5 + 4) = period;
                        columnId++;
                        cNr++;
                    }
                    break;
                }
            }
        }
        pointId++;
    }
    for (int i = 0; i < j; i++) {
        ampl = fabs(a[i] - a[i + 1]) / 2;
        mean = (a[i] + a[i + 1]) / 2;
        period = (t[i + 1] - t[i]) * 2;
        atime = t[i];
        if (ampl > 0) {
            *(rfy + columnId*5 + 0) = ampl;
            *(rfy + columnId*5 + 1) = mean;
            *(rfy + columnId*5 + 2) = 0.50;
            *(rfy + columnId*5 + 3) = atime;
            *(rfy + columnId*5 + 4) = period;
            columnId++;
        }
    }

    //create 2D rfyResult(5 * (lenOfSig2ext - cNr))
    rfyResult = (float*)malloc(5 * (lenOfSig2ext - cNr) * sizeof(float));

    *lenOfRainflow = lenOfSig2ext - cNr;

    for (int i = 0; i < 5 * (lenOfSig2ext - cNr); i++) {
            rfyResult[i] = rfy[i];
    }

    free(rfy);
    return rfyResult;
}

// //rfhist in rainflow
__host__ __device__ float *rfhist(float *rfy, int lenOfRainflow, int *lenOfRfhist) {
    float *noy = NULL, *xoy = NULL;
    int x = 32;
    *lenOfRfhist = x;

    xoy = (float*)malloc(x * sizeof(float));
    noy = (float*)malloc(x * sizeof(float));
    memset(noy, 0, x * sizeof(float));

    //halfc=find(rfy(3,:)==0.5)
    int halfcNum = 0;
    for (int i = 0; i < lenOfRainflow; i++) {
        if (rfy[i * 5 + 2] == 0.50)
            halfcNum++;
    }

    int *halfc = NULL;
    halfc = (int*)malloc(halfcNum * sizeof(int));
    for (int i = 0, j = 0; i < lenOfRainflow && j < halfcNum;) {
        if (rfy[i * 5 + 2] == 0.50) {
            halfc[j] = i;
            j++;
        }
        i++;
    }

    float min = rfy[0], max = rfy[0];
    for (int i = 0; i < lenOfRainflow; i++) {
        if (rfy[i * 5] >= max) {
            max = rfy[i * 5];
        } else if (rfy[i * 5] <= min) {
            min = rfy[i * 5];
        }
    }

    float wid = (max - min) / x;
    for (int i = 0; i < x; i++) {
        xoy[i] = min + (float) (i + 0.50) * wid;
    }

    for (int i = 0; i < lenOfRainflow; i++) {
        int j;
        j = (int) floor((rfy[i * 5] - min) / wid);
        if (j != 0 && fabs((rfy[i * 5] - min) - wid * j) < 1e-10) {
            noy[j - 1] += 1;
        } else {
            noy[j] += 1;
        }
    }

    //if ~isempty(halfc) {
    //[N2 x]=hist(rf(r,halfc),x)  N2 = noy2, x = *xoy
    //N1=N1-0.5*N2  N1 = noy
    // }
    if (halfcNum != 0) {
        float *noy2 = (float*)malloc(x * sizeof(float));
        memset(noy2, 0, x * sizeof(float));
        float *rf = (float*)malloc(halfcNum * sizeof(float));
        for (int i = 0; i < halfcNum; i++) {
            int j = halfc[i];
            rf[i] = rfy[j * 5];
        }

        for (int i = 0; i < halfcNum; i++) {
            int j;
            j = (int) floor((rf[i] - min) / wid);
            if (j != 0 && fabs((rf[i] - min) - wid * j) < 1e-10) {
                noy2[j - 1] += 1;
            } else {
                noy2[j] += 1;
            }
        }

        for (int i = 0; i < x; i++) {
            noy[i] -= noy2[i] * 0.5;
        }

        free(noy2);
        free(rf);
    }

    float *rfhist = (float*)malloc(2 * x * sizeof(float));
    for (int i = 0; i < x; i++) {
        rfhist[i] = noy[i];
        rfhist[i + x] = xoy[i];
    }

    free(halfc);
    free(xoy);
    free(noy);

    return rfhist;
}

// __global__ void c_testPreData(float *c_aaa, int aaaRow, int aaaCol, float *c_Dysum) {
//     memset(c_Dysum, 0, sizeof(float) * 9);
    
//     int idx = threadIdx.x;
//     float *sigy = (float*)malloc(aaaRow * sizeof(float));
//     float *dty = (float*)malloc(aaaRow * sizeof(float));

//     for (int i = 0; i < aaaRow; i++) {
//         sigy[i] = c_aaa[i * aaaCol + idx + 2];
//         dty[i] = c_aaa[i * aaaCol + 1];
//     }
//     __syncthreads();

//     float *ext = NULL, *exttime = NULL;
//     int lenOfSig2ext;
//     ext = sig2ext(sigy, dty, aaaRow, &lenOfSig2ext);
//     //__syncthreads();
//     exttime = ext + lenOfSig2ext;
//     __syncthreads();

//     float *rfy = NULL;
//     int lenOfRainflow;
//     rfy = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);
//     __syncthreads();

//     float *noy = NULL, *xoy = NULL;
//     int lenOfRfhist;
//     noy = rfhist(rfy, lenOfRainflow, &lenOfRfhist);
//     //__syncthreads();
//     xoy = noy + lenOfRfhist;
//     __syncthreads();

//     for (int i = 0; i < lenOfRfhist; i++) {
//         c_Dysum[idx] += noy[i] * pow(xoy[i] * 0.21 / 70, 3.5);
//     }
//     __syncthreads();
//     printf("%e\n", c_Dysum[idx]);

//     free(sigy);
//     free(dty);
//     free(ext);
//     free(rfy);
//     free(noy);
// }

void testPreData() {
    for (int kk = 0; kk < 9; kk++) {
        float *sigy = (float*)malloc(aaaRow * sizeof(float));
        float *dty = (float*)malloc(aaaRow * sizeof(float));

        for (int i = 0; i < aaaRow; i++) {
            sigy[i] = aaa[i * 11 + kk + 2];
            dty[i] = aaa[i * 11 + 1];
        }

        float *ext = NULL, *exttime = NULL;
        int lenOfSig2ext;
        ext = sig2ext(sigy, dty, aaaRow, &lenOfSig2ext);
        exttime = ext + lenOfSig2ext;

        float *rfy = NULL;
        int lenOfRainflow;
        rfy = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);

        float *noy = NULL, *xoy = NULL;
        int lenOfRfhist;
        noy = rfhist(rfy, lenOfRainflow, &lenOfRfhist);
        xoy = noy + lenOfRfhist;

        for (int i = 0; i < lenOfRfhist; i++) {
            Dysum[kk] += noy[i] * pow(xoy[i] * 0.21 / 70, 3.5);
        }
        //printf("%e\n", Dysum[kk]);

        free(sigy);
        free(dty);
        free(ext);
        free(rfy);
        free(noy);
    }
}


// __device__ float sum(float *x) {
//     float sum = 0;
//     for (int i = 0; i < populationSize; i++) {
//         sum += x[i];
//     }
//     return sum;
// }

 //fitness Function
float HfitnessFcn(float *x) {
    //initial Dzsum in every generation
    float *Dzsum = (float*)malloc(9 * sizeof(float));
    memset(Dzsum, 0, sizeof(float) * 9);

    float *Tzb = (float*)malloc(aRow * 9 * sizeof(float));

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < aRow; j++) {
            Tzb[j * 9 + i] = x[0] * aa[0 * 9 + i] * a[j * 16 + 2] + x[1] * aa[1 * 9 + i] * a[j * 16 + 3] + x[2] * aa[2 *0 + i] * a[j * 16 + 4] +
                            x[3] * aa[3 * 9 + i] * a[j * 16 + 5] + x[4] * aa[4 * 9 + i] * a[j * 16 + 6] + x[5] * aa[5 * 9 + i] * a[j * 16 + 7] +
                            x[6] * aa[6 * 9 + i] * a[j * 16 + 8] + x[7] * aa[7 * 9 + i] * a[j * 16 + 9] + x[8] * aa[8 * 9 + i] * a[j * 16 + 10] +
                            x[9] * aa[9 * 9 + i] * a[j * 16 + 11];
        }
    }

    for (int k = 0; k < 9; k++) {
        float *sig = (float*)malloc(aRow * sizeof(float));
        float *dt = (float*)malloc(aRow * sizeof(float));
        for (int i = 0; i < aRow; i++) {
            sig[i] = Tzb[i * 9 + k];
            dt[i] = a[i * 9 + 1];
        }

        float *ext = NULL, *exttime = NULL;
        int lenOfSig2ext;
        ext = sig2ext(sig, dt, aRow, &lenOfSig2ext);
        exttime = ext + lenOfSig2ext;

        float *rf = NULL;
        int lenOfRainflow;
        rf = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);

        float *no = NULL, *xo = NULL;
        int lenOfRfhist;
        no = rfhist(rf, lenOfRainflow, &lenOfRfhist);
        xo = no + lenOfRfhist;

        for (int i = 0; i < lenOfRfhist; i++) {
            Dzsum[k] += no[i] * pow(xo[i] * 0.21 / 70, 3.5);
        }
        //printf("%e\n", Dzsum[k]);

        free(sig);
        free(dt);
        free(ext);
        free(rf);
        free(no);
    }

    float y = 0;
    for (int i = 0; i < 9; i++) {
        //constraint : c =(Dysum[i]-Dzsum[i]) <= 0
        float c = Dysum[i] - Dzsum[i];
        if (c <= 0) {
            y += pow(c, 2);
        } else {
            y = 100;
        }
    }
    //printf("%e\n", y);

    free(Dzsum);
    free(Tzb);

    return y;
}

__device__ float DfitnessFcn(float *x) {
    //initial Dzsum in every generation
    float *Dzsum = (float*)malloc(9 * sizeof(float));
    memset(Dzsum, 0, sizeof(float) * 9);

    float *Tzb = (float*)malloc(c_aRow * 9 * sizeof(float));

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < c_aRow; j++) {
            Tzb[j * 9 + i] = x[0] * c_aa[0 * 9 + i] * c_a[j * 16 + 2] + x[1] * c_aa[1 * 9 + i] * c_a[j * 16 + 3] + x[2] * c_aa[2 *0 + i] * c_a[j * 16 + 4] +
                            x[3] * c_aa[3 * 9 + i] * c_a[j * 16 + 5] + x[4] * c_aa[4 * 9 + i] * c_a[j * 16 + 6] + x[5] * c_aa[5 * 9 + i] * c_a[j * 16 + 7] +
                            x[6] * c_aa[6 * 9 + i] * c_a[j * 16 + 8] + x[7] * c_aa[7 * 9 + i] * c_a[j * 16 + 9] + x[8] * c_aa[8 * 9 + i] * c_a[j * 16 + 10] +
                            x[9] * c_aa[9 * 9 + i] * c_a[j * 16 + 11];
        }
    }

    for (int k = 0; k < 9; k++) {
        float *sig = (float*)malloc(c_aRow * sizeof(float));
        float *dt = (float*)malloc(c_aRow * sizeof(float));
        for (int i = 0; i < c_aRow; i++) {
            sig[i] = Tzb[i * 9 + k];
            dt[i] = c_a[i * 9 + 1];
        }

        float *ext = NULL, *exttime = NULL;
        int lenOfSig2ext;
        ext = sig2ext(sig, dt, c_aRow, &lenOfSig2ext);
        exttime = ext + lenOfSig2ext;

        float *rf = NULL;
        int lenOfRainflow;
        rf = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);

        float *no = NULL, *xo = NULL;
        int lenOfRfhist;
        no = rfhist(rf, lenOfRainflow, &lenOfRfhist);
        xo = no + lenOfRfhist;

        for (int i = 0; i < lenOfRfhist; i++) {
            Dzsum[k] += no[i] * pow(xo[i] * 0.21 / 70, 3.5);
        }
        printf("%e\n", Dzsum[k]);

        free(sig);
        free(dt);
        free(ext);
        free(rf);
        free(no);
    }

    float y = 0;
    for (int i = 0; i < 9; i++) {
        //constraint : c =(Dysum[i]-Dzsum[i]) <= 0
        float c = c_Dysum[i] - Dzsum[i];
        if (c <= 0) {
            y += pow(c, 2);
        } else {
            y = 100;
        }
    }
    printf("%e\n", y);

    free(Dzsum);
    //free(Tzb);

    return y;
}


__global__ void GfitnessFcn(float *populationArray, float *fitness){
    int idx = threadIdx.x;
    float *x = (float*)malloc(chromosomeSize * sizeof(float));
    memset(x, 0, 10 * sizeof(float));
    for (int j = 0; j < chromosomeSize; j++) {
        x[j] = populationArray[idx * chromosomeSize + j];
        //printf("%f\n", x[j]);
    }
    __syncthreads();
    fitness[idx] = DfitnessFcn(Tzb);
    free(x);
}



void initial(float *populationArray){
    for (int i = 0; i < populationSize; i++) {
        float *x = (float*)malloc(chromosomeSize * sizeof(float));
        for (int j = 0; j < chromosomeSize; j++) {
            int high_pos = rand();
            int low_pos = (rand() & ((1 << 16) - 1));
            high_pos = (high_pos & ((1 << 15) - 1));
            int value = low_pos + (high_pos << 16);
            populationArray[i * chromosomeSize + j] = (UB[j] - LB[j]) * ((float) value / ((1U << 31) - 1)) + LB[j];
            x[j] = populationArray[i * chromosomeSize + j];
        }
        float tmp_fit = HfitnessFcn(x);
        if (tmp_fit > 99) {
            i--;
        }
        free(x);
    }
}

// __device__ void *bestFitness() {

//     //bestRes[bestFitness][bestIndex]
//     float c_bestFitness = c_fitness[0];
//     int c_bestIndex = 0;
//     float *bestRes = (float*)malloc(2 * sizeof(float));
//     for (int i = 0; i < populationSize; i++) {
//         if (c_fitness[i] < c_bestFitness) {
//             c_bestFitness = c_fitness[i];
//             c_bestIndex = i;
//         }
//     }
//     bestRes[0] = c_bestFitness;
//     bestRes[1] = c_bestIndex;

//     return bestRes;
// }

// //select function 轮盘选择
// __device__ void selectFcn() {

//     float tmpPopulationArray[populationSize][chromosomeSize];
//     float tmpFitness[populationSize];
//     //每个个体被选择的概率
//     float *Fitness = malloc(populationSize * sizeof(float));
//     float sumFitness = 0;

//     for (int i = 0; i < populationSize; i++) {
//         Fitness[i] = 1 / fitness[i];
//     }

//     sumFitness = sum(Fitness);
//     for (int i = 0; i < populationSize; i++) {
//         populationPro[i] = Fitness[i] / sumFitness;
//     }
//     free(Fitness);

//     //轮盘选择
//     int *index = malloc(populationSize * sizeof(int));
//     for (int i = 0; i < populationSize; i++) {
//         float pick = ((float) rand()) / RAND_MAX;
//         while (pick < 0.0001)
//             pick = ((float) rand()) / RAND_MAX;

//         for (int j = 0; j < populationSize; j++) {
//             pick -= populationPro[j];
//             if (pick <= 0) {
//                 index[i] = j;
//                 break;
//             }
//         }
//     }

//     //是否精英选择
//     int elitismSize = populationSize;
//     if (elitism == TRUE) {
//         int *rank;
//         rank = rankForElitism();
//         elitismSize = (int) (populationSize - eliteCount);

//         //在新种群的最后保留eliteCount个个体
//         for (int i = elitismSize, k = 0; i < populationSize && k < eliteCount; i++, k++) {
//             for (int j = 0; j < chromosomeSize; j++) {
//                 tmpPopulationArray[i][j] = populationArray[rank[k]][j];
//             }
//             tmpFitness[i] = fitness[rank[k]];
//         }
//     }
//     for (int i = 0; i < elitismSize; i++) {
//         for (int j = 0; j < chromosomeSize; j++) {
//             tmpPopulationArray[i][j] = populationArray[index[i]][j];
//         }
//         tmpFitness[i] = fitness[index[i]];
//     }
//     free(index);

//     //产生新种群
//     for (int i = 0; i < populationSize; i++) {
//         for (int j = 0; j < chromosomeSize; j++) {
//             populationArray[i][j] = tmpPopulationArray[i][j];
//         }
//         fitness[i] = tmpFitness[i];
//     }
// }

// //cross function 每两个个体做判断
// __device__ void crossFcn() {
//     for (int i = 0; i < populationSize; i += 2) {
//         //判断当前两个个体是否做交叉
//         float pick1 = ((float) rand()) / RAND_MAX;
//         if (pick1 > crossRate)
//             continue;

//         for (int j = 0; j < chromosomeSize; j++) {
//             //判断两个个体中的染色体是否做交叉
//             int pick2 = rand();
//             if (pick2 & 1) {
//                 float tmp = populationArray[i][j];
//                 populationArray[i][j] = populationArray[i + 1][j];
//                 populationArray[i + 1][j] = tmp;
//             }
//         }
//     }
// }

// //mutation function
//  __device__ void mutationFcn() {
//     float scale = 0.5, shrink = 0.75;
//     for (int i = 0; i < populationSize; i++) {
//         scale -= scale * shrink * i / maxGeneration;

//         //判断当前个体是否变异
//         float pick1 = ((float) rand()) / RAND_MAX;
//         if (pick1 > mutationRate)
//             continue;

//         for (int j = 0; j < chromosomeSize; j++) {
//             //判断当前染色体是否变异
//             int pick2 = rand();
//             if (pick2 & 1) {
//                 float tmpChromosome;
//                 do {
//                     float pick3 = ((float) rand()) / RAND_MAX * 2 - 1;
//                     tmpChromosome = populationArray[i][j] + scale * (UB[j] - LB[j]) * pick3;
//                     //判断是否越界
//                 } while (tmpChromosome > UB[j] || tmpChromosome < LB[j]);
//                 populationArray[i][j] = tmpChromosome;
//             }
//         }
//     }
// }




int main(int argc, char *argv[]){
	time_t start = clock();
	srand(time(NULL));

	if(argc != 6){
		printf("ERROR\n");
		return 0;
	}

    BOOL success = TRUE;

    aRow = atoi(argv[2]);
    a = createMatrix(aRow, 16);
    success = getData(argv[1], a, aRow, 16);
    if (!success) {
        return 0;
    }
    aa = createMatrix(10, 9);
    success = getData(argv[3], aa, 10, 9);
    if (!success) {
        return 0;
    }

    aaaRow = atoi(argv[5]);
    aaa = createMatrix(aaaRow, 11);
    success = getData(argv[4], aaa, aaaRow, 11);
    if (!success) {
        return 0;
    }
    testPreData();


    cudaMemcpyToSymbol(c_a, a, aRow * 16 * sizeof(float));
    cudaMemcpyToSymbol(c_aa, aa, 10 * 9 * sizeof(float));
    cudaMemcpyToSymbol(c_aaa, aaa, aaaRow * 11 * sizeof(float));
    cudaMemcpyToSymbol(c_aRow, &aRow, sizeof(float));
    cudaMemcpyToSymbol(c_aaaRow, &aaaRow, sizeof(float));
    cudaMemcpyToSymbol(c_LB, LB, 10 * sizeof(float));
    cudaMemcpyToSymbol(c_UB, UB, 10 * sizeof(float));
    cudaMemcpyToSymbol(c_Dysum, Dysum, 9 * sizeof(float));

    
    float *populationArray;
    float *fitness;
    float *populationPro;
    float *X_10;
    fval = 100;
    //BOOL elitism = TRUE;

    cudaMallocManaged(&populationArray, populationSize * chromosomeSize * sizeof(float));
    cudaMallocManaged(&fitness, populationSize * sizeof(float));
    cudaMallocManaged(&populationPro, populationSize * sizeof(float));
    cudaMallocManaged(&X_10, 10 * sizeof(float));


    //cudaMemset(Dzsum, 0, 9 * sizeof(float));
    cudaMemset(populationArray, 0, populationSize * chromosomeSize * sizeof(float));
    cudaMemset(fitness, 0, populationSize * sizeof(float));
    cudaMemset(populationPro, 0, populationSize * sizeof(float));
    cudaMemset(X_10, 0, 10 * sizeof(float));

    //initial population
    initial(populationArray);

    //fitness function
    GfitnessFcn<<<1, 100>>>(populationArray, fitness);
    cudaDeviceSynchronize();
    // for (int i = 0; i < populationSize; i++) {
    //     float *x = (float*)malloc(chromosomeSize * sizeof(float));
    //     for (int j = 0; j < chromosomeSize; j++) {
    //         x[j] = populationArray[i * chromosomeSize + j];
    //     }
    //     fitness[i] = HfitnessFcn(x);
    //     free(x);
    // }

    

    cudaFree(c_Dysum);
    cudaFree(c_LB);
    cudaFree(c_UB);
    cudaFree(populationArray);
    cudaFree(fitness);
    cudaFree(populationPro);
    cudaFree(X_10);

   
    time_t stop = clock();
    printf("time:%e\n", ((float) (stop - start)) / CLOCKS_PER_SEC);


    free(a);
    free(aa);
    free(aaa);

    return 0;
}