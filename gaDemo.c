#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h> 
#include<stdbool.h>

#define populationSize 100
#define chromosomeSize 10
#define maxGeneration 500
#define crossRate 0.8
#define mutationRate 0.05
float LB[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; //lower bound
float UB[10] = {5, 4, 5, 4, 5, 5, 5, 5, 5, 4}; //upper bound
float **a;
float **aa;
float **aaa;
int aRow;
int aaaRow;


float **createMatrix(int rows, int cols) {
    float **matrix = malloc(rows * sizeof(float *));
    for (int i = 0; i < rows; ++i) {
        *(matrix + i) = malloc(cols * sizeof(float));
    }
    return matrix;
}

void freeMatrix(float **matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        free(*(matrix + i));
    }
    free(matrix);
}


//Get data from files
bool getData(const char *fileName, float **x, int rows, int cols) {
    // open file to read data
    FILE *fp;
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Open file %s error!!\n", fileName);
        return false;
    }

    // read data
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fscanf(fp, "%f", &*(*(x + i) + j));
        }
    }
    return true;
}

//Fitness function
// float fitnessFcn(float *x){
// 	float x1 = *x;
// 	float x2 = *(x + 1);
// 	float x3 = *(x + 2);
// 	float x4 = *(x + 3);
// 	float x5 = *(x + 4);
// 	float x6 = *(x + 5);
// 	float x7 = *(x + 6);
// 	float x8 = *(x + 7);
// 	float x9 = *(x + 8);
// 	float x10 = *(x + 9);

// 	float y=pow((Dysum1-Dzsum1),2)+pow((Dysum2-Dzsum2),2)+pow((Dysum3-Dzsum3),2)+pow((Dysum4-Dzsum4),2)+pow((Dysum5-Dzsum5),2)+pow((Dysum6-Dzsum6),2)+pow((Dysum7-Dzsum7),2)+pow((Dysum8-Dzsum8),2)+pow((Dysum9-Dzsum9),2);
// }

//sig2ext in rainflow
float *sig2ext(float *sigy, float *dty, int rows, int *lenOfArray) {

    //w1=diff(sig);
    //w2=logical([1;(w1(1:end-1).*w1(2:end))<=0;1]);
    //tmp1=(w1(1:end-1).*w1(2:end))<=0
    float w1[rows - 1], tmp1[rows - 2], w2[rows];
    int n = 0, m = 0, k = 0, l = 0;

    w2[0] = 1;
    w2[rows - 1] = 1;

    for (int i = 1; i < rows; i++) {
        w1[i - 1] = sigy[i] - sigy[i - 1];
    }
    for (int i = 1; i < rows - 1; i++) {
        tmp1[i - 1] = w1[i - 1] * w1[i];
        if (tmp1[i - 1] <= 0) {
            tmp1[i - 1] = 1;
        } else {
            tmp1[i - 1] = 0;
            n++;
        }
        w2[i] = tmp1[i - 1];
    }

    //ext1=sigy(w2);	exttime1=dty(w2);
    float ext1[rows - n], exttime1[rows - n];
    for (int i = 0, j = 0; i < rows - n && j < rows;) {
        if (w2[j] == 0) {
            j++;
        } else {
            ext1[i] = sigy[j];
            exttime1[i] = dty[j];
            i++;
            j++;
        }
    }

    //w3=diff(ext1); tmp2 = w3(1:end-1)==0 & w3(2:end)==0; w4=~logical([0; w3(1:end-1)==0 & w3(2:end)==0; 0]);
    float w3[rows - n - 1], tmp2[rows - n - 2], w4[rows - n];

    w4[0] = 1;
    w4[rows - n - 1] = 1;

    for (int i = 1; i < rows - n; i++) {
        w3[i - 1] = ext1[i] - ext1[i - 1];
    }
    for (int i = 1; i < rows - n - 1; i++) {
        if (w3[i - 1] == 0 && w3[i] == 0) {
            tmp2[i - 1] = 1;
            w4[i] = 0;
            m++;
        } else {
            tmp2[i - 1] = 0;
            w4[i] = 1;
        }
    }

    //ext2=ext1(w4); exttime2=exttime1(w4);
    float ext2[rows - n - m], exttime2[rows - n - m];
    for (int i = 0, j = 0; i < rows - n - m && j < rows - n;) {
        if (w4[j] == 0) {
            j++;
        } else {
            ext2[i] = ext1[j];
            exttime2[i] = exttime1[j];
            i++;
            j++;
        }
    }

    //w5=~logical([0; ext2(1:end-1)==ext2(2:end)]);
    float w5[rows - n - m];
    w5[0] = 1;
    for (int i = 1; i < rows - n - m; i++) {
        if (ext2[i - 1] == ext2[i]) {
            w5[i] = 0;
            k++;
        } else {
            w5[i] = 1;
        }
    }

    //ext3=ext2(w5);
    float ext3[rows - n - m - k];
    for (int i = 0, j = 0; i < rows - n - m - k && j < rows - n - m;) {
        if (w5[j] == 0) {
            j++;
        } else {
            ext3[i] = ext2[j];
            i++;
            j++;
        }
    }
    //w6=(exttime2(2:end)-exttime2(1:end-1))./2
    //exttime3=[exttime2(1:end-1)+w6.*~w5(2:end); exttime2(end)];
    //exttime4=exttime3(w5);
    float w6[rows - n - m - 1], exttime3[rows - n - m], exttime4[rows - n - m - k];
    exttime3[rows - n - m - 1] = exttime2[rows - n - m - 1];

    for (int i = 1; i < rows - n - m; i++) {
        w6[i - 1] = (exttime2[i] - exttime2[i - 1]) / 2.00;
    }

    for (int i = 0, j = 1; i < rows - n - m - 1 && j < rows - n - m;) {
        if (w5[j] == 0) {
            exttime3[i] = w6[i] * 1.00 + exttime2[i];
            i++;
            j++;
        } else {
            exttime3[i] = w6[i] * 0.00 + exttime2[i];
            i++;
            j++;
        }
    }

    for (int i = 0, j = 0; i < rows - n - m - k && j < rows - n - m;) {
        if (w5[j] == 0) {
            j++;
        } else {
            exttime4[i] = exttime3[j];
            i++;
            j++;
        }
    }

    //length(ext3)>2,  w7=diff(ext3); w8=logical([1; w7(1:end-1).*w7(2:end)<0; 1]); ext4=ext3(w8); exttime5=exttime4(w8);
    float *ext4 = NULL, *exttime5 = NULL;
    if (rows - n - m - k > 2) {
        float w7[rows - n - m - k - 1], w8[rows - n - m - k];
        w8[0] = 1;
        w8[rows - n - m - k - 1] = 1;
        for (int i = 1; i < rows - n - m - k; i++) {
            w7[i - 1] = ext3[i] - ext3[i - 1];
        }
        for (int i = 1; i < rows - n - m - k - 1; i++) {
            if (w7[i - 1] * w7[i] < 0) {
                w8[i] = 1;
            } else {
                w8[i] = 0;
                l++;
            }
        }

        *lenOfArray = rows - n - m - k - l;
        ext4 = malloc(2 * (*lenOfArray) * sizeof(float));

        for (int i = 0, j = 0; i < rows - n - m - k - l && j < rows - n - m - k;) {
            if (w8[j] == 0) {
                j++;
            } else {
                ext4[i] = ext3[j];
                // exttime5[i] = exttime4[j];
                ext4[i + (*lenOfArray)] = exttime4[j];
                i++;
                j++;
            }
        }
    }
    return ext4;
}

//rainFlow in rainflow
float **rainFlow(float *ext, float *exttime, int lenOfSig2ext, int *lenOfRainflow) {
    float **rfy = NULL, **rfyResult = NULL;
    float a[16384], t[16384], ampl, mean, period, atime;
    int cNr = 1;
    int j = -1;

    //create 2D rfy(5 * (lenOfSig2ext -1))
    rfy = (float **) malloc(5 * sizeof(float *));
    for (int i = 0; i < 5; i++) {
        rfy[i] = (float *) malloc((lenOfSig2ext - 1) * sizeof(float));
    }

    int columnId = 0;
    int pointId = 0;

    for (int i = 0; i < lenOfSig2ext; i++) {
        a[++j] = *(ext + pointId);
        //printf("%f\n", *(ext + pointId));
        t[j] = *(exttime + pointId);
        //printf("%f\n", *(exttime + pointId));
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
                        *(*(rfy + 0) + columnId) = ampl;
                        *(*(rfy + 1) + columnId) = mean;
                        *(*(rfy + 2) + columnId) = 0.50;
                        *(*(rfy + 3) + columnId) = atime;
                        *(*(rfy + 4) + columnId) = period;
                        columnId++;
                    }
//                    printf("%f\n", *(*(rfy + 0) + columnId));
//                    printf("%f\n", *(*(rfy + 1) + columnId));
//                    printf("%f\n", *(*(rfy + 2) + columnId));
//                    printf("%f\n", *(*(rfy + 3) + columnId));
//                    printf("%f\n", *(*(rfy + 4) + columnId));
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
                        *(*(rfy + 0) + columnId) = ampl;
                        *(*(rfy + 1) + columnId) = mean;
                        *(*(rfy + 2) + columnId) = 1.00;
                        *(*(rfy + 3) + columnId) = atime;
                        *(*(rfy + 4) + columnId) = period;
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
            *(*(rfy + 0) + columnId) = ampl;
            *(*(rfy + 1) + columnId) = mean;
            *(*(rfy + 2) + columnId) = 0.50;
            *(*(rfy + 3) + columnId) = atime;
            *(*(rfy + 4) + columnId) = period;
            columnId++;
        }
    }

    //create 2D rfyResult(5 * (lenOfSig2ext - cNr))
    rfyResult = (float **) malloc(5 * sizeof(float *));
    for (int i = 0; i < 5; i++) {
        rfyResult[i] = (float *) malloc((lenOfSig2ext - cNr) * sizeof(float));
    }
    *lenOfRainflow = lenOfSig2ext - cNr;

    for (int i = 0; i < 5; i++) {
        for (int k = 0; k < lenOfSig2ext - cNr; k++) {
            rfyResult[i][k] = rfy[i][k];
        }
    }
    return rfyResult;
}

//hist in rfhist
float *hist(float *rf, int lenOfRfhist, int x){
    float *noy, *xoy;

    xoy = malloc(x * sizeof(float));
    noy = malloc(x * sizeof(float));
    memset(noy, 0, sizeof(noy));

    float min = rf[0], max = rf[0];
    for (int i = 0; i < lenOfRfhist; i++) {
        if (rf[i] >= max) {
            max = rf[i];
        } else if (rf[i] <= min) {
            min = rf[i];
        }
    }
    float wid = (max - min) / x;
    for (int i = 0; i < x; i++) {
        xoy[i] = min + (float) (i + 0.50) * wid;
    }

    for (int i = 0; i < lenOfRfhist; i++) {
        int j;
        j = (int) floor((rf[i] - min) / wid);
        if (j == x) {
            noy[j - 1] += 1.00;
        } else {
            noy[j] += 1.00;
        }
    }

    float *hist = malloc(2 * x * sizeof(float));
    for(int i = 0; i < x; i++){
        hist[i] = noy[i];
        hist[i + x] = xoy[i];
    }

    return hist;
}

//rfhist in rainflow
float *rfhist(float **rfy, int lenOfRainflow, int *lenOfRfhist) {
    float *noy = NULL, *xoy = NULL;
    int x = 32;
    *lenOfRfhist = x;

    //halfc=find(rfy(3,:)==0.5)
    int halfcNum = 0;
    for (int i = 0; i < lenOfRainflow; i++) {
        if (rfy[2][i] == 0.50)
            halfcNum++;
    }

    int *halfc = NULL;
    halfc = malloc(halfcNum * sizeof(int));
    for (int i = 0, j = 0; i < lenOfRainflow && j < halfcNum;) {
        if (rfy[2][i] == 0.50) {
            halfc[j] = i;
            j++;
        }
        i++;
    }

    //[N1 x]=hist(rfy(r,:),x)     r = 1, x = *xoy, N1 = *noy
    noy = hist(*rfy, lenOfRainflow, x);
    xoy = noy + x;

    //if ~isempty(halfc) {
    //[N2 x]=hist(rf(r,halfc),x)  N2 = noy2, x = *xoy
    //N1=N1-0.5*N2  N1 = noy
    // }
    if(halfcNum != 0){
        float *noy2 = NULL;
        float *rf = malloc(halfcNum * sizeof(float));
        for(int i = 0; i < halfcNum; i++){
            int j = halfc[i];
            rf[i] = rfy[0][j];
        }
        noy2 = hist(rf, halfcNum, x);
        for(int i = 0; i < x; i++){
            noy[i] -= noy2[i] * 0.5;
        }
    }

    float *rfhist = malloc(2 * x * sizeof(float));
    for(int i = 0; i < x; i++){
        rfhist[i] = noy[i];
        rfhist[i + x] = xoy[i];
    }

    return rfhist;
}

//test preData
void testPreData() {
    int kk = 0;
    float sigy[aaaRow], dty[aaaRow];
    for (int i = 0; i < aaaRow; i++) {
        sigy[i] = aaa[i][kk + 2];
        dty[i] = aaa[i][1];
    }

    float *ext = NULL, *exttime = NULL;
    int lenOfSig2ext;
    ext = sig2ext(sigy, dty, aaaRow, &lenOfSig2ext);
    exttime = ext + lenOfSig2ext;

    float **rfy = NULL;
    int lenOfRainflow;
    rfy = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);


    float *noy = NULL, *xoy = NULL;
    int lenOfRfhist;
    noy = rfhist(rfy, lenOfRainflow, &lenOfRfhist);
    xoy = noy + lenOfRfhist;
}


int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Error!\n");
        return 0;
    }

    bool success = true;

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



    // for(int i = 0; i < aRow; i++){
    // 	for(int j = 0; j < 16; j++){
    // 		printf("%f ", a[i][j]);
    // 	}
    // 	printf("\n");
    // }
    // for(int i = 0; i < 10; i++){
    // 	for(int j = 0; j < 9; j++){
    // 		printf("%f ", aa[i][j]);
    // 	}
    // 	printf("\n");
    // }
    // for(int i = 0; i < aaaRow; i++){
    // 	for(int j = 0; j < 11; j++){
    // 		printf("%f ", aaa[i][j]);
    // 	}
    // 	printf("\n");
    // }

    freeMatrix(a, aRow);
    freeMatrix(aa, 10);
    freeMatrix(aaa, aaaRow);


    return 0;
}