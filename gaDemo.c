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
#define eliteCount 0.05*populationSize

typedef double FLOAT;

FLOAT LB[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; //lower bound
FLOAT UB[10] = {5, 4, 5, 4, 5, 5, 5, 5, 5, 4}; //upper bound
FLOAT **a;
FLOAT **aa;
FLOAT **aaa;
int aRow;
int aaaRow;
FLOAT Dysum[9];
FLOAT Dzsum[9];
FLOAT populationArray[populationSize][chromosomeSize];  //种群数组
FLOAT fitness[populationSize]; //每个种群的适应度
FLOAT bestFitnessOfGen; //每一代的最优适应度
int bestIndexOfGen; //每一代的最优适应度位置
FLOAT aveFitnessOfGen[maxGeneration]; //每一代的平均最优适应度
FLOAT X_10[chromosomeSize]; //最优适应度对应的x值
FLOAT fval; //最终最优适应度
int G; //取得最终最优适应度的迭代次数





FLOAT **createMatrix(int rows, int cols) {
    FLOAT **matrix = malloc(rows * sizeof(FLOAT *));
    for (int i = 0; i < rows; ++i) {
        *(matrix + i) = malloc(cols * sizeof(FLOAT));
    }
    return matrix;
}

void freeMatrix(FLOAT **matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        free(*(matrix + i));
    }
    free(matrix);
}

//Get data from files
bool getData(const char *fileName, FLOAT **x, int rows, int cols) {
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
            fscanf(fp, "%lf", &*(*(x + i) + j));
        }
    }
    return true;
}

//sig2ext in rainflow
FLOAT *sig2ext(FLOAT *sigy, FLOAT *dty, int rows, int *lenOfArray) {

    //w1=diff(sig);
    //w2=logical([1;(w1(1:end-1).*w1(2:end))<=0;1]);
    //tmp1=(w1(1:end-1).*w1(2:end))<=0

    //FLOAT w1[rows - 1], tmp1[rows - 2], w2[rows];
    FLOAT *w = malloc((rows) * sizeof(FLOAT));
    int n = 0, m = 0, k = 0, l = 0;

//    w2[0] = 1;
//    w2[rows - 1] = 1;

    for (int i = 1; i < rows; i++) {
        w[i - 1] = sigy[i] - sigy[i - 1];
    }
//    for (int i = 1; i < rows - 1; i++) {
//        tmp1[i - 1] = w[i - 1] * w[i];
//        if (tmp1[i - 1] <= 0) {
//            tmp1[i - 1] = 1;
//        } else {
//            tmp1[i - 1] = 0;
//            n++;
//        }
//        w2[i] = tmp1[i - 1];
//    }


    for (int i = rows - 2; i > 0; i--) {
        FLOAT tmp = w[i] * w[i - 1];
        if (tmp <= 0) {
            w[i] = 1;
        } else {
            w[i] = 0;
            n++;
        }
    }
    w[0] = w[rows - 1] = 1;

    //ext1=sigy(w2);	exttime1=dty(w2);
    //FLOAT ext1[rows - n], exttime1[rows - n];

    FLOAT *ext = malloc((rows - n) * sizeof(FLOAT));
    FLOAT *exttime = malloc((rows - n) * sizeof(FLOAT));
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

    //w3=diff(ext1); tmp2 = w3(1:end-1)==0 & w3(2:end)==0; w4=~logical([0; w3(1:end-1)==0 & w3(2:end)==0; 0]);
    //FLOAT w3[rows - n - 1], tmp2[rows - n - 2], w4[rows - n];

//    w4[0] = 1;
//    w4[rows - n - 1] = 1;

    for (int i = 1; i < rows - n; i++) {
        w[i - 1] = ext[i] - ext[i - 1];
    }
//    for (int i = 1; i < rows - n - 1; i++) {
//        if (w3[i - 1] == 0 && w3[i] == 0) {
//            tmp2[i - 1] = 1;
//            w4[i] = 0;
//            m++;
//        } else {
//            tmp2[i - 1] = 0;
//            w4[i] = 1;
//        }
//    }

    for (int i = rows - n - 2; i > 0; i--) {
        if (w[i - 1] == 0 && w[i] == 0) {
            w[i] = 0;
            m++;
        } else {
            w[i] = 1;
        }
    }
    w[0] = w[rows - n - 1] = 1;



    //ext2=ext1(w4); exttime2=exttime1(w4);
    //FLOAT ext2[rows - n - m], exttime2[rows - n - m];
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

    //w5=~logical([0; ext2(1:end-1)==ext2(2:end)]);

    //FLOAT w5[rows - n - m];
    //w5[0] = 1;
    for (int i = 1; i < rows - n - m; i++) {
        if (ext[i - 1] == ext[i]) {
            w[i] = 0;
            k++;
        } else {
            w[i] = 1;
        }
    }
    w[0] = 1;

    //ext3=ext2(w5);

    //FLOAT ext3[rows - n - m - k];
    for (int i = 0, j = 0; i < rows - n - m - k && j < rows - n - m;) {
        if (w[j] == 0) {
            j++;
        } else {
            ext[i] = ext[j];
            i++;
            j++;
        }
    }
    //w6=(exttime2(2:end)-exttime2(1:end-1))./2
    //exttime3=[exttime2(1:end-1)+w6.*~w5(2:end); exttime2(end)];
    //exttime4=exttime3(w5);

    //FLOAT w6[rows - n - m - 1], exttime3[rows - n - m], exttime4[rows - n - m - k];
    //exttime3[rows - n - m - 1] = exttime2[rows - n - m - 1];

    FLOAT *w2 = malloc((rows - n - m - 1) * sizeof(FLOAT));
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

    //length(ext3)>2,  w7=diff(ext3); w8=logical([1; w7(1:end-1).*w7(2:end)<0; 1]); ext4=ext3(w8); exttime5=exttime4(w8);
    FLOAT *ext4 = NULL;
    *lenOfArray = 0;
    if (rows - n - m - k > 2) {

        //FLOAT w7[rows - n - m - k - 1], w8[rows - n - m - k];
        //w8[0] = 1;
        //w8[rows - n - m - k - 1] = 1;

        for (int i = 1; i < rows - n - m - k; i++) {
            w[i - 1] = ext[i] - ext[i - 1];
        }
//        for (int i = 1; i < rows - n - m - k - 1; i++) {
//            if (w7[i - 1] * w7[i] < 0) {
//                w8[i] = 1;
//            } else {
//                w8[i] = 0;
//                l++;
//            }
//        }

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
        ext4 = malloc(2 * (*lenOfArray) * sizeof(FLOAT));

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

//rainFlow in rainflow
FLOAT **rainFlow(FLOAT *ext, FLOAT *exttime, int lenOfSig2ext, int *lenOfRainflow) {
    FLOAT **rfy = NULL, **rfyResult = NULL;

    //function rfy5
    FLOAT a[16384], t[16384], ampl, mean, period, atime;
    int cNr = 1;
    int j = -1;

    //create 2D rfy(5 * (lenOfSig2ext -1))
    rfy = (FLOAT **) malloc(5 * sizeof(FLOAT *));
    for (int i = 0; i < 5; i++) {
        rfy[i] = (FLOAT *) malloc((lenOfSig2ext - 1) * sizeof(FLOAT));
    }

    int columnId = 0;
    int pointId = 0;

    for (int i = 0; i < lenOfSig2ext; i++) {
        a[++j] = *(ext + pointId);
        //printf("%.10f\n", *(ext + pointId));
        t[j] = *(exttime + pointId);
        //printf("%.10f\n", *(exttime + pointId));
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
    rfyResult = (FLOAT **) malloc(5 * sizeof(FLOAT *));
    for (int i = 0; i < 5; i++) {
        rfyResult[i] = (FLOAT *) malloc((lenOfSig2ext - cNr) * sizeof(FLOAT));
    }
    *lenOfRainflow = lenOfSig2ext - cNr;

    for (int i = 0; i < 5; i++) {
        for (int k = 0; k < lenOfSig2ext - cNr; k++) {
            rfyResult[i][k] = rfy[i][k];
        }
    }

    //free 2D rfy(5 * (lenOfSig2ext -1))
    for (int i = 0; i < 5; i++) {
        free(rfy[i]);
    }
    free(rfy);

    return rfyResult;
}

//rfhist in rainflow
FLOAT *rfhist(FLOAT **rfy, int lenOfRainflow, int *lenOfRfhist) {
    FLOAT *noy = NULL, *xoy = NULL;
    int x = 32;
    *lenOfRfhist = x;

    xoy = malloc(x * sizeof(FLOAT));
    noy = malloc(x * sizeof(FLOAT));
    memset(noy, 0, x * sizeof(FLOAT));

    //halfc=find(rfy(3,:)==0.5)
    int halfcNum = 0;
    for (int i = 0; i < lenOfRainflow; i++) {
        if (rfy[2][i] == 0.50)
            halfcNum++;
    }

//    printf("halfcNUM:\n");
//    printf("%d\n", halfcNum);

    int *halfc = NULL;
    halfc = malloc(halfcNum * sizeof(int));
    for (int i = 0, j = 0; i < lenOfRainflow && j < halfcNum;) {
        if (rfy[2][i] == 0.50) {
            halfc[j] = i;
            j++;
        }
        i++;
    }

//    printf("halfc:\n");
//    for (int i = 0; i < halfcNum; i++) {
//        printf("%d\n", halfc[i]);
//    }

    FLOAT min = rfy[0][0], max = rfy[0][0];
    for (int i = 0; i < lenOfRainflow; i++) {
        if (rfy[0][i] >= max) {
            max = rfy[0][i];
        } else if (rfy[0][i] <= min) {
            min = rfy[0][i];
        }
    }

    FLOAT wid = (max - min) / x;
    for (int i = 0; i < x; i++) {
        xoy[i] = min + (FLOAT) (i + 0.50) * wid;
    }


    //printf("min data %.10f wid %.10f\n", min, wid);
    for (int i = 0; i < lenOfRainflow; i++) {
        //printf("data rfy %.10f rfy-min %.10f rfy-min/wid %.10f floor %d\n",
        //     rfy[0][i], rfy[0][i] - min, (rfy[0][i] - min) / wid, (int) floor((rfy[0][i] - min) / wid));
        int j;
        j = (int) floor((rfy[0][i] - min) / wid);
        if (j != 0 && fabs((rfy[0][i] - min) - wid * j) < 1e-10) {
            noy[j - 1] += 1;
        } else {
            noy[j] += 1;
        }
    }

//    printf("noy:\n");
//    for (int i = 0; i < x; i++) {
//        printf("%f\n", noy[i]);
//    }
//    printf("xoy:\n");
//    for (int i = 0; i < x; i++) {
//        printf("%f\n", xoy[i]);
//    }

    //if ~isempty(halfc) {
    //[N2 x]=hist(rf(r,halfc),x)  N2 = noy2, x = *xoy
    //N1=N1-0.5*N2  N1 = noy
    // }
    if (halfcNum != 0) {
        FLOAT *noy2 = malloc(x * sizeof(FLOAT));
        memset(noy2, 0, x * sizeof(FLOAT));
        FLOAT *rf = malloc(halfcNum * sizeof(FLOAT));
        //printf("rf:\n");
        for (int i = 0; i < halfcNum; i++) {
            int j = halfc[i];
            rf[i] = rfy[0][j];
            //  printf("%f\n", rf[i]);
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

//        printf("noy2:\n");
//        for (int i = 0; i < x; i++) {
//            printf("%f\n", noy2[i]);
//        }

        //printf("noy:\n");
        for (int i = 0; i < x; i++) {
            noy[i] -= noy2[i] * 0.5;
            //    printf("%f\n", noy[i]);
        }
        free(rf);
    }

    FLOAT *rfhist = malloc(2 * x * sizeof(FLOAT));
    for (int i = 0; i < x; i++) {
        rfhist[i] = noy[i];
        rfhist[i + x] = xoy[i];
    }

    free(halfc);
    free(noy);

    return rfhist;
}

//test preData
void testPreData() {

    //printf("*************************Dysum:***************************************\n");
    for (int kk = 0; kk < 9; kk++) {
        //FLOAT sigy[aaaRow], dty[aaaRow];
        FLOAT *sigy = malloc(aaaRow * sizeof(FLOAT));
        FLOAT *dty = malloc(aaaRow * sizeof(FLOAT));
        for (int i = 0; i < aaaRow; i++) {
            sigy[i] = aaa[i][kk + 2];
            dty[i] = aaa[i][1];
        }

        FLOAT *ext = NULL, *exttime = NULL;
        int lenOfSig2ext;
        ext = sig2ext(sigy, dty, aaaRow, &lenOfSig2ext);
        exttime = ext + lenOfSig2ext;

//    for (int i = 0; i < lenOfSig2ext; i++) {
//        printf("%.10f\n", ext[i]);
//    }
//    printf("\n");
//    for (int i = 0; i < lenOfSig2ext; i++) {
//        printf("%.10f\n", exttime[i]);
//    }
//    printf("\n");

        FLOAT **rfy = NULL;
        int lenOfRainflow;
        rfy = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);

//    for (int i = 0; i < 5; i++) {
//        for (int j = 0; j < lenOfRainflow; j++) {
//            printf("%.10f\n", rfy[i][j]);
//        }
//    }
//    printf("\n");


        FLOAT *noy = NULL, *xoy = NULL;
        int lenOfRfhist;
        noy = rfhist(rfy, lenOfRainflow, &lenOfRfhist);
        xoy = noy + lenOfRfhist;

//
//    for (int i = 0; i < lenOfRfhist; i++) {
//        printf("%.10f\n", noy[i]);
//    }
//    printf("\n");
//    for (int i = 0; i < lenOfRfhist; i++) {
//        printf("%.10f\n", xoy[i]);
//    }
//    printf("\n");

        for (int i = 0; i < lenOfRfhist; i++) {
            Dysum[kk] += noy[i] * pow(xoy[i] * 0.21 / 70, 3.5);
        }
        //printf("%e\n", Dysum[kk]);
    }


    //test calculate Dzsum & y

    /*
    FLOAT **Tzb= NULL;
    //create 2D Tzb(aRow * 9)
    Tzb = (FLOAT **) malloc(aRow * sizeof(FLOAT *));
    for (int i = 0; i < aRow; i++) {
        Tzb[i] = (FLOAT *) malloc(9 * sizeof(FLOAT));
    }

    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 150; j++) {
            //0.7812*aa(1,i)*a(j,3)+4.0000*aa(2,i)*a(j,4)+0.5029*aa(3,i)*a(j,5)+4.0000*aa(4,i)*a(j,6)+3.1156*aa(5,i)*a(j,7)+0.5011*aa(6,i)*a(j,8)+1.4676*aa(7,i)*a(j,9)+0.6966*aa(8,i)*a(j,10)+4.9997*aa(9,i)*a(j,11)+4.0000*aa(10,i)*a(j,12)
            Tzb[j][i] = 0.7812 * aa[0][i] * a[j][2] + 4.0000 * aa[1][i] * a[j][3] + 0.5029 * aa[2][i] * a[j][4] +
                        4.0000 * aa[3][i] * a[j][5] + 3.1156 * aa[4][i] * a[j][6] + 0.5011 * aa[5][i] * a[j][7] +
                        1.4676 * aa[6][i] * a[j][8] + 0.6966 * aa[7][i] * a[j][9] + 4.9997 * aa[8][i] * a[j][10] +
                        4.0000 * aa[9][i] * a[j][11];
        }
    }


    FLOAT Dzsum[9] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
    printf("*************************Dzsum:***************************************\n");

    for (int k = 0; k < 9; k++) {

        FLOAT sig[aRow], dt[aRow];
        for (int i = 0; i < aRow; i++) {
            sig[i] = Tzb[i][k];
            dt[i] = a[i][1];
        }

        FLOAT *ext2 = NULL, *exttime2 = NULL;
        int lenOfSig2ext;
        ext2 = sig2ext(sig, dt, aRow, &lenOfSig2ext);
        exttime2 = ext2 + lenOfSig2ext;

        FLOAT **rf = NULL;
        int lenOfRainflow;
        rf = rainFlow(ext2, exttime2, lenOfSig2ext, &lenOfRainflow);

        FLOAT *no = NULL, *xo = NULL;
        int lenOfRfhist;
        no = rfhist(rf, lenOfRainflow, &lenOfRfhist);
        xo = no + lenOfRfhist;

        for (int i = 0; i < lenOfRfhist; i++) {
            Dzsum[k] += no[i] * pow(xo[i] * 0.21 / 70, 3.5);
        }
        printf("%e\n", Dzsum[k]);
    }

    FLOAT y = 0;
    for(int i = 0; i < 9; i++){
        y += pow((Dysum[i] - Dzsum[i]), 2);
    }
    printf("FitnessFCN:%e\n", y);
    */

}

//fitness Function
FLOAT fitnessFcn(FLOAT *x) {

    FLOAT **Tzb = NULL;
    //create 2D Tzb(aRow * 9)
    Tzb = (FLOAT **) malloc(aRow * sizeof(FLOAT *));
    for (int i = 0; i < aRow; i++) {
        Tzb[i] = (FLOAT *) malloc(9 * sizeof(FLOAT));
    }

    FLOAT y = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 150; j++) {
            //Tzb(j,i)=x(1)*aa(1,i)*a(j,3)+x(2)*aa(2,i)*a(j,4)+x(3)*aa(3,i)*a(j,5)+x(4)*aa(4,i)*a(j,6)+x(5)*aa(5,i)*a(j,7)+x(6)*aa(6,i)*a(j,8)+x(7)*aa(7,i)*a(j,9)+x(8)*aa(8,i)*a(j,10)+x(9)*aa(9,i)*a(j,11)+x(10)*aa(10,i)*a(j,12);
            Tzb[j][i] = x[0] * aa[0][i] * a[j][2] + x[1] * aa[1][i] * a[j][3] + x[2] * aa[2][i] * a[j][4] +
                        x[3] * aa[3][i] * a[j][5] + x[4] * aa[4][i] * a[j][6] + x[5] * aa[5][i] * a[j][7] +
                        x[6] * aa[6][i] * a[j][8] + x[7] * aa[7][i] * a[j][9] + x[8] * aa[8][i] * a[j][10] +
                        x[9] * aa[9][i] * a[j][11];
        }
    }

    for (int k = 0; k < 9; k++) {

        //FLOAT sig[aRow], dt[aRow];
        FLOAT *sig = malloc(aRow * sizeof(FLOAT));
        FLOAT *dt = malloc(aRow * sizeof(FLOAT));
        for (int i = 0; i < aRow; i++) {
            sig[i] = Tzb[i][k];
            dt[i] = a[i][1];
        }

        FLOAT *ext = NULL, *exttime = NULL;
        int lenOfSig2ext;
        ext = sig2ext(sig, dt, aRow, &lenOfSig2ext);
        exttime = ext + lenOfSig2ext;

        FLOAT **rf = NULL;
        int lenOfRainflow;
        rf = rainFlow(ext, exttime, lenOfSig2ext, &lenOfRainflow);

        FLOAT *no = NULL, *xo = NULL;
        int lenOfRfhist;
        no = rfhist(rf, lenOfRainflow, &lenOfRfhist);
        xo = no + lenOfRfhist;

        for (int i = 0; i < lenOfRfhist; i++) {
            Dzsum[k] += no[i] * pow(xo[i] * 0.21 / 70, 3.5);
        }
        //printf("%e\n", Dzsum[k]);
    }

    for (int i = 0; i < 9; i++) {
        //constraint : c =(Dysum[i]-Dzsum[i]) <= 0
        FLOAT c = Dysum[i] - Dzsum[i];
        if (c <= 0) {
            y += pow(c, 2);
        } else {
            y = 100;
        }
    }

    return y;
}

//initial population
void initialPopulation() {
    for (int i = 0; i < populationSize; i++) {
        for (int j = 0; j < chromosomeSize; j++) {
            populationArray[i][j] = (UB[j] - LB[j]) * ((FLOAT) rand() / RAND_MAX) + LB[j];
        }
    }
}

//sum fitness
FLOAT sum(FLOAT *fitness) {
    FLOAT sum = 0;
    for (int i = 0; i < populationSize; i++) {
        sum += fitness[i];
    }
    return sum;
}

//best fitness position
FLOAT *bestFitness(FLOAT *fitness) {

    //bestRes[bestFitness][bestIndex]
    FLOAT bestFitness = fitness[0];
    int bestIndex = 0;
    FLOAT *bestRes = malloc(2 * sizeof(FLOAT));
    for(int i = 0; i < populationSize; i++){
        if(fitness[i] < bestFitness){
            bestFitness = fitness[i];
            bestIndex = i;
        }
    }
    bestRes[0] = bestFitness;
    bestRes[1] = bestIndex;

    return bestRes;
}

//rank fitness
int *rankForElitism(FLOAT *fitness) {

    // initialize rank array
    int *rank = malloc(populationSize * sizeof(int));
    for (int i = 0; i < populationSize; i++) {
        rank[i] = i;
    }

    // rank fitness in increase order
    for (int i = populationSize - 1; i > 0; i--) {
        for (int j = 0; j < i; j++) {
            if (fitness[rank[j]] > fitness[rank[j + 1]]) {
                int tmp_rank = rank[j];
                rank[j] = rank[j + 1];
                rank[j + 1] = tmp_rank;
            }
        }
    }

    return rank;
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

    //calculate Dysum
    testPreData();

    //initial population
    initialPopulation();


    //start GA
    //for maxGeneration



    //calculate fitness of every population
    for (int i = 0; i < populationSize; i++) {
        fitness[i] = fitnessFcn(populationArray[i]);
    }

    //每一代平均最优适应度
    FLOAT sumFit = sum(fitness);
    aveFitnessOfGen[00000] = sumFit / populationSize;

    //每一代最优适应度及其位置
    //bestRes[bestFitness][bestIndex]
    FLOAT *bestRes = bestFitness(fitness);
    bestFitnessOfGen = bestRes[0];
    bestIndexOfGen = (int)bestRes[1];

    //找到最新最优适应度及其对应位置包含的x值
    if(bestFitnessOfGen < fval){
        fval = bestFitnessOfGen;
        for(int i = 0; i < chromosomeSize; i++){
            X_10[i] = populationArray[bestIndexOfGen][i];
        }
        G = 00000;
    }

    //if(i = maxGeneration - 1) break;

    //select
    selectFcn(populationArray);

    //cross
    crossFcn(populationArray);

    //mutation
    mutationFcn(populationArray);

    //stop GA


    freeMatrix(a, aRow);
    freeMatrix(aa, 10);
    freeMatrix(aaa, aaaRow);

    return 0;
}