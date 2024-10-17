#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "coco.h"


#define N 100 //集団の個体数、推奨値 5*Dから10*Dです
#define D 20//個体の次元
#define F 0.5 // 差分突然変異におけるスケール係数
#define Cr 0.9 //交叉率
//0,sphere,1,rastrigin-cos,2,rosenbrock-100
#define METHOD 0
//0,basic,1,lamarckian,2,baldwinian
#define REPAIR 1
#define G 2000
#define ROOP 31
#define TITLE "bbob-mixint-de-sphere-lamarckian-20D"

typedef struct{
    int label;
    double value;
}RESULT;

void print2Ddouble(double[][D]);
void initGroup(double[][D],int[]);
void selection(double[][D],int[]);
void printDdouble(double[]);
void update(double[], double[]);
double f(double[], int[]);
double sphere(double[], int[]);
double rastrigin(double[]);
double rosenbrock(double[]);
void crossover(double[],double[],double[]);
void evaluate(double[][D], RESULT[],int []);
void generateV(double[][D], int[], double[],int[],int);
void print_in_chunks(double[],int);
int cmpDscValue(const void * n1, const void * n2)
{
	if (((RESULT *)n1)->value > ((RESULT *)n2)->value)
	{
		return 1;
	}
	else if (((RESULT *)n1)->value < ((RESULT *)n2)->value)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

int main(void){
    double x[N][D] = {};//集団
    double x_best[D] = {};// 最良の個体
    double v[D] = {}; //変異ベクトル
    double u[N][D]  = {};//トライアルベクトル
    double all_result[G] = {};
    RESULT result[N] = {};
    int vector[3] = {};
    int g = 0;
    int cnt = 100;
    double x_score = 0.0;
    double u_score = 0.0;
    int l[5] = {2,4,8,16,INT_MAX};
    int lcnt = 1;
    double seido = 0.0;

    FILE *fp;
    char filename[256];
    sprintf(filename, "./data/%s.dat", TITLE);
    fp = fopen(filename, "w");
    if(REPAIR == 0){
        for(int i = 0; i < 4; i++){
            l[i] = INT_MAX;
        }
    }
    srand((unsigned int) time(NULL));

    int out[7] = {1,0,-1,-2,-3,-5,-7};
    for(int r = 0; r < ROOP; r++){
        g = 1;
        cnt = 100;
        lcnt = 0;
        //STEP1 初期集団の個体を生成
        initGroup(x,l);
        // printf("x0\n");
        // print_in_chunks(x[0],5);
        evaluate(x,result,l);//初期集団の集団の評価値計算
        // for(int i = 0; i < N; i++){
        //     if(result[i].label == 0){
        //         printf("score:%lf\n",result[i].value);
        //         print_in_chunks(x[result[i].label],5);
        //         printf("\n");
        //     }
        // }
        seido = pow(10.0, out[lcnt]); // 1e^lcntの計算
        fprintf(fp,"0 %.15lf\n",result[0].value);
        while(1){
            if(10000*D < cnt)
                break;
            // if(g == G)
            //     break;
            fprintf(fp,"%.15lf %.100lf\n",log10((double)cnt/(double)D),result[0].value);
            //STEP2 終了条件
            if(result[0].value < seido && lcnt < 7){
                printf("x_1e%d:\n",out[lcnt]);
                print_in_chunks(x[result[0].label],5);
                printf("score:");
                printf("%.30lf ",result[0].value);
                printf("evaluate_cnt%d\n",cnt);
                lcnt++;
                seido = pow(10.0, out[lcnt]); // 1e^lcntの計算
            }

            for(int i = 0; i < N; i++){//STEP3　集団のi番目の個体をtargetベクトルx_iとする
                selection(x,vector);//STEP4　集団から3つの個体を選ぶ
                generateV(x,vector,v,l,i);//STEP5 baseベクトルx_bに摂動を加え、mutateベクトルvを生成する
                crossover(v,x[i],u[i]);//STEP6 mutateベクトルvとtargetベクトルx_iを交叉しtrialベクトルuを生成する
                // printf("x\n");
                // print_in_chunks(x[i],5);
                // printf("v\n");
                // print_in_chunks(v,5);
                // printf("u\n");
                // print_in_chunks(u[i],5);
            }

            //STEP7 解を更新
            for(int i = 0; i < N; i++){
                x_score = f(x[i],l);
                u_score = f(u[i],l);
                result[i].value = x_score;
                result[i].label = i;
                if(u_score <= x_score){
                    update(u[i],x[i]);
                    result[i].value = u_score;
                }
            }
            cnt += N;//評価回数の加算
            qsort(result, N, sizeof(RESULT), cmpDscValue);

            //STEP8 世代を進める
            g++;
        }
        
        printf("x_best:\n");
        print_in_chunks(x[result[0].label],5);
        printf("score:");
        printf("%.30lf ",result[0].value);
        printf("evaluate_cnt%d\n",cnt);
        break;
    }
    return 0;
}

void print_in_chunks(double arr[], int chunk_size) {
    for (int i = 0; i < D; i += chunk_size) {
        for (int j = i; j < i + chunk_size && j < D; j++) {
            printf("%lf ", arr[j]);
        }
        printf("\n");
    }
}

void print2Ddouble(double x[][D]){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            printf("%lf ",x[i][j]);
        }
        printf("\n");
    }
}

void printDdouble(double x[]){
    for(int i = 0; i < D; i++){
        printf("%.15lf ",x[i]);
    }
}

void initGroup(double x[][D],int l[]){
    int cnt = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < D; j++){
            if(REPAIR == 0 || cnt == 4){ 
                x[i][j] = 10.0*((double)rand() / RAND_MAX) - 5.0;
            }
            else{
                x[i][j] =  ((double)rand() / RAND_MAX) * (l[cnt] - 1);
            }
            cnt++;
            if(cnt == 5){
                cnt = 0;
            }
        }
    }
}

void selection(double x[][D],int vector[]){
    vector[0] = rand() % N;
    do {
        vector[1] = rand() % N;
    } while (vector[1] == vector[0]);

    do {
        vector[2] = rand() % N;
    } while (vector[2] == vector[0] || vector[2] == vector[1]);
}

void update(double x1[], double x2[]){
    for(int i = 0; i < D; i++){
        x2[i] = x1[i];
    }
}

double f(double x[],int l[]){
    double score = 0.0;
    switch(METHOD){
        case 0:
            score = sphere(x,l);
            break;
        case 1:
            score = rastrigin(x);
            break;
        case 2:
            score = rosenbrock(x);
            break; 
    }
    return score;
}

double sphere(double x[],int l[]){
    // double y[16] = {};
    // double y_star;
    // double z[16] = {};
    // int  cnt = 0;
    // double x_tmp[D] = {};
    // double x_opt = 0.0;
    // double score = 0.0;

    // if(REPAIR != 0){
    //     for(int i = 0; i < D; i++){
    //         x_tmp[i] = x[i];
    //     }

    //     for(int i = 0; i < D; i++){
    //         if(cnt != 4){
    //             // 補助値を計算
    //             for(int  j = 0; j < l[cnt]; j++){
    //                 y[j] = (j+1) * (double)(l[cnt] - 1) / (double)(l[cnt] + 1);
    //             }

    //             // 最適値に最も近い補助値 y^* を求める
    //             double min_dist = fabs(y[0] - x_opt);
    //             y_star = y[0];
    //             for (int j = 1; j < l[cnt]; j++) {
    //                 double dist = fabs(y[j] - x_opt);
    //                 if (dist < min_dist) {
    //                     min_dist = dist;
    //                     y_star = y[j];
    //                 }
    //             }
    //             // 離散化された z を計
    //             for (int j = 0; j < l[cnt]; j++) {
    //                 z[j] = y[j] - (y_star - x_opt);
    //                 //printf("%d %lf %lf %lf\n",j,y[j],z[j],z[j] + (double)((l[cnt] - 1)) / (double)(l[cnt]+1));
    //             }

    //             // x を整数に変換
    //             for (int j = 0; j < l[cnt]; j++) {
    //                 if(x[i] < z[j] + (double)((l[cnt] - 1)) / (double)(l[cnt]+1)){
    //                     x_tmp[i] = j;
    //                     if(REPAIR == 1){
    //                         x[i] = j;
    //                     }
    //                     break;
    //                 }
    //                 else if(j == l[cnt] - 1){
    //                     x_tmp[i] = l[cnt] - 1;
    //                     if(REPAIR == 1){
    //                         x[i] = l[cnt] - 1;
    //                     }
    //                     break;
    //                 }
    //             }
    //         }
    //         cnt++;
    //         if(cnt == 5){
    //             cnt = 0;
    //         }
    //     }
    //     for(int i = 0; i < D; i++){
    //         score += x_tmp[i]*x_tmp[i];
    //     }
    //     return score;
    // }

    // for(int i = 0; i < D; i++){
    //     score += x[i]*x[i];
    // }
    // return score;


    double y[16] = {};
    double y_star;
    double z[16] = {};
    int  cnt = 0;
    double x_tmp[D] = {};
    double x_opt = 0.0;
    double score = 0.0;

    if(REPAIR != 0){
        for(int i = 0; i < D; i++){
            x_tmp[i] = x[i];
        }

        for(int i = 0; i < D; i++){
            if(cnt != 4){
                // 補助値を計算
                for(int  j = 0; j < l[cnt]; j++){
                    y[j] = j;
                }

                //連続値に最も近い補助値 y^* を求める
                double min_dist = fabs(y[0] - x_tmp[i]);
                y_star = y[0];
                for (int j = 1; j < l[cnt]; j++) {
                    double dist = fabs(y[j] - x_tmp[i]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        y_star = y[j];
                    }
                }

                // x を整数に変換
                x_tmp[i] = y_star;
                //printf("%lf %lf\n",x[i],x_tmp[i]);
                if(REPAIR == 1){
                    x[i] = y_star;
                }
            }
            cnt++;
            if(cnt == 5){
                cnt = 0;
            }
        }
        for(int i = 0; i < D; i++){
            score += x_tmp[i]*x_tmp[i];
        }
        return score;
    }

    for(int i = 0; i < D; i++){
        score += x[i]*x[i];
    }
    return score;


    // double  x_tmp[D] = {};
    // int cnt = 0;
    // int l_cnt = 0;
    // double min_dist;
    // double dist;
    // int y_star = 0;
    // double score = 0.0;
    // // printf("x\n");
    // // print_in_chunks(x,5);
    // // printf("y\n");
    // if(REPAIR != 0){
    //     for(int i = 0; i < D; i++){
    //         x_tmp[i] = x[i];
    //     }
    //     for(int i = 0; i < D; i++){
    //         if(cnt != 4){
    //             min_dist = fabs(x[i] - l_cnt);
    //             y_star = l_cnt;
    //             l_cnt++;
    //             while(l_cnt < l[cnt]){
    //                 //printf("%d %lf %d\n",cnt,x[i],l_cnt);
    //                 dist = fabs(x[i] - l_cnt);
    //                 if(dist < min_dist){
    //                     min_dist = dist;
    //                     y_star = l_cnt;
    //                 }
    //                 l_cnt++;
    //             }
    //             x_tmp[i] = y_star;
    //             if(REPAIR == 1){
    //                 x[i] = y_star;
    //             }
    //             l_cnt = 0; 
    //             cnt++;
    //         }
    //         else{
    //             cnt = 0;
    //         }
    //     }
    //     //print_in_chunks(x_tmp,5);
    //     for(int i = 0; i < D; i++){
    //         score += x_tmp[i]*x_tmp[i];
    //     }
    //     return score;
    // }
    // for(int i = 0; i < D; i++){
    //     score += x[i]*x[i];
    // }
    // return score;


    // double  x_tmp[D] = {};
    // int cnt = 0;
    // int l_cnt = 0;
    // double min_dist;
    // double dist;
    // int y_star = 0;
    // double score = 0.0;
    // // printf("x\n");
    // // print_in_chunks(x,5);
    // // printf("y\n");
    // if(REPAIR != 0){
    //     for(int i = 0; i < D; i++){
    //         x_tmp[i] = round(x[i]);
    //     }

    //     for(int i = 0; i < D; i++){
    //         if(REPAIR == 1){
    //             x[i] = round(x[i]);
    //         }
    //     }
    //     //print_in_chunks(x_tmp,5);
    //     for(int i = 0; i < D; i++){
    //         score += x_tmp[i]*x_tmp[i];
    //     }
    //     return score;
    // }
    // for(int i = 0; i < D; i++){
    //     score += x[i]*x[i];
    // }
    // return score;
}

double rastrigin(double x[]){
    double score = 10*D;
    for(int i = 0; i < D; i++){
        score += x[i]*x[i] - 10*cos(2*M_PI*x[i]);
    }
    return score;
}

double rosenbrock(double x[]){
    double score = 0;
    for(int i = 0; i < D-1; i++){
        score += 100*(x[i+1] - x[i]*x[i])*(x[i+1] - x[i]*x[i]) + (x[i] - 1)*(x[i] - 1);
    }
    return score;
}
 
void crossover(double v[],double x[],double u[]){
    double rnd_vals[D];
    int j_rand = rand() % D;  // j_rand is a random index between 0 and dim-1

    // Generate random values between 0 and 1
    for (int i = 0; i < D; i++) {
        rnd_vals[i] = (double)rand() / RAND_MAX;
    }

    // Set rnd_vals[j_rand] to 0.0
    rnd_vals[j_rand] = 0.0;

    // Perform binomial crossover
    for (int i = 0; i < D; i++) {
        if (rnd_vals[i] <= Cr) {
            u[i] = v[i];
        } else {
            u[i] = x[i];
        }
    }
}

void evaluate(double x[][D], RESULT result[],int l[]){
    for(int i = 0; i < N; i++){
        result[i].label = i;
        result[i].value = f(x[i],l);
    }
    qsort(result, N, sizeof(RESULT), cmpDscValue);
}

void generateV(double x[][D],int vector[], double v[],int l[],int j){
    int cnt = 0;
    for (int i = 0; i < D; i++) {
        v[i] = x[vector[0]][i] + F * (x[vector[1]][i] - x[vector[2]][i]);

        if(REPAIR == 0 || cnt == 4){
            if (v[i] < -5.0){
                //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(-5.0 - x[vector[0]][i]);  
                v[i] = (-5.0 + x[j][i]) / 2.0;
            }
            else if(v[i] > 5.0){
                //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(5.0 - x[vector[0]][i]);
                v[i] = (5.0 + x[j][i]) / 2.0;
            }
        }
        else{
            if (v[i] < 0){
                //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(0 - x[vector[0]][i]);  
                v[i] = (0.0 + x[j][i]) / 2.0;
            }
            else if(v[i] > l[cnt]-1){
                //v[i] = x[vector[0]][i] + ((double)rand() / RAND_MAX)*(l[cnt]-1 - x[vector[0]][i]);
                v[i] = (l[cnt] - 1.0 + x[j][i]) / 2.0;
            }
        }
        cnt++;
        if(cnt == 5){
            cnt = 0;
        }
    }
}