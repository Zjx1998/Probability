#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define PI acos(-1)
//#define _DEBUG
void Torus_Cov(int* T_cov, int n, int iter, int debug){
// Ce script calcule une sequence de temp couverture en une torus annulus T_N par MCMC

    srand(7);
    int L_t[n][n];                      // Local time matrix
    for(int i = 0;i < iter; i++){
        if(i % 10000 == 0)
            printf("%d iteration processed!\n", i);
        int X_t     = 0;
        int Y_t     = 0;
        for(int j = 0; j < n; j++)
            for(int k = 0; k < n; k++)
                L_t[j][k] = 0;
        int n_V     = n * n;                // # of not-v node
        int t       = 0;
        while(n_V != 0){
            if(L_t[X_t][Y_t] == 0)
                n_V = n_V - 1;
            L_t[X_t][Y_t]   += 1;
            t   = t + 1;
            int rnd = rand() % 4;
            if(rnd == 3){
                #ifdef _DEBUG
                    printf("right\n");
                #endif
                X_t = (X_t + 1) % n;
            }
            else if(rnd == 2){
                #ifdef _DEBUG
                    printf("left\n");
                #endif
                X_t = (X_t - 1 + n) % n;
            }
            else if(rnd == 1){
                #ifdef _DEBUG
                    printf("up\n");
                #endif
                Y_t = (Y_t + 1) % n;
            }
            else{
                #ifdef _DEBUG
                    printf("down\n");
                #endif
                Y_t = (Y_t - 1 + n) % n;
            }
        }
        T_cov[i]    = t;
    }
}

int main(){
    time_t start,end;
    start = time(NULL);
    int n = 10, iter = 10000000;
    int* T_cov;
    T_cov = (int *)malloc(sizeof(int) * (iter + 1));
    Torus_Cov(T_cov, n, iter, 0);
    end = time(NULL);		 
    printf("time = %lf\n",difftime(end,start));

    // write the data
    int write = 0;
    if(write == 1){
        FILE* fp = fopen("C50.txt","ab+");
        for(int i = 0; i < iter; i++)
            fprintf(fp, "%d\n", T_cov[i]);
        fclose(fp);
    }
    free(T_cov);
    return 0;
}