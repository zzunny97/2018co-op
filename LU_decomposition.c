#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>
#define N 8000


typedef double  Data;

void dim1mat_filler(Data* mat);
void dim2mat_filler(Data** mat);
void swap(Data* a, Data* b);


Data original_A[N][N];
Data L[N][N] = {0};
Data U[N][N] = {0};


int main(int argc, char const *argv[]){
    struct timeval start_point, end_point;
    double operating_time;
    gettimeofday(&start_point, NULL);  
    srand(time(NULL));

    Data X[N], Z[N];
    
    Data** A = (Data**)malloc(sizeof(Data*)*N);
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        A[i] = (Data*)malloc(sizeof(Data)*N);
    }
    dim2mat_filler(A);

    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            original_A[i][j] = A[i][j];
        }
    }

    Data* B = (Data*)malloc(sizeof(Data)*N);
    dim1mat_filler(B);

    printf("A * X = B\n");
    printf("L * U = A\n");
    printf("L * U * X = B\n");
    printf("L * Z = B\n");
    printf("U * X = Z\n");

    /*
    printf("------------ matrix A ------------\n");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%.1f\t", A[i][j]);
        }
        printf("\n");
    }
    */
    printf("row partitioning\n");
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        if(A[i][i]==0){
            int i_index = i;
            while(A[i_index][i]==0){
                i_index++;
            }
            for(int j=0; j<N; j++){
                swap(&A[i][j], &A[i_index][j]);
            }
        }
    }
    /*
    printf("------------ matrix A ------------\n");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%.1f\t", A[i][j]);
        }
        printf("\n");
    }
    */

    printf("gaussian elimination\n");
    //#pragma omp parallel for
    for(int k=0; k<N; k++){
        for(int i=k+1; i<N; i++)
            A[i][k] /= A[k][k];
        for(int i=k+1; i<N; i++){
            for(int j=k+1; j<N; j++){
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
    /*
    printf("------------ matrix A ------------\n");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%.1f\t", A[i][j]);
        }
        printf("\n");
    }
    */
    
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        for(int j=0; j<=i; j++){
            if(i==j)    L[i][j] = 1;
            else        L[i][j] = A[i][j];
        }
    }
    /*
    printf("------------ matrix L ------------\n");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%.1f\t", L[i][j]);
        }
        printf("\n");
    }
*/
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        for(int j=i; j<N; j++){
            U[i][j] = A[i][j];
        }
    }
    /*
    printf("------------ matrix U ------------\n");
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%.1f\t", U[i][j]);
        }
        printf("\n");
    }
*/
    
    Data sum;
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        sum = B[i];
        for(int j=0; j<i; j++){
            sum -= Z[j]*L[i][j];
        }
        Z[i] = sum;
    }
/*
    printf("------------ matrix Z ------------\n");
    for(int i=0; i<N; i++){
        printf("%.1f\t", Z[i]);
    }
    printf("\n");

*/
    //#pragma omp parallel for
    for(int i=N-1; i>=0; i--){
      sum = Z[i];
      for(int j=N-1;j>i;j--)
         sum -= X[j] * U[i][j];
      X[i] = sum / U[i][i];
   }
/*
    #printf("------------ matrix X ------------\n");
    for(int i=0; i<N; i++){
        printf("%.1f\t", X[i]);
    }
    printf("\n");
    */

    printf("Is this LU decomposition correct?\n");
    Data check[N];
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        check[i] = 0;
        for(int j=0; j<N; j++){
            //printf("%f * %f = %f\n", A[i][j], X[j], A[i][j]*X[j]);
            check[i] += original_A[i][j] * X[j];
        }
    }
    //#pragma omp parallel for
    for(int i=0; i<N; i++){
        if(check[i] != B[i]){
            printf("incorrect in check[%d] = %lf, B[%d] = %lf\n",i,check[i], i, B[i]);
            continue;
        }
        else{
            continue;
        }
    }
/*
    printf("------------ matrix B ------------\n");
    for(int i=0; i<N; i++){
        printf("%.1f\t", B[i]);
    }
    printf("\n");

    printf("------------ matrix check ------------\n");
    for(int i=0; i<N; i++){
        printf("%.1f\t", check[i]);
    }
    printf("\n");
*/  
    gettimeofday(&end_point, NULL);
    operating_time = ((double)((end_point.tv_sec*1000000+end_point.tv_usec) - (start_point.tv_sec*1000000+start_point.tv_usec)))/(double)1000000;
    printf("operating time = %f\n", operating_time);
    return 0;
}

void swap(Data* a, Data* b){
    Data temp;
    temp = *a;
    *a = *b;
    *b = temp;
}


void dim1mat_filler(Data* mat){
    for(int i=0; i<N; i++){
        mat[i] = rand()%5 + 1;
    }
}


void dim2mat_filler(Data** mat){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            mat[i][j] = rand()%5 + 1;
        }
    }
}




