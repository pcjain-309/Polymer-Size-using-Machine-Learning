#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_FILE_NAME 100
 

float mean(int k,float B[])
{
 int i;
 float sum = 0.0, average;
 for (i = 0; i < k; i++)
 {
  sum = sum + B[i];
 }
 average = (sum/k);
 return average;
}

float std(int k1,float B1[])
{
 int i1;
 float sum1 = 0.0, sd;
 for (i1 = 0; i1 < k1; i1++)
 {
  sum1 = sum1 + (B1[i1] - mean(k1,B1))*(B1[i1] - mean(k1,B1));
 }
 sd = sqrt(sum1/(k1-1));
 return sd;
}

// cc -O3 Rg2_lammps.c -lm (Running Command)
main()

{
    FILE * fp;
    int col, a, b, s, i;  // Line counter (result)
    char filename[MAX_FILE_NAME], str[MAX_FILE_NAME];
    char c;  // To store a character read from file
    int M, N, N_steps;
    float skin, rho;
    int iii = 0;
    FILE *fp11 = fopen("def.chain","r");
    char line[1000];
    char first[3];
    int line_number = 0;
    int flag = 0;
    while (fgets(line,1000,fp11)!=NULL)
    {
      while(line[iii]!=' ' && line_number == 8 ){
        first[iii] = line[iii];
        iii++;
        flag = 1;
      }
      line_number++;
      if(flag) break;
      // break;
    }
    char ansss[iii];
    for(int jjj = 0; jjj<iii; jjj++){
      ansss[jjj] = first[jjj];
    }
    N = atoi(ansss);
    printf("%d\n",N);
    // return 0;
    M = 1,skin = 1.0, rho = 0.85;
    a = 4000000, b = 4000000+a, s = 10000, col = 7;
    char *file = "test";
    // Get file name from user. The file should be
    // either in current folder or complete path should be provided
    N_steps = (b - a)/s + 1;
    int j, k, P, ii, jj, kk, NN, step;
    float A[M*N][col], Rg[M], Rgxx[M], Rgyy[M], Rgzz[M],  Re[M], Rexx[M], Reyy[M], Rezz[M], Rb[M], Rbxx[M], Rbyy[M], Rbzz[M], RR[M]; 
    float x[N], y[N], z[N], xx[N], yy[N], zz[N], xb[N-1], yb[N-1], zb[N-1], C[N_steps][14], lx, ly, lz, xo, xhi, yo, yhi, zo, zhi;
    float xmean, ymean, zmean;
    kk = 0;
    for (i = a; i <= b; i = i + s)
    {
    sprintf(filename, "%s.%d.txt", file, i);
     k = 0;
     fp = fopen(filename, "r");
     if (fp == NULL)
     {
        printf("Could not open file %s", filename);
        return 0;
     }
     for (j = 0; j < M*N+9; j++)
     {
      if (j < 9)
       {
         fgets (str, 200, fp);
       } 	       
      else if (j > 8 && col == 8)
      {
       fscanf(fp, "%f %f %f %f %f %f %f %f", &A[k][0], &A[k][1], &A[k][2], &A[k][3], &A[k][4], &A[k][5], &A[k][6], &A[k][7]);
       k++;
      }
      else if (j > 8 && col == 7)
      {
       fscanf(fp, "%f %f %f %f %f %f %f", &A[k][0], &A[k][1], &A[k][2], &A[k][3], &A[k][4], &A[k][5], &A[k][6]);
       printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", A[k][0], A[k][1], A[k][2], A[k][3], A[k][4], A[k][5], A[k][6]);
       k++;
      }
     }
     for (ii = 0; ii < M; ii++)
     {
      for (jj = 0; jj < N; jj++)
      {
       x[jj] = A[ii*N+jj][col-6];
       y[jj] = A[ii*N+jj][col-5];
       z[jj] = A[ii*N+jj][col-4];
      }
       xmean = mean(N,x);
       ymean = mean(N,y);
       zmean = mean(N,z);
      for (jj = 0; jj < N; jj++)
      {
       xx[jj] = (x[jj] - xmean)*(x[jj] - xmean);
       yy[jj] = (y[jj] - ymean)*(y[jj] - ymean);
       zz[jj] = (z[jj] - zmean)*(z[jj] - zmean);
      }
      for (jj = 0; jj < N-1; jj++)
      {
       xb[jj] = (x[jj] - x[jj+1])*(x[jj] - x[jj+1]);
       yb[jj] = (y[jj] - y[jj+1])*(y[jj] - y[jj+1]);
       zb[jj] = (z[jj] - z[jj+1])*(z[jj] - z[jj+1]);
      }
     Rbxx[ii]   = mean((N-1),xb);
     Rbyy[ii]   = mean((N-1),yb);
     Rbzz[ii]   = mean((N-1),zb);
     Rb[ii]     = Rbxx[ii] + Rbyy[ii] + Rbzz[ii];
     Rexx[ii]     = (x[0] - x[N-1])*(x[0] - x[N-1]);
     Reyy[ii]     = (y[0] - y[N-1])*(y[0] - y[N-1]);

     Rezz[ii]     = (z[0] - z[N-1])*(z[0] - z[N-1]);
     Re[ii]     = Rexx[ii] + Reyy[ii] + Rezz[ii];
     Rgxx[ii]   = mean(N,xx);
     Rgyy[ii]   = mean(N,yy);
     Rgzz[ii]   = mean(N,zz);
//     Rg[ii]     = mean(N,xx) + mean(N,yy) + mean(N,zz);
     Rg[ii]     = Rgxx[ii] + Rgyy[ii] + Rgzz[ii];
     RR[ii]     = Re[ii]/Rg[ii];
     }
     C[kk][0]  = i;
     C[kk][1]  = sqrt(mean(M,Rbxx)); 
     C[kk][2]  = sqrt(mean(M,Rbyy)); 
     C[kk][3]  = sqrt(mean(M,Rbzz)); 
     C[kk][4]  = sqrt(mean(M,Rb)); 
     C[kk][5]  = mean(M,Rgxx); 
     C[kk][6]  = mean(M,Rgyy); 
     C[kk][7]  = mean(M,Rgzz); 
     C[kk][8]  = mean(M,Rg); 
     C[kk][9]  = mean(M,Rexx); 
     C[kk][10] = mean(M,Reyy); 
     C[kk][11] = mean(M,Rezz); 
     C[kk][12] = mean(M,Re);
     C[kk][13] = mean(M,RR);
 //   printf("\n");
    fclose(fp);
    kk = kk + 1;
   }
    float tau = 0.005 ;
    int ts =  1;
    FILE * fw;
    char fname[MAX_FILE_NAME];
    sprintf(fname, "Rg2_Globule.txt");
    fw = fopen(fname,"a");
    float ans = 0;
    for(int a = 0; a<N_steps;a++){
      ans+= sqrt(C[a][8]);
    }
    ans = ans/N_steps;
    float B = N;
    fprintf(fw,"%f\t%f\n",ans, B );
    sprintf(fname, "test_Rg2.txt");
    fw = fopen(fname, "w");
    fputs("time\t\tRbxx\t\tRbyy\t\tRbzz\t\tBond-length\tRgxx\t\tRgyy\t\tRgzz\t\tRg2\t\t\tRexx\t\tReyy\t\tRezz\t\tRee2\t\tRatio\n",fw);
    for (i = 0; i < N_steps; i++)
    {
     fprintf(fw, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tau*ts*(C[i][0]-a),C[i][1],C[i][2],C[i][3],C[i][4],C[i][5],C[i][6],C[i][7],C[i][8],C[i][9],C[i][10],C[i][11],C[i][12],C[i][13]);
    } 
    fclose(fw);
    return 0;
}
