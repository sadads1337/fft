#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define OPENMPTH 16

const int align = 32;

double PI = 3.14159265359;

//const double Lx = 5.115; //физические размеры среды
//const double Lz = 5.115;
//const double T = 1.0; // время
const double WN7 = 5.0; // частота источника

const double tau = 0.001;
const double dx = 0.005; // 40 узлов на длину волны
const double dz = dx;

//tau = T/Nt; //шаг по времени
//dx = Lx/(Nx-1); //шаги по пространству, должны быть одинаковыми
//dz = Lz/(Nz-1);

#define Nt 1000
#define Nx 1024 //кол-во узлов сетки (от 0 до N-1)
#define Nz 1024
#define x_source 400 // координаты иточника
#define z_source 512

double *S4, *AK2, *AK4, *ROX, *ROZ, *F; //массивы коэффицентов
double *u, *w, *u1, *w1, *sxx, *szz, *sxz, *sxx1, *szz1, *sxz1; //массивы скоростей и смещений

int size, size1; //размер массивов для сеточных функций в байтах

int shift = Nz+1 - align*((Nz+1)/align); //выранивание для Nz+1 эл-та

#define u2(i,k) u2[k+(i)*Nz]
#define u(i,k) u[k+(i)*Nz]
#define w(i,k) w[k+(i)*Nz]
#define u1(i,k) u1[k+(i)*Nz]
#define w1(i,k) w1[k+(i)*Nz]
#define sxx(i,k) sxx[k+(i)*Nz]
#define szz(i,k) szz[k+(i)*Nz]
#define sxz(i,k) sxz[k+(i)*Nz]
#define sxx1(i,k) sxx1[k+(i)*Nz]
#define szz1(i,k) szz1[k+(i)*Nz]
#define sxz1(i,k) sxz1[k+(i)*Nz]
#define l(i,k) l[k+(i)*Nz]
#define m(i,k) m[k+(i)*Nz]
#define ro(i,k) ro[k+(i)*Nz]
#define S4(i,k) S4[k+(i)*Nz]
#define AK2(i,k) AK2[k+(i)*Nz]
#define AK4(i,k) AK4[k+(i)*Nz]
#define ROX(i,k) ROX[k+(i)*Nz]
#define ROZ(i,k) ROZ[k+(i)*Nz]

void getero_media() {
  double *l, *m, *ro;

  l = (double*)malloc(size);
  m = (double*)malloc(size);
  ro = (double*)malloc(size);

  for (int i = 0; i < Nx/2; i++)
    for (int k = 0; k < Nz; k++){
      ro(i, k) = 1.0;
      m(i, k) = 1.0;// Vs = 1
      l(i, k) = 2.0;// Vp = 2
    }

  for (int i = Nx/2; i < Nx; i++)
    for (int k = 0; k < Nz; k++) {
      ro(i, k) = 2.0;
      m(i, k) = 2.0; // Vs = 1
      l(i, k) = 4.0; // Vp = 2
    }

  for (int i = 0; i != Nx - 1; ++i)
    for (int k = 0; k != Nz - 1; ++k) {
      AK4(i, k) = (l(i, k) + 2 * m(i, k))*tau / dx;
      AK2(i, k) = l(i, k)*tau / dx;
      ROX(i, k) = 2 * tau / ((ro(i + 1, k) + ro(i, k))*dx);
      ROZ(i, k) = 2 * tau / ((ro(i, k + 1) + ro(i, k))*dx);
      double b = 1 / m(i + 1, k + 1) + 1 / m(i, k + 1) + 1 / m(i + 1, k) + 1 / m(i, k);
      //a=b=c=4.0;
      S4(i, k) = tau*4.0 / (dx*b);
    }

  free(l);
  free(m);
  free(ro);
}

void  gomo_media() {
  double l, m, ro;

  ro = 1.0; //Vp = ((l+2*m)/ro)^(0.5) = 2.0
  m = 1.0;  //Vs = (m/ro)^(0.5) = 1.0
  l = 2.0; //Vs/(WN7*dx) = 40 (20)

  for (int i = 0; i < size1; i++) {
    S4[i] = tau*m / dx;
    AK2[i] = l*tau / dx;
    AK4[i] = (l + 2 * m)*tau / dx;
    ROX[i] = 2 * tau / ((ro + ro)*dx);
    ROZ[i] = 2 * tau / ((ro + ro)*dx);
  }
}

void source(int IG, double WN7, double DT, double DZ, int K8){
  F=(double*)malloc((K8)*sizeof(double));
  double C,B,T,T1,T2,T3,T4,PII,X,A,S;
  int J,N,GI;
  //double t0;
  //WRITE(*,*)'BEGIN CALC F(T)'
  PII=acos(-1.0);
  int NPT1=K8;
  GI=IG;
  C=DT/DZ;
  B=0.00;
  T1=1.500;
  T2=3.500;
  T3=5.500;
  T=2.00*PII*WN7;

  if(GI<=2.0100) T4=T1;
  else {
    if((GI<=4.0100)&&(GI>=2.0100)) T4=T2;
    else T4=T3;
  }
  N=int(T4/(DT*WN7));
  X=T4/(2.00*WN7);
  A=(T/GI)*(T/GI);
  S=-X;
  for(J=0; J<NPT1; J++){
    if(J+1<=N+1) F[J]=C*exp((-A)*(S*S))*cos(T*S+B);
    else F[J]=0.00;
    S=S+DT;
  }

//!	F1(:)=F1(:)*DT/DZ*100000.0
//!	F1(:)=F1(:)*DT/DZ

  // printf("   %f ",F[2]);
  for(int i=0; i<NPT1; i++) F[i]=F[i]*DT/(DZ*DZ);

}

void allocation(){
  posix_memalign((void **)&u, align, size);
  posix_memalign((void **)&w, align, size);
  posix_memalign((void **)&u1, align, size);
  posix_memalign((void **)&w1, align, size);
  posix_memalign((void **)&sxx, align, size);
  posix_memalign((void **)&szz, align, size);
  posix_memalign((void **)&sxx1, align, size);
  posix_memalign((void **)&szz1, align, size);
  posix_memalign((void **)&sxz, align, size);
  posix_memalign((void **)&sxz1, align, size);

  //выравнивание для [1+Nz]
  u = u - shift;
  w=w - shift;
  u1=u1 - shift;
  w1=w1 - shift;
  sxx=sxx - shift;
  szz=szz - shift;
  sxz=sxz - shift;
  sxx1=sxx1 - shift;
  szz1=szz1 - shift;
  sxz1=sxz1 - shift;

  for (int i = 0; i < size1; i++){
    u[i] = 0;
    w[i] = 0;
    u1[i] = 0;
    w1[i] = 0;
    sxx[i] = 0;
    szz[i] = 0;
    sxz[i] = 0;
    sxx1[i] = 0;
    szz1[i] = 0;
    sxz1[i] = 0;
  }
}

void sigm (double *u, double *w,
           double *sxx, double *szz, double *sxz,
           double *sxx1, double *szz1, double *sxz1){
  double ux, wz;
  int i, k;

//схема 2-го порядка
#pragma omp parallel
  {
#pragma omp for nowait private(k)
    for (i = 1; i<Nx-1; i++)
#pragma omp simd
        for(k = 1; k<Nz-1; k++) {
          ux = u(i,k)-u(i-1,k);
          wz = w(i,k)-w(i,k-1);

          sxx1(i,k) = sxx(i,k) + AK4(i,k)*ux + AK2(i,k)*wz;
          szz1(i,k) = szz(i,k) + AK2(i,k)*ux + AK4(i,k)*wz;
        }
#pragma omp for private(k)
    for (i = 1; i<Nx-1; i++)
#pragma omp simd
        for(k = 1; k<Nz-1; k++) {
          sxz1(i, k) = sxz(i, k) + S4(i, k)*(u(i,k + 1) - u(i, k) + w(i + 1,k) - w(i,k));
        }
  }

}

void vel(double* u,double* w,
         double* u1,double* w1,
         double* sxx,double* szz,double* sxz){

  //double dsx, dsz;
  int i, k;

//схема 2-го порядка
#pragma omp parallel
  {
#pragma omp for nowait private(k)
    for (i = 1; i < Nx - 1; i++)
#pragma omp simd
        for (k = 1; k < Nz - 1; k++) {
          u1(i - 1, k) = u(i - 1, k) + ROX(i - 1, k)*(sxx(i, k) - sxx(i - 1, k) + sxz(i - 1, k) - sxz(i - 1, k - 1));
          w1(i, k - 1) = w(i, k - 1) + ROZ(i, k - 1)*(sxz(i, k - 1) - sxz(i - 1, k - 1) + szz(i, k) - szz(i, k - 1));
        }

    //граничные условия
#pragma omp for
    for (i = 1; i < Nx - 1; i++)
      w1(i, 0) += ROZ(i, 0)*(3 * szz(i, 0) - szz(i, 1));
  }
}

void source_force(double *sxx1, double *szz1, int t)
{
  int src_point = z_source + Nz * x_source;
  double d = F[t];
  sxx1[src_point] = sxx1[src_point] + d;
  szz1[src_point] = szz1[src_point] + d;
}

//функция вывода
void write (double *u2, double *w2, int l){
  FILE *fo;
  char f_u[20];
  sprintf(f_u,"uxz%i.sct",l);
  fo = fopen( f_u, "wb+" );
  //printf ("start");
  for (int i=0; i!=Nx; ++i){
    for (int j=0; j!=Nz; ++j){
      fwrite(&u2(i,j), 1, sizeof(double), fo);
    }
  }
  fclose(fo);
  //printf("\n");
}

int main (void){

  size = Nx*Nz*sizeof(double);
  size1 = Nx*Nz;
  //tau = T/Nt; //шаг по времени
  //dx = Lx/(Nx-1); //шаги по пространству, должны быть одинаковыми
  //dz = Lz/(Nz-1);

  if(tau>(dx/(2*sqrt(2.0)))){
    printf ("uslovie ustoichivosty ne vipolneno");
    return (0);
  }

  int max_num, num;
  max_num = omp_get_num_procs();
  num = OPENMPTH;
  omp_set_num_threads(num);
  printf("\n max_num_th = %i, set_num_th = %i \n", max_num, num);

  //коэфиценты среды
  posix_memalign( (void **)&S4, align, size);
  posix_memalign( (void **)&AK2, align, size);
  posix_memalign( (void **)&AK4, align, size);
  posix_memalign( (void **)&ROX, align, size);
  posix_memalign( (void **)&ROZ, align, size);

  S4=S4 - shift;
  AK2=AK2 - shift;
  AK4=AK4 - shift;
  ROX=ROX - shift;
  ROZ=ROZ - shift;

  for (int i = 0; i < size1; i++) {
    S4[i] = 0.0;
    AK2[i] = 0.0;
    AK4[i] = 0.0;
    ROX[i] = 0.0;
    ROZ[i] = 0.0;
  }

  // неоднородная среда:
  //getero_media();

  // однородная среда:
  gomo_media();

  // источник:
  source(4, WN7, tau, dz, Nt);

  //выделение памяти для основных массивов
  allocation();

  double t1;

  t1 = omp_get_wtime();

  //основные вычисления
  for(int t=1; t<=Nt; t+=2){
    //printf("%i: \n", t);
    sigm (u1, w1, sxx1, szz1, sxz1, sxx, szz, sxz);
    source_force(sxx, szz, t);
    vel(u1, w1, u, w, sxx, szz, sxz);

    //printf("%i:  \n", t + 1);
    sigm (u, w, sxx, szz, sxz, sxx1, szz1, sxz1);
    source_force(sxx1, szz1, t+1);
    vel(u, w, u1, w1, sxx1, szz1, sxz1);
  }

  t1 = omp_get_wtime()-t1;
  printf("\n It took me %f seconds.\n", t1);

  write(u, w, Nt);
  printf("  print %i \n", Nt);

  //deallocation();
}
