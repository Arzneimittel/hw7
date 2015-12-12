#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

void K_Vektoren(double* k,double x1, double x2,double x3,double x4,double mu){
k[0] =  x3;
k[1] =  x4;
k[2] =  x1 + 2.0*x4 - ((1.0-mu)*(x1+mu)/pow(sqrt(pow(x1+mu,2)+pow(x2,2)),3)) - (mu*(x1-1+mu)/pow(sqrt(pow(x1-1+mu,2)+pow(x2,2)),3));
k[3] =  x2 - 2.0*x3 - ((1.0-mu)*x2/pow(sqrt(pow(x1+mu,2)+pow(x2,2)),3)) - (mu*x2/pow(sqrt(pow(x1-1+mu,2)+pow(x2,2)),3));
}
void K_Step(double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* x,double dt,double mu){
  double x1=x[0],x2=x[1],x3=x[2],x4=x[3];
  const double a21=1.0/5.0;
  const double a31=3.0/40.0, a32=9.0/40.0; 
  const double a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0;
  const double a51=19372.0/6561.0, a52=-25360.0/2187, a53=64448.0/6561.0, a54=-212.0/729.0;
  const double a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64= 49.0/176.0, a65=-5103.0/18656.0;
  const double a71=35.0/384.0, a72= 0.0, a73=500.0/1113.0, a74= 125.0/192.0, a75= -2187.0/6784.0, a76=11.0/84.0;
  K_Vektoren( k1, x1,  x2, x3, x4, mu);
  K_Vektoren( k2, x1+dt*(a21*k1[0]),  x2+dt*(a21*k1[1]), x3+dt*(a21*k1[2]), x4+dt*(a21*k1[3]), mu);
  K_Vektoren( k3, x1+dt*(a31*k1[0]+a32*k2[0]),  x2+dt*(a31*k1[1]+a32*k2[1]), x3+dt*(a31*k1[2]+a32*k2[2]), x4+dt*(a31*k1[3]+a32*k2[3]), mu);
  K_Vektoren( k4, x1+dt*(a41*k1[0]+a42*k2[0]+a43*k3[0]),  x2+dt*(a41*k1[1]+a42*k2[1]+a43*k3[1]), x3+dt*(a41*k1[2]+a42*k2[2]+a43*k3[2]), x4+dt*(a41*k1[3]+a42*k2[3]+a43*k3[3]), mu);
  K_Vektoren( k5, x1+dt*(a51*k1[0]+a52*k2[0]+a53*k3[0]+a54*k4[0]),  x2+dt*(a51*k1[1]+a52*k2[1]+a53*k3[1]+a54*k4[1]), x3+dt*(a51*k1[2]+a52*k2[2]+a53*k3[2]+a54*k4[2]), x4+dt*(a51*k1[3]+a52*k2[3]+a53*k3[3]+a54*k4[3]), mu);
  K_Vektoren( k6, x1+dt*(a61*k1[0]+a62*k2[0]+a63*k3[0]+a64*k4[0]+a65*k5[0]),  x2+dt*(a61*k1[1]+a62*k2[1]+a63*k3[1]+a64*k4[1]+a65*k5[1]), x3+dt*(a61*k1[2]+a62*k2[2]+a63*k3[2]+a64*k4[2]+a65*k5[2]), x4+dt*(a61*k1[3]+a62*k2[3]+a63*k3[3]+a64*k4[3]+a65*k5[3]), mu);
  K_Vektoren( k7, x1+dt*(a71*k1[0]+a72*k2[0]+a73*k3[0]+a74*k4[0]+a75*k5[0]+a76*k6[0]),  x2+dt*(a71*k1[1]+a72*k2[1]+a73*k3[1]+a74*k4[1]+a75*k5[1]+a76*k6[1]), x3+dt*(a71*k1[2]+a72*k2[2]+a73*k3[2]+a74*k4[2]+a75*k5[2]+a76*k6[2]), x4+dt*(a71*k1[3]+a72*k2[3]+a73*k3[3]+a74*k4[3]+a75*k5[3]+a76*k6[3]), mu); 
}
void Kutta_5(double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* x,double dt,double mu){
     const double b15=35.0/384.0, b25= 0.0, b35=500.0/1113.0, b45= 125.0/192.0, b55= -2187.0/6784.0, b65=11.0/84.0, b75=0.0;
     K_Step(k1,k2,k3,k4,k5,k6,k7,x,dt,mu);
     for(int i=0;i<4;i++){
     x[i]+= dt*(b15*k1[i]+b25*k2[i]+b35*k3[i]+b45*k4[i]+b55*k5[i]+b65*k6[i]+b75*k7[i]);
     }
    }
void Kutta_4(double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* x,double dt,double mu){
     const double b14=5179.0/57600.0, b24= 0.0, b34=7571.0/16695.0, b44= 393.0/640.0, b54= -92097.0/339200.0, b64=187.0/2100.0, b74=1.0/40.0;
     K_Step(k1,k2,k3,k4,k5,k6,k7,x,dt,mu);
     for(int i=0;i<4;i++){
     x[i]+= dt*(b14*k1[i]+b24*k2[i]+b34*k3[i]+b44*k4[i]+b54*k5[i]+b64*k6[i]+b74*k7[i]);
     }
}
void norm_func(double* norm, double* x, double* y,double& max1, double& max2, double& MAX){
  for(int i=0;i<4;i++){
    norm[0] = abs(y[0] - x[0]);
    norm[1] = abs(y[1] - x[1]);
    norm[2] = abs(y[2] - x[2]);
    norm[3] = abs(y[3] - x[3]);
    max1 = max(norm[0],norm[1]);
    max2 = max(norm[2],norm[3]);
    MAX =  max(max1,max2);
  }
}
int main(void){
  double k1[4],k2[4],k3[4],k4[4],k5[4],k6[4],k7[4],x[4],y[4],xold[4],yold[4],norm[4],max1,max2,MAX,dt=1e-5,tol=1e-5;
  x[0]=0.994;x[1]=0.0;x[2]=0.0;x[3]=-2.00158510637908;
  y[0]=0.994;y[1]=0.0;y[2]=0.0;y[3]=-2.00158510637908;
  const double mu=0.012277471;
  double t=0.0,L=17.065216560157;
  ofstream out ("kutta_über_dt_1e-5.txt");
  //const double c2=1.0/5.0, c3=3.0/10.0, c4=4.0/5.0, c5=8.0/9.0, c6=1.0, c7=1.0;
  out << t << "\t" << y[0]<< "\t" << y[1] << "\t" << 0 <<endl;
 while(t<L){
  for(int i=0;i<4;i++){yold[0]=y[0];yold[1]=y[1];yold[2]=y[2];yold[3]=y[3];}
  Kutta_5(k1,k2,k3,k4,k5,k6,k7,x,dt,mu);
  Kutta_4(k1,k2,k3,k4,k5,k6,k7,y,dt,mu);
  norm_func(norm,x,y,max1,max2,MAX);
  for(int i=0;i<4;i++){y[0]=yold[0];y[1]=yold[1];y[2]=yold[2];y[3]=yold[3];}
  dt = dt * pow((tol/MAX),0.20);
  Kutta_4(k1,k2,k3,k4,k5,k6,k7,y,dt,mu);
  for(int i=0;i<4;i++){x[0]=y[0];x[1]=y[1];x[2]=y[2];x[3]=y[3];}
  t+=dt;
    out << t << "\t" << y[0]<< "\t" << y[1] <<"\t" << dt <<endl;
}
return 0;
out.close();
}