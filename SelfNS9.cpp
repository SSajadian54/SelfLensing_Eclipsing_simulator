#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include <cmath>
#include "VBBinaryLensingLibrary.h"
using std::cout;
using std::endl;
using std::cin;

///=============================================================================

const double RA=180.0/M_PI;
const double KBol=1.380649*pow(10.0,-23.0);
const double Hplank=6.62607015*pow(10.0,-34.0);
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity=299792458.0;//velosity of light  m/s
const double Z_sun=float(0.0152);//#solar metalicity
const double Msun=1.98892*pow(10.,30.0); //in [kg].
const double LSun= 3.846*pow(10.0,33.0);//erg/s
const double Mjupiter=1.898*pow(10,27.0); 
const double Mearth=  5.9722*pow(10.0,24.0);
const double AU=1.495978707*pow(10.0,11.0);
const double Rsun=6.9634*pow(10.0,8.0); ///solar radius [meter]
const double Avks=double(8.20922);

///=============================================================================

const double wave[4]= {0.673,0.7865,0.532,0.797};//https://en.wikipedia.org/wiki/Photometric_system  G, T, BP, RP
const double AlAv[4]= {0.791986524645539, 0.617155245862836  ,1.0386670640894,0.601810722049874};
const double sigma[4]={0.017, 0.02,0.02, 0.02};// G, T, BP, RP  Table (2) Cardeli

//x=1.0/wave    ///https://heasarc.gsfc.nasa.gov/docs/tess/the-tess-space-telescope.html
//y=x-1.82
//a=1.0+0.17699*y-0.50447*y*y-0.02427*y*y*y+0.72085*y*y*y*y+ 0.01979*y*y*y*y*y-0.7753*y*y*y*y*y*y+0.32999*y*y*y*y*y*y*y
//b=1.41338*y+2.28305*y*y+1.07233*y*y*y-5.38434*y*y*y*y-0.62251*y*y*y*y*y+5.3026*y*y*y*y*y*y-2.09002*y*y*y*y*y*y*y
//AlAv=a+b/3.1

const double thre=0.001; 
const int    NB=int(1000); 
const int    nbh=int(30);  
const int    Nc=int(9586);///limb-darkening
const int    nw=1772;///rows in WDCat.dat
const int    Np=106; ///rows in CDPPMAG.txt
const int    Nl=110; /// wavelength TESS throughput 
const int    Ng=9575;///row in table29Claret.dat
const int nns=90;//Age-Luminosity of Netrun stars   Thermal luminosities of cooling neutron stars 

///=============================================================================
const int no=1;//Tobs is one 2*13.7 days
const int sector=0;///CTL
const int nmmax=156778 ;// 162681;
const double cadence=double(2.0/60.0/24.0);//in days 

///=============================================================================
struct source{
    double Ds;
    double nsbl, blend, magb, Ai, Mab, Map, lumi;
    double Tstar, Rstar, mass, loggs;
    double ros, limb, lon, lat, grav;
    double Logg[Nc],Tef[Nc], Limb[Nc];  
    double logg[Ng],teff[Ng],grad[Ng];  
    double cdpp;  
    double prior;
};
struct lens{//WD
    int num;  
    double ecen, inc, tet, tp, period, a; 
    double phi, RE, stepb;   
    double ratio, q, MBH, RBH, ext, Map, Mab, age, lumi;
    double lage[nns], llum[nns];
};
struct detection{
    int    numt; 
    double dchi, t, Tobs;
    double gap[26], dur;
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
   int    flag; 
};
struct tessfile{
   double lumi;  
   double Map, error, blend, teff, radius, dist;
   double elong, elat;
   double glong, glat;
   double prior; 
   double mass;
   double tmag[Np], cdpp[Np];  
};
struct doppler{
  double wave[Nl], throu[Nl]; 
  double waven, Fpl0, Fpl1;  
};
///=============================================================================
int    Extinction(extinc & ex,source & s);
void   func_lens(lens & l, source & s, extinc & ex);
double ErrorTESS(double maga); 
double Interpol(double ds, extinc & ex);
double RandN(double sigma,double);
double RandR(double down, double up);
double Fluxlimb(double limb, double rstar);  
double Kepler(double phi, double ecen); 
double Bessel(int n,double x); 
double LimbF(source & s); 
double GravF(source & s);     
double CDPPM(tessfile & t, double); 
double DopplerF(doppler & dp, double , double); 
double Fplank(double , double);  
double Ellipsoid(double, double, double, double, double, double);
time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;

///===========================================================================//
///                                                                           //
///                  Main program                                             //
///                                                                           //
///===========================================================================//	

int main()
{
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   time( &_timeNow);
   printf("START time:   %s",ctime(&_timeNow));
      
   
   VBBinaryLensing vbb;
   vbb.Tol=1.e-3;
   vbb.a1 =0.0;  
   vbb.LoadESPLTable("./ESPL.tbl");
  
  
   source s;
   lens l;
   extinc ex;
   detection d; 
   tessfile t;
   doppler dp; 
    
   FILE* film;  
   FILE* fild; 
   FILE* limbt;  
   FILE* distr;
   FILE* magdi; 
   FILE* wdcat; 
   FILE* cdppf;
   FILE* dopp;
   FILE* gradt; 
        
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   char   filenam0[40],filenam1[40], filenam2[40], filenam7[40];
   int    num, nsim, ndet, ffd, FlagD, signD;
   long int ID; 
   double deltaA, ndtran, ntran, Delt1, Delt2, disp;   
   double chi2,chi1, emt, timee, Av; 
   double ksi, x0, y0, x1, y1, z1, dis, yc, zc, dec, ra, xp, vx, yp; 
   double phase, u, us, Astar, As, frac, frac0, flux, maga, Occul;   
   double yb,zb,dt, zlim, rstar, maga2;  
   double Amin, Amax=0.0;  
   double depth, snr; 
   double error, paral;
   double impactL, impactO, Rho, u2, u3, RE, proj;



   d.Tobs=double(13.7*no*2.0);
   for(int i=0;i<26;++i) d.gap[i]=-1.0; 
   for(int i=1;i<=int(2.0*no); ++i) d.gap[i-1]=double(i-1.0+i*12.7);
   d.dur=1.0;
  
///=============================================================================      
    
    
   //double lage[nns], llum[nns]; 
   cdppf=fopen("./files/AGELNS.dat","r");
   for(int i=0; i<nns; ++i){
   fscanf(cdppf,"%lf  %lf\n",&l.lage[i], &l.llum[i]);
   if((l.lage[i]< l.lage[i-1] and i>0) or (i>0 and l.llum[i]>l.llum[i-1])){ 
   cout<<"Error(Age-lum NS) i:  "<<i<<"\t lage: "<<l.lage[i]<<"\t llum:  "<<l.llum[i]<<endl; int uue; cin>>uue;}}
   fclose(cdppf); 
   cout<<"********** File AGELNS.dat was read ************"<<endl;    
    
        
    
    
   cdppf=fopen("./files/CDPPMAG.txt","r");
   for(int i=0; i<Np; ++i){
   fscanf(cdppf,"%lf  %lf\n",&t.tmag[i], &t.cdpp[i]);
   if((t.tmag[i]<t.tmag[i-1] and i>0) or (i>0 and t.cdpp[i]<t.cdpp[i-1])){ 
   cout<<"Error counter:  "<<i<<"\t tmag:  "<<t.tmag[i]<<"\t cdpp:  "<<t.cdpp[i]<<endl;   int uue; cin>>uue;}}
   fclose(cdppf); 
   cout<<"********** File CDPPMAG.dat was read ************"<<endl;    
    
    
    
   double smet;   
   limbt=fopen("./files/table24bClaret.dat","r");
   for(int i=0; i<Nc; ++i){
   fscanf(limbt,"%lf  %lf  %lf   %lf\n",&s.Tef[i], &s.Logg[i], &smet, &s.Limb[i]);}
   fclose(limbt); 
   cout<<"********** File table24bClaret.dat was read ************"<<endl;   
   
   
   
   
   gradt=fopen("./files/table29Claret.dat","r");//Gravity-Darkening coefficient 
   for(int i=0; i<Ng; ++i){
   fscanf(gradt,"%lf  %lf  %lf   %lf\n",&s.teff[i],&s.logg[i],&smet,&s.grad[i]);}
   fclose(gradt); 
   cout<<"********** File table24bClaret.dat was read ************"<<endl;
   
   
   
   
   dopp=fopen("./files/TESS_Throught.txt","r");
   for(int i=0; i<Nl; ++i){
   fscanf(dopp,"%lf  %lf\n",&dp.wave[i], &dp.throu[i]);}
   fclose(dopp); 
   cout<<"********** File TESS_Throught.txt was read ************"<<endl;    
   
      
     
   sprintf(filenam0,"./files/light/lcF2/%c%d%c%d.dat",'N',sector,'_',9);
   distr=fopen(filenam0,"a+"); 
   fclose(distr);              
      
///=============================================================================
   
   for(nsim=1;  nsim<=50000; ++nsim){
   
   
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
  // srand(time(0)); 

   

   num=int(RandR(0.0,nmmax-1.0));
   sprintf(filenam7,"./files/%c%c%c%c.dat",'C','T','L','n'); //maded on 8th May 2024
   magdi=fopen(filenam7,"r"); 
   if(!magdi){cout<<"cannot read CTLn.dat "<<sector<<endl; exit(0);}
   for(int i=0; i<num; ++i){
   fscanf(magdi,"%d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf  %lf %lf %lf\n",
   &ID, &ra, &dec, &t.glong, &t.glat, &t.elong, &t.elat, &t.Map, &t.error, &t.teff, &t.radius, &t.dist, &t.blend, &t.prior,&t.mass,&t.lumi);
   t.blend=1.0/(1.0+t.blend);}   
   fclose(magdi); 
   cout<<"file Information was read !!!!"<<endl;



   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   // srand(time(0)); 




   s.lumi =t.lumi; 
   s.prior=t.prior;
   s.mass= t.mass; 
   s.Map=  t.Map;   
   s.blend=t.blend;  
   s.Tstar=t.teff;  
   s.Rstar=t.radius; 
   s.loggs=log10(G*s.mass*Msun*100.0/(s.Rstar*s.Rstar*Rsun*Rsun));  
   s.Ds=   t.dist/1000.0;//kpc
   s.lat=double(t.glat);  
   s.lon= double(t.glong);
   if(s.lon<=0.0) s.lon +=360.0;
   s.magb=double(s.Map+2.5*log10(s.blend));
   s.nsbl=double(1.0/s.blend);   
   s.cdpp=CDPPM(t, s.magb);    
  
   
   
   ex.flag=-1;
   ex.flag=Extinction(ex,s);
   if(ex.flag>0) Av=double(Interpol(s.Ds,ex)*Avks);
   else          Av=double(0.7*s.Ds); 
   if(Av<0.0)    Av=0.0;
   s.Ai=fabs(Av*AlAv[1])+RandN(sigma[1],1.5);
   if(s.Ai<0.0) s.Ai=0.0;   
   s.Mab=s.Map-s.Ai-5.0*log10(s.Ds*100.0);  
   s.limb=LimbF(s);
   //s.grav=GravF(s);     
   func_lens(l,s,ex);//WD
   dt=double(cadence/1.0); 
   ntran= double(d.Tobs/l.period);    
      
      
      
   cout<<"**********************************************************"<<endl;
   cout<<"One source star is chosen !!!!! "<<endl; 
   cout<<"mass:  "<<s.mass<<"\t Mab:  "<<s.Mab<<"\t blend:  "<<s.blend<<endl;
   cout<<"Tstar: "<<s.Tstar<<"\tRstar:  "<<s.Rstar<<"\t Ds:  "<<s.Ds<<endl;
   cout<<"Tobs:  "<<d.Tobs<<"\t period:  "<<l.period<<"\t T/P: "<<d.Tobs/l.period<<endl;   
   cout<<"inc:   "<<l.inc<<"\t  teta:  "<<l.tet<<"\t dt:  "<<dt<<endl;
   cout<<"**********************************************************"<<endl;
   
   
   
   FlagD=-1;  
   if(l.period>float(1.0*cadence) and double(360.0/l.period)>1.5){
   
   
   
   if(nsim<0){
   ffd=1;
   sprintf(filenam1,"./files/light/lcF2/%c%c%d%c%d.dat",'L','_',int(nsim),'_',2);
   fild=fopen(filenam1,"w");
   sprintf(filenam2,"./files/light/lcF2/%c%c%d%c%d.dat",'M','_',int(nsim),'_',2);
   film=fopen(filenam2,"w");}
   else   ffd=0; 
   
   
   
   FlagD=1; 
   Amin=1000.0;  Amax=0.0;  error=0.0; xp=-1.0; x1=0.0;  yp=-1.0;   signD=0; 
   d.numt=0;     chi1=0.0;  chi2=0.0;  timee=0.0;   Rho=0.0;  
   impactL=100000000.0; impactO=100000000.0;   Rho=0.0;  RE=0.0; 
   
   
   
   for(d.t=0.0;  d.t<=l.period;  d.t+=dt){
   l.phi=double((d.t-l.tp)*2.0*M_PI/l.period); 
   if(l.ecen<0.01) ksi=l.phi;
   else            ksi=Kepler(l.phi , l.ecen);
   xp=x1; 
   //yp=y1; 
   x0=l.a*(cos(ksi)-l.ecen);//major axis [m]
   y0=l.a*sin(ksi)*sqrt(1.0-l.ecen*l.ecen);//minor axis [m]
   y1=              y0*cos(l.tet)+x0*sin(l.tet);//[m]
   x1= cos(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));//[m]
   z1=-sin(l.inc)*(-y0*sin(l.tet)+x0*cos(l.tet));//[m]
   if(d.t>0.0)  vx=double(x1-xp)/(dt*24.0*3600.0);//[m/s] 
  
    
   proj=double(s.Ds*KP/(s.Ds*KP-x1));///without dimention 
   dis= sqrt(x1*x1 + y1*y1 + z1*z1)+1.0e-50;//[m] 
   disp=sqrt(y1*y1 + z1*z1)+1.0e-50;///meter
   phase=acos(-x1/dis);
   l.RE= sqrt(4.0*G*Msun*l.MBH)*sqrt(fabs(x1)*proj)/velocity+1.0e-50;//meter
   s.ros=fabs(s.Rstar*Rsun*proj/l.RE); 
   u= double(disp/l.RE);//without dimention
   us=double(disp-l.RBH-s.Rstar*Rsun*proj);//[m]
   
   
   

   u2=double(disp/(s.Rstar*Rsun*proj));//[u/rho*]
   u3=double(disp/(l.RBH+s.Rstar*Rsun*proj));//[u/rho*+rbh]
   if(x1<0.0 and u2<impactL){impactL=u2;  Rho=s.ros;   RE=l.RE/l.a;}//Lensing
   if(x1>0.0 and u3<impactO){impactO=u3;  }//Occultation
   //if(d.t>0.0) Delt1=DopplerF(dp, vx, s.Tstar);//Doppler boosting 
   //else        Delt1=0.0;  
   //Delt2=Ellipsoid(s.limb,double(s.Rstar*Rsun/l.a),l.q,l.inc,phase,s.grav);//Ellipsoidal variations 
   Delt1=0.0;
   Delt2=0.0; 
   
   
   

   Astar=1.0;     
   if(x1<-0.001){//Self-lensing
   if(s.ros>100.0){  
   if(u<s.ros) Astar=double(1.0+2.0/s.ros/s.ros);
   else        Astar=double(u*u+2.0)/sqrt(u*u*(u*u+4.0));}
   else{       vbb.a1=s.limb;  Astar=vbb.ESPLMag2(u, s.ros);}}
   
  
   Occul=1.0; 
   if(x1>0.0 ){//Occultation  
   if(us<=0.0){
   if(disp<=fabs(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)-l.RBH)){frac=1.0;  frac0=1.0;}
   else{   
   frac=0.0;  frac0=0.0;    
   for(int i=0; i<nbh; ++i){
   for(int j=0; j<nbh; ++j){
   yb=double(-l.RBH + i*l.stepb);   
   zb=double(-l.RBH + j*l.stepb); 
   zlim=sqrt(l.RBH*l.RBH - yb*yb); 
   if(fabs(zb)<=zlim){
   yc=y1-yb;    
   zc=z1-zb;  
   rstar=sqrt(yc*yc+zc*zc)/(s.Rstar*Rsun*s.Ds*KP/(s.Ds*KP-x1)); 
   frac0+=1.0; 
   if(rstar<=1.0) frac+=1.0;
   }}}}   
   Occul=double(1.0-frac/frac0);}}
   
   
   if((x1<-0.01 and fabs(phase)>M_PI/2.0) or l.RE<-0.001 or s.ros<-0.01 or dis<-0.01 or x1>dis or disp<-0.01 or 
   phase<-0.00001 or phase>float(M_PI*1.01) or disp>dis or fabs(Delt1)>1.0 or fabs(Delt2)>1.0 or s.limb>1.01 or s.grav>1.01 
   or Astar<0.9  or Occul<-0.01 or Occul>1.01 or vx>velocity){ 
   cout<<"ERROx1:  "<<x1<<"\t phase: "<<phase<<"\t ros:  "<<s.ros<<endl;
   cout<<"RE:    "<<RE<<"\t dis:  "<<dis<<"disp:  "<<disp<<endl;
   cout<<"limb: "<<s.limb<<"\t gravityD:  "<<s.grav<<endl; 
   cout<<"Delt1: "<<Delt1<<"\t Delt2:  "<<Delt2<<endl; 
   cout<<"frac:  "<<frac<<"\t frac0:  "<<frac0<<"\t vx:  "<<vx<<endl;
   cout<<"u:  "<<u<<"\t Astar:  "<<Astar<<"\t us:  "<<us<<"\t Occul:  "<<Occul<<endl; 
   int uue; cin>>uue;}
   
   
   
   As=(Astar*(1.0 + Delt1 + Delt2 + Delt1*Delt2)+l.ratio*Occul)/(1.0+l.ratio);  
   As=double(As*s.blend+1.0-s.blend);
   if(As<Amin)  Amin=As;   
   if(As>Amax)  Amax=As;
   maga= s.magb-2.5*log10(As);
   
   
   if(ffd>0) fprintf(film,"%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf %.4lf  %.4lf\n",
   d.t,maga,Astar,Delt1,Delt2,Occul,As,double(x1/l.a),double(phase*RA),double(l.RE/Rsun),s.ros,u,us);
   
   
   
   if(d.t>0 and d.t<=d.Tobs and 
       float((d.t-d.gap[0])*(d.t-d.gap[0]-d.dur))>=0.0 and float((d.t-d.gap[1])*(d.t-d.gap[1]-d.dur))>=0.0 
   and float((d.t-d.gap[2])*(d.t-d.gap[2]-d.dur))>=0.0 and float((d.t-d.gap[3])*(d.t-d.gap[3]-d.dur))>=0.0
   and float((d.t-d.gap[4])*(d.t-d.gap[4]-d.dur))>=0.0 and float((d.t-d.gap[5])*(d.t-d.gap[5]-d.dur))>=0.0
   and float((d.t-d.gap[6])*(d.t-d.gap[6]-d.dur))>=0.0 and float((d.t-d.gap[7])*(d.t-d.gap[7]-d.dur))>=0.0
   and float((d.t-d.gap[8])*(d.t-d.gap[8]-d.dur))>=0.0 and float((d.t-d.gap[9])*(d.t-d.gap[9]-d.dur))>=0.0
   and float((d.t-d.gap[10])*(d.t-d.gap[10]-d.dur))>=0.0 and float((d.t-d.gap[11])*(d.t-d.gap[11]-d.dur))>=0.0
   and float((d.t-d.gap[12])*(d.t-d.gap[12]-d.dur))>=0.0 and float((d.t-d.gap[13])*(d.t-d.gap[13]-d.dur))>=0.0
   and float((d.t-d.gap[14])*(d.t-d.gap[14]-d.dur))>=0.0 and float((d.t-d.gap[15])*(d.t-d.gap[15]-d.dur))>=0.0
   and float((d.t-d.gap[16])*(d.t-d.gap[16]-d.dur))>=0.0 and float((d.t-d.gap[17])*(d.t-d.gap[17]-d.dur))>=0.0
   and float((d.t-d.gap[18])*(d.t-d.gap[18]-d.dur))>=0.0 and float((d.t-d.gap[19])*(d.t-d.gap[19]-d.dur))>=0.0
   and float((d.t-d.gap[20])*(d.t-d.gap[20]-d.dur))>=0.0 and float((d.t-d.gap[21])*(d.t-d.gap[21]-d.dur))>=0.0
   and float((d.t-d.gap[22])*(d.t-d.gap[22]-d.dur))>=0.0 and float((d.t-d.gap[23])*(d.t-d.gap[23]-d.dur))>=0.0
   and float((d.t-d.gap[24])*(d.t-d.gap[24]-d.dur))>=0.0 and float((d.t-d.gap[25])*(d.t-d.gap[25]-d.dur))>=0.0){
   timee+=dt; 
   if(timee>=cadence){
   timee-=cadence; 
   emt= fabs(ErrorTESS(maga));
   ///deltaA=fabs(pow(10.0,-0.4*emt)-1.0)*Astar; 
   maga2=maga+RandN(emt,2.5);
   error+= emt;  
   chi1 += pow((maga2-s.magb)/emt,2.0);
   chi2 += pow((maga2-  maga)/emt,2.0);
   if(ffd>0) fprintf(fild,"%.4lf %.5lf %.7lf\n", d.t, maga2, emt);//, Astar+RandN(deltaA,1.5), deltaA);
   d.numt+=1; }}
   if(int(d.numt)%1000==0) cout<<"dt:  "<<d.t<<"\t numt:  "<<d.numt<<"\t period:  "<<l.period<<endl;}//time loop
   
   
   
   
   if(ffd>0){fclose(film); fclose(fild);} 
   if(fabs(Amax-1.0)>fabs(1.0-Amin)){depth=fabs(Amax-1.0); signD=+1;}//Lensing
   else                             {depth=fabs(1.0-Amin); signD=-1;}//Occultation
   d.dchi=fabs(chi2-chi1)/d.numt;
   error= double(error/d.numt);}
   
   else{d.dchi=0.0;  error=fabs(ErrorTESS(s.magb));  depth=fabs(1.0-1.0/(1.0+l.ratio));   FlagD=0;   impactL=-10.0; 
   impactO=-10.0; RE=-10.0; Rho=-10.0; } 
      
     
   ndtran=double(l.period/cadence); 
   snr=   double(sqrt(ntran*1.0)*depth*1000000.0/s.cdpp); 
   if(snr>5.0 and ntran>0.5 and ndtran>1.0) ndet+=1;  
   

   distr=fopen(filenam0,"a+"); 
   fprintf(distr,
   "%d       %d     %d  %.5lf  %.5lf  "///5
   "%.6lf  %.6lf  %.2lf  %.6lf  %.6lf  %.4lf  %.4lf  %.4lf   %.6lf  %.3lf  %.5lf  %.7lf "//17
   "%.7lf  %.9lf  %e  %.7lf  %.5lf  %.5lf  %.8lf  %.5lf  %.8lf  %.5lf  %e  "//28
   "%.6lf  %.5lf  %.5lf  %e  %e   %e   %e  %.8lf  %.6lf %.8lf  %.8lf  %e  %e\n", //41
   FlagD, nsim, signD, s.lat, s.lon, //5   
   s.Ds, s.mass, s.Tstar, s.Rstar, s.limb, s.Mab, s.Map, s.magb, s.blend, s.nsbl, s.Ai, s.prior,//17
   l.MBH, l.RBH/1000.0, l.lumi, l.age, l.inc*RA, l.tet*RA, l.ecen, l.period,  log10(l.a/(s.Rstar*Rsun)),l.tp, l.ratio,//28
   error, d.dchi, snr,s.cdpp, depth,double(1.0-Amin),double(Amax-1.0),Rho,s.loggs,impactL,impactO,RE,s.lumi);//41
   
   
   fclose(distr);      
   cout<<"=============================================================="<<endl;
   cout<<"FlagD:  "<<FlagD<<"\t nsim:  "<<nsim<<"\t ndet:    "<<ndet<<endl;
   cout<<"l.Mab:  "<<l.Mab<<"\t l.Map: "<<l.Map<<"\t ratio:  "<<l.ratio<<endl;
   cout<<"latit:  "<<s.lat<<"\t longt: "<<s.lon<<endl;
   cout<<"nsim:   "<<nsim<<"\t l.MBH:  "<<l.MBH<<"\t RBH: "<<l.RBH<<endl;
   cout<<"Map:    "<<s.Map<<"\t magb:  "<<s.magb<<"\t l.period:  "<<l.period<<endl;
   cout<<"dchi:   "<<d.dchi<<"\t d.numt: "<<d.numt<<endl;
   cout<<"SNR:    "<<snr<<"\t N_transit:   "<<ntran<<"\t error:  "<<error<<endl;
   cout<<"Amax:   "<<Amax<<"\t Amin:  "<<Amin<<"\t depth :    "<<depth<<endl;
   cout<<"cdpp:   "<<s.cdpp<<"\t period/Tobs: "<<l.period/d.Tobs<<endl;
   cout<<"ros*:  "<<Rho<<"\t impact_Lensing:  "<<impactL<<"\t impact_O:  "<<impactO<<endl; 
   cout<<"RE:  "<<RE<<"\t depth:  "<<depth<<"\t Amax-1.0:  "<<Amax-1.0<<endl;
   cout<<"1-Amin:  "<<1.0-Amin<<"\t SignD:  "<<signD<<endl;
   cout<<"impact:  "<<tan(l.inc)*l.a/(s.Rstar*Rsun)<<endl;
   cout<<"=============================================================="<<endl;
   }
   fclose(distr); 
   fclose(_randStream);
   return(0);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double LimbF(source & s){

    int num=-1,num2=-1; int i;
    double distm, dist;
    if(s.Tstar>s.Tef[Nc-1])  num=int(Nc-1); 
    else if(s.Tstar<s.Tef[0])num=0;  
    else{
    for(int j=1; j<Nc; ++j){
    if(float((s.Tstar-s.Tef[j])*(s.Tstar-s.Tef[j-1]))<=0.0){ num=int(j-1); break;}}}

    distm=100000000.0; 
    for(int j=0; j<20; ++j){
    i=int(num+j);
    if(i>=0 and i<Nc){
    dist=sqrt(pow(s.Tstar-s.Tef[i],2.0)+pow(s.loggs-s.Logg[i],2.0)); 
    if(dist<distm) {distm=dist;  num2=i;} }}
    
    if(num<0 or num2<0 or num>=Nc or num2>=Nc or s.Limb[num2]<0.0 or s.Limb[num2]>1.0 or 
    fabs(s.Tstar-s.Tef[num2])>2500.0){// or fabs(s.Logg[num2]-s.loggs)>2.0){
    cout<<"Error num: "<<num<<"\t  num2:  "<<num2<<endl;  
    cout<<"TeffC:    "<<s.Tef[num2]<<"\tLoggC:  "<<s.Logg[num2]<<endl;
    cout<<"Tstar:   "<<s.Tstar<<"\t logg:  "<<s.loggs<<endl;  int yye;  cin>>yye;}
    
    cout<<"s.Tstar:   "<<s.Tstar<<"\t s.logg:  "<<s.loggs<<endl;
    cout<<"TeffC:     "<<s.Tef[num2]<<"\tLoggC:  "<<s.Logg[num2]<<endl;
    cout<<"limb_Drakening:   "<<s.Limb[num2]<<"\tnum2:  "<<num2<<endl;
    cout<<"***************************"<<endl;
    return(s.Limb[num2]);  
}    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double GravF(source & s){
    int num=-1,num2=-1; int i;
    double distm, dist;  
    if(s.Tstar>s.teff[Ng-1])  num=int(Ng-1); 
    else if(s.Tstar<s.teff[0])num=0;  
    else{
    for(int j=1; j<Ng; ++j){
    if(float((s.Tstar-s.teff[j])*(s.Tstar-s.teff[j-1]))<=0.0){num=int(j-1); break;}}}

    distm=100000000.0; 
    for(int j=0; j<20; ++j){
    i=int(num+j);
    if(i>=0 and i<Ng){
    dist=sqrt(pow(s.Tstar-s.teff[i],2.0)+pow(s.loggs-s.logg[i],2.0)); 
    if(dist<distm){distm=dist;  num2=i;}}}
    
    if(num<0 or num2<0 or num>=Ng or num2>=Ng or s.grad[num2]>1.0 or s.grad[num2]<0.0 or 
    fabs(s.Tstar-s.teff[num2])>2500.0 ){//or fabs(s.logg[num2]-s.loggs)>5.0){
    cout<<"Error num: "<<num<<"\t  num2:  "<<num2<<endl; 
    cout<<"TeffC:    "<<s.teff[num2]<<"\tLoggC:  "<<s.logg[num2]<<endl;
    cout<<"Tstar:   "<<s.Tstar<<"\t logg:  "<<s.loggs<<endl;  int yye;  cin>>yye;}
    
    cout<<"TeffC:    "<<s.teff[num2]<<"\tLoggC: "<<s.logg[num2]<<endl;
    cout<<"s.Tstar:  "<<s.Tstar<<"\t s.logg:    "<<s.loggs<<endl;
    cout<<"Gravity_Dakening:   "<<s.grad[num2]<<"\tnum2:  "<<num2<<endl;
    cout<<"***************************"<<endl;
    return(s.grad[num2]);  
}    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double DopplerF(doppler & dp,  double vx, double Tstar){
   double df0=0.0, df1=0.0, dw;
   for(int i=1; i<Nl;  ++i ){
   dp.waven= dp.wave[i]*(1.0 + vx/velocity);
   dp.Fpl0=Fplank(dp.wave[i], Tstar);
   dp.Fpl1=Fplank(dp.waven  , Tstar);
   dw= double(dp.wave[i]-dp.wave[i-1]); 
   df0 +=double(dp.Fpl0*dp.throu[i]*dw); 
   df1 +=double(dp.Fpl1*dp.throu[i]*dw);}
   return( double(df1-df0)/df0); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double Fplank(double wav, double Ts) {
   wav=wav*pow(10.0,-9.0);
   double con= Hplank*velocity/(KBol*Ts*wav); 
   return(fabs(2.0*Hplank*velocity*velocity*pow(wav,-5.0)/(exp(con)-1.0)) );
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double Ellipsoid(double limb, double ampl, double qm, double inc, double phase, double gc){
   double cosi2=double(cos(inc)*cos(inc)); 
   double li3=double(3.0-limb);   
   double L0, Delf2;  
   L0=fabs(1.0+(15.0+limb)*(1.0+gc)*pow(ampl,3.0)*(2.0+5.0*qm)*(2.0-3.0*cosi2)/(60.0*(3.0-limb))+
       9.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*qm*(8.0-40.0*cosi2+35.0*cosi2*cosi2)/(256.0*(3.0-limb)))+1.0e-50; 
   Delf2=double(15*limb*(2.0+gc)*pow(ampl,4.0)*qm*(4.0*cos(inc)-5.0*pow(cos(inc),3.0))/(32.0*li3))*cos(phase) +
         (-3.0*(15.0+limb)*(1.0+gc)*pow(ampl,3.0)*qm*cosi2/(20.0*li3)-
         15.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*qm*(6.0*cosi2-7.0*cosi2*cosi2)/(64.0*li3))*cos(2.0*phase)+
         (-25.0*limb*(2.0+gc)*pow(ampl,4.0)*qm*pow(cos(inc),3.0)/(32.0*li3))*cos(3.0*phase) +
         (105.0*(1.0-limb)*(3.0+gc)*pow(ampl,5.0)*qm*pow(cos(inc),4.0)/(256.0*li3) )*cos(4.0*phase);
  return(Delf2/L0);      
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double CDPPM(tessfile & t, double Mags){

   double cdpp=-1.0, shib;  

   if(Mags<t.tmag[0])  cdpp=t.cdpp[0];  
   else if(Mags>=t.tmag[Np-1]){
   shib=double(t.cdpp[Np-1]-t.cdpp[Np-2])/(t.tmag[Np-1]-t.tmag[Np-2]);  
   cdpp=double(t.cdpp[Np-1]+shib*(Mags-t.tmag[Np-1]));}
   else{
   for(int i=1; i<Np; ++i){
   if((Mags-t.tmag[i])*(Mags-t.tmag[i-1])<=0.0) {
   shib=double(t.cdpp[i]-t.cdpp[i-1])/(t.tmag[i]-t.tmag[i-1]);  
   cdpp=t.cdpp[i-1]+shib*(Mags-t.tmag[i-1]); 
   break;}}}
   if(cdpp<0.0 or cdpp>100000.0 or cdpp<1.0){
   cout<<"Error   cdpp:  "<<cdpp<<"\t Mags:  "<<Mags<<endl; int uue; cin>>uue;} 
   return(cdpp); 
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double Kepler(double phi, double ecen){  
    double ksi=0;   
    double term, term0;  
    phi=double(phi*RA); 
    while(phi>360.0) phi=phi-360.0; 
    while(phi<0.0)   phi=phi+360.0;      
    if(phi>180)      phi=double(phi-360.0);
    if(phi<-181.0 or phi>181.0){ 
    cout<<"Error :  Phi:  "<<phi<<"\t ecent:  "<<ecen<<endl;   int yye;  cin>>yye;}
    phi=double(phi/RA);
    ksi=phi; 
    for(int i=1; i<NB; ++i){
    term= Bessel(i,i*ecen)*sin(i*phi)*2.0/i;  
    ksi+=term; 
    if(i==1) term0=fabs(term); 
    if(fabs(term)<double(thre*term0) and i>5)  break;}        
    return(ksi); 
}    
///#############################################################################
double ErrorTESS(double maga){
   double emt=-1.0, m;     
   
   if(maga<7.5)       emt=double(0.22*maga-5.850); 
   else if(maga<12.5) emt=double(0.27*maga-6.225);  
   else               emt=double(0.31*maga-6.725);    
   emt=emt+RandN(0.1,3.0);
   if(emt<-5.0) emt=-5.0;  
   emt=pow(10.0,emt);
   if(emt<0.00001 or emt>0.5 or maga<0.0){
   cout<<"Error emt:  "<<emt<<"\t maga:  "<<maga<<endl;}
   return(emt); 
}
///#############################################################################
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] and ds<ex.dis[i+1]){
  F = ex.Extks[i]+(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0 or F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; exit(0); }
  return(F);
}
///#############################################################################
double RandN(double sigma, double nn){
   double rr,f,frand;
   do{
   rr=RandR(-sigma*nn , sigma*nn); ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/(sigma*sigma));
   frand=RandR(0.0 , 1.0);
   }while(frand>f);
   return(rr);
}
///#############################################################################
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
   double sig,Lon,Lat;
   int uue, flag=0;
   if(s.lon<0.0){ sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
   else sig=1.0;
   double delt=fabs(s.lon)-floor(fabs(s.lon));

     
   if(delt>1.0 or delt<0.0){cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
   else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
   else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
   else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
   else               Lon=(floor(fabs(s.lon))+0.75)*sig;
   if(fabs(s.lon)<0.24999999)      Lon=360.00;
   if(fabs(s.lon-360.0)<0.2499999) Lon=360.00;
   //cout<<"s.lat:  "<<s.lon<<"\t s.lat:  "<<s.lat<<endl;
   //cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;


   if(s.lat<0.0) sig=-1.0;
   else sig=1.0;
   delt=fabs(s.lat)-floor(fabs(s.lat));
   if(delt>1.0 or delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
   else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
   else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
   else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
   else                Lat=(floor(fabs(s.lat))+0.75)*sig;
    
   if(fabs(Lon)<0.2499999) Lon=360.00; 
   if(Lat==-0.00)  Lat=0.00;
   if(fabs(Lon)<0.24999999)    Lon=360.00;
   cout<<"Lon:    "<<Lon<<"\t     Lat:  "<<Lat<<endl;
     
   char filename[40];
   FILE *fpd;
   sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',float(Lat),'_',float(Lon) );
   fpd=fopen(filename,"r");
   if(!fpd){
   cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
   //FILE *SD;
   //SD=fopen("./files/Ext/saved_direction.txt","r");
   //for(int i=0; i<64881; ++i) {
   //fscanf(SD,"%lf %lf \n",&latit,&lonti);
   //if(fabs(Lat-latit)<0.1 and fabs(Lon-lonti)<0.1){
   //cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
   //cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
   flag=-1;}
   else{
   flag=1;
   for(int i=0; i<100; ++i){
   fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
   if(ex.dis[i]<0.2  or ex.dis[i]>50.0 or ex.Extks[i]<0.0){
   cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
   cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0;}}
   fclose(fpd);}
   //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl
   return(flag);
}
///==============================================================//
///                                                              //
///                  func lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, extinc & ex)
{    
  int nk; 
  double ftest, fmbh;
  double emax, Roche, color, Av; 
  double extG, extBP, extRP;  


  /*
  double cosi=RandR(cos(20.0*M_PI/180.0) , cos(0.0) );  ///[ cos(20), cos(0)]
  l.inc=acos(cosi);
  if(l.inc<-0.001 or l.inc>double(20.0*M_PI/180.0) ){
  cout<<"Error inc: "<<l.inc<<""<<cosi<<endl;    int uus; cin>>uus;}
  */
  
  l.inc=RandR(0.0, 20.0)*M_PI/180.0;//rad
  l.tet=RandR(0.1,359.9)*M_PI/180.0;//rad
  
 
  do{//non contacting binaries  
  do{
  l.MBH=RandR(0.5,2.9);
  fmbh=exp(-0.5*(l.MBH-1.37)*(l.MBH-1.37)/(0.15*0.15))+0.5*exp(-0.5*(l.MBH-1.73)*(l.MBH-1.73)/(0.3*0.3));
  ftest=RandR(0.0,1.5);   
  }while(ftest>fmbh);
  l.a=RandR(log10(3.0*s.Rstar),log10(1000000.0*s.Rstar));
  l.a=pow(10.0,l.a)*Rsun;//meter 
  l.period=sqrt(4.0*M_PI*M_PI/(G*Msun*(s.mass+l.MBH)))*pow(l.a,1.5)/(3600.0*24.0);//days
  l.q=double(l.MBH/s.mass);
  Roche=double(l.a-l.a*0.49*pow(l.q,2.0/3.0)/(0.6*pow(l.q,2.0/3.0)+log(1.0+pow(l.q,1.0/3.0))));
  if(Roche<-0.01 or l.period<-0.0001){
  cout<<"Error Roche:  "<<Roche/s.Rstar/Rsun<<"\t q:  "<<l.q<<endl; 
  cout<<"l.a:  "<<l.a<<"\t MBH:  "<<l.MBH<<endl; int uue; cin>>uue; }
  }while(Roche<double(s.Rstar*Rsun) );  
  
 

  emax=double(0.8-8.0*exp(-pow(6.0*l.period,0.35)));//Fig 7
  if(emax<0.0)  emax=0.0; 
  if(emax>1.0)  emax=1.0;
  if(emax<0.0 or emax>1.0){cout<<"Error emax:  "<<emax<<"\t period:  "<<l.period<<endl; int yyw; cin>> yyw; } 
  l.ecen=RandR(0.0, emax);
  l.tp=  RandR(0.0,l.period); 
 
 
  if(l.MBH>=0.3 and l.MBH<=1.7) l.RBH=RandR(11.5, 12.5)*1000.0;//in meter
  else if(l.MBH<0.3)            l.RBH=RandR(12.5, 15.0)*1000.0;//in meter 
  else if(l.MBH>1.7)            l.RBH=RandR(10.0, 11.5)*1000.0;// in meter
  l.stepb=double(2.0*l.RBH/nbh/1.0); 


  l.lumi=-1.0; 
  l.age=RandR(2.5,8.0);///Figure 
  if(l.age<5.5){
  if(l.MBH<1.2)      l.lumi=-0.43*(l.age-5.5)+32.0;
  else if(l.MBH<2.0) l.lumi=-0.37*(l.age-5.5)+30.8;
  else               l.lumi=-0.57*(l.age-5.5)+29.5;}
  
  else if(l.age<=6.0){
  if(l.MBH<1.2)      l.lumi=-9.0*(l.age-5.5)+32.0;
  else if(l.MBH<2.0) l.lumi=-6.6*(l.age-5.5)+30.8;
  else               l.lumi=-4.0*(l.age-5.5)+29.5;}
  
   else {
  l.lumi=-6.6*(l.age-5.5)+30.8; }
  
  /*
  nk=-1;
  if(l.age<=l.lage[0])    nk=0;
  else if(l.age>l.lage[nns-1]) nk=nns-2;
  else{
  for(int i=0;i<int(nns-1); ++i){
  if((l.age-l.lage[i])*(l.age-l.lage[i+1])<=0.0){nk=i;  break;}}}
  l.lumi=double(l.llum[nk]+double(l.llum[nk+1]-l.llum[nk])*(l.age-l.lage[nk])/(l.lage[nk+1]-l.lage[nk])); */
  l.lumi= pow(10.0,l.lumi)/LSun; 
  l.ratio=double(l.lumi/s.lumi); 
  
  
  if(l.lumi<0.0 or l.age<2.2 or l.age>9.0 or l.MBH<0.5 or l.MBH>2.9 or l.ecen<0.0 or l.ecen>1.0 
  or l.RBH<10000.0 or l.RBH>15000.0){
  cout<<"Error lum:  "<<l.lumi<<"\t age:  "<<l.age<<"\t MBH:  "<<l.MBH<<endl;
  cout<<"ecentrcity:  "<<l.ecen<<"\t RBH:   "<<l.RBH<<endl;   int uue;  cin>>uue; }
  
  cout<<"ratio:  "<<l.ratio<<"\t lumi_NS:  "<<l.lumi<<"\t lumi_Star:  "<<s.lumi<<endl;
  cout<<"inc:  "<<l.inc<<"\t tet:  "<<l.tet<<"\t MBH:  "<<l.MBH<<endl; 
  cout<<"period:  "<<l.period<<"\t semi:  "<<l.a/Rsun/s.Rstar<<"\t ecen_max:  "<<emax<<endl;
  cout<<"ecen:  "<<l.ecen<<"\t RBH(km):  "<<l.RBH/1000<<"\t tp:   "<<l.tp<<endl;
  cout<<"age:  "<<l.age<<"\t  Roche: "<<Roche/double(s.Rstar*Rsun)<<endl;
  cout<<"The end of func_lens    ***************************"<<endl;
  //int uue; cin>>uue; 
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
   double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
   return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double Fluxlimb(double limb, double rstar){
    return ( double(1.0-limb*(1.0-sqrt(fabs(1.0-rstar*rstar)))) );
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double Bessel(int n,double x)
{
    double j1=0.00000001,tet;
    int kmax=10000;
    for(int k=0; k<kmax; ++k){
    tet=double(k*M_PI/kmax);
    j1+=double(M_PI/kmax/1.0)*cos(n*tet-x*sin(tet)); }
    return(j1/M_PI);
}    
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&










