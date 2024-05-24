from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl
import pylab
rcParams["font.size"] = 18
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from matplotlib.ticker import FormatStrFormatter
from matplotlib.gridspec import GridSpec
################################################################################
sector=int(0)
nc=int(49) ## No.  parameters

f1=open('./files_Figs/D{0:d}_1.dat'.format(sector) )
N0=int(len(f1.readlines()))
Tobs=27.4;  


################################################################################

for ll in range(1): 
    f1=open("./files_Figs/D{0:d}_1.dat".format(sector),"r")
    nf= sum(1 for line in f1)  
    par=np.zeros((nf,nc))
    par=np.loadtxt("./files_Figs/D{0:d}_1.dat".format(sector)) 
    for i in range(nf):
        #par[i,36]=np.log10(par[i,36])## log10(F)
        #par[i,44]=np.log10(par[i,44])## log10(rho*)
        #par[i,46]=np.log10(par[i,46])## log10(u0/rho*) Log(impactL)
        #par[i,47]=np.log10(par[i,47])## log10(dp/Rstar) Log(impactO)
        FlagD, nsim, signD, lat, lon  =int(par[i,0]),int(par[i,1]),par[i,2],par[i,3],par[i,4]
        Ds, mass, Tstar, Rstar, limb  =par[i,5],par[i,6],par[i,7],par[i,8],par[i,9]
        Mab, Map,magb,fb,Nbl,Ext,prior= par[i,10],par[i,11],par[i,12],par[i,13],par[i,14],par[i,15],par[i,16]
        Ml, tefl, Rl, loggl, Mapl, Mabl=par[i,17],par[i,18], par[i,19],par[i,20],par[i,21],par[i,22]
        magG,magBP,magRP,MG,MBP,MRP,numl =par[i,23],par[i,24], par[i,25],par[i,26],par[i,27],par[i,28], par[i,29]
        inc,teta,ecen,period,lsemi,tp,ratio=par[i,30],par[i,31],par[i,32],par[i,33],par[i,34],par[i,35],par[i,36]
        error, dchi, snr,cdpp       =par[i,37], par[i,38], par[i,39], par[i,40]
        depth,depthO,depthL,rho=par[i,41],par[i,42],par[i,43],par[i,44]
        loggs, impactL, impactO,RE= par[i,45],par[i,46],par[i,47],par[i,48]  
        if(par[i,2]<0): par[i,2]=0 ## Occultation
        par[i,2]=int(par[i,2])
        sigA= abs(1.0-pow(10.0,-0.4*error))
        ######################################################################## 
        Tobs= float(13.7*2.0)
        ntran=float(Tobs/period)
        snr=  float(np.sqrt(ntran*1.0)*depth*1000000.0/cdpp)+0.01
        stat=r"$\rm{Not}-\rm{Detected}$"
        col='r'
        if(snr>5.0 and ntran>2.0): 
            stat=r"$\rm{Detected}$"
            col='g'
        ######################################################################## 
        if(FlagD>0): 
            nd=-1;  nm=-1;  
            try: 
                f1=open('./files_Figs/L_{0:d}_0.dat'.format(nsim) )
                nd=int(len(f1.readlines()))
                print(nd)
            except: 
                print("file does not exist",  nsim)    
            try:
                f2=open('./files_Figs/M_{0:d}_0.dat'.format(nsim) )
                nm=int(len(f2.readlines()))  
                print(nm)
            except: 
                print("file does not exist",  nsim)        
            print("nsim,   nd,  nm:  ", nsim, nd, nm)   
             #if(j==0):  stat=r"$\rm{Detected}$"; col='g'
                #elif(j==0):stat=r"$\rm{Not}-\rm{detected}$"; col='r'  
            if(nd>1 and nm>1): 
                dat=np.zeros((nd,3))
                dat=np.loadtxt('./files_Figs/L_{0:d}_0.dat'.format(nsim)) 
                for ii in range(nd): 
                    dat[ii,1]=dat[ii,1]+np.random.normal(0.0,dat[ii,2],1)*0.3
                
                mod=np.zeros((nm,13))
                mod=np.loadtxt('./files_Figs/M_{0:d}_0.dat'.format(nsim)) 
                ymax1,ymin1= np.max(mod[:,1]), np.min(mod[:,1])
                ymax2,ymin2= np.max(dat[:,1]), np.min(dat[:,1])
                ymax ,ymin= np.max(np.array([ymax1,ymax2])) , np.min(np.array([ymin1, ymin2]))
                plt.clf()
                plt.cla()
                fig=plt.figure(figsize=(8,6))
                ax1=fig.add_subplot(111)
                plt.errorbar(dat[:,0],dat[:,1],yerr=dat[:,2],fmt=".",markersize=2.,color='m',ecolor='gray',elinewidth=0.02, capsize=0,alpha=0.42,label=r"$\rm{TESS}~\rm{Data}$")
                plt.plot(mod[:,0],mod[:,1],'k-',label=r"$\rm{Model}~\rm{Light}~\rm{Curve}$", lw=1.7,alpha=1.0)
                plt.ylabel(r"$\rm{TESS}-\rm{magnitude}$",fontsize=19)
                plt.xlabel(r"$\rm{time}(\rm{days})$",fontsize=19)
                plt.title(
                r"$M_{\rm{WD}}(M_{\odot})=$"+'{0:.1f}'.format(Ml)+
                r"$,~T(\rm{days})=$"+'{0:.2f}'.format(period)+
                r"$,~\log_{10}[\mathcal{F}]=$"+'{0:.2f}'.format(np.log10(ratio))+
                r"$,~\rho_{\star}=$"+'{0:.1f}'.format(rho)+ "\n"+ 
                r"$b(R_{\star})=$"+'{0:.2f}'.format(impactL)+   
                r"$,~i(\rm{deg})=$"+'{0:.2f}'.format(inc)+    
                r"$,~\log_{10}[\Delta F]=$"+'{0:.2f}'.format(np.log10(depth))+     
                r"$,~\rm{SNR}=$"+'{0:.2f}'.format(snr),fontsize=17.0,color='k')
                pylab.ylim([ymin, ymax])
                plt.xlim([0.0,Tobs-1.0])
                plt.xticks(fontsize=18, rotation=0)
                plt.yticks(fontsize=18, rotation=0)
                plt.gca().invert_yaxis()
                ax1.grid("True")
                ax1.grid(linestyle='dashed')
                ax1.legend(title=stat,prop={"size":15.})
                fig=plt.gcf()
                fig.tight_layout(pad=0.2)  
                fig.savefig("./Figs/lightW/LightW{0:d}.jpg".format(nsim),dpi=200)
                ################################################################
                '''
                #d.t,maga,Astar,Delt1,Delt2,Occul,As,double(x1/l.a),double(phase*RA),double(l.RE/Rsun),s.ros,u,us);
                # 0    1   2      3     4     5    6      7               8                   9         10   11 12
                plt.cla() 
                plt.clf()
                plt.figure(figsize=(8, 6))
                plt.plot(mod[:,0]/period,mod[:, 7]/np.max(mod[:,7]), color="r",label="x1(LoS)",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:, 8]/np.max(mod[:,8]), color="b",label="Phase",  lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:, 9]/np.max(mod[:,9]), color="g",label="RE",     lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,10]/np.max(mod[:,10]),color="m",label="Rostar", lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,11]/np.max(mod[:,11]),color="k",label="u",      lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,12]/np.max(mod[:,12]),color="c",label="us",     lw=1.2, alpha=0.95)
                plt.xticks(fontsize=17, rotation=0)
                plt.yticks(fontsize=17, rotation=0)
                plt.xlabel(r"$time/Period$", fontsize=18)
                plt.grid("True")
                plt.legend()
                fig=plt.gcf()
                fig.savefig("./Figs/xyz{0:d}.jpg".format(nsim), dpi=200)
                ################################################################
                plt.cla()
                plt.clf()
                plt.figure(figsize=(8, 6))
                plt.plot(mod[:,0]/period,mod[:,2],    color="r",label="Magnification",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,3]+1.0,color="b",label="Doppler Boost",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,4]+1.0,color="g",label="Ellipsoidal  ",lw=1.2, alpha=0.95)
                plt.plot(mod[:,0]/period,(mod[:,5]*ratio+1.0)/(1.0+ratio),color="m",label="Occultation",lw=1.2,alpha=0.95)
                plt.plot(mod[:,0]/period,mod[:,6],    color="k",label="All effects",  lw=1.2, alpha=0.95)
                plt.xticks(fontsize=17, rotation=0)
                plt.yticks(fontsize=17, rotation=0)
                plt.xlabel(r"$time/Period$", fontsize=18)
                plt.grid("True")
                plt.legend()
                fig=plt.gcf()
                fig.savefig("./Figs/Astar{0:d}.jpg".format(nsim), dpi=200)
                ''' 

























