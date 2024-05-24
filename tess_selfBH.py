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


Rsun=6.9634*np.power(10.0,8.0)###; ///solar radius [meter]
G=6.67384*np.power(10.,-11.0);##// in [m^3/s^2*kg].
Msun=1.98892*np.power(10.,30.0);## //in [kg].
velocity=299792458.0;##//velosity of light  m/s
AU=1.495978707*np.power(10.0,11.0);

nam0=['Flagd', 'nsim', r"$SignD$",r"$\rm{Galactic}~\rm{Latitude}(\rm{deg})$", r"$\rm{Galactic}~\rm{Longitude}(\rm{deg})$",r"$D_{\rm{s}}(\rm{kpc})$",r"$M_{\star}(M_{\odot})$", r"$T_{\star}(\rm{K})$",r"$R_{\star}(R_{\odot})$",r"$\Gamma$", r"$M_{\star,T}(\rm{mag})$", r"$m_{\star,T}(\rm{mag})$", r"$m_{\rm{base},T}(\rm{mag})$", r"$f_{\rm{b}}$", r"$N_{\rm{b}}$", r"$\rm{Extinction}_{T}$", r"$\rm{Periority}$", r"$M_{\rm{c}}(M_{\odot})$", r"$R_{\rm{c}}(R_{\odot})$",r"$i(\rm{deg})$", r"$\theta(\rm{deg})$", r"$\epsilon$", r"$T(\rm{days})$", r"$\log_{10}[a/R_{\star}]$", r"$t_{\rm{p}}$", r"$\sigma_{\rm{m}}(\rm{mag})$", r"$\Delta \chi^{2}$", r"$\rm{SNR}$", r"$\rm{CDPP}$",r"$\rm{Depth}$",r"$\rm{Depth}_{\rm{O}}$",r"$\rm{Depth}_{\rm{L}}$",r"$\log_{10}[\rho_{\star}]$", r"$\log_{10}[g_{\star}(cm/s)]$", r"$\log_{10}[b_{\rm{l}}/R_{\star}]$",r"$\log_{10}[b_{\rm{o}}/(R_{\star}+R_{\rm{c}})]$", r"$R_{\rm{E}}$"]
################################################################################
sector=int(0)
ns=int(13)
nv=int(11) 
nc=int(37)
k1=int(0)  
k2=np.zeros((ns))



f1=open('./files/B{0:d}_tot.dat'.format(sector) )
N0=int(len(f1.readlines()))
fij=open("./Figs/HistB/RESBH.dat","w")
fij.close(); 


deep=np.zeros((ns))
pars=np.zeros((N0,nc))
pard=np.zeros((N0,nc,ns))
Effi=np.zeros((nv,12,ns+1))
xar=np.zeros((nv,12))
Dep=[]


xar[:,0]=np.array([0.0, 0.04,  0.08, 0.13, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.65])## Ds[0.0,1.0] kpc
xar[:,1]=np.array([0.1, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2., 2.2, 2.9, 4.9])## Rs[0.1,3.2]  Rsun
xar[:,2]=np.array([0.002, 0.004, 0.006,  0.008, 0.01,0.015, 0.02, 0.025, 0.03,0.04, 0.05])##log[Periority] [-2.8, -1.3]
xar[:,3]=np.array([2.0,7.0,12.0,17.0,22.,26.0,30.0,34.0,37.,40.0, 48.0])## M_{c}(Msun) [0.2, 1.4]
xar[:,4]=np.array([0.0,2.0, 4.0, 6.0,8.0,10.0,12.0,14.0,16.0,19.0,20.0]) # i(deg) [0.0,22.0]
xar[:,5]=np.array([-4.0,-3.7,-3.1,-2.4,-2.0,-1.5,-1.0, -0.5, -0.2, -0.1, -0.01])## log epsilon [0.0,1.0]
xar[:,6]=np.array([-1.5, -1.0, -0.7,-0.4, -0.1, 0.2, 0.4,0.6,1.0,1.2,1.5])## log10[Period (days)]
xar[:,7]=np.array([0.5,0.6, 0.7,0.8, 0.9, 1.05, 1.15, 1.35, 1.5, 1.75, 2.9]) ### log10[a/Rstar]
xar[:,8]=np.array([-1.5, -1.2, -0.75, -0.45, -0.25, 0.05, 0.5, 1.0 ,1.5, 2.0 ,2.5])##\log10(ImpactL)[-2.0,2.0]
xar[:,9]=np.array([-1.5, -1.2, -0.75, -0.45, -0.25, 0.05, 0.5, 1.0 ,1.5, 2.0 ,2.5])##\log10(ImpactO)[-2.0,2.0]
xar[:,10]=np.array([-1.7,-1.,-0.5,0.0,0.3,0.8,1.1,1.5, 1.7, 1.9, 2.5])## log10(rho*) [-1.0,2.5]
xar[:,11]=np.array([0.1, 0.4, 0.6, 0.8, 1.2, 1.5, 1.7, 2.2, 2.4, 2.6, 3.1])## log10(rho*) [-1.0,2.5]
################################################################################

def func(para,cc): 
    i1=-1
    for j in range(nv-1):  
        if(float((para-xar[j,cc])*(para-xar[j+1,cc]))<0.0  or  para==xar[j,cc]):
            i1=j;  
            break
    if(i1<0 and para<xar[0,cc]): i1=0
    elif(i1<0 and para>=xar[nv-1,cc] ):  i1=int(nv-1) 
    elif(i1<0): 
        print("There is a problem i1<0:   ", i1,   para,   xar[:,cc]   )
    return(i1)        

################################################################################

for ll in range(1): 
    f1=open("./files/B{0:d}_tot.dat".format(sector),"r")
    nf= sum(1 for line in f1)  
    par=np.zeros((nf,nc))
    par=np.loadtxt("./files/B{0:d}_tot.dat".format(sector)) 
    for i in range(nf):
        par[i,32]=np.log10(par[i,32])## log10(rho*)
        par[i,34]=np.log10(par[i,34])## log10(bl/Rstar) Log(impactL)
        par[i,35]=np.log10(par[i,35])## log10(bo/Rstar) Log(impactO)
        FlagD, nsim, signD, lat, lon   =int(par[i,0]),int(par[i,1]),par[i,2],par[i,3], par[i,4]
        Ds, mass, Tstar, Rstar, limb  =par[i,5],par[i,6],par[i,7],par[i,8],par[i,9]
        Mab,Map,magb,fb,Nbl,Ext,prior= par[i,10],par[i,11],par[i,12],par[i,13],par[i,14],par[i,15],par[i,16]
        Ml,Rl,inc,teta,ecen,period,lsemi,tp =par[i,17],par[i,18], par[i,19],par[i,20],par[i,21],par[i,22],par[i,23],par[i,24]
        error, dchi,snr, cdpp= par[i,25],par[i,26],par[i,27],par[i,28]
        depth,depthO, depthL,lrho=par[i,29],par[i,30],par[i,31],par[i,32]
        logg,limpactL, limpactO, RE= par[i,33],par[i,34],par[i,35],par[i,36]
        if(par[i,2]<0): par[i,2]=0 ## Occultation
        par[i,2]=int(par[i,2])
        sigA= abs(1.0-pow(10.0,-0.4*error))
        ########################################################################
        ini=np.array([Ds,Rstar,prior,Ml,inc,np.log10(ecen),np.log10(period),lsemi,limpactL,limpactO,lrho , mass ])
        if(FlagD>0.5 and period<float(360.0/2.0) ): 
            Dep.append(np.log10(depth/sigA))## depth normalized to Error       
            pars[k1,:]=par[i,:]
            k1+=1 
            for j in range(ns): 
                Tobs= float(13.7*2.0*(j+1.0))
                ntran=float(Tobs/period)
                snr=  float(np.sqrt(ntran*1.0)*depth*1000000.0/cdpp)+0.01
                par[i,27]=snr
                det=0
                if(snr>5.0 and ntran>2.0): 
                    det=1
                    pard[int(k2[j]),:,j]=par[i,:]  
                    k2[j]+=1
                    if(signD>0): 
                        deep[j]+=1.0## number of lensing-induced impacts
                for k in range(12): 
                    kk=int(func(ini[k],k))
                    Effi[kk,k,0]+=1.0
                    if(det>0):  
                        Effi[kk,k,j+1]+=1.0 
                       
                      
print("Probability of FlaD>0",    float(k1*100.0/nf) )                   
print("Range of log10(Rho*):  ",     np.min(pars[:k1,32]),   np.max(pars[:k1,32])  )                
print("Range of log10(ImpactL):  ",  np.min(pars[:k1,34]),   np.max(pars[:k1,34])  )      
print("Range of log10(ImpactO):  ",  np.min(pars[:k1,35]),   np.max(pars[:k1,35])  )     
print("Range of priority:  ",        np.min(pars[:k1,16]),   np.max(pars[:k1,16])  )    
print("Range of Rstar:  ",           np.min(pars[:k1,8]),   np.max(pars[:k1,8])  )      
################################################################################ 
fij=open("./Figs/HistB/RESBH.dat","a+")  
for i in range(ns): 
    k=int(k2[i])  
    efL=float(deep[i]*100.0/k)
    efO=float(100.0-efL)     
    
    b1= np.mean( np.abs(np.power(10.0,pard[:k,34,i])) + np.abs(np.power(10.0,pard[:k,34,i]))   )*0.5
    
    row=np.array([ i+1,float(27.4*(i+1)),np.mean(pard[:k,17,i]), np.mean(pard[:k,22,i]), np.mean(pard[:k,23,i]),np.mean(pard[:k,21,i]), np.mean(pard[:k,19,i]), np.mean(pard[:k,5,i]), np.mean(pard[:k,27,i]),np.mean(pard[:k,32,i]),b1,float(k*100.0/k1) ])
  
    np.savetxt(fij,row.reshape((-1,12)),fmt="$%d$ & $%.1f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$")  
    print("No_sector: ",  i ,   "fraction of detection:  ", round(k*100.0/k1,2) ,  k ,  k1)         
fij.write("****************************************************"+"\n")    
fij.close() 
## Ns, Tobs, Ml, T, log(a/Rstar), eccentricity, i, Ds, snr, Log(rho*), Log(u0), efficiency, (12)    

################################################################################    
w=int(k2[0])
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
plt.plot( pars[:k1,12],pars[:k1,28], "go",  ms=3.5, label=r"$\rm{Apparent}~\rm{Magnitude}$")
plt.xlabel(r"$\rm{m}_{T}(\rm{mag})$",  fontsize=18)
plt.ylabel(r"$CDPP$",fontsize=18)
plt.yscale('log')
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
plt.legend(prop={"size":15}, loc='best')
fig.tight_layout()
fig=plt.gcf()
fig.savefig("./Figs/HistB/cdppMag.jpg")

    
################################################################################ 
col=["k", "b", "r", "g", "c", "y", "m"]
for i in range(nc):
    plt.clf()
    plt.cla()
    fig=plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(pars[:k1,i],30,histtype='bar',ec='darkgreen',facecolor='green', alpha=0.4,rwidth=1.0)
    for j in range(int(ns/2)): 
        k=0
        if(j>6):  k=int(j%7)
        else:     k=j 
        plt.hist(pard[:int(k2[int(j*2)]),i,int(j*2)],30,histtype='step',color=col[k], alpha=1.0,lw=float(0.9+j*2*0.08),ls='--') 
    y_vals =ax.get_yticks()
    ax.set_yticks(y_vals)
    ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/k1))) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
    ax.set_xlabel(str(nam0[i]),fontsize=18,labelpad=0.1)
    plt.xticks(fontsize=18, rotation=0)
    plt.yticks(fontsize=18, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.tight_layout()
    fig.savefig("./Figs/HistB/HistoB{0:d}.jpg".format(i),dpi=200)
print("****  All histo_simulated_parameters are plotted ********************", i)

    
################################################################################
plt.cla()
plt.clf()
fig = plt.figure(figsize=(8,6))
gs = GridSpec(4,4)
ax_sc = fig.add_subplot(gs[1:4,0:3])
ax_hy = fig.add_subplot(gs[0:1,0:3], sharex=ax_sc)
ax_hx = fig.add_subplot(gs[1:4,3], sharey=ax_sc)

ax_sc.scatter(np.log10(pars[:k1,22]), pars[:k1,34], marker='o',c='k',s=1.0,label=r"$\rm{All}~\rm{Events}$")
ax_sc.scatter(np.log10(pard[:w,22,0]),pard[:w,34,0],marker='*',c='r',s=3.0,label=r"$\rm{Detectable}~\rm{Events}$")
ax_sc.set_xlabel(r"$\log_{10}[T(\rm{days})]$",  fontsize=21)
ax_sc.set_ylabel(r"$\log_{10}[b/R_{\star}]$",fontsize=21)
plt.subplots_adjust(wspace=-7.0)
plt.subplots_adjust(hspace=-7.0)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
ax_sc.set_xlim([-2.0,2.5])
ax_sc.set_ylim([-2.5,2.9])
#ax_sc.set_xticks([-2.0,2.4])
ax_sc.legend(prop={"size":15}, loc=3, fancybox=True, shadow=True)
ax_sc.legend(title=r"$\rm{BHMS}~\rm{Binary}$",prop={"size":15}, loc=3)



ax_hy.hist(np.log10(pars[:k1,22]),30,density=False,histtype='bar',ec='k',facecolor='k', alpha=0.4,rwidth=1.0)
ax_hy.hist(np.log10(pard[:w,22,0]),30,density=False,histtype='bar',ec='r',facecolor='r', alpha=0.4,rwidth=1.0)
y_vals =ax_hy.get_yticks()
ax_hy.set_yticks(y_vals)
ax_hy.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/k1))) for x in y_vals]) 
ax_hy.set_yscale('log')
#ax_hy.set_xticks([])
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)




ax_hx.hist(pars[:k1,34],30,density=False,histtype='bar',ec='k',facecolor='k', alpha=0.4,rwidth=1,orientation ='horizontal')
ax_hx.hist(pard[:w,34,0],30,density=False,histtype='bar',ec='r',facecolor='r',alpha=0.4,rwidth=1,orientation ='horizontal')
y_vals =ax_hx.get_xticks()
ax_hx.set_xticks(y_vals)
ax_hx.set_xticklabels(['{:.2f}'.format(float(1.0*x*(1.0/k1))) for x in y_vals]) 
ax_hx.set_xscale('log')
#ax_hx.set_yticks([])
plt.subplots_adjust(wspace=-2.0)
plt.subplots_adjust(hspace=-2.0)


fig.tight_layout(pad=0.1)
fig=plt.gcf()
fig.savefig("./Figs/HistB/MC_Impact.jpg")

################################################################################

for i in range(nc):
    plt.cla()
    plt.clf()
    fig5=plt.figure(figsize=(8,6))
    if(i==22):
        plt.scatter(np.log10(pars[:k1,i]),pars[:k1,29],  marker='o',c='k',s=1.5,label=r"$\rm{All}~\rm{Events}$")
        plt.scatter(np.log10(pard[:w,i,0]),pard[:w,29,0],marker='o',c='g',s=1.5,label=r"$\rm{Detectable}~\rm{Events}$")
        plt.xlim(-2.0,2.5)
    else: 
        plt.scatter(pars[:k1,i],pars[:k1,29],  marker='o',c='k',s=1.5,label=r"$\rm{All}~\rm{Events}$")
        plt.scatter(pard[:w,i,0],pard[:w,29,0],marker='o',c='g',s=1.5,label=r"$\rm{Detectable}~\rm{Events}$")
    plt.xlabel(str(nam0[i]),   fontsize=18)
    plt.ylabel(r"$\rm{Depth}$",fontsize=18)
    plt.yscale('log')
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    plt.legend(prop={"size":15}, loc="lower left",fancybox=True, shadow=True)
    fig5.tight_layout()
    fig5=plt.gcf()
    fig5.savefig("./Figs/HistB/Depth{0:d}.jpg".format(i) )
    
################################################################################
for i in range(12):## parameters
    for j in range(ns):## sector number 
        Effi[:,i,j+1]=(Effi[:,i,j+1]*100.0/(Effi[:,i,0]+0.00000016324512736) )
   
nam3=[r"$D_{\rm{l}}(\rm{kpc})$", r"$R_{\star}(R_{\odot})$", r"$\rm{Periority}$", r"$M_{\rm{BH}}(M_{\odot})$", r"$i(\rm{deg})$", r"$\log_{10}[\epsilon]$", r"$\log_{10}[T (\rm{days})]$", r"$\log_{10}[a/R_{\star}]$",r"$\log_{10}[b/R_{\star}]$",r"$\log_{10}[b_{\rm{o}}/(R_{\star}+R_{\rm{c}})]$",r"$\log_{10}[\rho_{\star}]$", r"$M_{\star}(M_{\odot})$"]##11

for i in range(12): 
    plt.clf()
    plt.cla()
    plt.figure(figsize=(8,6))
    fig, ax = plt.subplots()
    plt.step(xar[:,i],Effi[:,i,1],where='mid',lw=2.1,linestyle='-',color="gray",label=r"$T_{\rm{obs}}(\rm{days})=27.4$")
    plt.step(xar[:,i],Effi[:,i,6],where='mid',lw=2.1,linestyle='--',color="r",  label=r"$T_{\rm{obs}}(\rm{days})=164$")
    plt.step(xar[:,i],Effi[:,i,13],where='mid',lw=2.1,linestyle='-.',color="g", label=r"$T_{\rm{obs}}(\rm{days})=356.2$")
    plt.scatter(xar[:,i],Effi[:,i,1],marker="o",facecolors='k',        edgecolors='gray',s=36.0)
    plt.scatter(xar[:,i],Effi[:,i,6],marker="o",facecolors='darkred',  edgecolors='r',   s=36.0)
    plt.scatter(xar[:,i],Effi[:,i,13],marker="o",facecolors='darkgreen',edgecolors='g',   s=36.0)
    #plt.xlim([ xar[0,i], xar[nv-1,i] ])
    #plt.ylim([-2.0, np.log10(4.0)])
    #plt.yscale('log')
    plt.ylabel(r"$\log_{10}[\varepsilon(\%)]$",fontsize=21,labelpad=0)
    plt.xlabel(str(nam3[i]),fontsize=22,labelpad=0.1)
    plt.xticks(fontsize=21, rotation=0)
    plt.yticks(fontsize=21, rotation=0)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.grid("True")
    plt.grid(linestyle='dashed')
    if(i==0): 
        plt.legend()
        plt.legend(loc='best',fancybox=True, shadow=True)
        plt.legend(prop={"size":21})
    fig.tight_layout()    
    fig=plt.gcf()
    fig.savefig("./Figs/HistB/EffiB{0:d}.jpg".format(i),dpi=200)
print("**************************** parameter Efficiency plot", i)

################################################################################


plt.clf()
plt.cla()
fig=plt.figure(figsize=(8,6))
ax= plt.gca()              
plt.hist(Dep,30,histtype='bar',ec='darkgreen',facecolor='green', alpha=0.4,rwidth=1.0)
y_vals =ax.get_yticks()
ax.set_yticks(y_vals)
ax.set_yticklabels(['{:.2f}'.format(float(1.0*x*(1.0/k1))) for x in y_vals]) 
y_vals = ax.get_yticks()
plt.ylim([np.min(y_vals), np.max(y_vals)])
ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=18,labelpad=0.1)
ax.set_xlabel(r"$\log_{10}[\Delta F/\sigma_{\rm{m}}]$",fontsize=18,labelpad=0.1)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.tight_layout()
fig.savefig("./Figs/HistB/HistDepthlog.jpg",dpi=200)


################################################################################





























