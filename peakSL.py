from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl
import pylab
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import warnings
warnings.filterwarnings("ignore")
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
################################################################################
G= 6.67430*pow(10.0,-11.0)
velocity= 299792458.0
Msun= 1.989 *pow(10.0,30.0)
Rsun= 6.96340*pow(10.0,8)
KPC=  1000.0 * 3.0857*pow(10.0,16.0)
AU=   1.49*pow(10.0,11)
Mstar=2.0    ##red dwarfs [0.3 Msun, 0.35 Rsun]; Sun [1.0,1.0];  B-star [2.0,1.8];   Subginat [2.0,  3.6]
Rstar=3.6*Rsun
Dist= 1.0*KPC



def Thirdlaw(mblak,semi):
    return(2.0*np.pi*pow(semi,1.5)/np.sqrt(G*Msun*(Mstar+mblak))/(3600.0*24.0) )##days   



nd=int(500)
delf=np.zeros((nd, nd,4))
x=[]
y=[]
xo=[]
yo=[]
left=-1.0
bott=np.log10(5.0)
dx=float(2.7/nd/1.0)
dy=float(3.0/nd/1.0) 
for i in range(nd):##semi major axis [y-axis]
    for j in range(nd):## compact object mass[x-axis] 
        mbh =pow(10.0,left+j*dx)##[Msun]
        semi=pow(10.0,bott+i*dy)*Rsun##[m]
        period=Thirdlaw(mbh , semi)##days
        RE=np.sqrt(4.0*G*Msun*mbh)*np.sqrt(semi*Dist/(Dist+semi))/velocity+1.0e-50;## meter
        rhos= Rstar*Dist/(RE*(Dist+semi))
        delf[i,j,0]=np.log10(rhos)
        delf[i,j,1]=np.log10(float(2.0/rhos/rhos)) 
        delf[i,j,2]=np.log10(period)
        if(abs(period-13.7)<0.1): 
            xo.append(float(2.0/rhos/rhos))
            yo.append(i)
        q=float(mbh/Mstar)
        Roche=semi-semi*0.49*pow(q,2.0/3.0)/(0.6*pow(q,2.0/3.0)+np.log(1.0+pow(q,1.0/3.0)));## meter
        if(Roche<Rstar):  
            delf[i,j,3]=0
            x.append(j)
            y.append(i)
            
        else:  
            delf[i,j,3]=1    
print ("Delta F: ", np.min(xo),   np.max(xo),    np.mean(xo))
plt.cla()
plt.clf()
fig=plt.figure(figsize=(8,6))
plt.plot( delf[:,:,2], delf[:,:,1], "go",  ms=1.5)##, label=r"$\rm{Apparent}~\rm{Magnitude}$")
plt.xlabel(r"$\log_{10}[period]$",  fontsize=18)
plt.ylabel(r"$\log_{10}[Delta F]$",fontsize=18)
#plt.yscale('log')
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
#plt.grid("True")
#plt.grid(linestyle='dashed')
#plt.legend(prop={"size":15}, loc='best')
fig.tight_layout()
fig=plt.gcf()
fig.savefig("./Figs/PDgiant.jpg")        
        
                    
################################################################################
tt=int(9)
v=np.zeros((tt))
ticc=np.array([50.0 , 150.0 , 250.0 , 350.0 , 450.0])
def tickfun(x, start, dd0):
    return(start+x*dd0)     
nam=[r"$\log_{10}[\rho_{\star}]$",  r"$\log_{10}[\Delta F_{\rm{L}}]$",  r"$\log_{10}[T(\rm{days})]$",  r"$Flag$"]       

     
for i in range(3):          
    plt.cla()
    plt.clf()
    fig=plt.figure(figsize=(8,6))
    ax= plt.gca()  
    plt.imshow(delf[:,:,i],cmap='viridis',interpolation='nearest',aspect='equal', origin='lower')
    plt.scatter(x,y,  marker='o',c='w', s=4.0)
    #plt.scatter(xo,yo,marker='o',c='k', s=4.0)
    plt.clim()
    plt.title(str(nam[i])+r"$,~\rm{Sub}-\rm{giant}~\rm{star}$", fontsize=21)
    minn=np.min(delf[:,:,i])
    maxx=np.max(delf[:,:,i])
    step=float((maxx-minn)/(tt-1.0));
    for m in range(tt):
        v[m]=round(float(minn+m*step),1)
    cbar=plt.colorbar(orientation='vertical',shrink=0.9,pad=0.05,ticks=v)
    cbar.ax.tick_params(labelsize=20)
    plt.clim(v[0]-0.005*step,v[tt-1]+0.005*step)
    contours = plt.contour(delf[:,:,2], 10,  colors='black')
    plt.clabel(contours, inline=1, fontsize=16)
    plt.xticks(fontsize=20, rotation=0)
    plt.yticks(fontsize=20, rotation=0)
    plt.xlim(0.0,nd)
    plt.ylim(0.0,nd)
    ax.set_xticks(ticc,labels=[round(j,1) for j in tickfun(ticc, left, dx) ]  )
    ax.set_yticks(ticc,labels=[round(j,1) for j in tickfun(ticc, bott, dy) ]  )
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel(r"$\log_{10}[M_{\rm{c}}(M_{\odot})]$",fontsize=21,labelpad=0.05)
    plt.ylabel(r"$\log_{10}[a(R_{\odot})]$",         fontsize=21,labelpad=0.05)
    fig=plt.gcf()
    fig.tight_layout(pad=0.1)
    fig.savefig("./Figs/MapSG{0:d}.jpg".format(i),dpi=200)        

        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
