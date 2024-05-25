import numpy as np 
#import VBBinaryLensingLibrary as vb
#VBB=vb.VBBinaryLensing()
#VBB.Tol=1.0e-3; ### accuracy
#VBB.SetLDprofile(VBB.LDlinear);
#VBB.LoadESPLTable("./ESPL.tbl"); 
G= 6.67430*pow(10.0,-11.0)
c= 299792458.0
Msun= 1.989 *pow(10.0,30.0)
Rsun= 6.96340*pow(10.0,8)
KPC= 1000.0 * 3.0857*pow(10.0,16.0)
AU=1.49*pow(10.0,11)

print("ralative area:  ",  pow(2.0*G*Msun*45.0/c/c/Rsun/3.9,2.0) )
#input("Enter a number ")
################################################################################
Period=100.0##11.393423143561519# [days]  from Gaia
def Thirdlaw( mstar, mblak):
    return(np.power(G*Msun*(mstar+mblak)*Period*Period*24.0*24.0*3600.0*3600.0/(4.0*np.pi*np.pi),1.0/3.0) )    

print("Period:  ",    np.sqrt(4.0*np.pi*np.pi/(G*Msun*(1.0+1.4)))*pow(100.0*Rsun,1.5)/(3600.0*24.0)/365.0  )



fij=open("./table4.dat","w")
fij.close(); 

Nt=np.array([1031.0,258.0,38.0,9.0,5.0,5.0,3.0,2.0,1.0,1.0,1.0,16.0,13.0])*1000.0
ef=np.zeros((13,3))
ef[:,0]=np.array([14.64, 16.71, 18.00, 19.01, 19.75, 20.29, 20.76, 21.20, 21.54, 21.84, 22.08, 22.31, 22.55])*0.01##WD
ef[:,1]=np.array([13.63, 15.36, 16.36, 16.96, 17.43, 17.75, 18.04, 18.19, 18.41, 18.63, 18.86, 19.02, 19.14])*0.01##NS
ef[:,2]=np.array([10.57, 11.17, 11.54, 11.86, 12.13, 12.38, 12.57, 12.77, 12.92, 13.08, 13.23, 13.36, 13.47])*0.01##BH

f0=np.array([0.321,  0.332,  0.392])
f1=np.array([0.001, 0.0004,  0.000004])
f2=float(20.0/90.0)
Nc=np.zeros((14))

for i in range (3): 
    fij=open("./table4.dat","a+")
    Nc[:13]=ef[:,i]*Nt*f0[i]*f1[i]*f2
    Nc[13]=np.sum(Nc[:13])
    np.savetxt(fij,Nc.reshape((-1, 14)), fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$")
    fij.close()




input("Enter a number ")
RE= np.sqrt(4.0*G*Msun*Thirdlaw(1.0,0.6))/c
print ( Rsun/RE,   2.0*RE*RE/Rsun/Rsun )
input("Enter a number ")


##TIC=  153123772
ecn=0.83
mstar=2.08
mbh=np.arange(30.0, 40.0, 0.5)
Dl=float(1.0/0.978796586561136)*KPC  ##[m]

Rstar=3.8

semi=Thirdlaw(mstar,43.5)*(1.0+ecn)
xls=Dl/(Dl+semi)
RE=np.sqrt(4.0*G*43.5*Msun*(Dl+semi)*xls*(1.0-xls))/c
ros=Rstar*Rsun*Dl/RE/(Dl+semi)
print  ("parameter barkhord: ",semi, RE, ros,   semi*np.tan(1.5*np.pi/180.0)/RE/ros)
input("Enter a number ")

#VBB.ESPLMag2(u[k],rho[k]) 

#ros0=np.sqrt(2.0/0.0055)
#tan1=np.tan(1.5*np.pi/180.0)# tan(1 degree)
#print("ros0:  ",  ros0)
for i in range( len(mbh) ):  
    aa= Thirdlaw( mstar,mbh[i])
    q=mbh[i]/mstar
    Roche=0.49*pow(q,2.0/3.0)/(0.6*pow(q,2./3.)+np.log(1.0+pow(q,1.0/3.0)))
    print("MBH, ROCHE: ",  round(mbh[i],2),  round(q, 2) ,    round(Roche*aa*(0.2)/Rstar/Rsun,2),      round(Roche,2) )
    
    semi=aa*(1.0+ecn)
    Ds=Dl+semi
    xls=Dl/Ds
    RE=np.sqrt(4.0*G*mbh[i]*Msun*Ds*xls*(1.0-xls))/c
    ros=Rstar*Rsun*xls/RE 
    vel= 2.0*np.pi*aa/Period/24.0
    tme=2.0*Rstar*Rsun*(1.0+ecn)/vel
    #if(abs(ros-19.0)<0.5): 
    #    print("Ros, MBH,  Rtar:  ",   ros,   mbh[i],     Rstar,    tme)
    #vel=np.sqrt(G*(mstar+mbh[i])*Msun*(2.0/semi-1.0/aa))


    
input("Enter a num")
    
        
        
'''

b= tan1*semi/Rsun
x= np.sqrt(Rstar**2.0- b**2.0)*2.0*Rsun
tim=x/vel/3600.0
if( abs(ros-ros0)<0.01): 
    print(mbh[i],    ros,    tim,   b)
            
semi= np.power(G*(M1+M2)*T*T/(4.0*np.pi*np.pi),1.0/3.0) ##[m]        
v= 2.0*np.pi*semi/T
Dl=Ds-semi
Dls=Ds-Dl
rhos1= 2.0*Rsun*(Dl/Ds)/np.sqrt(4.0*G*M2*Dls*Dl/Ds/(c*c)) 
A1=1.0+2.0/rhos1/rhos1
Dl=Ds-1.3*semi
Dls=Ds-Dl
rhos2= 2.0*Rsun*(Dl/Ds)/np.sqrt(4.0*G*M2*Dls*Dl/Ds/(c*c)) 
A2=1.0+2.0/rhos2/rhos2


print("Semi,   v:  ",  semi/AU , v, rhos1 , rhos2,  abs(rhos2-rhos1), A1, A2, abs(A1-A2) )
input("Enter a number ")

print (2.0*Rsun/v/(3600.0), v*0.001, semi)       
        

input("Enter a number")


semi= np.power(G*(M1+M2)*T*T/(4.0*np.pi*np.pi),1.0/3.0) ##[m]
P= semi*(1.0-ecn*ecn)#m
Pmax=P/(1.0-ecn)#m
Pmin=P/(1.0+ecn)#m
Rs= 2.0*Rsun*Ds/(Ds+Pmax)
RE= np.sqrt(4.0*G*M2)*np.sqrt(Ds*Pmax/(Ds+Pmax))/c
rho=Rs/RE
Amax=1.0+2.0/rho/rho






print ("Semi(Rs):  ",  2.0*43.77*1000.0*T/Rsun/2.0/(2.0*np.pi*(1.0+ecn)) )
print ("Semi(AU):  ",  2.0*43.77*1000.0*T/AU/(2.0*np.pi*(1.0+ecn))  )
ss=2.0*43.77*1000.0*T/(2.0*np.pi)
print ("Sum of mass:  ",  (ss*ss*ss)*4.0*(np.pi*np.pi)/(T*T * G *Msun) )
print ("**********************************************************************")
print ("Mass of the second lens(Msun):  ",  M2/Msun)
print ("**********************************************************************")
print ("Semi major axis(AU):  ",     np.round(semi/AU, 3) )
print ("**********************************************************************")
print ("Projected displacement:  ",np.round( (Pmin/AU)*1000.0/(Ds/KPC) , 3) )
print ("**********************************************************************")
print ("RE(Rsun):  ",    np.round(RE/Rsun, 3)  )
print ("**********************************************************************")
print ("Rho:  ",   np.round(rho, 2)  )
print ("**********************************************************************")
print ("Amax:  ",  np.round(Amax, 5) )
print ("**********************************************************************")
print ("velocity : ",  2.0*np.pi*semi*0.001/T)

print ("Velocity_Agol: ",  102*0.96*Rsun*0.001/88.0/24.0/3600.0)
input("Enter a number   ")

#xls=Ds/(Ds+Pmax)
#Da= 1.00412043/0.999940072862261 -  0.999940072862261
#rho=np.sqrt(2.0/Da)
#Rl= 2.0*G*M2/c/c ## Shoartzshild radius 
#ecen=0.75
#vs= 2.0*np.pi*semi*np.sqrt(1.0-ecen)/T
#Ts= float(2.0*Rs/vs/3600.0) 
#TE= float(RE/vs/3600.0)


#print ("Time_scale,    TS,  TE (hours), rho:  ",   Times,  Ts,   TE,   Rs/RE, semi/AU )
#print ("Velocity:  ",   vs*0.001)

#Rl= 0.0108*Rsun*np.sqrt(pow(M2/Mch, -2.0/3.0)- pow(M2/Mch, 2./3.)) ## radius white dwarf
#Rl=10.0*1000.0;
#mass=(Da*Rs*Rs*0.8*0.8+Rl*Rl)*0.5/RE/RE
#mass= (rho*RE/Rs)**-2.0
#print ("Ds,   xls:  ", Ds, xls)
#print ("Delta magnification:  ",  Da)
#print ("Rs, Rl, RE:  ", Rs,  Rl,  RE)
#print ("mass(M_sun): ",   mass)
#print ("*******************************************************")
input("Enter a number ")
################################################################################
##TIC=  402786070

M1=2.03*Msun ##[kg]
M2=19.0*Msun### mass of compact object ##[kg]
T=8.4*24.0*3600.0 ## period  [s]
Ds=float(1.0/1.379349035082779)*KPC  ##[m]
semi= pow(G*(M1+M2)*T*T/(4.0*np.pi*np.pi),1.0/3.0) ##[m]
xls=1.0-semi/Ds
Da= 10140.0/10090.0-1.0
rho=np.sqrt(2.0/Da)
Rs= 1.8*Rsun*xls
RE= np.sqrt(4.0*G*Ds*Msun) * np.sqrt(xls*(1.0-xls))/c
mass= (rho*RE/Rs)**-2.0 
print ("Ds,   xls:  ", Ds, xls)
print ("Delta magnification:  ",  Da)
print ("rho:  ", rho)
print ("mass(M_sun): ",   mass)

print ("*******************************************************")

################################################################################
##TIC=  303109427

M1=2.2*Msun  ##[kg]
M2=19.0*Msun### mass of compact object [kg] 
T=11.02*24.0*3600.0 ## period  [s]
Ds=float(1.0/0.7333686124629026)*KPC ## [m]
semi= pow(G*(M1+M2)*T*T/(4.0*np.pi*np.pi),1.0/3.0)
xls=1.0-semi/Ds
Da= 4796.0/4768.0-1.0
rho=np.sqrt(2.0/Da)
Rs= 1.8*Rsun*xls
RE= np.sqrt(4.0*G*Ds*Msun) * np.sqrt(xls*(1.0-xls))/c
mass= (rho*RE/Rs)**-2.0
print ("Ds,   xls:  ", Ds, xls)
print ("Delta magnification:  ",  Da)
print ("rho:  ", rho)
print ("mass(M_sun): ",   mass)
print ("*******************************************************")

'''
