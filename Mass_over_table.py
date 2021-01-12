import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

show=4.4        #Variable for choosing which graph to show (apart from the trajectory one)

# show=1.1 -> R_m vs time graph 
# show=1.2 -> Radial position (|Rm|) vs time graph 
# show=1.3 -> Phi_m vs time graph
# show=1.4 -> Angular position vs time graph 
# show=1.5 -> Vertical position vs time graph 
# show=2.1 -> dR_m/dt vs time graph 
# show=2.2 -> Radial velocity vs time graph
# show=2.3 -> Angular velocity vs time graph 
# show=3.1 -> d2R_m/dt2 vs time graph 
# show=3.2 -> Radial acceleration vs time graph 
# show=3.3 -> Angular acceleration vs time graph
# show=4.1 -> Kinetic energy of mass "m" vs time graph 
# show=4.2 -> Kinetic energy of mass "M" vs time graph
# show=4.3 -> Potential energy of mass "M" vs time graph
# show=4.4 -> Total energy of the system vs time graph 
# show=5.1 -> Problematic Term 1 vs time graph 
# show=5.2 -> Problematic Term 2 vs time graph 
# show=5.3 -> Problematic Term 3 vs time graph 


h=0.01                             #Variable for changing the xlim and ylim   
rr=0.0001                          #Threshold radius [m]

simtime=8                          #Total simulation time [s]
dt=0.00001                         #Time delta for the finite difference scheme [s] 
M=1                                #Mass of the hanging body [Kg]
m=0.5                              #Mass of the surface body [Kg]
miu=0.3                            #Friction coefficient of the surface [-]
g=9.8                              #Acceleration due to gravity [m/s^2]
L=6                                #Cord's total length                               

Phi=np.zeros(int(simtime/dt))       #Inicialization of the Angular Position array 
Rad=np.zeros(int(simtime/dt))       #Inicialization of the Radial Position array
VelPhi=np.zeros(int(simtime/dt))    #Inicialization of the Angular Velocity array
VelRad=np.zeros(int(simtime/dt))    #Inicialization of the Radial Velocity array
AccelRad=np.zeros(int(simtime/dt))  #Inicialization of the Radial Acceleration array
AccelPhi=np.zeros(int(simtime/dt))  #Inicialization of the Angular Acceleration array
KE_m=np.zeros(int(simtime/dt))      #Inicialization of the Kinetic Energy array for mass "m"
KE_M=np.zeros(int(simtime/dt))      #Inicialization of the Kinetic Energy array for mass "M"
U_M=np.zeros(int(simtime/dt))       #Inicialization of the Potential Energy array for mass "M"


Phi0=0                             #Initial angular position [rad]
Rad0=2                             #Initial radial position [m]
VelPhi0=3                          #Initial angular velocity [rad/s]
VelRad0=2                          #Initial radial velocity [m/s] 

#When miu=0, VelRad0=0 and VelPhi0=np.sqrt((M*g/(m*Rad0))) we obtain uniform circular motion.

Phi[1]=Phi0       
Rad[1]=Rad0
Phi[0]=Phi0-dt*VelPhi0
Rad[0]=Rad0-dt*VelRad0
VelRad[0]=VelRad0
VelPhi[0]=VelPhi0
 

q=0 #Origin crossings counter

#Finite difference scheme--------------------------------------------------------------------------
for i in range(1,int(simtime/dt)-1):
    
    VelRad[i]=(Rad[i]-Rad[i-1])/dt
    VelPhi[i]=(Phi[i]-Phi[i-1])/dt
              
    a=VelRad[i]
    b=VelPhi[i]
    
    if a==0 and b==0 and abs(Rad[i])>rr:
        Rad[i+1]= 2*Rad[i] - Rad[i-1] + (dt**2)*(abs(Rad[i])*b**2 - np.sign(Rad[i])*M*g/m)*(1 + M/m)**(-1)
        Phi[i+1]= 2*Phi[i] - Phi[i-1] - (dt**2)*(2*a*b*np.sign(Rad[i])/abs(Rad[i]))               
    elif (a!=0 or b!=0) and abs(Rad[i])>rr:  
        Rad[i+1]= 2*Rad[i] - Rad[i-1] + (dt**2)*(abs(Rad[i])*b**2 - np.sign(Rad[i])*M*g/m - (miu*g*a)/np.sqrt(a**2 + (Rad[i]*b)**2))*(1 + M/m)**(-1)
        Phi[i+1]= 2*Phi[i] - Phi[i-1] - (dt**2)*(2*a*b*np.sign(Rad[i])/abs(Rad[i]) + (miu*g*b)/np.sqrt(a**2 + (Rad[i]*b)**2)) 
    elif (a!=0 or b!=0) and abs(Rad[i])<=rr:  
        Rad[i+1]= 2*Rad[i] - Rad[i-1] + (dt**2)*(abs(Rad[i])*b**2 - np.sign(Rad[i])*M*g/m - (miu*g*a)/np.sqrt(a**2 + (Rad[i]*b)**2))*(1 + M/m)**(-1)
        Phi[i+1]= 2*Phi[i] - Phi[i-1] - (dt**2)*(miu*g*b)/np.sqrt(a**2 + (Rad[i]*b)**2) 
    elif a==0 and b==0 and abs(Rad[i])<=rr:  
        Rad[i+1]= 2*Rad[i] - Rad[i-1] + (dt**2)*(abs(Rad[i])*b**2 - np.sign(Rad[i])*M*g/m )*(1 + M/m)**(-1)
        Phi[i+1]= 2*Phi[i] - Phi[i-1] 
     
    AccelRad[i-1]=(VelRad[i]-VelRad[i-1])/dt   
    AccelPhi[i-1]=(VelPhi[i]-VelPhi[i-1])/dt   
        
    if np.sign(Rad[i])!=np.sign(Rad[i+1]):
        q=q+1
#------------------------------------------------------------------------------------------


#Angular acceleration correction
AccelRad[len(AccelRad)-1]=AccelRad[len(AccelRad)-3]
AccelRad[len(AccelRad)-2]=AccelRad[len(AccelRad)-3]
AccelPhi[len(AccelRad)-1]=AccelPhi[len(AccelPhi)-3]
AccelPhi[len(AccelRad)-2]=AccelPhi[len(AccelPhi)-3]


#Calculation of the real physical variables
flips=np.zeros(q+1) 
flips[0]=0   
qq=0    
for i in range(1,len(Rad)):
    if np.sign(Rad[i-1])!=np.sign(Rad[i]):
        flips[qq]=i
        qq=qq+1
flips[q]=len(Rad)     

VelRad2=np.sign(Rad)*VelRad.copy()  #Real Radial velocity Vector
AccelRad2=np.sign(Rad)*AccelRad.copy() #Real acceleration Vector
Phi2=Phi.copy()  #Real angular Vector
Rad2=abs(Rad.copy())  
qqq=1    
for j in range(0,len(flips)-1):
    Phi2[int(flips[j]):int(flips[j+1])]=Phi2[int(flips[j]):int(flips[j+1])]+3.14159265*qqq
    qqq=qqq+1
    
        
#WARNINGS
if L<max(Rad):
    print("WARNING: The maximum radial distance calculated in this simulation exceeds the value of L. The simulation is only accurate in so far the radial distance doesn't exceed this value. ")    
if (M/m)<=miu:
    print("WARNING: The coefficient of friction (miu) is greater than the cocient M/m. The simulation will be accurate only as long as (M/m)>miu ")    
        
#Calculation of the position of the surface body in cartesian coordinates
X=Rad*np.cos(Phi) 
Y=Rad*np.sin(Phi)
Z=abs(Rad)-L 

#Calculation of the energy vectors
KE_m=0.5*m*(VelRad**2 + (Rad*VelPhi)**2)
KE_M=0.5*M*VelRad**2 
U_M=(L+Z)*M*g
KT= KE_m + KE_M + U_M

#Calculation of the problematic terms
sqtermR=(miu*g*VelRad)/np.sqrt(VelRad**2 + (Rad[i]*VelPhi)**2)
sqtermPhi=(miu*g*VelPhi)/np.sqrt(VelRad**2 + (Rad[i]*VelPhi)**2)
cocitermPhi=2*VelRad*VelPhi/abs(Rad[i])   


#Grahphics
if show==1.1:

    #R_m vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), Rad)
    plt.xlim((0, simtime))
    plt.ylim((min(Rad)*(1+h), max(Rad)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('R_m [m]')  
    plt.title( 'R_m(t) vs. time')  
    plt.show() 
    
elif show==1.2:

    #Radial position (|Rm|) vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), Rad2)
    plt.xlim((0, simtime))
    plt.ylim((min(Rad2)*(1-h), max(Rad2)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('|R_m| [m]')  
    plt.title( 'Radial position of the body of mass "m" (|R_m(t)|) vs. time')  
    plt.show()
    
elif show==1.3:

    #Phi_m vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), Phi)
    plt.xlim((0, simtime))
    plt.ylim((min(Phi)*(1-h), max(Phi)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('Phi_m (without pi shift) [rad]')  
    plt.title('Phi_m (without pi shift) vs. time')  
    plt.show() 

elif show==1.4:

    #Angular position vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), Phi2)
    plt.xlim((0, simtime))
    plt.ylim((min(Phi2)*(1-h), max(Phi2)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('Phi_m (with pi shift) [rad]')  
    plt.title('Angular position of the body of mass "m" (Phi_m with pi shift) vs. time')  
    plt.show()     
    
elif show==1.5:

    #Vertical position vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), Z)
    plt.xlim((0, simtime))
    plt.ylim(( min(Z)*(1+h), max(Z)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('Z [m]') 
    plt.title('Vertical position of the body of mass "M" vs. time')   
    plt.show() 
    
elif show==2.1:

    #dR_m/dt vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), VelRad)
    plt.xlim((0, simtime))
    plt.ylim(( min(VelRad)*(1+h), max(VelRad)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('dR_m/dt [m/s]') 
    plt.title('dR_m/dt vs. time')   
    plt.show() 
    
elif show==2.2:

    #Radial velocity vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), VelRad2)
    plt.xlim((0, simtime))
    plt.ylim(( min(VelRad2)*(1+h), max(VelRad2)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('VelRad [sign(R_m)*dR_m/dt] [m/s]') 
    plt.title('Radial velocity of the body of mass "m" [sign(R_m)*dR_m/dt] vs. time')   
    plt.show() 
    
elif show==2.3:

    #Angular velocity vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), VelPhi)
    plt.xlim((0, simtime))
    plt.ylim(( min(VelPhi)*(1+h), max(VelPhi)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('VelPhi [rad/s]') 
    plt.title('Angular velocity of the body of mass "m" vs. time')   
    plt.show()  
    
elif show==3.1:

    #d2R_m/dt2 vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), AccelRad)
    plt.xlim((0, simtime))
    plt.ylim(( min(AccelRad)*(1+h), max(AccelRad)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('d2R_m/dt2 [m/s^2]') 
    plt.title('d2R_m/dt2 vs. time')   
    plt.show()  
    
elif show==3.2:

    #Radial acceleration vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), AccelRad2)
    plt.xlim((0, simtime))
    plt.ylim(( min(AccelRad2)*(1+h), max(AccelRad2)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('AccelRad [sign(R_m)*d2R_m/dt2] [m/s^2]') 
    plt.title('Radial acceleration of the body of mass "m" [sign(R_m)*d2R_m/dt2] vs. time')   
    plt.show() 
    
elif show==3.3:

    #Angular acceleration vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), AccelPhi)
    plt.xlim((0, simtime))
    plt.ylim(( min(AccelPhi)*(1+h), max(AccelPhi)*(1-h)))
    plt.xlabel('time [s]')
    plt.ylabel('AccelPhi [d2Phi_m/dt2] [m/s^2]') 
    plt.title('Angular acceleration of the body of mass "m" [d2Phi_m/dt2] vs. time')   
    plt.show()     
    
elif show==4.1:

    #Kinetic energy of mass "m" vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), KE_m)
    plt.xlim((0, simtime))
    plt.ylim(( min(KE_m)*(1-h), max(KE_m)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('KE_m [J]') 
    plt.title('Kinetic energy of the body of mass "m" vs. time')   
    plt.show() 
    
elif show==4.2:

    #Kinetic energy of mass "M" vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), KE_M)
    plt.xlim((0, simtime))
    plt.ylim(( min(KE_M)*(1-h), max(KE_M)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('KE_M [J]') 
    plt.title('Kinetic energy of the body of mass "M" vs. time')   
    plt.show() 

elif show==4.3:

    #Potential energy of mass "M" vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), U_M)
    plt.xlim((0, simtime))
    plt.ylim(( min(U_M)*(1-h), max(U_M)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('U_M [J]') 
    plt.title('Potential energy of the body of mass "M" vs. time')   
    plt.show()

elif show==4.4:

    #Total energy of the system vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), KT)
    plt.xlim((0, simtime))
    plt.ylim(( min(KT)*(1-h), max(KT)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('KT [J]') 
    plt.title('Total energy of the system vs. time')   
    plt.show()
    
elif show==5.1:

    #Problematic Term 1 vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), sqtermR)
    plt.xlim((0, simtime))
    plt.ylim(( min(sqtermR)*(1-h), max(sqtermR)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('sqtermR') 
    plt.title('sqtermR vs. time')   
    plt.show()    
    
elif show==5.2:

    #Problematic Term 2 vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), sqtermPhi)
    plt.xlim((0, simtime))
    plt.ylim(( min(sqtermPhi)*(1-h), max(sqtermPhi)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('sqtermPhi') 
    plt.title('sqtermPhi vs. time')   
    plt.show()     
    
elif show==5.3:

    #Problematic Term 3 vs time graph 
    plt.plot(np.linspace(0, simtime, int(simtime/dt)), cocitermPhi)
    plt.xlim((0, simtime))
    plt.ylim(( min(cocitermPhi)*(1-h), max(cocitermPhi)*(1+h)))
    plt.xlabel('time [s]')
    plt.ylabel('cocitermPhi') 
    plt.title('cocitermPhi vs. time')   
    plt.show() 
    

    

#Animation of the body moving
fig = plt.figure()
ax = plt.subplot(1, 1, 1) 


data_skip = int(0.004375*len(Phi)) #Data sampling rate for the animation

def init_func():    #Initialization function. It runs just once.
    ax.clear()
    plt.xlabel('x [m]')  #Axes labeling 
    plt.ylabel('y [m]')
    plt.xlim((-max(Rad), max(Rad))) #Setting of axes limits 
    plt.ylim((-max(Rad), max(Rad)))
    plt.title('Motion of the body of mass "m" in time')  
    ax.scatter(0, 0, marker='o', color='k') #Point at the origin is drawn

def update_plot(i):    #Graphing Function
    if (i%(2*data_skip)==0):
        ax.clear()
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.xlim((-max(Rad), max(Rad)))
        plt.ylim((-max(Rad), max(Rad)))
        plt.title('Motion of the body of mass "m" over time')  
        ax.scatter(0, 0, marker='o', color='k')
        ax.scatter(X[i], Y[i], marker='o', color='r')  
        ax.plot([0,X[i]],[0, Y[i]], color='#808080')            #Gray rope, liking the moving body with the center hole, is drawn
        ax.plot(X[0:i], Y[0:i], color='r', linestyle='dashed')  #The body's trajectory is drawn
        
anim = FuncAnimation(fig,                                          #Animating function
                     update_plot,
                     frames=np.arange(0, len(X), data_skip),
                     init_func=init_func,
                     interval=15)





