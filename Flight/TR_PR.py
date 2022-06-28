import math

#import numpy as np
import matplotlib.pyplot as plt

#Single Velocity Value or Multiple

singleValue = False

#Initial Dimensions and Constants

WingArea = 47 # Area of the wings [m^2]
Density = 1.225 #Air density [kg/m^3]
#VelocityRange = [0, 100] #[m/s] this is an example. Put this below. And switch single value to False.
Velocity = [50, 100] #Use a single int here for testing single velocity configurations [m/s]
Mass = 10000 # of the aircraft [kg]
AR = 6.5 #Aspect ratio
CD0 = 0.032 #Zero lift drag coefficient
efficiency = 0.87 

def LiftCoefficient(m, q, s):
    
    L = m * 9.81
    
    if isinstance(q, list):
        
        cl = []
        for pd in q:
            clValue = L/(pd*s)
            cl.append(clValue)
    
    else:
        
        cl = L/(q*s)
    
    return(cl)


def DynamicPressure(density, v):
    
    if isinstance(v, list):
        
        pd = []
        for velocity in range(v[0],v[1]+1):
            pdValue = 0.5 * density * velocity**2
            pd.append(pdValue)
    else:
        
        pd = 0.5 * density * v**2
    
    return(pd)

def InducedDragCoefficient(cl, efficiency, AR):
    
    if isinstance(cl, list):
        
        idc = []
        for col in cl:
            idcValue = col**2/(math.pi * efficiency * AR)
            idc.append(idcValue)
    
    else:
        
        idc = cl**2/(math.pi * efficiency * AR)
        
    return(idc)

def ThrustRequired(idc, cd0, q, s):
    
    if isinstance(idc, list):
        
        tr = []
        for i in range(0, len(idc)):
            trValue = (idc[i]+cd0) * s *  q[i]
            
            tr.append(trValue)
    else:
        
        tr = (idc + cd0) * s * q
    
    return(tr)

def PowerRequired(TR, V):
    
    if isinstance(TR, list):
        
        pr = []
        for i in range(0, len(TR)):
            prValue = TR[i] * V[i]
            
            pr.append((prValue))
    else:
        pr = TR * Velocity
    
    return(pr)

def CoefficientLiftMinPower(cd0, ar, e):
    
    CLMinP =  math.sqrt(3* cd0 * math.pi * ar * e)
    
    return(CLMinP)

def CoefficientLiftMinDrag(cd0, ar, e):
    
    CLMinD = math.sqrt(cd0 * math.pi * ar * e)
    
    return(CLMinD)

def CoefficientLiftToVelocity(cl, m, density, wingarea):
    
    Lift = m * 9.81
    
    velocity = math.sqrt(Lift/(0.5 * density * wingarea * cl))
    
    return(velocity)

def MinThrustRequired(q, wingarea, cd0):
    
    tr = 2 * q * wingarea * cd0
    
    return(tr)

def MinPowerRequired(cd0, mass, cl, v):
    
    L = mass * 9.81
    
    pr = (L/(cl/(4*cd0))) * v
    
    return(pr)

            

if singleValue == False:
    V = []
    
    for v in range(Velocity[0],Velocity[1] + 1):
        V.append(v)
else:
    V = Velocity;
    


DP = DynamicPressure(Density, Velocity)
CL = LiftCoefficient(Mass, DP, WingArea)
ID = InducedDragCoefficient(CL, efficiency, AR)
TR = ThrustRequired(ID, CD0, DP, WingArea)
PR = PowerRequired(TR, V)


#Min Thrust Required

CLminD = CoefficientLiftMinDrag(CD0, AR, efficiency)

VminD = CoefficientLiftToVelocity(CLminD, Mass, Density, WingArea)

DPminD = DynamicPressure(Density, VminD)

TRmin = MinThrustRequired(DPminD, WingArea, CD0)

#Min Power Required

ClminPR = CoefficientLiftMinPower(CD0, AR, efficiency)

VminPR = CoefficientLiftToVelocity(ClminPR, Mass, Density, WingArea)

PRmin = MinPowerRequired(CD0, Mass, ClminPR, VminPR)



if singleValue:
    print('INPUTS')
    
    print('Velocity: ' + str(Velocity) + ' m/s')
    print('Density: ' + str(Density) + ' kg/m^3')
    
    print('Mass: ' + str(Mass) + ' kg')
    
    print('WingArea: ' + str(WingArea) + ' m^2')
    print('AR: ' + str(AR))
    
    
    print('CD0: ' + str(CD0))
    print('Efficiency ' + str(efficiency * 100) + '%')
    
    print('')
    
    print('OUTPUTS')
    
    print('Dynamic Pressure: ' + str(DP) + ' Pa')
    print('Lift Coefficient: ' + str(CL))
    print('Drag Coefficient: ' + str(ID + CD0) )
    print('Induced Drag Coefficient: ' + str(ID))
    print('Thrust Required: ' + str(TR) + ' N')
    print('Thrust Required Minimum: ' + str(TRmin) + ' N')
    print('Velocity for Min Drag: ' + str(VminD) + ' m/s')
    print('Power Required: ' + str(PR) + ' Watts')
    print('Power Required Minimum: ' + str(PRmin) + ' Watts')
    print('Velocity for Min Power Required: ' + str(VminPR) + ' m/s')
    
    

else:
    print('Graphs')
    
    plt.plot(V,PR)
    plt.ylabel('Power Required')
    plt.xlabel('Velocity')
    plt.show()
    
    plt.plot(V,TR)
    plt.ylabel('Thrust Required')
    plt.xlabel('Velocity')
    plt.show()
    
    print('Thrust Required Minimum: ' + str(TRmin) + ' N')
    print('Velocity for Min Drag: ' + str(VminD) + ' m/s')
    
    print('Power Required Minimum: ' + str(PRmin) + ' Watts')
    print('Velocity for Min Power Required: ' + str(VminPR) + ' m/s')
