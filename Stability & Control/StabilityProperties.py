import math

#

a = float(input("Mass Value:")) #Mass
b = float(input("Damping Value:")) #Damping
c = float(input("Spring Coefficient Value:")) #Spring Coefficient

#Find Angular Frequency
natAngFreq = math.sqrt(c/a)

#Find Damping Ratio
dampingRatio = b/(2*natAngFreq)

#Find Damped Angular Frequency
dampAngFreq = natAngFreq * math.sqrt(1-dampingRatio**2)

#Find Natural Period
naturalPeriod = (2*math.pi)/natAngFreq

#Find Damped Period
dampedPeriod = (2*math.pi)/ dampAnqFreq

#Print Results
print("Mass: " + str(a))
print ("Damping Ratio: "+ str(dampingRatio))
print ("Natural Angular Frequency: " + str(natAngFreq))
print ("Damped Angular Frequency: "+ str(dampAngFreq))
print ("Damped Period: " + str(dampedPeriod))
print ("Natural Period: " + str(naturalPeriod))