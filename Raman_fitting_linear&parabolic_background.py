import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize


array = np.genfromtxt("S1_Raman_center_2700.txt", dtype='float')
x,y = array.T.copy()

plt.plot(x, y,'r') 
plt.title("spectra fitting", fontsize=25) 
plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)
plt.xlabel('Raman shift(1/cm)', loc ="center", fontsize=18)
plt.ylabel('Counts', loc ="center", fontsize=18)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()


def linear(x, m, b):
    return m*x + b

def parabolic(x, a, b, c):
    return a*x**2 + b*x + c

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

def _1Lorentzian(x, amp, cen, wid):
    return amp*wid**2/((x-cen)**2+wid**2)

def _slope_Lorentzian(x, m, b, amp, cen, wid):
    return m*x + b + amp*wid**2/((x-cen)**2+wid**2)

def _parabolic_Lorentzian(x, a, b, c, amp, cen, wid):
    return a*x**2 + b*x + c + amp*wid**2/((x-cen)**2+wid**2)


### Lorentzian fitting with parabolic background
p0=[0.001, -0.098, 1825, 90, 2712, 25]
popt, pcov = optimize.curve_fit(_parabolic_Lorentzian, x, y, p0)
perr = np.sqrt(np.diag(pcov))

print ("amplitude = %0.2f (+/-) %0.2f" % (popt[3], perr[3]))
print ("center = %0.2f (+/-) %0.2f" % (popt[4], perr[4]))
print ("sigma = %0.2f (+/-) %0.2f" % (popt[5], perr[5]))


plt.plot(x, y,'r')
plt.plot(x, _parabolic_Lorentzian(x, *popt),'b')
plt.title("spectra fitting", fontsize=25) 
plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)
plt.xlabel('Raman shift(1/cm)', loc ="center", fontsize=18)
plt.ylabel('Counts', loc ="center", fontsize=18)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()


def FWHM(X,Y):
    half_max = _1Lorentzian(x, *popt[3:]).max()*0.5
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(_1Lorentzian(x, *popt[3:])[0:-1])) - np.sign(half_max - np.array(_1Lorentzian(x, *popt[3:])[1:]))        #plot(X[0:len(d)],d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    return X[right_idx] - X[left_idx] #return the difference (full width)


print('FWHM', FWHM(x,_1Lorentzian(x, *popt[3:])))


plt.plot(x, _1Lorentzian(x, *popt[3:]),'b')
plt.plot(x, np.ones(y.shape)*_1Lorentzian(x, *popt[3:]).max()*0.5,'g')
plt.title("spectra fitting", fontsize=25) 
plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)
plt.xlabel('Raman shift(1/cm)', loc ="center", fontsize=18)
plt.ylabel('Counts', loc ="center", fontsize=18)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()




##############################################################

### Lorentzian fitting with linear background

p0=[-0.098, 1825, 90, 2712, 25]
popt, pcov = optimize.curve_fit(_slope_Lorentzian, x, y, p0)
perr = np.sqrt(np.diag(pcov))

print ("amplitude = %0.2f (+/-) %0.2f" % (popt[2], perr[2]))
print ("center = %0.2f (+/-) %0.2f" % (popt[3], perr[3]))
print ("sigma = %0.2f (+/-) %0.2f" % (popt[4], perr[4]))


plt.plot(x, y,'r')
plt.plot(x, _slope_Lorentzian(x, *popt),'b')
plt.title("spectra fitting", fontsize=25) 
plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)
plt.xlabel('Raman shift(1/cm)', loc ="center", fontsize=18)
plt.ylabel('Counts', loc ="center", fontsize=18)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()


def FWHM(X,Y):
    half_max = _1Lorentzian(x, *popt[2:]).max()*0.5
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(_1Lorentzian(x, *popt[2:])[0:-1])) - np.sign(half_max - np.array(_1Lorentzian(x, *popt[2:])[1:]))        #plot(X[0:len(d)],d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    return X[right_idx] - X[left_idx] #return the difference (full width)
    

print('FWHM', FWHM(x,_1Lorentzian(x, *popt[2:])))


plt.plot(x, _1Lorentzian(x, *popt[2:]),'b')
plt.plot(x, np.ones(y.shape)*_1Lorentzian(x, *popt[2:]).max()*0.5,'g')
plt.title("spectra fitting", fontsize=25) 
plt.xticks(fontsize= 16)
plt.yticks(fontsize= 16)
plt.xlabel('Raman shift(1/cm)', loc ="center", fontsize=18)
plt.ylabel('Counts', loc ="center", fontsize=18)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()

