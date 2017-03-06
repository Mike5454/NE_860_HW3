from __future__ import division
import scipy 
import matplotlib.pyplot as plt
import numpy as np    

def compute_pij(x, SigmaT, SigmaS, S, left, right):
    Pij=[]
    DeltI =[]
    for i in range(0,len(x)-1):
        DeltI.append(x[i+1] - x[i])
        Pij.append( 1- 1/(2*SigmaT[i]*(x[i+1] - x[i])) * (1-2*scipy.special.expn(3,SigmaT[i]*(x[i+1] - x[i]))))
#     Calculation of optical length tau
    tau = np.zeros((len(x)-1,len(x)-1))
    for i in range(0,len(x)-1):
        for j in range(0,len(x)-1):
            if i == j: 
                tau[i,j] = 0
            elif i <= j:
                tcalc = 0
                for l in range(i+1,j):
                    tcalc = DeltI[l]*SigmaT[l] + tcalc
                tau[i,j] = tcalc
            elif i >= j:
                tcalc = 0
                for l in range(j,i-1):
                    tcalc = DeltI[l]*SigmaT[l+1] + tcalc
                tau[i,j] = tcalc
#    print " This is the tau Matrix" 
#    print tau
#    print "\n"

    # Obtain Probabilities and Matrix 
    Mat = np.zeros((len(x)-1,len(x)-1))
    for i in range(0,len(x)-1):
        for j in range(0,len(x)-1):
            if i == j:
                Mat[i,j] = 1 - (1/(2*SigmaT[i]*DeltI[i]))*(1-2*scipy.special.expn(3,SigmaT[i]*DeltI[i]))
            else:
                Mat[i,j] = (1/(2*SigmaT[j]*DeltI[j]))*(scipy.special.expn(3,tau[i,j]) - scipy.special.expn(3,tau[i,j] + \
                    DeltI[i]*SigmaT[i]) - scipy.special.expn(3,tau[i,j] + DeltI[j]*SigmaT[j])+scipy.special.expn(3,tau[i,j] + DeltI[i]*SigmaT[i] + DeltI[j]*SigmaT[j]))
#    print " This is the Pij Matrix" 
#    print Mat
#    print "\n"
    return Mat
    
def solve_cpm(x, SigmaT, SigmaS, S):
    # Make the H Matrix to caluclate the values of phi
    print x
    HMat = np.zeros((len(x)-1,len(x)-1))
    for i in range(0,len(x)-1):
        for j in range(0,len(x)-1):
            if i == j: 
                HMat[i,j] = 1 - (SigmaS[i]/SigmaT[i]) * compute_pij(x, SigmaT, SigmaS, S, left, right)[i,j]
            else:
                HMat[i,j] = -(SigmaS[j]/SigmaT[j]) * compute_pij(x, SigmaT, SigmaS, S, left, right)[i,j] 
#    print " This is the H Matrix" 
#    print HMat
#    print "\n"
    
    # Define the Solution Matrix
    Sol = np.zeros((len(x)-1, 1))
    for i in range(0,len(x)-1):
        solve = 0
        for j in range(0, len(x)-1):
            solve = S[j]*compute_pij(x, SigmaT, SigmaS, S, left, right)[i,j]*(x[j+1]-x[j]) + solve
        Sol[i,0] = solve 
#    print " This is the Solution Matrix"
#    print Sol
#    print "\n"
    
    # Solve for the X-Matrix
    XMat = np.linalg.solve(HMat, Sol)
#    print "This is the X-Matrix"
#    print XMat
#    print "\n"
    # Divide components in X-Matrix to get flux solutions
    phi =[]
    for i in range(0,len(x)-1):
        phi.append(XMat[i,0]/((x[i+1]-x[i])*SigmaT[i]))
    return phi
    
    
def plot(xval, yval):
    plt.figure
    params = {'mathtext.default': 'regular' }          
    plt.figure(figsize = (15,10.5))
    plt.rcParams.update(params)
    
    plt.xlabel("Cell Number",  fontname="Arial", fontsize=30)
    plt.tick_params(which='major', length=15, labelsize=25)
    plt.tick_params(which='minor', length=7)
    #grid(b=True, which='major', color='light grey', linestyle='-')
    plt.grid(True, which='minor', color='lightgrey', linestyle='-')
    plt.grid(True, which='major', color='dimgrey', linestyle='-')
    
    plt.ylabel(r'Flux $\phi$', fontname="Arial", fontsize=30)
    plt.title ("Collision Probability Method",fontsize=30)
    plt.rc('font',family='Arial')
    #plt.errorbar(phits_eve,phits_tave,yerr=phits_FError, linestyle="None",capsize=8)
    #plt.errorbar(HZETRN_eave, HZETRN_Results,xerr=HZETRN_err, linestyle="None", ecolor='r', capsize=0)
    plt.plot(xval, yval, 'k-', label = 'Data', linewidth = 4)
    
    plt.legend(loc=1,prop={'size':20})
    plt.show()
    
left = 'vacuum'
right = 'vacuum'

x = map(float, [0.0, 1.0, 2.0, 3.0])
SigmaT = map(float, [1.0, 1.5, 1.0])
SigmaS = map(float, [1.0, 1.5, 1.0])
S = map(float, [0.0, 1.0, 0.0])

pij = compute_pij(x, SigmaT, SigmaS, S, left, right)
phi = solve_cpm(x, SigmaT, SigmaS, S)
x2 = x[1:]
print "These are the flux values"
print phi
print "\n"
s = np.linspace(1, len(x2), 1000)
v = np.polyfit(x2, phi, 2)
t = s**2*v[0] + s*v[1]+v[2]
plot(s,t)