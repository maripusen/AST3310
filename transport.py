import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as intp


mum = 1.66053904 * 10**(-27)
N_A=6.0221e23

"""
This bit is just convertion of the last project to fit the needs of this one better, thus making it easier
to calculate the luminosity as it varies with T and

"""

"Defining all number densities as a function of density"
def npr(rho):
    return rho* 0.7/(mum)

def nhe4(rho):
    return rho* 0.29/(4*mum)

def nhe3(rho):
    return (rho* 10**(-10))/(3*mum)

def nli(rho):
    return (rho * 10 ** (-7)) / (7 * mum)

def nbe(rho):
    return (rho * 10 ** (-7)) / (7 * mum)

def nni(rho):
    return (rho * 10 ** (-11)) / (14 * mum)

def ne(rho):
    return rho/(2*(mum)) * (1+ 0.7)


mev = 1.6*10**(-13)
"Defining Q values for PP chain in joules"
Qpp = 1.177 * mev
Qdh = 5.494 * mev
QHeHe = 12.86 * mev
Qhe3he4 = 1.586 * mev
Qbee = 0.049 * mev
QLih = 17.346 * mev
Qbeh = 0.137 * mev
QBbe = 8.367 * mev
QBeHe = 2.995 * mev


"Defining the CNO chain Q values in joules"
QCH = 1.944 * mev
QNCe = 1.513 * mev
QCHN = 7.551 * mev
QNHO = 7.297 * mev
QON = 1.757 * mev
QNHC = 4.966 * mev


"Implementing a function to calculate the reaction rates based on number densities"
def r_ik(ni,nk, rho, lam, chk):
    if chk == "1":
        rate = ni * nk /(rho * (2)) * lam
    else:
        rate = ni * nk / (rho) * lam
    return rate

"Nicking the main results from project 1 and modifying them to fit the needs of this project. "

def epsilon(T,rho):
    T_9 = float(T)/10**9
    T_9str1=T_9/(1+4.95e-2*T_9)
    T_9str2=T_9/(1+0.759*T_9)
    lambda_pp=(4.01e-15*T_9**(-2/3)*np.exp(-3.380*T_9**(-1/3))*(1+0.123*T_9**(1/3)+1.09*T_9**(2/3)+0.938*T_9))/N_A/1e6
    lambda_33=(6.04e10*T_9**(-2/3)*np.exp(-12.276*T_9**(-1/3))*(1+0.034*T_9**(1/3)-0.522*T_9**(2/3)-0.124*T_9+0.353*T_9**(4/3)+0.213*T_9**(5/3)))/N_A/1e6
    lambda_34=(5.61e6*T_9str1**(5/6)*T_9**(-3/2)*np.exp(-12.826*T_9str1**(-1/3)))/N_A/1e6
    lambda_e7=(1.34e-10*T_9**(-1/2)*(1-0.537*T_9**(1/3)+3.86*T_9**(2/3)+0.0027*T_9**(-1)*np.exp(2.515e-3*T_9**(-1))))/N_A/1e6
    lambda_17mrk=(1.096e9*T_9**(-2/3)*np.exp(-8.472*T_9**(-1/3))-4.83e8*T_9str2**(5/6)*T_9**(-3/2)*np.exp(-8.472*T_9str2**(-1/3))+1.06e10*T_9**(-3/2)*np.exp(-30.442*T_9**(-1)))/N_A/1e6
    lambda_17=(3.11e5*T_9**(-2/3)*np.exp(-10.262*T_9**(-1/3))+2.53e3*T_9**(-3/2)*np.exp(-7.306*T_9**(-1)))/N_A/1e6
    lambda_p14=(4.9e7*T_9**(-2/3)*np.exp(-15.228*T_9**(-1/3)-0.092*T_9**2)*(1+0.027*T_9**(1/3)-0.778*T_9**(2/3)-0.149*T_9+0.261*T_9**(4/3)+0.127*T_9**(5/3))+2.37e3*T_9**(-3/2)*np.exp(-3.011*T_9**(-1))+2.19e4*np.exp(-12.53*T_9**(-1)))/N_A/1e6
    if T_9 < 10**(-3): #and N_A * lambda_e7 >  1.57*10**(-7) /ne(rho):
        lambda_e7 = (1.57*10**(-7))/(ne(rho)*N_A/10**6)

    RPP0 = r_ik(npr(rho), npr(rho), rho, lambda_pp, "1")
    RPP10 = r_ik(nhe3(rho),nhe3(rho), rho, lambda_33, "1")
    RPP20 = r_ik(nhe3(rho),nhe4(rho), rho, lambda_34, "n")
    #print(RPP0, RPP10,RPP20,r_ik(nbe(rho), ne(rho), rho, lambda_e7,"n") ,r_ik(nbe(rho), npr(rho), rho, lambda_17,"n"),r_ik(nli(rho), npr(rho), rho, lambda_17mrk,"n"))
    if RPP0 < 2 * RPP10 + RPP20:
        ratio = RPP0 / (2 * RPP10 + RPP20)
        RPP10 = ratio * RPP10
        RPP20 = ratio * RPP20

    RPP21 = r_ik(nbe(rho), ne(rho), rho, lambda_e7,"n")
    RPP22 = r_ik(nli(rho), npr(rho), rho, lambda_17mrk,"n")

    RPP31 = r_ik(nbe(rho), npr(rho), rho, lambda_17,"n")
    RPP32 = RPP31
    RPP33 = RPP32
    if RPP20 < RPP21 + RPP31:
        ratio = RPP20/(RPP21+RPP31)
        RPP21 = RPP21 * ratio
        RPP31 = RPP31 * ratio
        RPP32 = RPP31
        RPP33 = RPP32
    if RPP21 < RPP22:
        RPP22 = RPP21
    RCNO = r_ik(nni(rho), npr(rho),rho,lambda_p14,"n")
    CNO = RCNO*(QCH + QCHN + QON + QNHC + QNCe + QNHO)
    PP0 = RPP0 * (Qpp + Qdh)
    PP10 = RPP10  * QHeHe

    PP20 = RPP20  * Qhe3he4
    PP21 = RPP21  * Qbee
    PP22 = RPP22  * QLih

    PP31 = RPP31  * Qbeh
    PP32 = RPP32  * QBbe
    PP33 = RPP33  * QBeHe

    rates = [PP0, PP10, PP20, PP21, PP22, PP31 +PP32+ PP33,CNO]

    return sum(rates)



"Defining a few necessary constants"
k = 1.38065*10**(-23)
sigma = 5.670374419 * 10**(-8)
c = 299792458
m_u = 1.66053904*10**(-27)
G = 6.67430*10**(-11)
mu = 0.618237

"Defining initial stellar conditions"
L =  3.828*10**26
R = 696340000
M = 1.9891 * 10**30

L0 = 1 * L
R0 = 1 * R
M0 = 1 * M
rho0 = 1 *  1.42 * 10**(-7) * 1408
T0 =  1 * 5770

Cp = 5/2 * k/(mu*m_u)


"Firstly dealing with everything to do with the reading and handling of the opacity file"

def reading(file):
    "Function to read input file and format it so that it outputs a set of arrays"
    data = open(file ,'r')
    lines = data.readlines()
    n = np.int((np.shape(lines)[0]))
    tlog = np.zeros(n)
    opacity = [i for i in range(0,n)]
    errordata = np.zeros(n)
    "Running through the lines and extracting the values in the data file"
    for i, line in enumerate(lines):
        if i == 0:
            R = np.asarray(line.split()[1:])
        if i == 1:
            None
        if i != 1 and i != 0:
            tlog[i] = line.split()[0]
            opacity[i]  = line.split()[1:]
            None
    "Returning and slicing to appropriate dimensons"
    return R.astype(float), tlog[2:].astype(float),np.asarray(opacity[2:]).astype(float)

LogR, LogT, LogK = reading("data.txt")


def interpolate(R,T,k):
    R, T = np.meshgrid(R,T)
    "Returning a function that interpolates the values and returns a warning if out of bounds"
    return intp.interp2d(R,T,k, kind = "linear")

def kappaSI():
    interpolated = interpolate(LogR,LogT,LogK)
    return interpolated

def pressure(rho,T):
    "A simple function to determine the pressure"
    Pg = rho*k*T/(mu*m_u)
    Prad = 4*sigma/(c*3) * T**4
    return Pg + Prad

def g(M,R):
    'Simple function for calculating and updating g as the simulation evolves'
    return 6.67430*10**(-11) * M /(R**2)

def Hp(P,g,rho):
    return P/(g*rho)

def density(P,T):
    density = ((mu * m_u)/(k * T)) * (P -(4*sigma*T**4)/(3*c))
    return density

def u(T):
    u = 3/2 * 1/(mu*m_u)*K*T
    return u

"Defining the four differential equations "
def drdm(r,rho):
    return (1/(4*np.pi*r**2 *rho))

def dPdm(r,m):
    return -(6.67430*10**(-11) * m)/(4*np.pi*r**4)

def dLdm(epsilon):
    return epsilon

def dTdm(kap,L,r,T):
    return - (3*kap * L)/(256*np.pi**2*sigma*r**4*T**3)

def dynamicdm(step,derive):
    p = 0.01
    dm = step / derive * p
    return dm


"""
This section of the code is dedicated to the gradients and calculating them and the fluxes .
"""

def gradstar(nabs,nabad,rho,T,g,M,R,hp,kappa,alpha):
    #Simple function for calculating the * gradient
    lm = alpha * hp
    U = (64*sigma*T**3)/(3*kappa* rho**2 *Cp)*np.sqrt(hp/(g))
    coeffs = np.asarray((1,U/lm**2 ,4*U**2/lm**4 ,-U*(nabs-nabad)/lm**2))
    xi = np.roots(coeffs)
    xi = xi[np.isreal(xi)]
    xi = np.real(xi)
    #Extracting the relevant root of the polynomial

    return nabs - (lm**2)*(xi**3) /U

def gradstable(L, rho,T, kappa, R,hp):
    #Simple function for the stable gradient
    return 3*L*kappa*rho*hp/(64*np.pi*(R**2) * sigma * T**4)

def gradad():
    #Adiabatic gradient is constant
    return 2/5

def conflux(gradstable,gradstar,T,kappa,rho,hp):
    #Function for calculating the convective flux
    return (gradstable -gradstar)*(16*sigma*T**4)/(3*kappa*rho*hp)

def radflux(T,kappa,rho,hp,nabstar):
    #Function for calculating the radiative flux
    return (16*sigma*T**4)/(3*kappa*rho*hp)*nabstar

"""
This section of the code is dedicated to solving the differential equations using the forward euler approach.
"""

def polate(T, rho, Sanity = False):
    #interpolates or extrapolates in case of values outside boundaries
    if Sanity == False:
        x = np.log10(T)
        y = np.log10(rho / (1000 * (T * 10 ** (-6)) ** 3))
        if x < min(LogT) or x > max(LogT) or y < min(LogR) or y > max(LogR):
            print("Achtung, Warning, avertissement, the opacity value will be extrapolated, this will lead to inaccuracies, you have been warned")
    else:
        x = T
        y = rho

    Kappa = LogK

    b = np.zeros(2)
    c = np.zeros(2)

    if x > LogT[-1]:
        numb = -2
        b[0] = LogT[-2]
        b[1] = LogT[-1]
    elif x < LogT[0]:
        numb = 0
        b[0] = LogT[1]
        b[1] = LogT[0]
    else:
        a = LogT[LogT > x][0]
        numb = np.where(a == LogT)[0][0] - 1
        b[0] = a
        b[1] = LogT[numb]
    if y > LogR[-1]:
        numb2 = -2
        c[0] = LogR[-2]
        c[1] = LogR[-1]
    elif y < LogR[0]:
        numb2 = 0
        c[0] = LogR[1]
        c[1] = LogR[0]
    else:
        a = LogR[LogR > y][0]
        numb2 = np.where(a == LogR)[0][0] - 1
        c[0] = LogR[numb2]
        c[1] = a

        #Using a linear method for interpolation found on the internet(wikipedia)
    f = ( Kappa[numb + 1, numb2] * (c[1] - y) * (b[1] - x) \
        + Kappa[numb + 1, numb2 + 1] * (y - c[0]) * (b[1] - x) \
        + Kappa[numb, numb2] * (c[1] - y) * (x - b[0]) \
        + Kappa[numb, numb2 + 1] * (y - c[0]) * (x - b[0]) ) \
          / ( (c[1] - c[0]) * (b[1] - b[0]) )
    if Sanity == False:
        intp = 10 ** f * 0.1
    else:
        intp = f
    return intp


def solve():
    #Initialising the differential equations
    P0 = pressure(rho0,T0)
    M = [M0]
    rho = [rho0]
    R = [R0]
    P = [P0]
    T = [T0]
    L = [L0]
    kappa = [polate(T0,rho0)]

    gr0 = g(M0,R0)
    hp0 = Hp(P0,gr0,rho0)

    gradstab =[gradstable(L0,rho0,T0,kappa[0],R0,hp0)]
    gradst = [gradstar(gradstab[0],gradad(),rho0,T0,gr0,M0,R0,hp0,kappa[0],1)]
    conf = [conflux(gradstab[0],gradst[0],T0,kappa[0],rho0,hp0)]
    radf =[radflux(T0,kappa[0],rho[0],hp0,gradst[0])]
    i = 0
    while M[i] > 0:
        #Solving the differential equations using the euler method
        #Calculating the values for the tangent for this step
        kappa.append(polate(T[i],rho[i]))
        gr  = g(M[i],R[i])
        hpi = Hp(P[i],gr,rho[i])
        R1 = drdm(R[i],rho[i])
        P1 = dPdm(R[i],M[i])
        L1 = dLdm(epsilon(T[i],rho[i]))
        gradstab.append(gradstable(L[i],rho[i],float(T[i]),kappa[i],R[i],hpi))
        if gradad() < gradstab[-1]:
            gradst.append(gradstar(gradstab[-1],gradad(),rho[i],T[i],gr,M[i],R[i],hpi,kappa[i],1))
        else:
            gradst.append(gradstab[-1])

        coff =-3/16 * kappa[i]*rho[i]/(sigma*T[i]**3)* 1/(4*np.pi * R[i]**2*rho[i])
        T1 = coff* radflux(T[i],kappa[i],rho[i],hpi,gradst[i+1])
        radf.append(radflux(T[i],kappa[i],rho[i],hpi,gradst[i]))
        conf.append(conflux(gradstab[-1],gradst[-1],T0,kappa[i],rho[i],hpi))
        #Calculating dm
        dm1 = abs(dynamicdm(R[i-1],R1))
        dm2 = abs(dynamicdm(P[i-1],P1))
        dm3 = abs(dynamicdm(L[i-1],L1))
        dm4 = abs(dynamicdm(T[i-1],T1))
        dm = min([dm1,dm2,dm3,dm4])
        #This bit was used for troubleshooting and is left in for possible use and analysis
        """
        c = 1
        if i // c > (i - 1) // c:
            index = i
            print("|%5i|%11.4e|%11.4e|%11.4e|%11.4e|%11.4e|%12.4e|%12.4e|%12.8e|%6.5e|"%(index,np.asarray(R[index]), np.asarray(M[index]), np.asarray(P[index]), np.asarray(L[index]), \
                                                                                           np.asarray(T[index]), np.asarray(gradstab[index]), np.asarray(gradst[index]),np.asarray(rho[index]),dm),\
                                                                                           dLdm(epsilon(np.asarray(T[index]),np.asarray(rho[index]))),np.asarray(radflux(T[i],kappa[i],rho[i],hpi,gradst[-1])),hpi,np.asarray(kappa[index]))
        """
        #Adding the neew step to the solution
        R.append(R[i] - R1*dm)
        P.append(P[i] - P1*dm)
        L.append(L[i] - L1*dm)
        T.append(T[i] - T1*dm)

        rho.append(density(P[i+1],T[i+1]))
        M.append(M[i] - dm)

        i = i + 1

        #Setting a few conditions on the end values of the calculations and ending the loop for these conditions
        if M[i]-dm < 0:
            M[i] = 0

        elif len(M) >= 30000:
            M[i] = 0
    #Returning the solved values
    return R,P,L,T,M,rho,kappa, conf, gradst, gradstab,radf


#This bit of the code is dedicated to implementing the tests for
def intertest():
    "Defining a test to verify the interpolated values"
    logts =    np.asarray([3.750,3.755,3.755,3.755,3.755,3.770,3.780,3.795,3.770,3.775,3.780,3.795,3.800])
    logrs =    np.asarray([-6.00 ,-5.95 ,-5.80 ,-5.70 ,-5.55 ,-5.95 ,-5.95 ,-5.95 ,-5.80 ,-5.75 ,-5.70 ,-5.55 ,-5.50])
    compvals = np.asarray([-1.55,-1.51,-1.57,-1.61,-1.67,-1.33,-1.20,-1.02,-1.39,-1.35,-1.31,-1.16,-1.11])
    compvals2 =np.asarray([2.88,3.11,2.68,2.46,2.12,4.70,6.25,9.45,4.05,4.43,4.94,6.89,7.69]) * 10**(-3)
    logR, logT, logK = reading("data.txt")

    k = 0
    tol = 0.05 # It is practical to have a tolerance value, in this case it is set to a relative 5%

    for i in range(0, len(logts)):
        k = polate(logts[i],logrs[i],Sanity = True)
        "Testing for the likeness of the interpolated and comparison values"
        if 1- abs(k/compvals[i]) > tol:
            print(k, compvals[i],i,logts[i],logrs[i],k/compvals[i])
            print("Sanity check has failed, interpolated values are more than 2% different from comparison values")
        else:
            k+=1

    for j in range(0, len(logts)):
        l = polate(10**logts[j], 10**logrs[j]*1000)
        if 1-abs(l/(compvals2[j])) > tol:
            print(l,compvals2[j],j,logts[j],logrs[j],l/compvals2[j])

            print("Sanity check has failed, interpolated values are more than 2% different from comparison values")
        else:
            k+=1
    if k == 2*len(logrs):
        output = "Interpolation test has passed"
    else:
        output = "Interpolation has failed due to "
def gradtest():
    expected = [3.26,32.5*10**6,5.94*10**5,1.175*10**-3,0.4,65.62,0.88,0.12]
    #Initialising relevant initial conditions
    L =  3.828*10**26
    R = 696340000
    M = 1.9891 * 10**30

    T = 0.9e6
    rho = 55.9
    R = 0.84 *R
    M = 0.99 *M
    gr = g(M,R)
    lum = L
    kappa = 3.98
    nabad = 2/5
    alpha = 1
    #Calculating all the necessary bits  for the sanity check
    hp = Hp(pressure(rho,T),gr,rho)

    U = (64*sigma*T**3)/(3*kappa* rho**2 *Cp)*np.sqrt(hp/(gr))

    lm = alpha * hp
    Q = np.pi * lm
    S = 2*np.pi *lm**2
    d = lm
    nabs = gradstable(lum,rho,T,kappa,R, hp)
    coeffs = np.asarray((lm**2/U , 1 ,U*S/(Q*d*lm ),-(3.26-nabad)))
    xi = np.roots(coeffs)
    xi = xi[np.isreal(xi)]
    xi = np.real(xi)
    v = np.sqrt((gr*lm**2)/(4*hp))*xi
    nabstar = gradstar(3.26, gradad(),rho,T,gr,M,R,hp,kappa,1)
    con = conflux(nabs,nabstar,T,kappa,rho,hp)
    rad = radflux(T,kappa,rho,hp,nabstar)
    rcon = con/(con+rad)
    rrad = rad/(con + rad)
    #Leaving it to the user to verify the validity of the calculated values
    print("""All these values are calculated with a mu = 0.618 that implies a minor difference in calculated and expected values
Expected nabla_stable; %6.4e, calculated; %6.4e
Expected Hp; %6.4e, calculated; %6.4e
Expected U; %6.4e, calculated; %6.4e
Expected Xi; %6.4e, calculated; %6.4e
Expected nablastar; %9.8e, calculated; %9.8e
Expected v; %6.4e, calculated; %6.4e
Expected relative convective flux; %6.4e, calculated; %6.4e
Expected relative radiative flux; %6.4e, calculated; %6.4e""" %(expected[0],nabs,expected[1],hp,expected[2],U,expected[3],xi,expected[4],nabstar,expected[5],v,expected[6],rcon,expected[7],rrad))
    test = input("Does this satisfy your sanity?[y/n]")
    if test == "y":
        print("Sanity has received a positive user input, sanity is a free test! ")
    else:
        print("Sanity has received a negative user input and has thus not passed.")

def cross_section(R, L, F_C, show_every=20):
    """
    plot cross section of star
    :param R: radius, array
    :param L: luminosity, array
    :param F_C: convective flux, array
    :param show_every: plot every <show_every> steps
    """

    R_sun = 6.96E8      # [m]
    L_sun = 3.846E26    # [W]

    plt.figure(figsize=(800/100, 800/100))
    fig = plt.gcf()
    ax  = plt.gca()

    r_range = 1.2 * R[0] / R_sun
    rmax    = np.max(R)

    ax.set_xlim(-r_range, r_range)
    ax.set_ylim(-r_range, r_range)
    ax.set_aspect('equal')

    core_limit = 0.995 * L_sun

    j = 0
    for k in range(0,len(R)-1):
        j += 1
        # plot every <show_every> steps
        if j%show_every == 0:
            if L[k] >= core_limit:     # outside core
                if F_C[k] > 0.0:       # plot convection outside core
                    circle_red = plt.Circle((0, 0), R[k]/rmax, color='red', fill=False)
                    ax.add_artist(circle_red)
                else:                  # plot radiation outside core
                    circle_yellow = plt.Circle((0, 0), R[k]/rmax, color='yellow', fill=False)
                    ax.add_artist(circle_yellow)
            else:                      # inside core
                if F_C[k] > 0.0:       # plot convection inside core
                    circle_blue = plt.Circle((0, 0), R[k]/rmax, color='blue', fill=False)
                    ax.add_artist(circle_blue)
                else:                  # plot radiation inside core
                    circle_cyan = plt.Circle((0, 0), R[k]/rmax, color='cyan', fill=False)
                    ax.add_artist(circle_cyan)

    # create legends
    circle_red    = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color='red', fill=True)
    circle_yellow = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color='yellow', fill=True)
    circle_blue   = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color='blue', fill=True)
    circle_cyan   = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color='cyan', fill=True)

    ax.legend([circle_red,circle_yellow,circle_cyan,circle_blue],\
              ['Convection outside core','Radiation outside core','Radiation inside core','Convection inside core']\
              , fontsize=13)
    plt.xlabel(r'$R$', fontsize=13)
    plt.ylabel(r'$R$', fontsize=13)
    plt.title('Cross section of star', fontsize=15)
    plt.show()


if __name__ == "__main__":
    #cross_section(M, R,L,Flux_con)
    print("Do you wish to run the sanity checks?")
    inp1 = input("Enable sanity checks?[y/n]:")
    if inp1 == "y":
        inp2 = input("Do you wish to run the interpolation sanity check?[y/n]:")
        if inp2 == "y":
            print(intertest())
        else:
            print("Interpolation test has not been run")
        inp3 = input("Do you wish to run the gradient sanity check?[y/n]:")
        if inp3 == "y":
            print(gradtest())
        else:
            print("Gradient test has not been run")
            print("Proceeding to the main body of the program")
    else:
        print("Sanity check has not been run.")
        print("Proceeding to the main body of the program")
    print("Producing cross section of the star with initial paramaters %6.4e,[kg/m^3], %6.4e [W/m^2], %6.4e [kg], %6.4e [m]"%(rho0,L0,M0,R0))

    R, P, L, T, M, rho, kappa, conf, grads, gradstab,radf = solve()


    R = np.asarray(R)
    P = np.asarray(P)
    L = np.asarray(L)
    T = np.asarray(T)
    M = np.asarray(M)
    rho = np.asarray(rho)
    kappa = np.asarray(kappa)

    cross_section(R,L,conf,show_every=5)


    print("""Final values:
Final luminosity; %8.6e
Final Mass; %8.6e
Final Radius; %8.6e
Final temperature; %8.6e"""% (100*(1-L[-1]/L0),100*(1-M[-2]/M0),100*(1-R[-1]/R0),T[-1]))

    print(len(M))
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(R/R0,M, label="Mass")
    plt.xlabel("R/R0")
    plt.legend(fontsize=13)

    plt.subplot(1,2,2)
    plt.axhline(0.995*L0, label = r"$0.995 \cdot L_0$",c = "r")
    plt.axhline(0.05*L0, label = r"$0.05 \cdot L_0$",c = "r")
    plt.plot(R/R0,L,label = "luminosity")
    plt.xlabel("R/R0")
    plt.legend(fontsize=13)

    plt.figure()
    plt.subplot(1,2,1)
    plt.yscale("log")
    plt.plot(R/R0,P,label = "Pressure")
    plt.xlabel("R/R0")
    plt.legend(fontsize=13)


    plt.subplot(1,2,2)
    plt.yscale("log")
    plt.plot(R/R0,rho,label = "density")
    plt.xlabel("R/R0")
    plt.legend(fontsize=13)

    plt.figure()
    plt.subplot(1,2,1)
    plt.yscale("log")
    plt.plot(R/R0,np.asarray(conf)/(np.asarray(conf)+np.asarray(radf)),label = "Relative convective flux")
    plt.xlabel("R/R0")

    plt.legend(fontsize=20)
    plt.subplot(1,2,2)
    plt.yscale("log")
    plt.xlabel("R/R0")
    plt.axhline(gradad(),label = r"$\nabla_{ad}$",c ="r")
    plt.plot(R/R0,gradstab,label = r"$\nabla_{stable}$",c="g")
    plt.plot(R/R0,grads,label = r"$\nabla^*$",c="b")
    plt.legend(fontsize=20)
    plt.show()

    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(R/R0,kappa,label = "kappa")
    plt.xlabel("R/R0")
    plt.legend(fontsize=13)
    plt.show()

    plt.subplot(1,2,2)
    plt.plot(R/R0,T,label = "Temperature")
    plt.legend(fontsize=13)
    plt.xlabel("R/R0")

    plt.show()
