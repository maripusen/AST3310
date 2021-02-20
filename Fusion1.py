import numpy as np
import matplotlib.pyplot as plt

print("To run the sanity check the normal plotting is disabled, thus to achieve a check and plot, the script has to be run twice. \n Due to the sanity check only being run for one temperature the plotting and check is not being run at the same time.")
inp = input("This is to disable or enable the sanity check, do you wish to enable the check? input is [y/n]:")


Tvar = np.logspace(3.76,7.13,num = 10000, base = 10)
rho = 1.62 * 10**5
mu = 1.66053904 * 10**(-27)
Tsol = 1.57*10**10

if inp in ["y", "y "]:
    T_9 = 1.57*10**7 /10**9
else:
    T_9=Tvar / 10**9


N_A=6.0221e23

"Defining all number densities as a function of density"
def npr(rho):
    return rho* 0.7/(mu)

def nhe4(rho):
    return rho* 0.29/(4*mu)

def nhe3(rho):
    return (rho* 10**(-10))/(3*mu)

def nli(rho):
    return (rho * 10 ** (-7)) / (7 * mu)

def nbe(rho):
    return (rho * 10 ** (-7)) / (7 * mu)

def nni(rho):
    return (rho * 10 ** (-11)) / (14 * mu)

def ne(rho):
    return rho/(2*(mu)) * (1+ 0.7)


T_9str1=T_9/(1+4.95e-2*T_9)
T_9str2=T_9/(1+0.759*T_9)
lambda_pp=(4.01e-15*T_9**(-2/3)*np.exp(-3.380*T_9**(-1/3))*(1+0.123*T_9**(1/3)+1.09*T_9**(2/3)+0.938*T_9))/N_A/1e6
lambda_33=(6.04e10*T_9**(-2/3)*np.exp(-12.276*T_9**(-1/3))*(1+0.034*T_9**(1/3)-0.522*T_9**(2/3)-0.124*T_9+0.353*T_9**(4/3)+0.213*T_9**(5/3)))/N_A/1e6
lambda_34=(5.61e6*T_9str1**(5/6)*T_9**(-3/2)*np.exp(-12.826*T_9str1**(-1/3)))/N_A/1e6
lambda_e7=(1.34e-10*T_9**(-1/2)*(1-0.537*T_9**(1/3)+3.86*T_9**(2/3)+0.0027*T_9**(-1)*np.exp(2.515e-3*T_9**(-1))))/N_A/1e6
lambda_17mrk=(1.096e9*T_9**(-2/3)*np.exp(-8.472*T_9**(-1/3))-4.83e8*T_9str2**(5/6)*T_9**(-3/2)*np.exp(-8.472*T_9str2**(-1/3))+1.06e10*T_9**(-3/2)*np.exp(-30.442*T_9**(-1)))/N_A/1e6
lambda_17=(3.11e5*T_9**(-2/3)*np.exp(-10.262*T_9**(-1/3))+2.53e3*T_9**(-3/2)*np.exp(-7.306*T_9**(-1)))/N_A/1e6
lambda_p14=(4.9e7*T_9**(-2/3)*np.exp(-15.228*T_9**(-1/3)-0.092*T_9**2)*(1+0.027*T_9**(1/3)-0.778*T_9**(2/3)-0.149*T_9+0.261*T_9**(4/3)+0.127*T_9**(5/3))+2.37e3*T_9**(-3/2)*np.exp(-3.011*T_9**(-1))+2.19e4*np.exp(-12.53*T_9**(-1)))/N_A/1e6
if isinstance(lambda_e7,np.ndarray):
    for i in range(len((Tvar))):
        if Tvar[i] < 10**6: #and N_A * lambda_e7 >  1.57*10**(-7) /ne(rho):
            lambda_e7[i] = (1.57*10**(-7))/(ne(rho)*N_A/10**6)


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


def RPP(T):
    if isinstance(lambda_p14,float) == True:
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
    else:
        RPP21 = r_ik(nbe(rho), ne(rho), rho, lambda_e7,"n")
        RPP22 = r_ik(nli(rho), npr(rho), rho, lambda_17mrk,"n")
        RPP0 = r_ik(npr(rho), npr(rho), rho, lambda_pp, "1")
        RPP10 = r_ik(nhe3(rho),nhe3(rho), rho, lambda_33, "1")
        RPP20 = r_ik(nhe3(rho),nhe4(rho), rho, lambda_34, "n")
        RPP31 = r_ik(nbe(rho), npr(rho), rho, lambda_17,"n")
        RPP32 = RPP31
        RPP33 = RPP32
        for i in range(len(Tvar)):
            if RPP0[i] < 2 * RPP10[i] + RPP20[i]:
                ratio = RPP0[i] / (2 * RPP10[i] + RPP20[i])
                RPP10[i] = ratio * RPP10[i]
                RPP20[i] = ratio * RPP20[i]

            #print(RPP21 + RPP31 ,RPP20)

            if RPP20[i] < RPP21[i] + RPP31[i]:
                ratio2 = RPP20[i]/(RPP21[i]+RPP31[i])
                RPP21[i] = RPP21[i] * ratio2
                RPP31[i] = RPP31[i] * ratio2
                RPP32[i] = RPP31[i]
                RPP33[i] = RPP32[i]

                if RPP21[i] < RPP22[i]:
                    RPP22[i] = RPP21[i]

            if RPP21[i] < RPP22[i]:
                RPP22[i] = RPP21[i]
    rates = [RPP0, RPP10, RPP20, RPP21, RPP22, RPP31, RPP32, RPP33]

    return np.asarray(rates)

def PP(T):
    R = RPP(Tvar)
    RPP0 = R[0]
    RPP10 = R[1]
    RPP20 = R[2]
    RPP21 = R[3]
    RPP22 = R[4]
    RPP31 = R[5]
    RPP32 = R[6]
    RPP33 = R[7]
    #print(RPP20,RPP21,RPP22)
    #print(RPP10)
    #print("---------------------------------------")
    #print(RPP20)
    #print("---------------------------------------")
    #print(RPP21)
    #print("---------------------------------------")
    #print(RPP31)

    PP0 = RPP0 * rho * (Qpp + Qdh)
    PP10 = RPP10 * rho * QHeHe

    PP20 = RPP20 * rho * Qhe3he4
    PP21 = RPP21 * rho * Qbee
    PP22 = RPP22 * rho * QLih

    PP31 = RPP31 * rho * Qbeh
    PP32 = RPP32 * rho * QBbe
    PP33 = RPP33 * rho * QBeHe
    branch = [PP0, PP10, PP20, PP21, PP22, PP31 + PP32 + PP33]
    #print(RPP10*QHeHe+(Qpp + Qdh)* RPP10)
    #print((RPP31 * Qbeh + RPP32 * QBbe + RPP33 * QBeHe) + ((Qpp + Qdh) *  RPP(Tvar)[2]))
    if inp not in ["y", "y "]:
        PP0 = RPP0 * (Qpp + Qdh)
        PP1 = (RPP10 * QHeHe) + (Qpp + Qdh) * RPP10

        PP2 = (RPP21 * Qbee + RPP22 * QLih) + ((Qpp + Qdh) * RPP20) + Qhe3he4 * RPP20

        PP3 = (RPP31 * Qbeh + RPP32 * QBbe + RPP33 * QBeHe) + ((Qpp + Qdh) *  RPP20) + Qhe3he4 * RPP20
        neutrino = np.asarray([0.815 * mev * RPP21 + 0.265 * mev * RPP20, 6.711 * RPP32 * mev + 0.265 * mev * RPP20])
        branch = np.asarray([PP1, PP2, PP3])

    #print(RPP20,RPP21)
        return np.asarray((branch,neutrino))
    else:
        return np.asarray(branch)
def CNO(T):
    RCNO = r_ik(nni(rho), npr(rho),rho,lambda_p14,"n")
    #print(RCNO)
    CNO = RCNO*rho*(QCH + QCHN + QON + QNHC + QNCe + QNHO)
    #print(CNO/rho)

    if inp not in ["y", "y "]:
        CNO = RCNO * (QCH + QCHN + QON + QNHC + QNCe + QNHO)
        neutrino = RCNO*(0.707*mev + 0.997*mev)
        return np.asarray((CNO,neutrino))
    else:
        return np.asarray(CNO)



def epsilon(T):
    sumPP = np.zeros(len(CNO(Tvar)))
    sumPP2 = sumPP
    red = PP(Tvar)[:-1]
    red2 = PP(Tvar)[:-1]
    for i in range(len(sumPP)):
        sumPP[i] = sum(red[:,i])
    eps = sumPP + CNO(Tvar)
    neutrino = sumPP2 + CNO[1]
    return eps

def sanity(T,inp):
    tol = 0.01
    expCNO = 9.18*10**(-8)
    expPP = [404, 8.68*10**(-9), 4.86*10**(-5), 1.49*10**(-6), 5.29*10**(-4),1.63*10**(-6)]
    PPc = PP(Tsol)
    if inp in ["y", "y "]:
        for i in range(len(expPP)):
            if abs(expPP[i]/PPc[i] -1) >= tol:
                print(abs(expPP[i]/PPc[i] -1))
                print("Sanity failed")
                break

            assert abs(CNO(Tsol)/expCNO -1) < tol

        return("Sanity check has passed, all is well")
    else:
        return "Sanity check has not been run"


if __name__ == "__main__":
    if inp in ["y", "y "]:
        print(sanity(Tsol, inp))
    else:
        print("Sanity check has not been run, thus producing plot:")
        #print(np.shape(RPP(Tvar)), RPP(Tvar)[2])
        #print((PP(Tvar)[0] + PP(Tvar)[1] + PP(Tvar)[2] + CNO(Tvar))/epsilon(Tvar))
        #norm = epsilon(Tvar)
        norm = (PP(Tvar)[0][0] + PP(Tvar)[0][1] + PP(Tvar)[0][2] + CNO(Tvar)[0])
        print("Calculated relative neutrino energy loss is calculated to be:\n",(CNO(Tvar)[1][120]/(CNO(Tvar)[0][120] + CNO(Tvar)[1][120]), PP(Tvar)[1][1]/( PP(Tvar)[0][2] +  PP(Tvar)[1][1])))

        plt.title("Relative energy production ", fontsize = 20)
        plt.plot(Tvar, (PP(Tvar)[0][2])/norm, Label = "PP3 Branch")
        plt.plot(Tvar, (PP(Tvar)[0][1])/norm, Label = "PP2 Branch")

        plt.plot(Tvar, (PP(Tvar)[0][0])/norm, Label = "PP1 Branch")
        plt.plot(Tvar, CNO(Tvar)[0]/norm, Label = "CNO Branch")
        plt.xscale("log")
        plt.legend(fontsize = 20)
        plt.show()
