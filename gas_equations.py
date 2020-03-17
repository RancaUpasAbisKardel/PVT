def sg_gas(mwg):
    # Gas Specific Gravity from Molecular Weight
    y = mwg/28.96
    return y


def nat_gas_Ppc_Tpc(yg):
    # Pseudo-critical Pressure for Common Natural Gas Reservoirs
    Ppc = 677 + 15*yg - 37.5*pow(yg, 2)
    Tpc = 168 + 325*yg - 12.5*pow(yg, 2)
    return Ppc, Tpc


def gas_cond_Ppc_Tpc(yg):
    # Pseudo-critical Pressure for Common Gas Condensate Reservoirs
    ppc = 706 - 51.7*yg - 11.1*pow(yg, 2)
    tpc = 187 + 330*yg - 71.5*pow(yg, 2)
    return ppc, tpc


def Ppr(p, t, ppc, tpc):
    # Pseudo-reduced pressure
    # p in psia
    # t in F
    ppr = p/ppc
    tpr = (t+460)/tpc
    return ppr, tpr


def Ppc_Tpc_WichertAziz_correction(Ppc, Tpc, yh2s, yco2):
    # yh2s : mole fraction of H2S in gas mixture
    # yco2 : mole fraction of CO2 in gas mixture
    # applicable for natural gas and gas condensate
    a = yh2s + yco2
    b = yco2
    e = 120*(pow(a, 0.9)-pow(a, 1.6))+15*(pow(b, 0.5)-pow(b, 4))
    tpc2 = Tpc - e
    ppc2 = (Ppc*tpc2)/(Tpc+b*(1-b)*e)
    return ppc2, tpc2


def Ppc_Tpc_CarrKobayashiBurrows_Correction(Ppc, Tpc, yh2s, yco2, yn2):
    # only applicable to natural gas
    # yh2s, yco2, yn2 : mole fraction of H2S, CO2, N2
    tpc2 = Tpc - 80*yco2 + 130*yh2s - 250*yn2
    ppc2 = Ppc + 440*yco2 + 600*yh2s - 170*yn2
    return ppc2, tpc2


# Gas Compressibility Factor (z)
def z_HallYarborough(ppr, tpr):
    t = 1/tpr
    # initial guess
    y = 0.0125*ppr*t*pow(2.718281828, -1.2*pow(1-t, 2))
    err = 100
    # constants
    x1 = -0.06125*ppr*t*pow(2.718281828, -1.2*pow(1-t, 2))
    x2 = 14.76*t - 9.76*pow(t, 2) + 4.58*pow(t, 3)
    x3 = 90.7*t - 242.2*pow(t, 2) + 42.4*pow(t, 3)
    x4 = 2.18 + 2.82*t
    while err > pow(10, -12):
        fdy = ((1+4*y+4*y*y-4*pow(y, 3)+pow(y, 4))/pow(1-y, 4)) - 2*x2*y + x3*x4*pow(y, x4-1)
        fy = x1 - x2*pow(y, 2) + x3*pow(y, x4) + (y+y*y+pow(y, 3)+pow(y, 4))/pow(1-y, 3)
        yy = y - fy/fdy
        err = abs(y-yy)
        y = yy
    z = (0.06125*ppr*t/y)*pow(2.718281828, -1.2*pow(1-t, 2))
    return z


def z_DranchukAbouKassem(ppr, tpr):
    ec = 2.718281828    # constant e value
    a = [0, 0.3265, -1.07, -0.5339, 0.01569, -0.05165,
          0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.721]
    r1 = a[1] + (a[2]/pow(tpr, 1)) + (a[3]/pow(tpr, 3)) + (a[4]/pow(tpr, 4)) + (a[5]/pow(tpr, 5))
    r2 = 0.27*(ppr/tpr)
    r3 = a[6] + (a[7]/pow(tpr, 1)) + (a[8]/pow(tpr, 2))
    r4 = a[9]*((a[7]/pow(tpr, 1))+(a[8]/pow(tpr, 2)))
    r5 = (a[10]/pow(tpr, 3))
    # Initial Guess
    err = 100
    y = 0.27*(ppr/tpr)
    fy = r1*y - (r2/y) + r3*(pow(y, 2)) - r4*pow(y, 5) + \
         r5*(1+a[11]*pow(y, 2))*pow(ec, -a[11]*pow(y, 2)) + 1
    while err > pow(10, -12):
        fdy = r1 + r2*pow(y, -2) + 2*r3*y - 5*r4*pow(y, 4) + \
              2*r5*y*pow(ec, -a[11]*pow(y, 2))*(1+2*a[11]*pow(y, 3)) - \
              a[11]*pow(y, 2)*(1+a[11]*pow(y, 2))
        yy = y - fy/fdy
        err = abs(y - yy)
        y = yy
    z = 0.27*ppr/(tpr*y)
    return z