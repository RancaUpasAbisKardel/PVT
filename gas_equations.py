def sg_gas(mwg):
    # Gas Specific Gravity from Molecular Weight
    y = mwg/28.96
    return y


def nat_gas_Ppc(yg):
    # Pseudo-critical Pressure for Common Natural Gas Reservoirs
    Ppc = 677 + 15*yg - 37.5*pow(yg, 2)
    return Ppc


def nat_gas_Tpc(yg):
    # Pseudo-critical Temperature for Common Natural Gas Reservoirs
    Tpc = 168 + 325*yg - 12.5*pow(yg, 2)
    return Tpc


def gas_cond_Ppc(yg):
    # Pseudo-critical Pressure for Common Gas Condensate Reservoirs
    y = 706 - 51.7*yg - 11.1*pow(yg, 2)
    return y


def gas_cond_Tpc(yg):
    # Pseudo-critical Temperature for Common Gas Condensate Reservoirs
    y = 187 + 330*yg - 71.5*pow(yg, 2)
    return y


def Ppr(p, ppc):
    # Pseudo-reduced pressure
    # p in psia
    y = p/ppc
    return y


def Tpr(t, tpc):
    # Pseudo-reduced temperature
    # t in R
    y = t/tpc
    return y


def z_HallYarborough(tpr, ppr):
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