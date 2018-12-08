'''
Set of functions used for coordinate transformations.
'''




#LIBRARIES---------------------------------------------



import numpy as np




#CONSTANTS---------------------------------------------



a_Ea = 6378.137
f_Ea = 1/298.257223563
e_Ea = np.sqrt(1-(1-f_Ea)**2)




#FUNCTIONS---------------------------------------------



def rad(deg):                       # degrees to radians
    return deg*np.pi/180



def deg(rad):                       # radians to degrees
    return rad*180/np.pi


# Gets the Julian time
def t2JD (t):
    return 367*t.year - int(7*(t.year+int((t.month+9)/12))/4) + int(275*t.month/9) + t.day + 1721013.5 + t.hour/24 + t.minute/24/60 + t.second/24/60/60 + t.microsecond/24/60/60/1000000


# Gets the Julian 2000 time
def TJ2000 (jd):
    return (jd - 2451545)/36525


# Gets the GMST angle from the Julian 2000 time (the function takes UT1 in julian format but the difference with UTC is always less than a second which is less than the precision I need)
def thGMST(Tut1, t):
    th = 100.4606184 + 36000.77005361*Tut1 + 0.00038793*Tut1**2 - 2.6e-8*Tut1**3
    w = 1.002737909350795 * 2*np.pi / (24*3600)
    th = rad(th)
    gmst = th + w*(t.hour*3600 + t.minute*60 + t.second)
    while gmst < 0:
        gmst += 2*np.pi
    while gmst >= 2*np.pi:
        gmst -= 2*np.pi
    return gmst


# ECI -> ECEF
def ECI2ECEF(x, y, z, gmst):
    rot3 = np.array([[np.cos(gmst),np.sin(gmst),0],[-np.sin(gmst),np.cos(gmst),0],[0,0,1]])
    return np.dot(rot3, (x,y,z))


# LLH -> ECEF
def LLH2ECEF (gdlat, lon, hellp):
    N = a_Ea/np.sqrt(1-e_Ea**2*np.sin(gdlat)**2)
    x = (N+hellp)*np.cos(gdlat)*np.cos(lon)
    y = (N+hellp)*np.cos(gdlat)*np.sin(lon)
    z = ((1-f_Ea)**2*N+hellp)*np.sin(gdlat)
    return np.array([x,y,z])


# ECEF -> SEZ
def ECEF2SEZ(x, y, z, gdlat, lon, hellp): #vx,vy,vz, ):
    x,y,z = (x,y,z) - LLH2ECEF(gdlat,lon,hellp)
    rot3long = np.array([[np.cos(lon),np.sin(lon),0],[-np.sin(lon),np.cos(lon),0],[0,0,1]])
    aux = np.pi/2 - gdlat
    rot2lat = np.array([[np.cos(aux),0,-np.sin(aux)],[0,1,0],[np.sin(aux),0,np.cos(aux)]])
    S,E,Z = np.dot(np.dot(rot2lat,rot3long),np.array([x,y,z]))
    #vS,vE,vZ = np.dot(np.dot(rot2lat,rot3long),np.array([vx,vy,vz]))
    return S,E,Z #,vS,vE,vZ


# SEZ -> obs
def SEZ2obs (S, E, Z): #, velcalc, vS, vE, vZ):
    ran = np.sqrt(S**2 + E**2 + Z**2)
    el = np.arcsin(Z/ran)
    #if not(velcalc):
        #vS = -1
        #vE = 0
        #vZ = 0
    if abs(el - np.pi) > 0.00001:
        azi = np.arctan2(E/np.sqrt(S**2+E**2), -S/np.sqrt(S**2+E**2))
    else:
        azi = np.arctan2(vE/np.sqrt(vS**2+vE**2), -vS/np.sqrt(vS**2+vE**2))
    #vran = np.dot([S,E,Z],[vS,vE,vZ])/ran
    #vel = (vZ - vran*np.sin(el))/np.sqrt(S**2+E**2)
    #vazi = (vS*E-vE*S)/(S**2+E**2)
    return azi,el,ran #,vazi,vel,vran


# ECI -> obs
def ECI2obs(x, y, z, gdlat, lon, hellp, t):
    gdlat = rad(gdlat)
    lon = rad(lon)
    T = TJ2000(t2JD(t))
    gmst = thGMST(T, t)
    ecef = ECI2ECEF(x, y, z, gmst)
    s, e, z = ECEF2SEZ(ecef[0], ecef[1], ecef[2], gdlat, lon, hellp)
    azi, el, ran = SEZ2obs(s, e, z)
    azi = deg(azi)
    el = deg(el)
    return np.array([azi, el, ran])



#maplon = np.arctan2(ecef[1]/(ecef[0]**2+ecef[1]**2), ecef[0]/(ecef[0]**2+ecef[1]**2))
#maplat = np.arcsin(ecef[2]/ran)
