'''
Set of functions used for orbit propagation.
'''



#LIBRARIES---------------------------------------------


from sgp4.earth_gravity import wgs84
from sgp4.earth_gravity import wgs72
from sgp4.earth_gravity import wgs72old
from sgp4.io import twoline2rv
from datetime import datetime
from datetime import timedelta
import numpy as np


#FUNCTIONS---------------------------------------------


def sgp4Prop(line1, line2, t0, tf, dt, model):
    if model == '84':
        model = wgs84
    elif model == '72':
        model = wgs72
    else:
        model = wgs72old
    temp = t0
    time = [t0]
    X = np.array([[0,0,0]])
    sat = twoline2rv(line1, line2, model)
    for d in range(0, (tf - t0).days + 1):
        lim = datetime(temp.year, temp.month, temp.day, tf.hour, tf.minute, 0, 0)
        while temp < lim:
            time.append(temp)
            r, v = sat.propagate(temp.year, temp.month, temp.day, temp.hour, temp.minute, temp.second)
            X = np.vstack((X,r))
            temp += dt
        X = np.vstack((X,(0,0,0)))
        time.append(temp)
        temp = t0 + timedelta(days=1)
    return time, X


def KeplerProp(line1, line2, t0, tf, dt):
    a=0 # COMPLETE