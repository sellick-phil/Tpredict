'''
Final project for Astrodynamics
Astronautical Engineering, National University of San Martín
Professor: DÍAZ RAMOS, Manuel Francisco
Student: BURRONI, Tomás Ignacio
2019/12/13

This code was fully developed by Tom Burroni, using algorithms from 'Fundamentals of Astrodynamics and Applications' by David A. Vallado and the sgp4 package developed by Brandon Rhodes. It was built for use in the Miguelete Ground Station at UNSAM but it is free for anyone else to use.

To get the manuals, report issues, or anything else do not hesitate to contact me.

contact: tburroni@unsam.edu.ar
'''


#LIBRARIES---------------------------------------------


import numpy as np
import urllib.request as url
import os
import Propagator as pr
import Transformations as tr
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
from tabulate import tabulate




#FILES---------------------------------------------


folder = 'Defaults/'
fgenDFTs = 'Defaults/General.txt'
fsatDFTs = 'Defaults/Satellites.txt'
fGSDFTs = 'Defaults/Ground_Stations.txt'
celestrak = 'https://www.celestrak.com/NORAD/elements/'




#CONSTANTS---------------------------------------------


satStandard = np.array([' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\r','\n'])
coarse = 60
fine = 1




#FUNCTIONS---------------------------------------------


# Loads the General Defaults to a list
def loadDFTs():
    opt = []
    genDFTs = open(fgenDFTs,'r')
    read = genDFTs.readline().rstrip()
    # If the file is not empty reads the rest of the file and loads it into opt (name of the default, used only to show the user and access the file)
    while read != '':
        opt.append(read)
        read = genDFTs.readline().rstrip()      # rstrip removes the '\n' at the end, it has to be removed because it uses the name to search for the file
    genDFTs.close()
    return opt


# Creates a general default (includes N satellites and a ground station)
def create_genDFT():
    print ('\nInsert a name for the default (one line without spaces)')
    name = input()      # Name of the default, also used for the file
    
    temp = 'a'
    sats = []           # The sats are accumulated here
    locs = []           # The last bit of the celestrak site where the sat and TLE are found
    sat = 'a'
    while temp != 'stop' and temp != 'Stop':    # Gets sats until the user inputs 'stop'
        print ('\nInsert the satellite\'s name exactly as it appears on celestrak or enter \'stop\' to finish adding satellites')
        temp = input()          # Temporary variable with the sat's name for formatting purposes
        i = 0
        sat = np.copy(satStandard)  # It must have this specific format in amount of spaces, etc.
        if temp != 'stop' and temp != 'Stop':
            for let in temp:        # Replaces the first n places of sat while leaving the rest untouched to respect the format
                sat[i] = temp[i]
                i += 1
            loc = input('\nComplete the path:\n' + celestrak)  # Gets the specific document where the sat is found
            web = url.urlopen(celestrak + loc)
            i2 = 0
            verif = 0
            for line in web:        # Checks if the sat is really there (incorrect document or typing mistake), by comparing one line at a time the whole document until it's found
                a = 1
                i2 = 0
                for let in sat:
                    if sat[i2] != chr(line[i2]):    # Checks full coincidence char by char in each line
                        a = 0
                        break                       # If at least one char mismatches it jumps to the next line
                    i2 += 1
                if a:                               # If it's found it stops the search
                    verif = 1
                    break
            if verif:
                print ('Satellite found!')
                sats.append(sat)        # If the sat was found it's added to the list, if not it's ignored
                locs.append(loc)
            else:
                print ('Satellite not found! Check everything is written correctly')
                
    print ('\nInsert the location of the Ground Station:')  # Once no more sats are added it gets the info on the ground station
    gs_llh = np.array([float(input('Latitude [deg]: ')), float(input('Longitude [deg]: ')), float(input('Height above the ellipsoid [m] (height over sea level + height over geoid, <http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/intptW.html>): ')), int(input('Time zone: '))])
    
    genDFTs = open(fgenDFTs, 'a')
    genDFTs.write(name + '\n')      # Adds the default to the general list, this one is then used to show the available defaults
    genDFTs.close()
    
    new = open(folder + name + '.txt', 'w')     # Creates a new file with the satelites, and ground station info
    for coord in gs_llh:
        new.write(str(coord) + '\n')        # It starts with the ground station so that when the file is read it's easier to know how many sats there are (just reads until the end of the file)
    i = 0
    for a in sats:
        for b in a:
            new.write(b)
        new.write(locs[i] + '\n')
        i += 1
    new.close()
    print ('New general default successfully created\n')


# Loads a general default on the pack
def getDFT(pack, name):
    fdft = open(folder + name + '.txt', 'r')
    for i in range(0,4):            # Loads the coordinates for the ground station
        pack.gs_llh[i] = float(fdft.readline())
    temp = fdft.readline()
    while temp != '':               # Loads all the sats while checking 'end of file'
        pack.sats.append(temp)
        pack.locs.append(fdft.readline())
        temp = fdft.readline()
    fdft.close()
    for i in range(0,len(pack.sats)):
        pack.sats[i] = pack.sats[i].rstrip('\n')    # Removes the line breaks for printing
    return pack


# Creates a custom pack that is not saved, functionality not yet developed
def custom():
    a = 0




#CLASSES---------------------------------------------


class Package():     # class that contains all the satellites and the ground station for the analysis
    def __init__(self):
        self.sats = []                          # Satellites' names
        self.locs = []                          # Satellites's locations in celestrak
        self.gs_llh = np.array([0.0, 0.0, 0.0, 0])  # Ground Station's location in llh and time zone
        self.line1 = []                         # TLE's
        self.line2 = []
        self.X = 0                              # X, Y, Z coordinates in ECI
        self.time = []                          # Time vector (list of datetime values)
        self.obs = 0                            # Azimuth, elevation and range
        self.passes = 0
        
    def getTLEs(self):          # Gets the TLEs for all the sats
        i = 0
        for sat in self.sats:       # For each sat in the pack
            web = url.urlopen(celestrak + self.locs[i])    # Searches online for the latest TLEs
            i2 = 0
            l1 = l2 = 0
            tle1p = []
            tle2p = []
            for line in web:        # Searches for the TLE's that belong to each specific sat by comparing one line at a time the whole document until it's found
                if not l1:
                    a = 1
                    i2 = 0
                    for let in sat:
                        if sat[i2] != chr(line[i2]):   # Checks full coincidence char by char
                            a = 0
                            break
                        i2 += 1
                    if a:                   # Once the sat is found it waits for the next line to save the TLE
                        l1 = 1
                elif not l2:
                    tle1 = line             # Saves the first line and waits for the next line
                    l2 = 1
                else:
                    tle2 = line             # Saves the second line and breaks the loop
                    break
            for i3 in range(0,69):
                tle1p.append(chr(tle1[i3]))     # Turns the binary string to a char array
                tle2p.append(chr(tle2[i3]))
            self.line1.append(''.join(tle1p))   # Turns the array to a string
            self.line2.append(''.join(tle2p))
            del tle1p                       # Deletes the array so as not to store trash (if the next name is shorter the last letters would not be overwritten)
            del tle2p
            i += 1
    
    def propagate(self, grav, t0, tf, dt):      # Propagates using the sgp4 package
        i = 0
        for sat in self.sats:
            self.time, x = pr.sgp4Prop(self.line1[i], self.line2[i], t0, tf, dt, grav)      # This function is in Propagator.py
            x = np.expand_dims(x, axis=0)       # Adjusts the format of the array for concatenation
            if not i:               # If it's the first array it can't concatenate because of the dimension error
                self.X = x
            else:
                self.X = np.concatenate((self.X,x))
            i += 1
    
    def Kepler(self):           # Propagates a Keplerian orbit, functionality not yet developed
        a=0
    
    def transform(self):        # Tranforms ECI to obs
        self.obs = np.copy(self.X)      # Copies the values of X just to have the same dimensions
        for i in range(0,len(self.sats)):   # Goes through each sat
            for i2 in range(0,np.ma.size(self.X,axis=1)):   # For all time
                self.obs[i][i2] = tr.ECI2obs(self.X[i][i2][0], self.X[i][i2][1], self.X[i][i2][2], self.gs_llh[0], self.gs_llh[1], self.gs_llh[2], self.time[i2])   # This function is in Transformations.py
    
    def find_passes(self, grav, minel):
        i = 0
        first = 1
        counter = 0
        for i in range(0,len(self.sats)):
            tcomp = datetime(2000,1,1,0,0,0,0)
            for i2 in range(0,np.ma.size(self.X,axis=1)):
                if self.time[i2] - tcomp > timedelta(0,20*60):
                    if self.obs[i][i2][1] > 0:
                        counter += 1
                        if i2:
                            i2 -= 1
                        temp = self.time[i2]
                        start = 1
                        passx = np.array([i,datetime(2000,1,1,0,0,0,0),0,0,0])
                        el = 1
                        while el > 0 or start:
                            r = pr.sgp4Prop_fine(self.line1[i], self.line2[i], temp, grav)
                            azi, el, ran = tr.ECI2obs(r[0],r[1],r[2],self.gs_llh[0], self.gs_llh[1], self.gs_llh[2],temp)
                            passx = np.vstack((passx,(i,self.time[i2],azi,el,ran)))
                            if start:
                                if el > 0:
                                    start = 0
                            temp += timedelta(0,fine)
                        maxel = np.amax(passx,axis=0)[3]
                        if maxel > minel:
                            if first:
                                self.passes = np.copy(passx)
                                first = 0
                            else:
                                self.passes = np.vstack((self.passes,passx))
                        del passx
                        tcomp = temp
        return counter




#MAIN PROGRAM---------------------------------------------


print ('\n\n-------------------------------Tpredict--------------------------------')
print ('Developed by: Tomas Ignacio Burroni\n\n')

# Making sure the files exist, if an error appears then a folder must have been deleted. By opening in append mode nothing gets deleted if it exists and if it doesn't the file gets created
genDFTs = open(fgenDFTs,'a')
satDFTs = open(fsatDFTs,'a')
GSDFTs = open(fGSDFTs,'a')
genDFTs.close()
satDFTs.close()
GSDFTs.close()

pack = Package()

# Reads the first line to see if something exists
genDFTs = open(fgenDFTs,'r')
read = genDFTs.readline()
DFT = []
genDFTs.close()

# If the file is empty
c = 0
if read == '':
    sel = 0
    while sel < 1 or sel > 2:
        print ('There were no general defaults found, select an option:')
        print ('\t1) Create a default')
        print ('\t2) Make a custom one-time prediction')
        sel = int(input())
        if sel == 1:
            create_genDFT()                     # Creates a new general default
        elif sel == 2:
            c = 1
            custom()                            # Functionality not yet developed
        else:
            print ('Invalid selection, insert again')

opt = loadDFTs()    # Loads the defaults to a list

# If the file is not empty or the user created a default (ergo, didn't choose custom) it goes to the normal menu until the user selects a default or a custom prediction
while c == 0:
    sel = 0
    while sel < 1 or sel > 4:
        print ('Select an option:')
        print ('\t1) Load one of the following defaults:')
        i = 0
        for gen in opt:         # Prints out all of the defaults
            i += 1
            print ('\t\t', i, '-', gen)
        print ('\t2) Create a default')
        print ('\t3) Make a custom one-time prediction')
        print ('\t4) Delete a default')
        sel = int(input())
        if sel == 1:
            c = 1
            print ('Select the desired default')
            sel2 = int(input())
            while sel2 < 1 or sel2 > len(opt):
                print ('Invalid option, insert again')
                sel2 = int(input())
            pack = getDFT(pack, opt[sel2-1])        # Loads default pack
        elif sel == 2:
            create_genDFT()                         # Creates a new general default
            del opt
            opt = loadDFTs()                        # Reloads the defaults to the list
        elif sel == 3:
            c = 1
            custom()                                # Functionality not yet developed
        elif sel == 4:
            print ('Select the desired default to delete')
            sel2 = int(input())
            while sel2 < 1 or sel2 > len(opt):
                print ('Invalid option, insert again')
                sel2 = int(input())
            os.remove(folder + opt[sel2-1] + '.txt')    # Deletes the specific file
            del opt[sel2-1]                             # Deletes the name from the list
            genDFTs = open(fgenDFTs,'w')
            for name in opt:
                genDFTs.write(name + '\n')              # Rewrites the defaults file with the updated list
            genDFTs.close()
            print ('Default deleted\n')
        else:
            print ('Invalid selection, insert again\n')

# Gets the TLEs
pack.getTLEs()
print ('TLEs downloaded\n')

# Gets the time window
print ('Adjust your prediction window\nInsert the starting time')
hs1 = int(input('Hours: '))
mn1 = int(input('Minutes: '))
sec = 0
usec = 0
print ('\nInsert the stop time')
hs2 = int(input('Hours: '))
mn2 = int(input('Minutes: '))
while hs2 < hs1 or (hs2 == hs1 and mn2 <= mn1):
    print ('\nStop time cannot be earlier than starting time, insert again')
    hs2 = int(input('Hours: '))
    mn2 = int(input('Minutes: '))

utcloc = int(input('\nInsert the computer\'s time zone: '))

print ('\nInsert starting day for prediction ("day") or start now ("today")')
sel = input()
if sel == 'today' or sel == 'Today':
    t0 = datetime.now()         # Gets current local time
    t0 = datetime(t0.year, t0.month, t0.day, hs1, mn1, sec, usec)
else:
    yr = int(input('Year: '))
    mth = int(input('Month: '))
    day = int(input('Day: '))
    t0 = datetime(yr, mth, day, hs1, mn1, sec, usec)

print ('\nInsert last day for prediction')
yr = int(input('Year: '))
mth = int(input('Month: '))
day = int(input('Day: '))
tf = datetime(yr, mth, day, hs2, mn2, sec, usec)
while tf < t0:
    print ('\nLast day cannot be earlier than starting day, insert again')
    yr = int(input('Year: '))
    mth = int(input('Month: '))
    day = int(input('Day: '))
    tf = datetime(yr, mth, day, hs2, mn2, sec, usec)

t0 -= timedelta(hours = utcloc)     # Calculations use UTC time
tf -= timedelta(hours = utcloc)

# Gets the minimum accepted elevation
minel = float(input('\nInsert minimum elevation [deg]: '))

# Gets the gravity model to be used
print ('\nSelect the Earth gravity model:')
print ('\t1) WGS84')
print ('\t2) WGS72')
print ('\t3) WGS72old')
print ('\t4) Ellipsoidal Earth, Keplerian orbit')
sel = int(input())
while sel < 1 or sel > 4:
    print ('Invalid option, insert again')
    sel = int(input())
if sel == 4:
    #pack.Kepler(t0, tf)        # Functionality not yet developed
    a=0
else:
    if sel == 1:
        grav = '84'
    elif sel == 2:
        grav = '72'
    else:
        grav = '72old'
    pack.propagate(grav, t0, tf, timedelta(0,coarse))       # Propagates the orbit for the defined time intervals. It doesn't go from t0 to t1 but rather from time in t0 to time in tf for each day between t0 and tf

# Transforms to obs
pack.transform()

# Searches for passes
counter = pack.find_passes(grav, minel)

def mapr(r):            # Remap the radial axis, the default goes from 0 inside to x outside
    return 90 - r


header = ['Date & time', 'Satellite']
table = np.array([datetime(2000,1,1,0,0,0,0), 'NULL'])

# Polar plot
p = 0
for i in range(0,len(pack.passes)):
    if pack.passes[i][1] == datetime(2000,1,1,0,0,0,0):
        for i2 in range(i+1,len(pack.passes)):
            if pack.passes[i2][1] == datetime(2000,1,1,0,0,0,0):
                break
        fig = plt.figure(p)
        ax = fig.add_subplot(1, 1, 1, polar=True)
        invazi = 2*np.pi - pack.passes[i+1:i2,2]*np.pi/180        # The angle also goes the other way and it has to be in radians. I still have to turn the labels around
        ax.set_theta_zero_location('N')
        ax.plot(invazi, mapr(pack.passes[i+1:i2,3]),'bo',markersize=5)
        ax.set_yticks(range(0, 90, 10))                     # Define the yticks
        ax.set_yticklabels(map(str, range(90, 0, -10)))     # Change the labels
        ax.set_ylim(0,90)
        xlabels = list(map(str, range(360, 0, -45)))
        xlabels[0] = 'N'
        xlabels[2] = 'W'
        xlabels[4] = 'S'
        xlabels[6] = 'E'
        ax.set_xticklabels(xlabels)                         # Change the labels
        fig.suptitle(pack.sats[pack.passes[i][0]],x=0.02,y=0.98,ha='left',va='top',size='xx-large')
        ax.set_title(pack.passes[i+1][1] + timedelta(hours = utcloc))
        table = np.vstack((table, np.array([pack.passes[i+1][1],pack.sats[pack.passes[i+1][0]].rstrip()])))
        p += 1

plt.show()

for t in range(1,len(table)):
    table[t][0] += timedelta(hours = utcloc)

table = table[table[:,0].argsort()]
print (tabulate(table[1:,:], header, tablefmt='fancy_grid'))






#DEBUGGING---------------------------------------------

'''
# Triple plot, azimuth, elevation and range vs time
plt.figure(1)

tmin = np.arange(np.ma.size(pack.X,axis=1))/2

plt.subplot(311)
plt.plot(tmin,pack.obs[0][:,0],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('azi [deg]')

plt.subplot(312)
plt.plot(tmin,pack.obs[0][:,1],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('el [deg]')

plt.subplot(313)
plt.plot(tmin,pack.obs[0][:,2],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('ran [km]')

plt.show()

# x, y and z vs time
plt.figure(1)

tmin = np.arange(np.ma.size(pack.X,axis=1))

plt.subplot(311)
plt.plot(tmin,pack.X[0][:,0],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('X [km]')

plt.subplot(312)
plt.plot(tmin,pack.X[0][:,1],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('Y [km]')

plt.subplot(313)
plt.plot(tmin,pack.X[0][:,2],'bo',markersize=1)
plt.xlabel('t [mn]')
plt.ylabel('Z [km]')

plt.show()

# World map, sat's path
plt.figure()

plt.plot(pack.obs[0][:,0],pack.obs[0][:,1],'bo',markersize=1)
plt.grid(which='both',axis='both')
plt.ylabel('Latitude [deg]')
plt.xlabel('Longitude [deg]')

plt.show()
'''
