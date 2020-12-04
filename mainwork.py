#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 19:49:09 2020

@author: user
"""

#%%
import matplotlib as mlp
import numpy as np
import matplotlib.pyplot as plt
import seis_utils
import scipy.stats # for linregress of PL
import datetime
import calendar
import math
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


file_in = "solomon.csv"

'''
                       LOAD DATA
'''
#load dates
                                            #(4,1,2,1,2) to get only dates 
m_dates = np.genfromtxt(file_in, delimiter = (4,1,2,1,2,1,2,1,2,1,2), 
                       usecols=(0,2,4,6,8,10), dtype = int, skip_header=1).T #transpose
                       #(0,2,4, etc.) to get only dates and skip the dashes, 
                       #colons, and the T between dates and times 

#m_data is a matrix whose column 1 are the dates and times of each event, from
#the most recent (column 1) in the dataset to the oldest (last column). If the 
#matrix were not transposed, the dates and times would be arranged by rows  


nEq= len(m_dates[0])  #the number of elements in a row (here I took the first)
                      #is equal to the number of earthquakes in the dataset 

#load times
a_decYr = np.zeros(nEq) #initialization of the vector of decimal years (to 0)

for i in range(nEq):
    a_decYr[i] = seis_utils.dateTime2decYr( m_dates[:,i]) 
    #the function dateTime2dcYr takes as input a date-time column, converts the 
    #days that have passed in the given year (till the given date) into secs, 
    #as well as the hours and minutes in the given day, and computes the 
    #decimal part of the year by dividing by the total number of seconds in 
    #a year: 24*3600*356 (or 366 for leap years)
    
#location and magnitude 
m_LocMag = np.genfromtxt (file_in, delimiter = ",", skip_header=1, 
                          usecols=(1,2,3,4), 
                          dtype = float).T  #transpose
                          #here I collect the other columns, from the second 
                          #to the fifth, from the file_in, i.e., the excel 
                          #file with all the earthquakes 

# here I put together the dec years with location & magnitude, column-wise                      
m_Eq = np.vstack ((a_decYr, m_LocMag)) 


#sort in ascending time order  
sel = np.argsort( m_Eq[0])
m_Eq = m_Eq.T[sel].T

# density of events over time 
a_tbin, a_rate = seis_utils.compRate (m_Eq[0])

plt.figure(1)
ax=plt.subplot (211)
ax.scatter( a_decYr,m_LocMag[3], s=np.exp(m_LocMag[3]), marker ='+', c= m_LocMag[3]) #s is the size of the marker 
ax.set_xlabel( 'Year')
ax.set_ylabel('Magnitude')
plt.title('Magnitude observations')
plt.xscale('linear')
plt.yscale('linear')
plt.subplots_adjust(bottom=-3.5, right= 3.0, top=2.9)

#compute the earthquake rate using the function from seis_utils.py
a_tbin, a_rate = seis_utils.compRate (m_Eq[0])

a_cumul = np.cumsum( np.ones ( nEq ))


plt.figure(2)
ax1=plt.subplot(311)
ax1.plot(a_tbin, a_rate, 'b-')
ax1.set_xlabel( 'Time (Yr)')
ax1.set_ylabel('Eq Rate')
ax1.set_title('Eq-rate', fontsize=18)

ax2=plt.subplot(312)
ax2.plot(m_Eq[0], a_cumul, 'r-')
ax2.set_xlabel( 'Time (Yr)')
ax2.set_ylabel('Cumulative EQs')
ax2.set_title('Cumulative EQks over time', fontsize=18)
plt.subplots_adjust(bottom=-3.5, right= 3.0, top=2.9)

#==========Tofind  when the main shock is occurred and its magnitude
# MAIN SHOCK IS TAKEN TO BE DOUBLETS WHICH HAPPENED IN 2014 WITH M7.5
max_mag = 0 #initialization
ms_lat = 0 #initialization 
ms_log = 0 #initialization 

ms_decYr = m_Eq[0,0]  #initialization

for i in range(nEq):
   if m_Eq[4, i]> max_mag:
       max_mag = m_Eq[4, i]
       ms_decYr = m_Eq[0, i] 
       ms_lat =  m_Eq[1, i] 
       ms_log =  m_Eq[2, i] 

# To find Rmax
Rmax = 10**(0.1238*max_mag+0.983)

''' To look for the events within the area of radius Rmax centered at the 
main shock. First I consider the full time series, then only the events after
the main shock. To compute the distance of each earthquake from the location 
of the MS, the harversine formula is used, from the seis_utils.py set of 
functions, which computes the distance betweeen 2 points over the Earth's 
surface (approximated as a sphere of radius 6371 km) as a function of the 
angles of those points with respect to the center of the Earth, i.e., their 
latitudes and longitudes.  
'''
#dimesnions of m_Eq: the columns are the events, the rows the info about each 
d1, d2 = np.size(m_Eq, axis=0) , np.size(m_Eq, axis=1)

# initialization of the matrix where I save all the events within the circle 
# of radius Rmax around the location of the main shock, both before and after 
# shock. Since I do not know in advance the number of events, I use the full 
#size of the matrix that contains the full dataset to initialize the new matrix
m_EqR = np.zeros((d1, d2)) 

# initialization of the matrix where I save all the events within the circle 
# of radius Rmax around the location of the main shock, only after shock 
# Since I do not know in advance the number of events, I use the full size
# of the matrix that contains the full dataset to initialize the new matrix
m_EqRas = np.zeros((d1, d2)) 

dist = np.zeros(d2)

countR = 0 # number of events before, and aftershock within the radius Rmax
countRas = 0 # number of events aftershock within the radius Rmax

for i in range(nEq):
    #distance of event i from the location of the mainshock 
    dist[i] = seis_utils.haversine( ms_log, ms_lat, m_Eq[2, i], m_Eq[1, i])
    
    
    if dist[i] < Rmax:
        m_EqR[:,countR] = m_Eq[:,i]
        countR = countR+1 #update the number of events within the circle with 
                          #radius = Rmax   
        if m_Eq[0,i]>ms_decYr:   #to detect only aftershock events
            m_EqRas[:,countRas] = m_Eq[:,i]
            countRas = countRas+1 
            
            
# redoing the same procedure above, over the events within the radius Rmax,to plot only those events

# density of events over time 
a_tbinR, a_rateR = seis_utils.compRate (m_EqR[0, 0:countR])

plt.figure(3)
ax1=plt.subplot (311)
ax1.scatter(m_EqR[0, 0:countR], m_EqR[4, 0:countR], s=np.exp(m_EqR[4, 0:countR]), marker ='o', color ='k' ) #s is the size of the marker 

ax1.set_xlabel( 'Time (Yr)')
ax1.set_ylabel('Magnitude')

# ax2=plt.subplot(312)
# ax2.plot( a_tbinR, a_rateR, 'b-')
# ax2.set_xlabel( 'Time (Yr)')
# ax2.set_ylabel('Eq Rate')
plt.subplots_adjust(bottom=-3.5, right= 3.0, top=2.9)


'''
To convert the time of occurrence of the aftershock events into decimal  
years from the main shock 
'''     
m_EqRas2 = m_EqRas.copy()
m_EqRas3 = m_EqRas2.copy()

for i in range(countRas):
    m_EqRas2[0,i] = m_EqRas[0,i]- ms_decYr #time from the MS in decimal years

    #to compute the number of days, the difference in decimal years is suitably multipied by 365 or 366
    #depending on if the year is leap or not. 
    if m_EqRas2[0,i]<(2015-ms_decYr): #i.e., if we are within the same year of the mainshock, 2014
        m_EqRas3[0,i] = m_EqRas2[0,i]*(365+calendar.isleap( int(math.floor(m_EqRas[0,i])))) #num of days since the MS
    else: #i.e., if we are in 2015, the next year after the mainshock in 2014 
        DaysWithin2014 = (2015-ms_decYr)*(365+calendar.isleap( int(math.floor(m_EqRas[0,i])))-1) # -1 because we have to check the previous year, not the current 
        DaysWithin2015 = (m_EqRas[0,i]-2015)* (365+calendar.isleap( int(math.floor(m_EqRas[0,i]))))
        
        m_EqRas3[0,i] = DaysWithin2014 + DaysWithin2015
        

'''
Aftershock rates 
'''
# density estimate of aftershock events over time (expressed in years and then 
# in days after the main shock)
# the chosen value of k is sqrt(nEq), rounded to the closest integer (computed 
# within the function seis_utils.compRate)


# HERE EXPRESSED IN YEARS
a_tbinRas, a_rateRas = seis_utils.compRate (m_EqRas[0, 0:(countRas)])

plt.figure(4)
ax3=plt.subplot (411)
ax3.scatter(m_EqRas[0, 0:(countRas)], m_EqRas[4, 0:(countRas)], s=np.exp(m_EqRas[4, 0:(countRas)]), marker ='o', color ='k' ) #s is the size of the marker 

ax3.set_xlabel( 'Time (Yr)')
ax3.set_ylabel('Magnitude')

ax4=plt.subplot(412)
ax4.plot( a_tbinRas, a_rateRas, 'b-')
ax4.set_xlabel( 'Time (Yrs after MS)')
ax4.set_ylabel('Eq Rate')
ax4.set_xlim( ax3.get_xlim()) 
#plt.savefig('Aftershocks(Yrs).png')
plt.subplots_adjust(bottom=-3.5, right= 3.0, top=2.9)

# HERE EXPRESSED IN DAYS AFTER THE MAIN SHOCK 
a_tbinRas3, a_rateRas3 = seis_utils.compRate (m_EqRas3[0, 0:(countRas)])

plt.figure(5)
ax5=plt.subplot (511)
ax5.scatter(m_EqRas3[0, 0:(countRas)], m_EqRas3[4, 0:(countRas)], s=np.exp(m_EqRas3[4, 0:(countRas)]), marker ='o', color ='k' ) #s is the size of the marker 

ax5.set_xlabel( 'Time (Days after the MS)')
ax5.set_ylabel('Magnitude')

ax6=plt.subplot(512)
ax6.plot( a_tbinRas3, a_rateRas3, 'b-')
ax6.set_xlabel( 'Time (Days after the MS)')
ax6.set_ylabel('Eq Rate')
ax6.set_xlim( ax5.get_xlim()) 
#plt.savefig('Aftershock(days).png')
plt.subplots_adjust(bottom=-3.5, right= 3.0, top=2.9)

  
'''
Rate-decay fitting (time expressed in days)
'''

# m_EqRas3 contains the events information organized in terms of days after the 
# main shock, and it is limited to the aftershock events, and so a_tbinRas3, 
# which has been derived from m_EqRas3. 
# Further, some events after the first 650 days from the MS are neglected in order to fit the  
# power-law to points that clearly show such behavior

sel = a_tbinRas3 <= 1000

a_t_AS, a_rate_as = a_tbinRas3[sel], a_rateRas3[sel] 

#fit PL with linregress in the log-log plane 

p_slope, intercept, R2, sign, stdErr = scipy.stats.linregress(np.log10(a_t_AS), np.log10(a_rate_as))
K = 10**intercept

a_y_hat = K * a_t_AS**p_slope

plt.figure(6)
ax7= plt.subplot (111)
ax7.loglog( a_t_AS, a_rate_as, 'ko')
ax7.loglog( a_t_AS, a_y_hat, 'r--')
ax7.set_title( f"p={round(p_slope, 2)}")
ax7.set_xlabel('Time after MS (days)')
ax7.set_ylabel( 'Aftershock rate')



'''
# Histogram of Magnitude on 
# '''
# def plot_loghist(m_LocMag, bins):
#   hist, bins = np.histogram(m_LocMag, bins=bins)
#   logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#   plt.figure(7)
#   ax8=plt.hist(m_LocMag, bins=logbins)
#   ax8=plt.yscale('log')
#   #ax8=plt.xscale('log')

# plot_loghist(np.random.rand(200), 10)
plt.show  
