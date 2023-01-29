# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 19:29:32 2022

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import NonlinearFitting
import os
import time


#%% example data loading
def SelectPeaks(folder='data',datetime_start='0-0',datetime_end='np.inf-np.inf',fork_name='auxfork-filterbox'):
    """
    For data of format used in Helsinky
    
    """
    months = glob(os.path.join(folder,'*'))
    fqs = []
    xs = []
    ys = []
    drives = []
    filenames = []
    for month in months:

        days = glob(os.path.join(month,'*'))
        for day in days:

            peaks_w_metadata = glob(os.path.join(day,'*'+fork_name+'.dat'))
            peaks = []
            for file in peaks_w_metadata:
                if 'T' in file:
                    peaks.append(file)
                else:
                    metadata_path = file
            metadata = np.loadtxt(metadata_path,dtype='str')
            

            for i,peak in enumerate(peaks):
                temp_string_file = peak[::-1].split('T')
                datetime = int(temp_string_file[0][:-7:-1])+int(temp_string_file[1][7::-1])*10**6
                
                temp_string_input1 = datetime_start.replace(',',':').replace('/',':').replace(' ',':').replace('-',':').split(':')
                temp_string_input2 = datetime_end.replace(',',':').replace('/',':').replace(' ',':').replace('-',':').split(':')
                
                t_start = t_end = 0
                j=10
                for part1,part2 in zip(temp_string_input1,temp_string_input2):
                    t_start += int(part1)*10**j
                    t_end += int(part2)*10**j
                    j-=2

                if j!=-2:
                    raise Exception('Incorect time format')
                    
                if datetime >= t_start and datetime <= t_end:
                    data = np.loadtxt(peak,comments='#')
                    fqs.append(np.copy(data[:,1]))
                    xs.append(np.copy(data[:,2]))
                    ys.append(np.copy(-data[:,3]))

                    filenames.append(peak[10:])
                    drives.append(float(metadata[i,5]))
    return fqs,xs,ys,drives,filenames


plt.close('all')

start_time = r'2022/08/18-01:24:21'
end_time = r'2022/08/18-05:04:24'
"""
Here is where data loading happens. Feel free to supply your own, fqs, xs, ys are lists arrays.
Drives, filenames are lists of numbers, str, respectively, both are optional
"""
fqs,xs,ys,drives,filenames = SelectPeaks(datetime_start = start_time ,datetime_end=end_time)
"""
Peak loading happens here ^^

start_time, end_time : string
year, month, day, hour, min, sec
each part must be separated by any of the following:
'/', '-', ':', ' ', ',' 

otherwise selection wont work properly
example: '2022-8,3:10 22/6'
important for string parsing

works on data with Helsinky format only
"""

#voltage velocity conversion
xs = [x*7.415e-3 for x in xs]
ys = [y*7.415e-3 for y in ys]

#%% important part, how loading data and fitting works
t_start = time.time()

Helsinky_test = NonlinearFitting.HelsinkyFit(fqs=fqs,xs = xs, ys=ys,drives=drives, filenames = filenames)
Helsinky_test.remove_transient() #optional method for strongly duffing datapoints, removes data between top and bottom of peaks
Helsinky_test.fit(vary_pars={'g3':False},init_values={'g3':0})

t_end = time.time()
print('Fit took',t_end - t_start,'seconds')


#%% this shows how figures are created
t_start = time.time()

Helsinky_test.figures(path='figures',extension='.png',show=False,save=~Helsinky_test.forwards)

t_end = time.time()
print('Figures took',t_end - t_start,'seconds')

#%% and here data are saved
Helsinky_test.save_params(folder='output_textfiles',name='08_18',which_peaks = ~Helsinky_test.forwards)
Helsinky_test.save_reports(folder='output_results',which_peaks= ~Helsinky_test.forwards)
