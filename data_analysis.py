#!/usr/bin/python

'''
data_analysis.py : Graphical viewer of openQCM NEXT data curves
date created : 2021/10/26
@author : Luca PACCARD (paccard@insa-toulouse.fr)
'''


import matplotlib.pyplot as plt
import csv
import pandas as pd
import sys
from scipy.interpolate import make_interp_spline
import numpy as np

filename = str(sys.argv[1])

file = open(filename)
csvreader = csv.reader(file)
header = next(csvreader)
data = pd.read_csv(filename)

##################################################################################
# SINGLE SCAN
##################################################################################
if len(header) == 6:

    t = data.Relative_time

    # Plot of Resonance Frequency
    f = data.Resonance_Frequency

    plt.figure('Single_Scan_Resonance_Frequency')
    plt.plot(t, f)
    plt.xlabel("Relative Time (s)")
    plt.ylabel("Resonance Frequency (Hz)")
    plt.title('Resonance Frequency (Single Scan)')
    
    # Plot of Dissipation
    d = data.Dissipation
    plt.figure('Single_Scan_Dissipation')
    plt.plot(t, d)
    plt.xlabel("Relative Time (s)")
    plt.ylabel("Dissipation")
    plt.title('Dissipation (Single Scan)')

    plt.show()


##################################################################################
# MULTI SCAN
##################################################################################
elif len(header) == 14:

    t = data.Relative_time

    # Plot of Resonance Frequency
    f0 = data.Frequency_0
    f1 = data.Frequency_1
    f2 = data.Frequency_2
    f3 = data.Frequency_3
    f4 = data.Frequency_4

    plt.figure('Multi_Scan_Resonance_Frequency')
    
    plt.plot(t, f0, 'r', label='Frequency0 (Hz)')
    plt.plot(t, f1, 'b', label='Frequency1 (Hz)')
    plt.plot(t, f2, 'g', label='Frequency2 (Hz)')
    plt.plot(t, f3, 'c', label='Frequency3 (Hz)')
    plt.plot(t, f4, 'm', label='Frequency4 (Hz)')
    
    plt.xlabel("Relative Time (s)")
    plt.title('Resonance Frequency (Multi Scan)')
    plt.legend()
    
    # Plot of Dissipation
    d0 = data.Dissipation_0
    d1 = data.Dissipation_1
    d2 = data.Dissipation_2
    d3 = data.Dissipation_3
    d4 = data.Dissipation_4
    
    plt.figure('Multi_Scan_Dissipation')

    plt.plot(t, d0, 'r', label='Resonance0 (Hz)')
    plt.plot(t, d1, 'b', label='Resonance1 (Hz)')
    plt.plot(t, d2, 'g', label='Resonance2 (Hz)')
    plt.plot(t, d3, 'c', label='Resonance3 (Hz)')
    plt.plot(t, d4, 'm', label='Resonance4 (Hz)')
    
    plt.xlabel("Relative Time (s)")
    plt.title('Dissipation (Multi Scan)')
    plt.legend()

    plt.show()

else :
    print('Wrong file format')