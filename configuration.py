#!/usr/bin/env python
# coding: utf-8

# ## Demo introduction and preparition for single qubits characterization

# This demo file is for basic qubit charaterization. The including exerimental sequences are
# - Resonator spectroscopy
# - Qubit spectroscopy
# - Rabi flopping
# - T1 and T2
# - Mixer calibration
#
# Required instruments
# - HDAWG and UHFQA
# - LO MW sources for both qubit driving and qubit readout
# - Frequency up and down conversion units
# - Control PC, USB or LAN connection
#
# Recommanded connections for HDAWG and UHFQA
# - Ref clc out of UHFQA connected to Ref clc in of HDAWG
# - Trigger output 1 of HDAWG to Ref/Trigger 1 of UHFQA
# - Enable DIO connection if it is necessary
#

# ## Import Modules

# In[4]:


#  ---------- keyboard shortcuts in 'Help' tap ------------
# restart kernel when somthing changed in subfunctions

# from VISAdrivers import sa_api as sa
import time
import importlib.util
import json
import sys, os
import LO845m as LO
import numpy as np
import UHFQA as qa
import HDAWG as hd
import experiment_funcs as expf
import matplotlib.pyplot as plt
import csv
import glob
import scipy as scy
import plot_functions as pf
import h5py
from VISAdrivers.LabBrick_LMS_Wrapper import LabBrick_Synthesizer
pi = np.pi
# from VISAdrivers.vaunix_attenuator_wrapper import VaunixAttenuator

'''Instruments and connection'''
qa_id = 'dev2528'
awg_id = 'dev8233'
meas_device = "CandleQubit_5"

# Instrument Addresses
qubitLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0538::0::INSTR"
readoutLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0519::0::INSTR"
acStarkLO = 21841
# readout_attn = 26551

# initialize instruments as python objects
qubitLO = LO.LO(address=qubitLO_IP,reset=False)
readoutLO = LO.LO(readoutLO_IP,reset=False)
acStarkLO = LabBrick_Synthesizer()
acStarkLO.initDevice(21841)
# readout_attn = LabBrick_attn()
# readout_attn.initDevice(26551)

qubitLO.RF_ON()
readoutLO.RF_ON()
acStarkLO.setRFOn(bRFOn=True)

qubitLO.set_freq(3.21)
readoutLO.set_freq(7.2275)
acStarkLO.setPowerLevel(9.0)
acStarkLO.getUseInternalRef()
acStarkLO.setFrequency(7.3586e9)

'''Initialize connection with Zurich Instruments'''
daq, device_qa = qa.create_api_sessions_uhf('dev2528', use_discovery= 1, ip='127.0.0.1')
awg, device_awg = hd.create_api_sessions_hd('dev8233', use_discovery= 1, ip = '127.0.0.1')

# set clocks to 10 MHz reference
daq.setInt('/dev2528/system/extclk', 1)
awg.setInt('/dev8233/system/clocks/referenceclock/source', 1)

'''Channel offsets'''
# read the current channel offset values and store them for future reference in case the QA needs power cycling
# qubit mixer
offset_qubit_ch1 = awg.get('/dev8233/sigouts/0/offset')['dev8233']['sigouts']['0']['offset']['value']
offset_qubit_ch2 = awg.get('/dev8233/sigouts/2/offset')['dev8233']['sigouts']['2']['offset']['value']

# AC stark mixer
offset_ac_stark_ch1 = awg.get('/dev8233/sigouts/1/offset')['dev8233']['sigouts']['1']['offset']['value']
offset_ac_stark_ch2 = awg.get('/dev8233/sigouts/3/offset')['dev8233']['sigouts']['3']['offset']['value']

# readout mixer
offset_readout_ch1 = daq.get('/dev2528/sigouts/0/offset')['dev2528']['sigouts']['0']['offset']['value']
offset_readout_ch2 =  daq.get('/dev2528/sigouts/1/offset')['dev2528']['sigouts']['1']['offset']['value']

awg.setDouble('/dev8233/sigouts/0/offset',offset_qubit_ch1)
awg.setDouble('/dev8233/sigouts/2/offset',offset_qubit_ch2)
awg.setDouble('/dev8233/sigouts/1/offset',offset_ac_stark_ch1)
awg.setDouble('/dev8233/sigouts/3/offset',offset_ac_stark_ch2)

daq.setDouble('/dev2528/sigouts/0/offset',offset_readout_ch1)
daq.setDouble('/dev2528/sigouts/1/offset',offset_readout_ch2)

