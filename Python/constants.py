#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Constantes for Ganglion OpenBCI
Created on Sun Sep 30 23:09:22 2018

@author: felipe
"""
import numpy as np
import queue



###Streaming data constants###
DATA_PACKET_LENGTH = 20
#NB_CHANNELS = 4 #4 available, only 2 used
NB_CHANNELS = 2
PACKET_FIRST_WORD = 0
SCALE_FACTOR=1.2 * 8.3886070 * 1.5 * 51.0 #uV
START_STREAMING_MSG = 'b'
STOP_STREAMING_MSG = 's'
bits19_PACKET_COUNTER=101
bits18_PACKET_COUNTER=0
bits19_18=True #True 19 bits, False 18 bits
RESET_MSG = 'D'
FS=200
CH_FZ=1
CH_CZ=2
counter_all=0
IS_CONNECTED=False


###Stimulation Constansts###
STIM_TYPE=0 #0 Auditory, 1 tactile, 2 ambas.
STIM_SHAM=True; #True stimulation, False: Sham
#Pass Band filter, high 2, low 4, notch 2, 0.3-38 delay=84ms @1Hz
#250Hz
#b_pb=np.array([[0.0171126856325302, 0.0235650904067252, -0.0213205617166703, -0.0235650904067252, 0.00841575216828040, -0.0235650904067252, -0.0213205617166703, 0.0235650904067252, 0.0171126856325302]])
#a_pb=np.array([[1, -4.08109468848783, 8.02093995686550, -10.1293306984435, 8.81133830527304, -5.30905061055516, 2.14376138917718, -0.512883423854458, 0.0563412046278591]])
#200Hz high 2, low 4, notch 2, 0.5-40 Hz delay=17samples=85ms @1hz
b_pb=np.array([[0.0399465786719009,0.0798931573438018,-1.38777878078145e-17,-0.0798931573438017,-0.0798931573438018,-0.0798931573438017,-1.38777878078145e-17,0.0798931573438018,0.0399465786719009]])
a_pb=np.array([[1,-2.76876674457128,3.94704711582018,-4.31693656839793,3.40384449524308,-1.91499893836337,0.802769300485467,-0.174438025232564,0.0215929126928629]])
pb_order=np.shape(b_pb)[1]
previous_pb_filter=np.zeros((pb_order-1,1));
previous_pb=np.zeros((pb_order,1));
#Delta Band filter 0.5-4Hz delay=93samples=371.6ms @1hz
#250Hz
#b_delta=np.array([[4.16214652680027e-06, 0, -1.66485861072011e-05, 0, 2.49728791608016e-05, 0, -1.66485861072011e-05, 0, 4.16214652680027e-06]])
#a_delta=np.array([[1, -7.72571009713194, 26.1205925064617, -50.4804403199552, 60.9933246120536, -47.1808884342465, 22.8181453235623, -6.30830655127718, 0.763282960534981]])
#200Hz 0.5-4.5Hz Butterwoth 4 with Chebyshev 4 0.5-2.5 Apass=0.1, delay=66 samples=330ms @1Hz
#a_delta=np.array([[1,-7.66558215121082,25.7206243825777,-49.3403887680484,59.1883282386414,-45.4665279188275,21.8413615997012,-5.99917979292508,0.721364410097870]])
#b_delta=np.array([[1.08718641090292e-05,0,-4.34874564361169e-05,0,6.52311846541754e-05,0,-4.34874564361169e-05,0,1.08718641090292e-05]])
#200Hz 0.5-4Hz Butterwoth 4 with Chebyshev I 4 0.6-2 Apass=0.1, delay=66 samples=330ms @1Hz
a_delta=np.array([[1,-7.74086774596518,26.2266008150716,-50.7978912113993,61.5210173931342,-47.7067964696863,23.1324077317930,-6.41257023736723,0.778099724428693]])
b_delta=np.array([[3.76226851507215e-06,0,-1.50490740602886e-05,0,2.25736110904329e-05,0,-1.50490740602886e-05,0,3.76226851507215e-06]])
delta_order=np.shape(b_delta)[1]
previous_delta_filter=np.zeros((delta_order-1,1))
previous_delta=np.zeros((delta_order,1))

# Envelope filter, 2 seconds
#a_envolture=np.array([[1, -3.96716259594885, 5.90202586149088, -3.90255878482324, 0.967695543813138]])
#b_envolture=np.array([[1.53324552060196e-09, 6.13298208240784e-09, 9.19947312361176e-09, 6.13298208240784e-09, 1.53324552060196e-09]])
#Envelope filter Low pass 4 order, Buttherworth delay=156 samples=780 ms @0.7 Hz
a_envelope=np.array([[1,-3.93432582079874,5.80512542105514,-3.80723245722885,0.936433243152019]])
b_envelope=np.array([[2.41362231351615e-08,	9.65448925406460e-08,	1.44817338810969e-07,	9.65448925406460e-08,	2.41362231351615e-08]])
envelope_order=np.shape(b_envelope)[1]
previous_envelope_filter=np.zeros((envelope_order-1,1));
previous_envelope=np.zeros((envelope_order,1));
#envolture_90_previous=np.zeros((envolture_order-1,1));

#alpha filter Chebysev 8th order, 8-13 Apass=0.01
a_alpha=np.array([[1,	-7.15354063480188,	22.7664593754429,	-42.0784389223266,	49.3860894782566,	-37.6874282288573,	18.2638216857298,	-5.14077530323027,	0.643898799241065]])
b_alpha=np.array([[7.96551308373823e-05,	0,	-0.000318620523349529,	0,	0.000477930785024294,	0,	-0.000318620523349529,	0,	7.96551308373823e-05]])
alpha_order=np.shape(b_alpha)[1]
previous_alpha_filter=np.zeros((alpha_order-1,1));
previous_alpha=np.zeros((alpha_order,1));
#beta filter 16-30  Chebysev 8th order, Apass=0.01
a_beta=np.array([[1,	-5.10477314870557,	12.5595158089290,	-19.1726194329363,	19.7882955681030,	-14.1016221716494,	6.78924424539719,	-2.02849353188025,	0.294131897567411]])
b_beta=np.array([[0.00335604138270074,	0,	-0.0134241655308030,	0,	0.0201362482962045,	0,	-0.0134241655308030,	0,	0.00335604138270074]])
beta_order=np.shape(b_beta)[1]
previous_beta_filter=np.zeros((beta_order-1,1));
previous_beta=np.zeros((beta_order,1));

#Theta filter 4.5-7 Chebysev 8th order, Apass=0.01
a_theta=np.array([[1,	-7.65330041228977,	25.7499186004591,	-49.7443308646542,	60.3473492708559,	-47.0778797316863,	23.0635405077905,	-6.48758409532795,	0.802287553418907]])
b_theta=np.array([[5.54855578104519e-06,	0,	-2.21942231241808e-05,	0,	3.32913346862711e-05,	0,	-2.21942231241808e-05,	0,	5.54855578104519e-06]])
theta_order=np.shape(b_beta)[1]
previous_theta_filter=np.zeros((theta_order-1,1));
previous_theta=np.zeros((theta_order,1));



a_low=np.array([[1,	-2.99854927967837,	2.99709961146064,	-0.998550331400904]])
b_low=np.array([[4.76712100408731e-11,	1.43013630122619e-10,	1.43013630122619e-10,	4.76712100408731e-11]])
low_order=np.shape(b_low)[1]
previous_delta_power_filter=np.zeros((low_order-1,1))
previous_delta_power=np.zeros((low_order,1))
previous_theta_power_filter=np.zeros((low_order-1,1))
previous_theta_power=np.zeros((low_order,1))
previous_alpha_power_filter=np.zeros((low_order-1,1))
previous_alpha_power=np.zeros((low_order,1))
previous_beta_power_filter=np.zeros((low_order-1,1))
previous_beta_power=np.zeros((low_order,1))
previous_eeg_power_filter=np.zeros((low_order-1,1))
previous_eeg_power=np.zeros((low_order,1))

#SWS detection
max_bd=0.04
min_delta=0.4
min_alpha_bd_diff=0.004
previous_Psws=0.2
Psws=0

#Flags decoding

flag_packet0=False
flag_order=True
decoding_case=0



#Flags stimulating
flag_shot=False
flag_hold=False
flag_search=False
flag_pause=False
flag_sws=False
flag_closed_eyes=False
flag_forced=False
count_forced=0

pulses=0




#Phase detector and SWS detector
origen_eeg=np.zeros((2,1))
shiftphase_eeg=np.zeros((2,1))
ENVELOPE_OFFSET=10 #uV
SHIFTPHASE_GAIN=0.8/0.01885; #Mean magnitude, for normalization; 0.8 from 200/250
STIM_ANGLE=np.pi/4 #rad
EPSILON_ANGLE=5*np.pi/180 #rad
SO_THRESHOLD=-40 #uV
DELTA_THRESHOLD=10 #uV
PHASE_ADD=0
PHASE_MULTIPLIER=1 
MIN_EEG=-400 #uV
PAUSE_TIME=2.5 #s
COUNT_SO=0
time_so=0
time_pause=0
sws_counter=0
closed_counter=0
opened_counter=0




