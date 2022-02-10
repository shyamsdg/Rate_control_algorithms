# Overview
    Minstrel(nonHT) is a rate control algorithm implemented in MadWifi and Linux. 

    The basic principle is to probe the environment and adapt the rate based on statistics collected on the probability of successful transmission. The algorithm adapts the rate to the highest rate that it considers successful, and spends a fraction of its time doing 'look around' by trying other rates. 

    Minstrel(nonHT) is appropriate for non-HT configurations; for HT (i.e. 802.11n or higher), users should use Minstrel-HT. Minstrel(nonHT) will error exit if the user tries to configure it with a Wi-Fi MAC that supports 802.11n or higher.

# About this code
    This Matlab code executes the Minstrel-HT algorithm, where NS3-MinstrelHTWifimanager was used to decode, calibrate and implement in Matlab.

# More Information
    Attached a detailed Flowchart in PDF format in Github; The comments in the code explains the detailed implementation of the algorithm, use Live Script for easier understanding.