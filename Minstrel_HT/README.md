# Overview of Minstrel (nonHT) and Minstrel-HT
 Minstrel (nonHT) and **`Minstrel-HT`** is a rate control algorithm implemented in MadWifi and Linux. 

 The basic principle of  Minstrel (nonHT) is to probe the environment using Frame Success Ratio (FSR) and adapt the rate based on statistics collected about the probability of successful transmission. The algorithm adapts the rate to the highest rate that it considers successful, and spends a fraction of its time doing 'look around' by trying other rates (sampling). 

 Minstrel(nonHT) is appropriate for nonHT/legacy(i.e. 802.11/b/a/g) WiFi configurations. For HT WiFi configurations (i.e. 802.11n or higher), **`Minstrel-HT`** rate adaptation algorithm is to be used. 

 **`Minstrel-HT`** is designed for high-latency devices that implement a Multi Rate Retry (MRR) chain. This kind of device does not give feedback (i.e. open-loop system) for every frame/packet retransmission, but only when a frame/packet is correctly transmitted (an Ack is received) or a frame/packet transmission completely fails (all retransmission attempts fail). 
 
 The MRR chain is used to advise the hardware about which rate to use when re-transmitting a frame/packet. 
 
 In **`Minstrel-HT`**, the Sampling is done differently from nonHT/nonHT/legacy Minstrel. **`Minstrel-HT`** tries to sample all rates in all groups at least once and to avoid many  consecutive samplings. 

 For more details about MRR chain, sampling and the working of this algorithm and also its implementation in MATLAB please refer _About this code_ and _Source code and inspiration_.

# About this code
Minstrel-HT algorithm is implemented in MATLAB. The folder contains the following files:
1. `Minstrel HT_NS3_FV.pdf`: A detailed flowchart of the Minstrel-HT algorithm.
2. `MinstrelHT_Live_Script.mlx`: MATLAB live editor version of the Minstrel-HT algorithm. Implemented in this format for easier understanding of different blocks of the algorithm. 
3. `MinstrelHT_Matlab_Code.m`: Mintrel-HT algorithm implemented in MATLAB script file. This is used to run and evaluate the performance of the algorithm.
4. `helperFrequencyOffset.m`: This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.
5. `htNoiseEstimate.m`: This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.

# Source code and inspiration
Our Minstrel-HT algorithm implementation in MATLAB is inspired from NS3-MinstrelHTWifiManager. The NS3 version is used to debug and calibrate the performance of the MATLAB implementaion.

The NS3 source code can be found [here](https://www.nsnam.org/doxygen/classns3_1_1_minstrel_ht_wifi_manager.html). 

   