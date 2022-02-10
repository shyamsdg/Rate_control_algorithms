# Overview
    An Ideal RAA initially creates a table of SINR and MCS pairs. The SINR thresholds in this table ensure selecting an MCS that leads to a BER below a certain value. For example, the default value is 10^-5), the SINR is fed back from the receiver to the transmitter via a perfect out-of-band mechanism. 
    
    The main drawback of this mechanism is the use of an out-of-band channel for sending back the feedback which is not available in Industrial, Scientific, and Medical (ISM) bands used by IEEE 802.11. 

# About this code
    This Matlab code executes the Ideal Rate adaptation algorithm where it is been designed for a PER<=10%.

# More Information
    The comments in the code explains the detailed implementation of the algorithm, use Live Script for easier understanding.