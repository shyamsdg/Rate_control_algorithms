# Overview
An Ideal Rate adaptation algorithm (RAA) initially creates a table of Signal to Noise Ratio (SNR) and Modulation Coding Scheme (MCS) pairs. The SNR thresholds in this table ensure selecting an MCS that leads to a Packet Error Rate (PER) below a certain value. The SNR is fed back from the receiver to the transmitter via a perfect out-of-band mechanism. 

The main drawback of this mechanism is the use of an out-of-band channel for sending back the feedback which is not available in Industrial, Scientific, and Medical (ISM) bands used by IEEE 802.11. 

# About this code
Ideal RAA is implemented in MATLABT, where it is been designed for a PER<=10% condition.

The folder contains the following files:
1. `Ideal_RAA_Live_Script.mlx` : MATLAB live editor version of the Ideal RAA algorithm. Implemented in this format for easier understanding of different blocks of the algorithm.
2. `Ideal_RAA_MATLAB_Code.m` : Ideal RAA algorithm implemented in MATLAB script file. This is used to run and evaluate the performance of the algorithm.
3. `helperFrequencyOffset.m` : This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.
4. `helperNoiseEstimate.m` : This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.
5. `vhtNoiseEstimate.m` : This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.
6. `vhtSingleStreamChannelEstimate.m` : This is a helper function from MATLAB, required to run the MATLAB/Live_Script code.

# More Information
The comments in the code explains the detailed implementation of the algorithm, use Live Script for easier understanding.

## Contributing
    1. Clone repo and create a new branch
    2. Make changes and test
    3. Submit pull request with comprehensive description of changes