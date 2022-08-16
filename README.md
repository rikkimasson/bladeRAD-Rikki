# bladeRAD - Multi-functional Radar
This repository contains the necessary control scripts to Nuand bladeRF micro 2.0 Software Defined Radios as a multifunctional active and passive radar.  

# Introduction to Project 
SDRs provide highly flexible platforms that can be used to perform multiple functions at multiplefrequency bands through rapid reconfiguration of hardware via software. The cost-efficiency and ease of use make SDRs theoretically very suitable for experimental research in active and passive networked radar. The Nuand bladeRf micro 2.0 is a low cost highly capable SDR available on the COTS market. Of considerable interest, is the bladeRF’s high tuning range of up to 6.0 GHz, surpassing the LimeSDR, XTRX CS and Kerbos SDR. The bladeRF micro 2.0 A4 features 2x2 duplex MIMO architecture, each with a maximum bandwidth of 61.44 MHz. Control and programming of the SDR is most commonly achieved through the USB 3.0 interface between the on board Cypress FX3 micro-controller and a PC. The FX3 micro controller is used to administer control to all other components on the SDR board e.g. the VCXO and AD9361 RFIC chip. Support is available for controlling the board with GNU Radio via gr-osmosdr, Pothos via SoapySDR, SDRange, SDR Console, SDR via sdrsharp-bladeRF, MathWorks MATLAB.

## Top level description
This repository contains code for both the command and control of the bladeRF SDRs, in addition to code for conducting radar signal processing on the data captured by the SDRs. All processing is currently completed post-capture as the receive channel is saved to RAM for post-processing in MATLAB. 
There is currently code in this repository to control the radar as an:

1. FMCW active radar, using two synchronised SDRs - 1 for Tx, 1 for Rx - this is done to reduce the direct-path breakthrough between transmit and receive channels and reduce bandwidth requirement on the SDR to PC USB-3 conncection. 
2. Passive Radar, using a single SDR. 
3. Hybrid Radar, using three synchronised SDRs. Two for active FMCW radar and a single SDR for passive radar.  

## Installation and Usage

Requirements 
  1. Linux OS
  2. bladeRF CLI - https://github.com/Nuand/bladeRF/wiki/Getting-Started%3A-Linux
  3. Gnome-terminal - $ sudo apt install gnome-terminal
  4. MATLAB

Clone this repository in a folder named 'repos' in your home directory '~/home/repos/bladeRAD' - important to get this location correct. 

