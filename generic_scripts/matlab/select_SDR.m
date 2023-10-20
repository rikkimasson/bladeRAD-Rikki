function [serial_No] = select_SDR(SDR_No)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch SDR_No
    case 1
        serial_No = '4eedd4abca9b48a0a49d2a469f09b719';
    case 2
        serial_No = '6469e0a2ac0f4597937fa51de169675c';
    case 3
        serial_No = 'a0193abefb6f43b2aa7c24ed137ce1f5';
    case 4
        serial_No = '6aa6c4dc3dbb493fbc7ff4c2560c761d';
    case 10
        serial_No = 'cdb288010b364eb389ef1fa76db6e492';     
    otherwise 
        disp("SDR Number doesn't exist - select valid number between 1-3 ")
end
