% Extracting Effective/Equivalent Refractive Index Model of a Metamaterial
%
% Function for FDTD Spectra Extraction
%
% Author: Ekin Gunes Ozaktas, December 2022
%
% This function is called by effective_dsweep.m, and calculates the
% reflectance, transmittance, and absorptance for an inhomogeneous
% material from simulations. The inputs are the data file containing the
% simulated reflectance and transmittance. The outputs are the vector
% of frequencies and the transmittance, reflectance, and absorptance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq, transmitted, reflected, absorbed] = f_abs_spectra(datafilename)
load(datafilename)
freq = T.f;
reflected = squeeze(T.T); % Due to monitor naming conventions in the .fsp files, the naming here is counterintuitive.
transmitted = squeeze(R2.T);
absorbed = reflected - transmitted; % These definitions may change depending on the way the simulation is structured.
c = 3e8;

reflected = 1-reflected;


