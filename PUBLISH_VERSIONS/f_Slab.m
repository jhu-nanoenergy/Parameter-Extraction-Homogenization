% Extracting Effective/Equivalent Refractive Index Model of a Metamaterial
% 
% Function for Analytical Solution of Slab
% 
% Author: Ekin Gunes Ozaktas, December 2022
%
% This function is called by effective_dsweep.m, and calculates the
% reflectance, transmittance, and absorptance for a homogeneous slab. 
% The inputs are a matrix containing refractive index information, a vector
% of wavelengths, and the thickness. The outputs are the reflectance, 
% transmittance, and absorptance of the slab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TMM_transmitted, TMM_reflected, TMM_absorbed] = f_Slab(nk_name, lambda, d)
ul = 1e-9; % Unit of wavelength, 1e-9 for nm, etc.

% Raw data.
nkmatrix = load(nk_name);
nold = fliplr([ones(1,length(nkmatrix.exp_n)); (nkmatrix.exp_n)' + 1i*(nkmatrix.exp_k)']);
wvl_nm = flipud(nkmatrix.wvl_nm);
n(1,:) = interp1(wvl_nm, nold(1,:), lambda);
n(2,:) = interp1(wvl_nm, nold(2,:), lambda);

r12 = (n(2,:) - n(1,:))./(n(2,:) + n(1,:));
r21 = (n(1,:) - n(2,:))./(n(1,:) + n(2,:));
t12 = (2*n(1,:))./(n(2,:) + n(1,:));
t21 = (2*n(2,:))./(n(2,:) + n(1,:));

delta = 2 .* pi .* d.* (1./(lambda.*1e-9)) .* sqrt((n(2,:)).^2);

r123 = (r12 + r21.*exp(2*1i*delta))./(1 + r12.*r21.*exp(2*1i*delta));
t123 = (t12 .* t21 .* exp(1i*delta))./(1 + r12.*r21.*exp(2*1i*delta));

TMM_transmitted = (abs(t123)).^2;
TMM_reflected = (abs(r123)).^2;

TMM_absorbed = 1 - TMM_reflected - TMM_transmitted;














