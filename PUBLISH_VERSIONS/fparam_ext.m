% Extracting Effective/Equivalent Refractive Index Model of a Metamaterial
%
% Function to Extract Effective and Equivalent Parameters
%
% Author: Ekin Gunes Ozaktas, December 2022
%
% This function is called by effective_dsweep.m, and calculates the
% effective or equivalent refractive index model for an inhomogeneous
% material using continuity after starting from the branch with the lowest 
% absolute mean derivative. The inputs are the data file containing the
% simulated S parameters, the name to be given the the file the effective/
% equivalent model is to be stored in, the geometric thickness of the
% metamaterial, and optionally a file containing a true data file and 
% the binary variable to toggle this comparison. The outputs are the vector
% of wavelengths and the effective/equivalent refractive index.
%
% Acknowledgements: Parts of this code take as a basis the code by Zsolt 
% Szabó, to which a link is available in Szabó et al., "A Unique Extraction
% of Metamaterial Parameters Based on Kramers–Kronig Relationship", IEEE 
% Transactions on Microwave Theory and Techniques, November 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [wvl_nm, n_eff] = fparam_ext(datafilename, name, d);
load(datafilename); % File with simulation data.
uf = 1e15; % Unit of Frequency, eg. 1e15 for PHz.
uf_text = "PHz";
c = 2.998e8; % Speed of light in vacuum, m/s.

% Obtaining data from simulation data file. (This is currently adapted for
% Lumerical but can be easily reworked by only changing this portion of the code.) 
Frequency_vals = S.f/uf;
szfreq = length(Frequency_vals);
omega = 2*pi*Frequency_vals*uf; % Angular frequency, rad/s
k0 = omega./c; % Wavenumber in free space, rad/m
S11 = squeeze(S.S11_Gn);
S21 = squeeze(S.S21_Gn);

% Calculation of wave impedance.
Z_eff = sqrt( ((1.0 + S11).^2 - S21.^2)./((1.0 - S11).^2 - S21.^2) );
exp_ink0d = S21./( 1.0 - S11.*((Z_eff - 1.0)./(Z_eff + 1.0)));

% Choice of correct sign for the real part of the wave impedance
% and the imaginary part of the refractive index. This chunk was adapted
% from Szabo et al.
for i = 1:szfreq
    if (abs(real(Z_eff(i))) > 1e-2) % if large enough
        if (real(Z_eff(i)) < 0)
            Z_eff(i) = -Z_eff(i);
            exp_ink0d(i) = S21(i)./( 1.0 - S11(i).*((Z_eff(i) - 1.0)./(Z_eff(i) + 1.0)));
        end;
    else
        if (abs(exp_ink0d(i)) > 1.0)
            Z_eff(i) = -Z_eff(i);
            exp_ink0d(i) = S21(i)./( 1.0 - S11(i).*((Z_eff(i) - 1.0)./(Z_eff(i) + 1.0)));
        end;
    end;
end;

% Calculation of the primary branch of the effective refractive index.
n_eff_0_comp = (imag(log(exp_ink0d)) - 1i*real(log(exp_ink0d)))./(k0*d); % Complex refractive index
k_eff = imag(n_eff_0_comp); % Imaginary part of refractive index.
n_eff_0 = real(n_eff_0_comp); % Real part of refractive index.

% Finding the branch that starts off with the minimum absolute mean
% value of the derivative, as it is predicted to be the correct one. The
% bounds of where to search are based on what can be considered reasonable
% values.
lim = 5; % The index of the value until which to check the derivatives.
low_m = -5; % Lower bound from which to start searching.
high_m = round(2*(d/(200e-9/10))); % Upper bound for search, 
% taken as roughly 2 times the number of periods traversed in the material 
% at a wavelength of 200nm (a reasonable lower limit) and assuming n is 
% approximately no greater than 10. 
branch_num = 0; % Default
min_mean_diff = abs(mean(diff(n_eff_0(1:lim))));
for k = low_m:high_m
    test_n = n_eff_0 + 2.0*pi*k./(k0*d);
    mean_diff(k - low_m + 1) = abs(mean(diff(test_n(1:lim))));
    if mean_diff(k - low_m + 1) < min_mean_diff
        min_mean_diff = mean_diff(k - low_m + 1);
        branch_num = k;
    end
end

% Produce all the branches necessary
n_branches = zeros(length(branch_num:high_m), length(n_eff_0));
for k = low_m:high_m
    n_branches(k - low_m + 1, :) = n_eff_0 + 2.0*pi*(k)./(k0*d);
end

% Go through branches with continuity, changing the branch number to
% whichever of the current branch or the nearest neighbors are closest to
% the previous value of the refractive index.
re_n_eff = n_eff_0 + 2.0*pi*branch_num./(k0*d); % Initializing real refractive index.
m = branch_num;
val_n_old = re_n_eff(1);
for k = 1:length(re_n_eff)
    sz = size(n_branches);
    if m-low_m+1 == 1 % First branch, only check itself and above.
        [minimum, loc] = min(abs([n_branches(m - low_m + 1,k) - val_n_old, n_branches(m - low_m + 2,k) - val_n_old]));
        m_add = loc - 1;
    elseif m-low_m+1 == sz(1) % Last branch, only check itself and below.
        [minimum, loc] = min(abs([n_branches(m - low_m,k) - val_n_old, n_branches(m - low_m + 1,k) - val_n_old]));    
        m_add = loc - 2;
    else % Other branch, check itself, below, and above.
        [minimum, loc] = min(abs([n_branches(m - low_m,k) - val_n_old, n_branches(m - low_m + 1,k) - val_n_old, n_branches(m - low_m + 2,k) - val_n_old]));
        m_add = loc - 2;
    end
    m = m + m_add; % Change branch if necessary.
    branch_nums(k) = m;
    re_n_eff(k) = n_eff_0(k) + 2.0*pi*m./(k0(k)*d); % New value of real refractive index.
    val_n_old = re_n_eff(k); % Old value of real refractive index for new iteration.
end

n_eff = re_n_eff + 1j*k_eff; % Creation of complex refractive index.

% Saving data to file.
wvl_nm = (c./(omega./(2*pi)))*1e9;
exp_n = real(n_eff);
exp_k = imag(n_eff);
save(name, "wvl_nm", "exp_n", "exp_k");


