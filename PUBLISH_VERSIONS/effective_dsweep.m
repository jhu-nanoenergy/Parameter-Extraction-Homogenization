% Extracting Effective/Equivalent Refractive Index Model of a Metamaterial
% 
% Main Code to Run
% 
% Author: Ekin Gunes Ozaktas, December 2022
%
% This code takes in simulated data of the S parameters of a metamaterial
% and produces an effective or equivalent refractive index model. 
% The complex logarithm branching problem is resolved via choosing the 
% starting branch based on the starting derivative being minimized 
% (predicted homogeneity at long wavelengths) and then enforcing 
% continuity. A comparison is made between the simulated and modeled 
% slabs to determine the effective thickness.
% Plots of the final effective refractive index parameters varying with 
% wavelength are produced, alongside graphs of the reflectance, 
% transmittance, and absorptance of the simulated and modeled slabs.
%
% This code calls the f_abs_spectra.m, f_Slab.m, and fparam_ext.m
% functions.
%
% Acknowledgements: Parts of this code take as a basis the code by Zsolt 
% Szabó, to which a link is available in Szabó et al., "A Unique Extraction
% of Metamaterial Parameters Based on Kramers–Kronig Relationship", IEEE 
% Transactions on Microwave Theory and Techniques, November 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
datafilename = "Data_FDTD_Birefringent_Perpendicular.mat"; % File with simulation data
name = "Pattern-Model.mat"; % Name of file for effective index model
truefilename = "Si-Model.mat"; % Real data to compare with, if applicable (for testing purposes).
toggle_true_data = 0; % 0 or 1 depending on whether comparsion to true data for testing is desired.
ttoggle = 0; % 0 or 1 depending on whether figure titles are desired.
d_actual = 790e-9; % Thickness of material, in m.
% Range of thicknesses to test.
delta_d = d_actual/20; % Step size between thicknesses to test.
min_num = 6; % Minimum number of delta_d to go below the geometric thickness.
max_num = 30; % Maximum number of delta_d to go above the geometric thickness.
d_array = d_actual-min_num*delta_d:delta_d:d_actual+max_num*delta_d;
lambda = 400:1500; % Range of wavelengths, in nm.
c = 3e8; % Speed of light in vacuum, in m/s.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st iteration
d = d_array(1);
mse_total = zeros(1, length(d_array));
% Creation of effective/equivalent model.
[wvl_nm, n_eff] = fparam_ext(datafilename, name, d);
% Modelling the effective/equivalent slab.
[TMM_transmitted, TMM_reflected, TMM_absorbed] = f_Slab(name, lambda, d);
% Calculating the optical behavior (reflectance, transmittance, absorptance) for the inhomogeneous slab.
[freq, transmitted, reflected, absorbed] = f_abs_spectra(datafilename);
transmitted_fit = interp1(1e9*c./freq, transmitted, lambda);
reflected_fit = interp1(1e9*c./freq, reflected, lambda);
absorbed_fit = interp1(1e9*c./freq, absorbed, lambda);
% Calculation of MSE
mse_T = mean((TMM_transmitted - transmitted_fit).^2);
mse_R = mean((TMM_reflected - reflected_fit).^2);
mse_A = mean((TMM_absorbed - absorbed_fit).^2);
mse_total(1) = (mse_T + mse_R + mse_A)./3;
min_mse = mse_total(1);
min_loc = 1;

% Loop
for k = 2:length(d_array)
    d = d_array(k);
    % Creation of effective/equivalent model.
    [wvl_nm, n_eff] = fparam_ext(datafilename, name, d);
    % Modelling the effective/equivalent slab;
    [TMM_transmitted, TMM_reflected, TMM_absorbed] = f_Slab(name, lambda, d);
    % Calculating the optical behavior (reflectance, transmittance, absorptance) for the inhomogeneous slab.
    [freq, transmitted, reflected, absorbed] = f_abs_spectra(datafilename);
    transmitted_fit = interp1(1e9*c./freq, transmitted, lambda);
    reflected_fit = interp1(1e9*c./freq, reflected, lambda);
    absorbed_fit = interp1(1e9*c./freq, absorbed, lambda);
    % Calculation of MSE
    mse_T = mean((TMM_transmitted - transmitted_fit).^2);
    mse_R = mean((TMM_reflected - reflected_fit).^2);
    mse_A = mean((TMM_absorbed - absorbed_fit).^2);
    mse_total(k) = (mse_T + mse_R + mse_A)./3;
    if mse_total(k) < min_mse
        min_mse = mse_total(k);
        min_loc = k;
    end
end

% Reevaluate final case since effective thickness has been determined.
d = d_array(min_loc);
% Creation of effective/equivalent model.
[wvl_nm, n_eff] = fparam_ext(datafilename, name, d);
% Modelling the effective/equivalent slab;
[TMM_transmitted, TMM_reflected, TMM_absorbed] = f_Slab(name, lambda, d);
% Calculating the optical behavior (reflectance, transmittance, absorptance) for the inhomogeneous slab.
[freq, transmitted, reflected, absorbed] = f_abs_spectra(datafilename);
transmitted_fit = interp1(1e9*c./freq, transmitted, lambda);
reflected_fit = interp1(1e9*c./freq, reflected, lambda);
absorbed_fit = interp1(1e9*c./freq, absorbed, lambda);
% Calculation of MSE
mse_T = mean((TMM_transmitted - transmitted_fit).^2);
mse_R = mean((TMM_reflected - reflected_fit).^2);
mse_A = mean((TMM_absorbed - absorbed_fit).^2);
mse_total_end = (mse_T + mse_R + mse_A)./3;

% Kramers Kronig Transform of imaginary part of refractive index
omega = 2*pi*freq;
dw = omega(2) - omega(1); % Since frequencies are equally spaced
kappa = imag(n_eff);
n_KraKro = zeros(1,length(omega));
for k = 1:length(omega)
    sum_var = 0;
    for m = 1:length(omega)
        % Removing point that would result in infinity and replacing it with the previous point. 
        if omega(m) ~= omega(k)
            sum_var = sum_var + omega(m)*kappa(m)*dw ./ (omega(m)^2 - omega(k)^2);
        elseif m ~= 1
            sum_var = sum_var + omega(m-1)*kappa(m-1)*dw ./ (omega(m-1)^2 - omega(k)^2); % Use the previous value to avoid divergence
        elseif m ~= length(omega)
            sum_var = sum_var + omega(m+1)*kappa(m+1)*dw ./ (omega(m+1)^2 - omega(k)^2); % Use the previous value to avoid divergence
        end
    end
    n_KraKro(k) = 1 + (2/pi)*sum_var;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

font_size = 24;
figure() % Refractive index as a function of wavelength.
set(gcf,'Color','w');
hold on;
plot(wvl_nm, imag(n_eff), 'Color', 'r', 'LineWidth', 1.5);
plot(wvl_nm, real(n_eff), 'Color', 'b', 'LineWidth', 1.5);
legendtext = ["\kappa_{eff}", "{\itn}_{eff}"];
if toggle_true_data == 1
    plot(truefile.wvl_nm, truefile.exp_k, '--', 'Color', 'k', 'LineWidth', 1.5);
    plot(truefile.wvl_nm, truefile.exp_n, '--', 'Color', 'm', 'LineWidth', 1.5);
    legendtext = ["\kappa_{eff}", "{\itn}_{eff}", "\kappa_{ih}", "{\itn}_{ih}"];
end
hold off;
axis([400, 2000, -0.5, 4]);
legend(legendtext);
xlabel("Wavelength (nm)");
ylabel("Refractive Index");
if ttoggle == 1
    title("Extracted Refractive Index");
end
set(gca, 'fontname', 'Calibri', 'fontsize', font_size);
grid on;

colors = [0.5 0 1;
          1 0 0.5;
          0 0 0;
          0 1 0.5;
          0 0.5 1;
          1 0.5 0];
figure() % Reflectance, Transmittance, Absorptance.
set(gcf,'Color','w');
hold on
plot(1e9*c./freq, transmitted, 'Color', colors(1,:), 'LineWidth', 1.5);
plot(1e9*c./freq, reflected, 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1e9*c./freq, absorbed, 'Color', colors(3,:), 'LineWidth', 1.5);
plot(lambda, TMM_transmitted, "--", 'Color', colors(4,:), 'LineWidth', 1.5);
plot(lambda, TMM_reflected, "--", 'Color', colors(5,:), 'LineWidth', 1.5);
plot(lambda, TMM_absorbed, "--", 'Color', colors(6,:), 'LineWidth', 1.5);
hold off
legend("T_{inh}","R_{inh}","A_{inh}","T_{eff}", "R_{eff}", "A_{eff}")
axis([400 2000 -0.1 1])
xlabel("Wavelength (nm)")
if ttoggle == 1
    title("Reflectance, Transmittance, Absorptance Comparison")
end
set(gca, 'fontname', 'Calibri', 'fontsize', font_size);
grid on;

figure() % MSE from thickness optimization.
set(gcf,'Color','w');
hold on
plot(d_array*1e9, mse_total, 'Color', colors(1,:), 'LineWidth', 1.5);
plot(d_array(min_num + 1)*1e9, mse_total(min_num+1), 'diamond', 'Color', colors(3,:), 'LineWidth', 1.5);
plot(d_array(min_loc)*1e9, min_mse, 'o', 'Color', colors(2,:), 'LineWidth', 1.5);
hold off
legend("","Geometric Thickness","Effective Thickness")
xlabel("Effective Thickness (nm)")
ylabel("MSE")
if ttoggle == 1
    title("MSE for Different Effective Thicknesses")
end
text(0.55, 0.2, "Final MSE:", 'Units', 'normalized', 'fontname', 'Calibri', 'fontsize', font_size)
text(0.55, 0.1, "" + min_mse, 'Units', 'normalized', 'fontname', 'Calibri', 'fontsize', font_size)
set(gca, 'fontname', 'Calibri', 'fontsize', font_size);
grid on;

figure() % Refractive index as a function of wavelength with Kramers-Kronig Transform.
set(gcf,'Color','w');
hold on;
plot(wvl_nm, imag(n_eff), 'Color', 'r', 'LineWidth', 1.5);
plot(wvl_nm, real(n_eff), 'Color', 'b', 'LineWidth', 1.5);
plot(wvl_nm, n_KraKro, '--', 'Color', 'm', 'LineWidth', 1.5);
legendtext = ["\kappa_{eff}", "{\itn}_{eff}", "{\itn}_{eff, KK}"];
hold off;
axis([400, 2000, -0.5, 3]);
legend(legendtext);
xlabel("Wavelength (nm)");
ylabel("Refractive Index");
if ttoggle == 1
    title("Extracted Refractive Index");
end
set(gca, 'fontname', 'Calibri', 'fontsize', font_size);
grid on;


