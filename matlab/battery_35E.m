function bat = battery_35E(pack, SOC0, Qscale)
% =========================================================================
% BATTERY_35E  Battery model: Samsung INR18650-35E cells, 13s6p ("48 V") or 20s6p ("72 V").
%   OCV(SOC) + internal resistance. OCV from the 0.2 A discharge curve (lygte-info.dk,
%   "Discharge, capacity: Samsung INR18650-35E"), plus 0.2 A x 35 mOhm.
%   R_cell = 45 mOhm (between the 35 mOhm AC-IR of the datasheet and the ~57 mOhm effective
%   value derived from the 7 A discharge curve).
%     bat = battery_35E(48, 0.6)        % 48 V pack, initial SOC 60 %
%     bat = battery_35E(72, 0.95, 1e-3) % capacity reduced 1000 times (CV scenario)
% In Simulink: Controlled Voltage Source (Vocv) in series with the resistance bat.R;
%   SOC  = SOC0 + (1/(3600*bat.Q)) * integral(Ibat)   (Ibat > 0 = charging)
%   Vocv = bat.Ns * interp1(bat.SOCtab, bat.OCVtab, SOC)  (1-D Lookup Table block)
% =========================================================================
if nargin < 3, Qscale = 1; end
bat.SOCtab = [0.0000 0.0413 0.0855 0.1150 0.1593 0.1888 0.2625 0.3363 0.4100 0.4838 0.5575 0.6313 0.7050 0.7788 0.8525 0.9263 1.0000];
bat.OCVtab = [2.807 3.007 3.207 3.307 3.427 3.447 3.517 3.587 3.657 3.727 3.807 3.867 3.917 3.997 4.067 4.107 4.167];   % V/cell
bat.Np = 6;  bat.Qcell = 3.39;  bat.Rcell = 0.045;
switch pack
    case 48, bat.Ns = 13;  bat.Vcv = 54.6;  bat.Vcut = 13*2.65;
    case 72, bat.Ns = 20;  bat.Vcv = 84.0;  bat.Vcut = 20*2.65;
    otherwise, error('pack must be 48 or 72');
end
bat.Q    = bat.Np*bat.Qcell*Qscale;          % Ah (reduced for the CV / rising-voltage scenarios)
bat.R    = bat.Ns*bat.Rcell/bat.Np;          % Ohm (13s6p: 0.0975; 20s6p: 0.15)
bat.SOC0 = SOC0;
bat.V0   = bat.Ns*interp1(bat.SOCtab, bat.OCVtab, SOC0);   % initial OCV (initial voltage of Cin/Cd if needed)
end
