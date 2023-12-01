% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created October, 2023

function [Corrected] = CorrectionNonparalyzable (Uncorrected, Tau,TotalPulses, BinWidth)
%
% Inputs: Uncorrected: an array containing photon counts which have not
%                      been corrected
%         Tau:         the nonparalyzable dead time of a system       [sec]
%         TotalPulses: the total number of pulses including system
%                      integration time and extra integration 
%         BinWidth:    the altitude range bin width                     [m]
%         c:           the speed of light                             [m/s]
% Outputs: Corrected:  an array containing photon counts which have been
%                      corrected 
%
%% Constants
c = 299792458;
%% Correcting Saturation
Rate          = FindRates (Uncorrected, TotalPulses, BinWidth, c);
CorrectedRate = Rate  ./(1 - (Tau.*Rate ));
Corrected     = FindCounts (CorrectedRate, TotalPulses, BinWidth, c);
end

function [Rate] = FindRates (Data, TotalPulses, BinWidth, c)
%
% Inputs: Data:            the photon number data per time bin
%         TotalPulses:     the total number of laser pulses in the
%                          integration times
%         BinWidth:        the range resolving bin width                [m]
%         c:               Light speed                                [m/s]
%
% Outputs: Rate:           the count rate per range bin                [Hz]
% 
%%
tau    = (2*BinWidth)/c;                    % time to travel the bin width 
Rate = Data ./(tau.*TotalPulses);% Photon rate in Hz
end

function [Count] = FindCounts (Data, TotalPulses, BinWidth, c)
%
% Inputs: Data:            the photon number data per time bin
%         TotalPulses:     the total number of laser pulses in the
%                          integration times
%         BinWidth:        the range resolving bin width                [m]
%         c:               Light speed                                [m/s]
%
% Outputs: Rate:           the count rate per range bin                [Hz]
% 
%%
tau     = (2*BinWidth)/c;               % time to travel the bin width 
Count   = Data .*(tau.*TotalPulses);     % Photon rate in Hz
end