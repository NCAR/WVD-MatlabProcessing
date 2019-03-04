


function [Capabilities,Iterations,HardwareMap,Map] = ReadFakeJSondFile(Paths)


%% Defining type map  (Info that should end up in jsond files eventually) 
if strcmp(Paths.FigureType,'DIAL01') || strcmp(Paths.FigureType,'DIAL03') ||  ...
   strcmp(Paths.FigureType,'DIAL04') || strcmp(Paths.FigureType,'DIAL02')
    % Capabilities
    Capabilities.HSRL   = 0;
    Capabilities.O2DIAL = 0;
    Capabilities.O2HSRL = 0;
    Capabilities.WVDIAL = 1;
    % Number of iterations to perform on the retrievals 
    Iterations = 1;
    % Map in the stored cell arrays
    Map.Channels  = {'Online';'Offline'};
    Map.Offline   = 2;
    Map.Online    = 1;
    % Map to get to hardware location
    HardwareMap.ChannelName    = {'WVOffline';'WVOnline'};
    HardwareMap.Etalon         = [0;0];
    HardwareMap.Laser          = [0;1];
    HardwareMap.PhotonCounting = [0;8]; % Still need from file
    HardwareMap.Power          = [0;6]; % Still need from file
elseif strcmp(Paths.FigureType,'DIAL02')
%     % Capabilities
%     Capabilities.HSRL   = 0;
%     Capabilities.O2DIAL = 0;
%     Capabilities.O2HSRL = 0;
%     Capabilities.WVDIAL = 1;
%     % Number of iterations to perform on the retrievals 
%     Iterations = 1;
%     % Map in the stored cell arrays
%     Map.Channels  = {'Online';'Offline';'Molecular';'Combined'};
%     Map.Combined  = 4;
%     Map.Offline   = 2;
%     Map.Online    = 1;
%     Map.Molecular = 3;
%     % Map to get to hardware location
%     HardwareMap.ChannelName    = {'WVOnline';'WVOffline';'HSRLMolecular';'HSRLCombined'};
%     HardwareMap.PhotonCounting = [0;8;2;3];
%     HardwareMap.Power          = [0;6;1;1];
%     HardwareMap.Etalon         = [0;0;1;1];
%     HardwareMap.Laser          = [0;1;2;2];
elseif strcmp(Paths.FigureType,'DIAL05')
    % Capabilities
    Capabilities.HSRL   = 0;
    Capabilities.O2DIAL = 1;
    Capabilities.O2HSRL = 1;
    Capabilities.WVDIAL = 1;
    % Number of iterations to perform on the retrievals 
    Iterations = 3;
    % Map in the stored cell arrays
    Map.Channels  = {'Online';'Offline';'Molecular';'O2OnlineMol';'Combined';'O2OnlineComb'};
    Map.Combined     = 5;
    Map.Offline      = 2;
    Map.Online       = 1;
    Map.Molecular    = 3;
    Map.O2OnlineMol  = 4;
    Map.O2OnlineComb = 6;
    % Map to get to hardware location
    HardwareMap.ChannelName    = {'WVOnline';'WVOffline';'O2OfflineMol';'O2OnlineMol';'O2OfflineComb';'O2OnlineComb'};
    HardwareMap.PhotonCounting = [0;8;1;9;2;10];
    HardwareMap.Power          = [0;6;1;7;1;7];
    HardwareMap.Etalon         = [0;0;0;0;0;0];
    HardwareMap.Laser          = [0;1;3;4;3;4];
end
end