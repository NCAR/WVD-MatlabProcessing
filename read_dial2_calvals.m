cd JSon/

dat=loadjson(['dial',Options.System(6),'_calvals.json'],'SimplifyCell',1);
cd ..
%t_date = '11-Jun-2017'
t_date = datetime(num2str(Date),'InputFormat','yyMMdd');

for i=1:size(dat.Default_P,2)
    if (t_date >= datetime(dat.Default_P(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
        P0 = dat.Default_P(i).value;
    end
end

for i=1:size(dat.Molecular_Gain_Matlab,2)
    if (t_date >= datetime(dat.Molecular_Gain_Matlab(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
        receiver_scale_factor = dat.Molecular_Gain_Matlab(i).value;
        diff_geo_on = dat.Molecular_Gain_Matlab(i).diff_geo;
    end
end

for i=1:size(dat.switch_ratio,2)
    if (t_date >= datetime(dat.switch_ratio(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
        switch_ratio = dat.switch_ratio(i).value;
    end
end

for i=1:size(dat.Location,2)
    if (t_date >= datetime(dat.Location(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
        location = dat.Location(i).location;
    end
end
location; %write the location to the screen

% % % load('diff_geo_cor_170810.mat');
timing_range_correction = (1.25-0.2+0.5/2-1.0/2)*150;  % changed hardware timing to start after pulse through
blank_range = 450; % new pulse generator shifts gate timing so less outgoing pulse contamination

clear i

JSondeData.BasePressure        = P0;
JSondeData.BlankRange          = blank_range;
JSondeData.Data                = dat;
JSondeData.Date                = t_date;
JSondeData.DeadTime            = 37.25E-9;
% JSondeData.GeoCorrectionOn     = diff_geo_on;
% JSondeData.GeoCorrectionOff    = diff_geo_corr;
JSondeData.Location            = location;
JSondeData.RangeCorrection     = timing_range_correction;
% % % % JSondeData.ReceiverScaleFactor = receiver_scale_factor;
JSondeData.SwitchRatio         = switch_ratio;

clear blank_range diff_geo_corr diff_geo_on location MCS profiles2ave
clear switch_ratio time_per_column receiver_scale_factor dat
clear P0 timing_range_correction t_date