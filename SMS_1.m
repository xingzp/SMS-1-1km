%% Function used to retrieve Sentinel-1 SM; 

% Refer to paper "Retrieval of 1 km surface soil moisture from Sentinel-1 over bare soil and grassland on the Qinghai-Tibetan Plateau" for the information of equations
% A,C and D parameters can be obtained via: https://drive.google.com/drive/folders/1LHI_c_7ji5Ezvsd_GjgypQmLB4S-l2go?usp=drive_link 
% NDVI and SIG data can be downloaded from Google Earth Engine platform,
% and Please refer to my paper for more detailed processing steps.

function SM = SMS_1(NDVI,C,D,A,SIG)
    % theta is the normalized local incidence angle
    % in my study, SM was retrieved track by track following the Sentinal-1
    % orbits.
    theta = cos(deg2rad(40));
    B = 0.084;
    SIG_linear = 10.^(SIG/10);
    VWC = 0.5*NDVI.*NDVI-0.5*NDVI+0.15*((nanmax(NDVI)-nanmin(NDVI))./(1-nanmin(NDVI)));
    gama2 = exp(-2*VWC*B/theta);
    gama2(gama2>1|gama2<0) = nan;

    SIG_veg_linear = A.*theta.*(1-gama2);
    SIG_soil_linear1 = SIG_linear-SIG_veg_linear;
    SIG_soil_linear1(SIG_soil_linear1<=0)=nan;
    SIG_soil_linear2 = SIG_soil_linear1./gama2;
    SIG_soil = 10.*log10(SIG_soil_linear2);
    SM=(SIG_soil-C)./D;
    SM(SM>=1|SM<=0)=nan;
end