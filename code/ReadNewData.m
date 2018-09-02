% Organize data.

doMap=false;
PlacesYear;

% XX: matrix with data from year read.

nst=length(places);

year_reg=[31 28 31 30 31 30 31 31 30 31 30 31];
year_leap=0*year_reg; year_leap(2)=1;
c_reg=[0 cumsum(year_reg)];
c_leap=[0 cumsum(year_leap)];
leap_years=[2000 2004 2008 2012 2016];
Years=[2000:2017];
ndays=365+ismember(Years,leap_years);
c_days=[0 cumsum(ndays)];

for jst=45:nst,
    
    cd('../../CRNH0203-201508171551')
    cd 2016
    file_name=strcat('CRNH0203-2016-',places{jst}{5},'.txt');
    clearvars VarName*
    uiimport(file_name);
    
    pause
    
    cd('../../Current/Boulder')
    
    m=length(VarName2);
    
    XX=zeros(m,5);
    
    D=VarName2; % Date.
    T=VarName3; % Time
    
    year=floor(D/10000);
    month=floor((D-year*10000)/100);
    day=floor((D-year*10000-month*100));
    
    ItsLeap=ismember(year,leap_years);
    
    t=(c_days(year-1999)'+c_reg(month)'+ItsLeap.*c_leap(month)'+day)*24+T/100; % Time in hours from Jan 1st 2000 at 0am.
    
    XX(:,1)=t/24;      % Time in days from Jan. 1st 2000 at 0am.
    
    XX(:,2)=VarName10; % Temperature;
    
    XX(:,3)=VarName13; % Precipitation;
    
    XX(:,4)=VarName14; % Solar radiation;
    
    XX(:,5)=VarName15; % Flag for it (0: good, 3: bad);
    
    file_name=strcat('Data',places{jst}{1},'.mat');
    
    load(file_name);
    
    t_old=max(XXX(:,1));
    
    XXX=[XXX; XX(XX(:,1)>t_old,:)];
    
    %%%%%%%%%%
    
    cd('../../CRNH0203-201508171551')
    cd 2017
    file_name=strcat('CRNH0203-2017-',places{jst}{5},'.txt');
    clearvars VarName*
    uiimport(file_name);
    
    pause
    
    cd('../../Current/Boulder')
    
    m=length(VarName2);
    
    XX=zeros(m,5);
    
    D=VarName2; % Date.
    T=VarName3; % Time
    
    year=floor(D/10000);
    month=floor((D-year*10000)/100);
    day=floor((D-year*10000-month*100));
    
    ItsLeap=ismember(year,leap_years);
    
    t=(c_days(year-1999)'+c_reg(month)'+ItsLeap.*c_leap(month)'+day)*24+T/100; % Time in hours from Jan 1st 2000 at 0am.
    
    XX(:,1)=t/24;      % Time in days from Jan. 1st 2000 at 0am.
    
    XX(:,2)=VarName10; % Temperature;
    
    XX(:,3)=VarName13; % Precipitation;
    
    XX(:,4)=VarName14; % Solar radiation;
    
    XX(:,5)=VarName15; % Flag for it (0: good, 3: bad);
    
    file_name=strcat('Data',places{jst}{1},'.mat');
    
%     load(file_name);
    
%     t_old=max(XXX(:,1));
%     
%     XXX=[XXX; XX(XX(:,1)>t_old,:)];

    XXX=[XXX; XX];
    
    save(file_name,'XXX')
    
end

% XX(:,6)=VarName21; % Surface temperature;
% 
% XX(:,7)=VarName22; % Flag for it (0: good, 3: bad);
% 
% XX(:,6)=VarName27; % Relative humidity;
% 
% XX(:,7)=VarName28; % Flag for it (0: good, 3: bad);
% 
% XX(:,8)=VarName29; % Soil moisture at 5cm below the surface;
% 
% XX(:,9)=VarName30; % Soil moisture at 10cm below the surface;
% 
% XX(:,10)=VarName31; % Soil moisture at 20cm below the surface;
% 
% XX(:,11)=VarName32; % Soil moisture at 50cm below the surface;
% 
% XX(:,12)=VarName33; % Soil moisture at 100cm below the surface;
% 
% XX(:,13)=VarName34; % Soil temperature at 5cm below the surface;
% 
% XX(:,14)=VarName35; % Soil temperature at 10cm below the surface;
% 
% XX(:,15)=VarName36; % Soil temperature at 20cm below the surface;
% 
% XX(:,16)=VarName37; % Soil temperature at 50cm below the surface;
% 
% XX(:,17)=VarName38; % Soil temperature at 100cm below the surface;

% load DataBoulder;
% 
% XXX=[XXX;XX];
% 
% save('DataBoulder.mat','XXX');

% load DataMontrose;
% 
% XXX=[XXX;XX];
% 
% save('DataMontrose.mat','XXX');

% load DataDinosaur;
% 
% XXX=[XXX;XX];
% 
% save('DataDinosaur.mat','XXX');

% load DataNunn;
% 
% XXX=[XXX;XX];
% 
% save('DataNunn.mat','XXX');

% load DataLaJunta;
% 
% XXX=[XXX;XX];
% 
% save('DataLaJunta.mat','XXX');

% load DataLander;
% 
% XXX=[XXX;XX];
% 
% save('DataLander.mat','XXX');

% load DataManhattanKs;
% 
% XXX=[XXX;XX];
% 
% save('DataManhattanKs.mat','XXX');

% load DataSocorro;
% 
% XXX=[XXX;XX];
% 
% save('DataSocorro.mat','XXX');

% load DataStillwater;
% 
% XXX=[XXX;XX];
% 
% save('DataStillwater.mat','XXX');

% load DataMonahans;
% 
% XXX=[XXX;XX];
% 
% save('DataMonahans.mat','XXX');

% load DataEdinburg;
% 
% XXX=[XXX;XX];
% 
% save('DataEdinburg.mat','XXX');

% load DataLafayette;
% 
% XXX=[XXX;XX];
% 
% save('DataLafayette.mat','XXX');










