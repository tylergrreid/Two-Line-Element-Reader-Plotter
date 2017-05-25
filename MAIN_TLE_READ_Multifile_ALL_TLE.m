%% PLOT ALL OPERATIONAL ORBITS
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       Lab:                  Stanford GPS Lab
%       Project Start Date:   Feb 18, 2014
%       Last updated:         May 24, 2017
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% Given multiple NORAD Two-Line-Element (TLE) files, this matlab code plots
% the orbits of the satellites as well as the Earth. This is meant to be a
% simple orbit plotting / visulization tool. This particular version plots
% all the NORAD TLEs of the operational satellites. 
%
%% WORKSPACE SETUP

clc
clear
close all
format longG

%% LOAD PHYSICAL CONSTANTS INTO WORKSPACE

physical_constants

global mu

%% DECODE TLE DATA

% Plot the Earth. 
% If you want a color Earth, use 'neompa', 'BlueMarble'.
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'.
% A smaller sample step gives a finer resolution Earth.
h = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 1);

% Add path to data folder.
data_main_folder = [pwd, '/TLE_Files_All/'];

% Names of folders.
folders = {'Comm','Weather_Earth_Resources','Scientific','Nav','SBAS'};

% Make colors.
colors_default = lines(10);
colors(1,:) = colors_default(6,:); % Comm
colors(4,:) = colors_default(2,:); % Nav
colors(5,:) = colors_default(3,:); % SBAS
colors(3,:) = colors_default(4,:) * 1.7; % Science
colors(2,:) = colors_default(5,:); % Earth Observation

% Make dummy plot for legend purposes.
for i = 1:length(colors)
   h2(i) = plot([0.1,0.2],[0.1,0.2],'color',colors(i,:),'LineWidth',2);
   hold on
end

% Choose a representation date.
simStart = datenum('Mar 15 2016 00:00:00');

% Compute sidereal time.
GMST = utc2gmst(datevec(simStart));

% Go through all categories.
for cat = 1:length(folders)
    % Get data from all the folders.
    files_in_folder = dir([data_main_folder, folders{cat}]);
    files_in_folder.name;
    
    % Clear variables.
    clear filenames
    
    % Start count variable.
    file_count = 1;
    
    for j = 1:length(files_in_folder)
        filename_test = [data_main_folder,folders{cat},'/',...
            files_in_folder(j).name];
        if exist(filename_test,'file')==2 && ...
                strcmp(files_in_folder(j).name,'.DS_Store')==0
            filenames{file_count} = filename_test;
            file_count = file_count + 1;
        end
    end
    
    % Import the TLE data. 
    for k = 1:length(filenames)
        % Get the orbital elements.
        [coe] = two_line_elem_conv(horzcat(filenames{k}),'all');
        
        % Find latest epoch (all others can be run up from here).
        coeDateNums = datenum(coe.date);
        [val, ind] = min(coeDateNums);
        
        % Run through them all and compute trajectories for plotting.
        for i = 1:length(coeDateNums) % for each satellite
            % Define the max time from simStart.
            a = coe.a(i);
            n = sqrt(mu/a^3);
            tFinal = 2*pi/n*1.1/3600/24; % This gives us just over 1 orbit. 
            
            % Create a time vector. 
            tSim = linspace(simStart, simStart+tFinal,300);
            
            % Allocate space. 
            RSave = NaN(length(tSim), 3);
            
            for j = 1:length(tSim) % for each time step
                % Get the orbit data. 
                a = coe.a(i);
                n = sqrt(mu/a^3);
                e = coe.e(i);
                inc = coe.i(i)*pi/180;
                RAAN = coe.RAAN(i)*pi/180;
                omega = coe.omega(i)*pi/180;
                M = coe.M(i)*pi/180 + n*(tSim(j)-coeDateNums(i))*24*3600;
                
                % Adjust RAAN such that we are consisten with Earth's current
                % orientation. This is a conversion to Longitude of the
                % Ascending Node (LAN).
                RAAN = RAAN - GMST;
                
                % Convert to ECI and save the data.
                [X,~] = COE2RV(a,e,inc,RAAN,omega,M);
                RSave(j,:) = X';
            end
            
            % Plot the orbit.
            plot3(RSave(:,1)/R_e,RSave(:,2)/R_e,RSave(:,3)/R_e,'color',[colors(cat,:),1],'LineWidth',0.75)
            % plot3(RSave(:,1,i)/R_e,RSave(:,2,i)/R_e,RSave(:,3,i)/R_e,plot_style{colorI},'color',colors(colorI,:),'LineWidth',1.5)
            % plot3(RSave(1,1,i)/R_e,RSave(1,2,i)/R_e,RSave(1,3,i)/R_e,'.','color',colors(colorI,:),'MarkerSize',15)
            hold on
        end
    end
end

%% FORMAT FIGURE

% If you want a black background set to 'off', otherwise set to 'on' or
% just comment this out.
% set(gcf, 'InvertHardCopy', 'off'); % black background
set(gcf, 'InvertHardCopy', 'on'); % white background

% Set the view angle of the figure. 
view([-20, 9])

% Set zoom. 
% zoom reset
% zoom(10)
% xlim([-2000,2000])
% ylim([-2000,2000])
% zlim([-2000,2000])
% zoom(200) % Regular view.
% zoom(500) % Close up for cover photo.

% Turn clipping off. 
ax = gca;               
ax.Clipping = 'off';   

% Set legend.
% h3 = legend(h2,{'Communication','Earth Observation','Science','Navigation','SBAS'},'TextColor',[1,1,1],'location','southwest');
% legend('boxoff')

% Save the figure.
% To get a custom view, zoom and rotate the figure manually and re-run this
% line of code. 
exportfig(gcf,'TLE_ALL.jpeg','height',10,'width',14,'fontsize',24,'LineWidth',5,'resolution',300);
