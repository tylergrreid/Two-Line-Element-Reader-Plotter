%% TWO-LINE-ELEMENT READER / PLOTTER
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
% simple orbit plotting / visulization tool. 
% 
%% WORKSPACE SETUP

clc
clear all
close all
format longG

% Add path to the Earth plotting function. 
addpath([pwd, '/PlotEarth/PlotEarth']);

% Add path to the TLE data (these are just some examples with a focus on
% GNSS). 
addpath([pwd, '/TLE_Files']);

%% LOAD PHYSICAL CONSTANTS INTO THE GLOBAL WORKSPACE

physical_constants

global mu

%% SELECT THE TLE FILES TO PLOT

% All of GNSS. 
filenames = {'gps-ops-clean','galileo-full','sbas-all-clean',...
    'beidou-clean-complete','glo-ops-clean'};

% GPS + SBAS.
% filenames = {'gps-ops-clean','sbas-clean','qzss-full'};

% GLONASS. 
% filenames = {'glo-ops-clean'};

% Completed Galileo constellation. 
% filenames = {'galileo-full'};

% BeiDou
% filenames = {'beidou-05-25-2017'};

% GPS + Iridium
% filenames = {'gps-ops-clean','iridium'};

% GPS + SBAS GEOS. 
% filenames = {'gps-ops-clean','sbas-reduced-clean'};

% Transit. 
% filenames = {'nnss'};

% QZSS. 
% filenames = {'qzss-full'};

% Planet labs.
% filenames = {'rapideye', 'skybox', 'planet_mc'};

% Create colors for plotting. 
colors = lines(length(filenames));
% colors = hsv(length(filenames));

%% DECODE TLE DATA

% Plot the Earth. 
% If you want a color Earth, use 'neompa', 'BlueMarble'.
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'.
% A smaller sample step gives a finer resolution Earth.
h = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 2);

% Choose a date, by default this chooses the current date/time.
simStart = datenum(clock);
% simStart = datenum('Jan 09 2017 00:00:00');

% Compute sidereal time. 
GMST = utc2gmst(datevec(simStart)); % [rad]

% Import the TLE data. 
for k = 1:length(filenames)
    % Get the orbital elements. 
    [coe] = two_line_elem_conv(horzcat(filenames{k}, '.txt'), 'all');
    
    % Find latest epoch (all others can be run up from here).
    coeDateNums = datenum(coe.date);
    [val, ind] = min(coeDateNums);
    
    % Define the max time from simStart.
    a = coe.a(1);
    n = sqrt(mu/a^3);
    tFinal = 2*pi/n*1.1/3600/24; % This gives us just over 1 orbit. 
    
    % Create a time vector. 
    tSim = linspace(simStart, simStart + tFinal, 200);
    
    % Allocate space.
    RSave = NaN(length(tSim), 3, length(coeDateNums));
    
    % Run through all of the satellites in the TLE file 
    % and compute trajectories for plotting.
    for i = 1:length(coeDateNums) % For each satellite
        for j = 1:length(tSim) % For each time step
            % Get the orbit data. 
            a = coe.a(i);
            n = sqrt(mu / a^3);
            e = coe.e(i);
            inc = coe.i(i) * pi / 180;
            RAAN = coe.RAAN(i) * pi / 180;
            omega = coe.omega(i) * pi / 180;
            M = coe.M(i) * pi / 180 + ...
                n * (tSim(j) - coeDateNums(i)) * 24 * 3600; 
            
            % Adjust RAAN such that we are consisten with Earth's current
            % orientation. This is a conversion to Longitude of the
            % Ascending Node (LAN). 
            RAAN = RAAN - GMST;
            
            % Convert to ECI and save the data.
            [X,~] = COE2RV(a, e, inc, RAAN, omega, M);
            RSave(j,:,i) = X';
        end
    end
    
    % Plot the orbit.
    for i = 1:length(coeDateNums)
        colorI = k;
        plot3(RSave(:,1,i) / R_e, RSave(:,2,i) / R_e, RSave(:,3,i) / R_e,...
            'color', colors(colorI,:), 'LineWidth', 1)
        plot3(RSave(1,1,i) / R_e, RSave(1,2,i) / R_e, RSave(1,3,i) / R_e,...
            '.', 'color', colors(colorI,:), 'MarkerSize', 10)
        hold on
    end
end

%% SAVE FIGURE

% If you want a black background set to 'off', otherwise set to 'on' or
% just comment this out.
set(gcf, 'InvertHardCopy', 'on');

% Set the view angle of the figure. 
view([-20, 9])

% Reset the zoom. 
% zoom reset
zoom(1.25)

% Turn off axis clipping. 
ax = gca;               
ax.Clipping = 'off';    

% Export the figure. 
exportfig(gcf,horzcat(filenames{1},'Multi.tiff'),'height',6,'width',9,'fontsize',16,'LineWidth',10,'resolution',220);
