
%% Disclaimer
% This file is codes used in the paper "F.Fakhraei Roudasari, C. M. J.
% Tampère, Exploring sensitivity-based clustering of OD variables in
% dynamic demand calibration,  IT-MTS2017 Conference" developed by the
% KULeuven.
%
% Copyright (C) 2017  Farzad Fakhraei Roudsari, Leuven, Belgium
% contact: farzad.fr {@} kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% More information regarding dynamic traffic assignments at:
% http://www.mech.kuleuven.be/en/cib/traffic/downloads or contact:
% willem.himpe {@} kuleuven.be

%% Introduction
% The code here is organized into two parts corresponding to the sections
% of the paper. In Part III.A we considered all origins to one
% destination and we show by excluding some part of variables in the
% optimization direction, we can increase the calibration accuracy. in Part
% III.B we use visual inspection of approximated hessian matrix to define
% two separate clusters and we prove if we include all OD pairs in the
% calibration process, we stock in the local optima, whereas using
% the clustering approach by separating them and calibrate each cluster
% separately, the accuracy of the calibration will be maximized
%

% clear the command window
clear
clc
close all
close all hidden
clearvars -global
rng(0,'twister');

%add these folders to the search path;
addpath('Dynamic Traffic Assignment','Visualization Tools','Network Data')
javaclasspath('Dynamic Traffic Assignment');

% Network and demand data
load net5_old.mat                                       %Load network data


det =[1 4 7 8 9 10 26 27 28 33 34 35 36 38 39 40 41 44 45 46 53];                       %location of the detectors
dt = 5/60;                                              %simulation interval
totT = round(3/dt);                                     %simulation time horizon

% define the section and experiment
paper_experiment_part = 'B';                            % The section corresponds in the paper {A, B}

%----% for Part III.A %-----%
DATA_INPUT.partA.scenario = 3;
% Scenario 1: All OD pairs are calibrated simultaneously
% Scenario 2: First O5-D6 calibrated then all OD pairs
% Scenario 3: First all OD pairs except O5-D6 then all OD pairs including O5-D6.

%----% for Part III.B %-----%
DATA_INPUT.partB.scenario = 2;
% Scenario 1:   30% * Randomize (True Destination 6 value)  & 70% * Randomize (True Destination 16 value)
% Scenario 2:   70% * Randomize (True Destination 6 value)  & 30% * Randomize (True Destination 16 value)

DATA_INPUT.partB.method = 2;
%choose between calibrating ALL OD pairs at once or two seperate clusters
%possible input 1: All_at_once', 2: clusters'
%% switch between experiment results of part A and B of the paper
switch paper_experiment_part %% Scenario definition corresponding to section III.A
    case 'A'
        [DATA_OUTPUT,DATA_INPUT,DATA_STRUCT]= part_A(dt,totT,ODmatrices,timeSeries,...
            links,nodes,det,DATA_INPUT);
        
    case 'B'  %% Scenario definition corresponding to section III.B
        [DATA_OUTPUT,DATA_INPUT,DATA_STRUCT]= part_B(dt,totT,ODmatrices,timeSeries,...
            links,nodes,det,DATA_INPUT);
end

