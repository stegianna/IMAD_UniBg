% Author: Mirko Mazzoleni
% Date: 13/12/2017

%% Connect

clear 
% add utilities
addpath('..\StarterCode\')

startup()
showIP()

% The measured data will be that of a still phone which lies horizontally
% on the table

%% start server

import('com.liu.sensordata.*');  % Used to receive data.
try
  %% Create data link
  %3400 port
  server = StreamSensorDataReader(3400);
  % Makes sure to resources are returned.
  sentinel = onCleanup(@() server.stop());

  server.start();  % Start data reception.
catch e
  fprintf(['Unsuccessful connecting to client!\n' ...
    'Make sure to start streaming from the phone *after* '...
           'running this function!']);
  return;
end

%define measurements to be stored
meas = struct('t', zeros(1, 0),...
    'acc', zeros(3, 0),...
    'gyr', zeros(3, 0),...
    'mag', zeros(3, 0),...
    'orient', zeros(4, 0));

% Initial time
t0 = [];

while server.status()  % Repeat while data is available
  % Get the next measurement set, assume all measurements
  % within the next 5 ms are concurrent (suitable for sampling
  % in 100Hz).
  data = server.getNext(5);

  if isnan(data(1))  % No new data received
    continue;
  end
  
  t = data(1)/1000;  % Extract current time
  if isempty(t0)  % Initialize t0
      t0 = t;
  end

  % Process data
  % data(1)         Time [ms]
  % data(1, 2:4)    Accelerometer; x, y, z [m/s^2]
  % data(1, 5:7)    Gyroscope; x, y, z [rad/s]
  % data(1, 8:10)   Magnetometer; x, y, z [uT]
  % data(1, 11:13)  GPS [deg North, deg East, alt m]
  % data(1, 14)     Pressure [hPa]
  % data(1, 15)     Light [lux]
  % data(1, 16)     Proximity [cm]
  % data(1, 17)     Temperature [deg C]
  % data(1, 18:21)  Orientation [normalized quaternion]
  % data(1, 22:27)  Uncalibrated gyroscope; x, y, z, bx, by, bz [rad/s]
  % data(1, 28:33)  Uncalibrated magnetometer; x, y, z, bx, by, bz [uT]
  %
  % Nan values indicate missing data.
  
  % store measurements
  meas.t(end+1) = t - t0;
  meas.acc(:, end+1) = data(1, 2:4)';
  meas.gyr(:, end+1) = data(1, 5:7)';
  meas.mag(:, end+1) = data(1, 8:10)';
  
end

% save acquired data
save data_still meas


%% load and plot saved data
clear
% load saved data to analyze
load('data_still.mat')


% ========================
% ACCELEROMETER
% ========================
figure
% It should be 0 m/s^2 if the surface is perfectly plane
subplot 311; plot(meas.t, meas.acc(1,:), 'b', 'linewidth', 2); grid on; xlim([0, max(meas.t)])
% It should be 0 m/s^2 if the surface is perfectly plane
subplot 312; plot(meas.t, meas.acc(2,:), 'b', 'linewidth', 2); ylabel('Acceleration [m/s^2]'); grid on; xlim([0, max(meas.t)])
% It should be 9.822 m/s^2 if the surface is perfectly plane
subplot 313; plot(meas.t, meas.acc(3,:), 'b', 'linewidth', 2); 
xlabel('Time [s]'); grid on; xlim([0, max(meas.t)])

% From the plot it can be seen that there is an initial transient (probably
% due to the finger which presses the button to start streaming, or an 
% internal filter output): Let's take it out starting from 3 seconds up to 
% 36 seconds (do not keep the final piece in which pressing the stop button
% of the app induces small accelerations)
acc_x = meas.acc(1, find(meas.t>3, 1):find(meas.t>36, 1) );
acc_y = meas.acc(2, find(meas.t>3, 1):find(meas.t>36, 1) );
acc_z = meas.acc(3, find(meas.t>3, 1):find(meas.t>36, 1) );

% Create a vector of cleaned data
meas.acc_cleaned = [acc_x; acc_y; acc_z]; clearvars acc_x acc_y acc_z

% % Example to understand how to exclude the NaN numbers (columns: time
% % stamps, rows: sensor measurements)
% aaa=[NaN 1 NaN 4 6; 0 3 NaN 2 4; NaN 4 NaN 0 NaN  ]
% isnan(aaa)
% % If 1, one of the 3 axis gave a NaN at time instant t
% ~any(isnan(aaa),1)
% clearvars aaa

% Compute the mean acceleration discarding NaN values
acc_means = mean(meas.acc_cleaned(:, ~any(isnan(meas.acc_cleaned), 1)), 2);
% The difference between the mean value estimate and the true value is the
% sensor bias (supposing measurements taken with the phone on a plane
% table)
acc_biases = acc_means - [0; 0; 9.81];

% From the histogram it is possible to have an idea of the noise which
% affects the measurements
figure
subplot 311; hist(meas.acc_cleaned(1,:)); 
subplot 312; hist(meas.acc_cleaned(2,:)); 
subplot 313; hist(meas.acc_cleaned(3,:)); 

% Variance of the measurements
acc_vars=[ var(meas.acc_cleaned(1, ~any(isnan(meas.acc_cleaned), 1) ));
           var(meas.acc_cleaned(2, ~any(isnan(meas.acc_cleaned), 1) ));
           var(meas.acc_cleaned(3, ~any(isnan(meas.acc_cleaned), 1) ))];

       

% ========================
% GYROSCOPE
% ========================       

figure
% It should be 0 rad/s if the phone is still
subplot 311; plot(meas.t, meas.gyr(1,:), 'b', 'linewidth', 2); grid on; xlim([0, max(meas.t)])
% It should be 0 rad/s if the phone is still
subplot 312; plot(meas.t, meas.gyr(2,:), 'b', 'linewidth', 2); ylabel('Speed [rad/s]'); grid on; xlim([0, max(meas.t)])
% It should be 0 rad/s if the phone is still
subplot 313; plot(meas.t, meas.gyr(3,:), 'b', 'linewidth', 2); 
xlabel('Time [s]'); grid on; xlim([0, max(meas.t)])

gyr_raw_x = meas.gyr(1, find(meas.t>3, 1):find(meas.t>36, 1) );
gyr_raw_y = meas.gyr(2, find(meas.t>3, 1):find(meas.t>36, 1) );
gyr_raw_z = meas.gyr(3, find(meas.t>3, 1):find(meas.t>36, 1) );
meas.gyr_cleaned = [gyr_raw_x; gyr_raw_y; gyr_raw_z]; clearvars gyr_raw_x gyr_raw_y gyr_raw_z

gyr_means = mean(meas.gyr_cleaned(1:3, ~any(isnan(meas.gyr_cleaned), 1)), 2);
gyr_biases = gyr_means-[0; 0; 0];


figure
subplot 311; hist(meas.gyr_cleaned(1, :)); 
subplot 312; hist(meas.gyr_cleaned(2, :)); 
subplot 313; hist(meas.gyr_cleaned(3, :)); 

gyr_vars=[ var(meas.gyr_cleaned(1, ~any(isnan(meas.gyr_cleaned), 1) ));
           var(meas.gyr_cleaned(2, ~any(isnan(meas.gyr_cleaned), 1) ));
           var(meas.gyr_cleaned(3, ~any(isnan(meas.gyr_cleaned), 1) ))];




% ========================
% MAGNETOMETER 
% ========================    

% It is assumed that the magnetometer is already calibrated

figure
subplot 311; plot(meas.t,meas.mag(1, :), 'b', 'linewidth', 2); ; grid on; xlim([0, max(meas.t)])
subplot 312; plot(meas.t,meas.mag(2, :), 'b', 'linewidth', 2); ylabel('Magnetic field [ \muT ]'); ; grid on; xlim([0, max(meas.t)])
subplot 313; plot(meas.t,meas.mag(3, :), 'b', 'linewidth', 2); grid on; xlim([0, max(meas.t)])
xlabel('Time [s]');


figure
subplot 311; hist(meas.mag(1, :)); 
subplot 312; hist(meas.mag(2, :)); 
subplot 313; hist(meas.mag(3, :)); 

mag_vars=[ var(meas.mag(1, ~any(isnan(meas.mag), 1) ));
           var(meas.mag(2, ~any(isnan(meas.mag), 1) ));
           var(meas.mag(3, ~any(isnan(meas.mag), 1) ))];


% save biases and variances of the sensors
save calibration_data


