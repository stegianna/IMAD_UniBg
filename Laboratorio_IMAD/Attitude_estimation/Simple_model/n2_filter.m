% Author: Mirko Mazzoleni
% Date: 13/12/2017

%% Connect
clear
addpath('..\StarterCode\')

startup()
showIP()


%% Kalman Filter

% Load sensors biases and variances
load calibration_data

% create calibration struct for each sensor
calAcc.m = acc_means;
calAcc.R = acc_vars;
calAcc.b = acc_biases;

calGyr.m = gyr_means;
calGyr.R = gyr_vars;
calGyr.b = gyr_biases;

calMag.R = mag_vars;


% Call the filtering procedure
% In order to visualize the Google estimation, enable to send the
% "orientation" variable
[xhat, meas] = EKF_orientation(calAcc, calGyr, calMag);



% Quaternion estimated vs Google estimation
figure
subplot 411
plot(xhat.t, abs(xhat.x(1,:)), 'g','linewidth', 2);  hold on;
plot(meas.t, abs(meas.orient(1,:)), 'k--','linewidth', 2);
subplot 412
plot(xhat.t, abs(xhat.x(2,:)), 'g', 'linewidth', 2);  hold on;
plot(meas.t, abs(meas.orient(2,:)), 'k--', 'linewidth', 2);
subplot 413
plot(xhat.t, abs(xhat.x(3,:)), 'g', 'linewidth', 2);  hold on;
plot(meas.t, abs(meas.orient(3,:)), 'k--', 'linewidth', 2);
subplot 414
plot(xhat.t, abs(xhat.x(4,:)), 'g', 'linewidth', 2);  hold on;
plot(meas.t, abs(meas.orient(4,:)), 'k--', 'linewidth', 2);
legend('Estimate','True');


