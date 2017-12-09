% Author: Mirko Mazzoleni
% Date: 13/12/2017

%%
function [xhat, meas] = EKF_orientation(calAcc, calGyr, calMag)
% FILTERTEMPLATE  Filter template
%
% Calibration data for the accelerometer, gyroscope and magnetometer
% assumed available as structs with fields m (mean), R (variance) and b
% (biases).
%
% The function returns xhat as an array of structs comprising t
% (timestamp), x (state), and P (state covariance) for each
% timestamp, and meas an array of structs comprising t (timestamp),
% acc (accelerometer measurements), gyr (gyroscope measurements),
% mag (magnetometer measurements), and orint (orientation quaternions
% from the phone).  Measurements not available are marked with NaNs.
% Returns also the Kalman Gains K for each time stamp and Eigenvector of
% the dynamic matrix F-KH, for each time stamp.


%% Setup necessary infrastructure
import('com.liu.sensordata.*');  % Used to receive data.

%% Filter settings
% Initial time (initialize on first data received)
t0 = [];
%sampling time (s). Data are sampled on the smartphone at 100 Hz
Ts = 0.01; %10ms
% system order
nx = 4; ny=6;
% Process noise matrix
V1 = eye(nx)*1e-6;
% Output matrix
V2 = eye(ny).*repmat( [ calAcc.R; calMag.R ] ,1,ny);


% Filter state initialization
% Quaternion needs to have unit module
x = [0; 0; 0 ; 1];
% state estimation initialization x(1|0)
x_pred = [0; 0; 0 ; 1];
% State covariance matrix initialization P(1|0)
P = eye(nx, nx)*1e-5;

% Angular velocity measurements
omega = zeros(3,1);
% Acceleration measurements
acceleration = zeros(3,1);
% Magnetic field measurements
magnetic = zeros(3,1);

% Acceleration reference
g0 = [0 ; 0 ; 9.822];

mxy = sqrt( (-14.2)^2 + (-19.9)^2); mz=-33.9;
m0 = [0 ; mxy ; mz];


% Saved filter states x, covariance matrices P, innovation e
xhat = struct('t', zeros(1, 0),...
    'x', zeros(nx, 0),...
    'P', zeros(nx, nx, 0),...
    'e', zeros(ny,0));


% Saved measurents
meas = struct('t', zeros(1, 0),...
    'acc', zeros(3, 0),...
    'gyr', zeros(3, 0),...
    'mag', zeros(3, 0),...
    'orient', zeros(4, 0));
try
    %% Create data link
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

% Used for visualization.
figure(1);
subplot(1, 2, 1);
ownView = OrientationView('Own filter', gca);  % Used for visualization.
googleView = [];
counter = 0;  % Used to throttle the displayed frame rate.



%% Filter loop
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
    
    % The true acceleration signal is the measured ones minus the bias
    acc =  data(1, 2:4)' - calAcc.b;
    % Acc measurements are available.
    if ~any(isnan(acc))
        % Do not consider the contribute of external accelerations
        if(norm(acc)<=(9.822+0.5) && norm(acc)>=(9.822-0.5))
            acceleration(1:2) = acc(1:2);
            acceleration(3) = acc(3);
        end
    end
    
    gyr = data(1, 5:7)' - calGyr.b;
    % Gyro measurements are available.
    if ~any(isnan(gyr))
        omega(1:2) = gyr (1:2);
        omega(3) = gyr(3);
    end
    
    mag = data(1, 8:10)';
    % Mag measurements are available.
    if ~any(isnan(mag))
        % Do not consider disturbances on magnetic field
        if(norm(mag)<=70 && norm(mag)>=10)
            magnetic(1:2) = mag(1:2);
            magnetic(3) = mag(3);
        end
    end
    
    
    
    % ====================================================================
    % ======================= BEHOLD THE KALMAN FILTER ===================
    % ====================================================================
    
    % https://en.wikipedia.org/wiki/Extended_Kalman_filter
    % ========================== MEASUREMENTS AND UPDATE STEP
    
    % True output ==> y(k)
    y = [acceleration; magnetic];
    
    % Output predicted ==> y_hat( x(k|k-1) )
    Q = Qq( x_pred );
    y_hat = [Q'*g0; Q'*m0];
    
    % Output matrix H linearization ==> H(k) = H( x(k|k-1) )
    q0 = x_pred(1); q1 = x_pred(2); q2 = x_pred(3); q3 = x_pred(4);
    H = [ 9.822*(-2*q2)         9.822*2*q3         9.822*(-2*q0)           9.822*2*q1        ;
          9.822*2*q1            9.822*2*q0         9.822*2*q3              9.822*2*q2        ;
          9.822*4*q0            0                  0                       9.822*4*q3        ;
          mxy*2*q3+mz*(-2*q2)   mxy*2*q2+mz*2*q3   mxy*2*q1+mz*(-2*q0)     mxy*2*q0+mz*2*q1  ;
          mxy*4*q0+mz*2*q1      mz*2*q0            mxy*4*q2+mz*2*q3        mz*2*q2           ;
          mxy*(-2*q1)+mz*4*q0   mxy*(-2*q0)        mxy*2*q3                mxy*2*q2+mz*4*q3 ];
    
    
    % Kalman Gain ==> K0(k) = P(k|k-1)H(k) / H(k)P(k|k-1)H(k)' + V2
    K0 = (P*H')/(H*P*H'+V2);
    
    % Innovation ==> e(k) = y(k)-y_hat( x(k|k-1) )
    e = y-y_hat;
    % State update, filtering ==> x(k|k) = x(k|k-1) + K0(k)e(k)
    x = x_pred+K0*e;
    % Quaternion normalization
    x = quatnormalize_2(x')';
    
    % Covariance matrix update ==> P(k|k) = [I-K0(k)H(k)]*P(k|k-1)
    P = (eye(nx)-K0*H)*P;
    
    
    % ===================================  PREDICTION STEP
    
    % State prediction with non-linear system dynamics
    % x(k+1|k) = f(x(k|k),u(k))
    x_pred = (eye(4)+1/2*Somega(omega)*Ts)*x(1:4);
    % Quaternion normalization
    x_pred = quatnormalize_2(x_pred')';
    
    % Dynamic matrix F linearization
    wx = omega(1); wy = omega(2); wz = omega(3);
    % F(k) evaluated in x(k|k) and u(k)
    F = [     1       -1/2*Ts*wx  -1/2*Ts*wy   -1/2*Ts*wz  ;
          1/2*Ts*wx      1        1/2*Ts*wz    -1/2*Ts*wy  ;
          1/2*Ts*wy  -1/2*Ts*wz       1        1/2*Ts*wx   ;
          1/2*Ts*wz  1/2*Ts*wy    -1/2*Ts*wx         1     ];
    
    % Covariance prediction
    % P(k+1|k) = F(k)*P(k|k)*F(k)' + V1
    P = F * P * F' + V1;
    
    % ====================================================================
    % ============================ END KALMAN FILTER =====================
    % ====================================================================
    

    
    % Google's orientation estimate
    orientation = data(1, 18:21)';
    
    % Visualize result
    if rem(counter, 10) == 0
        setOrientation(ownView, x(1:4));
        title(ownView, 'OWN', 'FontSize', 16);
        if ~any(isnan(orientation))
            if isempty(googleView)
                subplot(1, 2, 2);
                % Used for visualization.
                googleView = OrientationView('Google filter', gca);
                
            end
            setOrientation(googleView, orientation);
            title(googleView, 'GOOGLE', 'FontSize', 16);
        end
    end
    counter = counter + 1;
    
    % Save estimates
    xhat.x(:, end+1) = x;
    xhat.P(:, :, end+1) = P;
    xhat.t(end+1) = t - t0;
    xhat.e(:,end+1) = e;
    
    meas.t(end+1) = t - t0;
    meas.acc(:, end+1) = acc;
    meas.gyr(:, end+1) = gyr;
    meas.mag(:, end+1) = mag;
    meas.orient(:, end+1) = orientation;
    
end
end
