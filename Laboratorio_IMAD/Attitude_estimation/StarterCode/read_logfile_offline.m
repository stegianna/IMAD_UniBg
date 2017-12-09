
addpath('StarterCode\')

startup()
showIP()



%% 

import('com.liu.sensordata.*');
dataFile = FileSensorDataReader('sensorLog.txt');
dataFile.start();

% Normal measurement data
dataFile.reset();
data = dataFile.getAll(5); % Typically group things separated less than half a period

% Process data
% One row per time stamp