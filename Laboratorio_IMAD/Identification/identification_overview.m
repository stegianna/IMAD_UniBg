% Author: Mirko Mazzoleni
% Date: 13/12/2017


%% Data simulation
clear
clc

rng('default');

me = 10; % white noise mean
a = [0.7 -0.2 -0.3]; % AR coefficients
N = 11000; % num data

e = randn(1, N) + me; % generate WN input
y = zeros(1, length(e)); % generate output

Ts = 0.02; % sample time [s]

for t = 4 : length(e)
    
    y(t) = a(1)*y(t-1) + a(2)*y(t-2) + a(3)*y(t-3) + e(t);
    
end
% Transfer function of the process
G = tf([1 0 0 0], [1 -a(1) -a(2) -a(3)] ,'Ts', Ts, 'variable', 'z^-1');
G

% Abs value of the poles of the transfer function
abs( pole(G) )

% Add stagionality to the process
s = 3 * cos(2 * pi/40*[1:N]);
y = y + s;

figure
plot(y, 'b', 'linewidth', 2); hold on
xlim([0, N]); ylim([-1, 21]); grid on; legend('Process', 'Location', 'best');
xlabel('Samples [-]'); ylabel('Values [-]'); grid on;

figure
plot(y, 'b', 'linewidth', 2); hold on
plot(s + mean(y), 'k--', 'linewidth', 4);
xlim([0, N]); ylim([-1, 21]); grid on; legend('Process', 'Stagionality', 'Location', 'best');
xlabel('Samples [-]'); ylabel('Values [-]')
xlim([100, 500]); grid on;

%% Pre-processing of the time series

% Remove initial condition
initial_condition = 1000;
y = y(initial_condition+1:end);
e = e(initial_condition+1:end);
N = N - initial_condition;

theoretical_mean = me/(1 - a(1) -a(2) -a(3)) % theoretical mean
my = mean(y); % sample estimate of the y mean
me = mean(e);

% ============================== TREND REMOVAL ============================
figure
subplot 211
plot(e(1:100), 'b', 'linewidth', 2); hold on;   
e_tilde = detrend(e); % detrend the data: remove mean and linear trends
plot(e_tilde(1:100), 'r', 'linewidth', 2); grid on;
legend('Raw input', 'Detrended input')

subplot 212
plot(y(1:100), 'b', 'linewidth', 2); 
hold on; 
y_tilde = detrend(y); 
plot(y_tilde(1:100), 'r', 'linewidth', 2)
legend('Raw output', 'Detrended output'); title('Output');
xlabel('Samples [-]'); ylabel('Values [-]'); grid on;
% =========================================================================



% =========================== STAGIONALITY REMOVAL ========================
% Identify the stagionality on one period, the periodicity is known = 40
for t = 1:40
    s_hat(t) = mean( y_tilde( t:40:(t+40*round(N/40)-1) )); % takes the values a multiples of 40
end

S_hat = s_hat;
% Extend the identified stagionality all over the period
for k=1:round(N/40)-1 
    S_hat = [S_hat s_hat];
end

figure
plot(y_tilde, 'b' ,'linewidth', 2); hold on;
plot(S_hat, 'g' ,'linewidth', 3);
plot([1:40], s_hat, 'k--' ,'linewidth', 3); grid on; xlabel('Samples [-]'); ylabel('Values [-]');
legend('Process with stagionality', 'Stagionality extended to all periods', 'Identified stagionality');
xlim([0, 500]); grid on;

% Remove stagionality
y_tilde_s = y_tilde - S_hat;
figure
plot(y_tilde, 'b--' ,'linewidth', 2); hold on;
plot(y_tilde_s, 'k' ,'linewidth', 3);
xlabel('Samples [-]'); ylabel('Values [-]'); grid on;
legend('Process with stagionality', 'Process without stagionality');
xlim([100, 500])
% =========================================================================

%% Dataset object creation

data = iddata(y_tilde_s', e_tilde', Ts); % create identification data object

N_ident = round(N/2);
N_valid = N - N_ident;

data_ident = iddata(y_tilde_s(1:N_ident)', e_tilde(1:N_ident)', Ts); % create identification data object
get(data_ident) % visualize properties

data_valid = iddata(y_tilde_s(N_valid + 1:end)', e_tilde(N_valid +1 :end)', Ts); % create identification data object



%% Compute the predictor

model = ar( data_ident.OutputData, 3, 'Ts', Ts); % uso un modello AR(3)
display(model)
G


% Prediction on identification data
y_hat = predict(model, data_valid, 1, 'estimate');


figure
plot(y_hat.OutputData, 'g', 'linewidth', 4), hold on;
plot(data_valid.OutputData, 'k--', 'linewidth', 2);
legend('1-step prediction', 'True data')
xlim([0, 200]); grid on; ylim([-5, 5]); xlabel('Samples [-]'); ylabel('Values [-]');


%% Nonparametric identification

% Sample covariance function
gamma_hat = covf(data , 1); % sample variance (tau = 0)
display(gamma_hat)


% Spectral density analytical way
figure
num_frequencies = 100;  % points onto which to evaluate the frequency response
[H, W] = freqz( cell2mat(G.Numerator), cell2mat(G.Denominator), num_frequencies);
plot(W, abs(H).^2, 'b*', 'linewidth', 2) %modulo della risposta in frequenza calcolato al quadrato
hold on;

w = -pi:pi/num_frequencies:pi; % frequencies where compute the spectrum
Gy = 1./abs(1 - a(1)*exp(1j*w) - a(2)*exp(2*1j*w) - a(3)*exp(3*1j*w) ).^2 * 1;   
plot(w, Gy, 'k--' ,'linewidth', 2); grid on; xlim([-pi, pi]);
xlabel('\omega [rad/sample]'); ylabel('\Gamma_y(\omega)');
legend('True spectrum from Freqz', 'True spectrum from definition')




% Estimation of the spectrum with formulas seen at lesson
figure
Gy_hat = (fft(y_tilde_s).*conj(fft(y_tilde_s))) / N; % module squared of the complex number returned by fft
w_hat = [1:round(N/2)]/(round(N/2))*pi;
plot(w_hat, Gy_hat(1:round(N/2)), 'b', 'linewidth', 2); hold on; 

% Regularization (10 folds with 1000 elements each)
n_folds  = 10;                       % folds number
elem_fold = round(N/n_folds);          % elements in each fold
for k = 1 : n_folds
    temp = fft(y_tilde_s(1+(k-1)*elem_fold:elem_fold+(k-1)*elem_fold));
    GG(k,:) = temp.*conj(temp)/elem_fold;
end
Gy_hat_hat = mean(GG);                          % estimate as the mean of the single spectra
w_hat_hat = [1:elem_fold/2]/(elem_fold/2)*pi;   % omega
plot(w_hat_hat, Gy_hat_hat(1:round(elem_fold/2)), 'g', 'linewidth', 4);

plot(w, Gy, 'k--' ,'linewidth', 3); grid on; xlim([-pi, pi]);
grid on; xlabel('\omega [rad/sample]'); ylabel('\Gamma_y(\omega)');
legend('Raw estimation', 'Regularized estimation', 'True spectrum'); 
xlim([-pi, pi]);



% Estimation of the spectrum with matlab functions
figure
% Spectral Density (PSD) estimate via periodogram method
[Pyy, w_hat_per] = periodogram(y_tilde_s); % use a power of 2 number of data in order to use the FFT more efficiently
plot(w_hat_per, Pyy, 'b', 'linewidth', 2); grid on; xlabel('\omega [rad/sample]'); ylabel('\Gamma_y(\omega)');
hold on; 
% Regularization
[Pyy , w_hat_welch] = pwelch(y_tilde_s, bartlett(1024), 0);
plot(w_hat_welch, Pyy, 'g', 'linewidth', 4); grid on; xlabel('\omega [rad/sample]'); ylabel('\Gamma_y(\omega)');
plot(w, Gy, 'k--' ,'linewidth', 3); grid on; xlim([-pi, pi]);
legend('Raw estimation', 'Regularized estimation', 'True spectrum'); 
xlim([-pi, pi]);



%% Estimate frequency response function from data

% Directly estimate frequency response from data. The frequency from which
% the response is evaluated are w = [1:128]/128*pi/Ts as default
G_hat = spa(data);
figure
bode(G, 'k--', G_hat, 'b'); grid on;
legend('True FRF', 'Estimated FRF'); 
set(findall(gcf,'type','line'), 'linewidth', 2)


%% Parametric estimation

m1 = ar(data.OutputData, 1, 'ls', 'Ts', Ts); % AR(1) model
m2 = ar(data.OutputData, 3, 'ls', 'Ts', Ts); % AR(3) model
% armax(Z, [na nb nc nk])
m3 = armax(data, [0 0 3 1], 'Ts', Ts); % MA(3) model
m4 = armax(data, [1 0 2 1], 'Ts', Ts); % ARMA(1, 2) model
m5 = ar(data.OutputData, 4, 'ls', 'Ts', Ts); % AR(4) model
m6 = ar(data.OutputData, 5, 'ls', 'Ts', Ts); % AR(5) model

% Compute the prediction error: if it is a white noise, we have chosen the
% rigth model
figure; resid(m1, data); grid on; set(findall(gcf,'type','line'), 'linewidth', 2)

figure; resid(m2, data); grid on; set(findall(gcf,'type','line'), 'linewidth', 2)% resid OK! 

figure; resid(m3, data); grid on; 
set(findall(gcf,'type','line'), 'linewidth', 2)

figure; resid(m4, data); grid on; set(findall(gcf,'type','line'), 'linewidth', 2)

figure; resid(m5, data); grid on; % resid OK! but overparametrized ==> overfitting risk
set(findall(gcf,'type','line'), 'linewidth', 2)

figure; resid(m6, data); grid on; % resid OK! but overparametrized ==> overfitting risk
set(findall(gcf,'type','line'), 'linewidth', 2)


% Chosen model family: AR(n). The  next step is to choose the optimal model
% Compute BIC/FPE criteria
max_ar_order = 20;
for i = 1:max_ar_order   
    fprintf('Computing model %i', i); fprintf('\n');
    model = ar(data.OutputData, i, 'ls');
    AICs(i) = aic(model, 'nAIC');
    FPEs(i) = fpe(model);
    BICs(i) = aic(model, 'BIC');
    MDLs(i) = mdl(model);
end

figure; plot(FPEs, 'bo-', 'linewidth', 2); 
ylim([min(FPEs), max(FPEs)]); xlabel('Model order'); ylabel('FPE(n)'); grid on; xticks([1:max_ar_order]); 
[~, best_model_FPE] = min(FPEs); % Overestimate model order

figure; plot(AICs, 'bo-', 'linewidth', 2); grid on; 
ylim([min(AICs), max(AICs)]); xticks([1:max_ar_order]); xlabel('Model order'); ylabel('AIC(n)');
[~, best_model_AIC] = min(AICs); % Overestimate model order

figure; plot(BICs, 'bo-', 'linewidth', 2); 
ylim([min(BICs), max(BICs)]); xlabel('Model order'); ylabel('BIC(n)'); grid on; xticks([1:max_ar_order]); 
[~, best_model_BIC] = min(BICs); % Correct model order

figure; plot(MDLs, 'bo-', 'linewidth', 2); 
ylim([min(MDLs), max(MDLs)]); xlabel('Model order'); ylabel('MDL(n)'); grid on; xticks([1:max_ar_order]); 
[~, best_model_MDL] = min(MDLs); % Correct model order



%% Compute final model and validate on validation data

final_model = ar(data_ident.OutputData, best_model_MDL, 'ls', 'Ts', Ts);

figure; resid(final_model, data_valid); set(findall(gcf,'type','line'), 'linewidth', 2);
figure; compare(final_model, data_valid, 1); % 1-step prediction comparison
set(findall(gcf,'type','line'), 'linewidth', 2); grid on;

y_hat = predict(final_model, data_valid, 1)'; % 1-step prediction
RSS = 1/length(data_valid.OutputData) * sum( (data_valid.OutputData - y_hat.OutputData).^2 );
RSS % the variance of the error


%% Final time series 1-step prediction

% recover signal mean and stagionality, adding the estimate of both  quantities
y_hat_s = y_hat.OutputData' + S_hat(1:length(data_valid.OutputData)) + my;
figure;
plot(y_hat_s, 'g', 'linewidth', 2); hold on;
plot(y, 'k--', 'linewidth', 2);
legend('1-step Prediction', 'True data'); grid on; xlim([100, 500]);
xlabel('Samples [-]'); ylabel('Values [-]');
