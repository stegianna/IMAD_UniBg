function x = mdl(m)
%MDL  Returns Rissanen's Minimum Description Length.

d = length(m.Report.Parameters.ParVector);   % d=number of model parameters
N = m.Report.DataUsed.Length;  % N=number of data points fitted
V = m.Report.Fit.LossFcn;     % V=loss function;
x = log(V) + log(N) * d/N; % Rissanen's Minimum Description Length