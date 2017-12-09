function [x,P]=tu_qw(x,P,omega,T,Rw)
% where omega is the measured angular rate, T the time since
% the last measurement, and Rw the process noise covariance matrix.

% state 10-dim : 4 quatern, 3 bias1, 3 bias2
q=x(1:4);
% bias_omega=x(5:7);
% bias_acc=x(8:10);

x(1:4)=(eye(4)+1/2*Somega(omega)*T)*q' + T/2*Sq(q)*Rw(1:3,1:3)*rand(3,1);

end
