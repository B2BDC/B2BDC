function [A,B,C] = msd1d(mVec,k,c)
% mVec: Vector of N masses
% Scalar (identical) stiffness
% Scalar (identical) damping
% Returns state-space matrices [A,B,C] such that x_dot = Ax+Bu, y = Cx
%     Input: u, the forces at each mass
%     Output: [all N positions (x(1:N)); all N velocities (x(N+1:2N))]
%
% You can obtain (for example), the frequency-response function (ie., the
% steady-state gain due to harmonic inputs) with
%
%   H = C*inv(jwI-A)*B
%
% where w is a frequency of interest. Specific elements of H can be
% accessed using normal array indexing.  Hence the input/output behavior
% from the 3rd input to the 2nd output is characterized by H(2,3).  So,
% the transfer function representing the Force at mass #1 and velocity response of
% mass #1 should be at H(N+1,1), since velocities are after position in the
% differential equation model, as written.

N = numel(mVec);
Kmat = zeros(N,N);
Cmat = zeros(N,N);
Kmat(1,1) = -k;
Cmat(1,1) = -c;
for i=2:N
    Kmat(i,i) = -2*k;
    Kmat(i,i-1) = k;
    Kmat(i-1,i) = k;
    Cmat(i,i) = -2*c;
    Cmat(i,i-1) = c;
    Cmat(i-1,i) = c;
end
A = [zeros(N) eye(N);inv(diag(mVec))*Kmat inv(diag(mVec))*Cmat];
B = [zeros(N);inv(diag(mVec))];
C = [eye(N) zeros(N);zeros(N) eye(N)];
%S = ss(A,B,C,0);


