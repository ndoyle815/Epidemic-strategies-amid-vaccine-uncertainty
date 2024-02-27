% function which uses the NextGeneration Matrix (NGM) approach to compute
% the basic reproduction number R_0 of the SEIR age model

function R = Get_R0(para)

% matrix elements ordered x = (E_{1-3}a, I_{1-3}a^A, I_{1-3}a^S, I_{1-3}a^PH, I_a^H)

T = zeros(5*para.n);
Sigma = zeros(5*para.n);

% transmission matrix T
T(1:para.n,para.n+1:2*para.n) = para.tau.*para.beta;
T(1:para.n,2*para.n+1:3*para.n) = para.beta;
T(1:para.n,3*para.n+1:4*para.n) = para.beta;
T(1:para.n,4*para.n+1:5*para.n) = para.rho.*para.beta;

% transition matrix Sigma
In = eye(para.n);

% E{1-3}a
Sigma(1:para.n,1:para.n) = -para.epsilon.*In;

% IA
Sigma(para.n+1:2*para.n,1:para.n) = para.epsilon.*(1 - para.da).*In;
Sigma(para.n+1:2*para.n,para.n+1:2*para.n) = -para.gamma.*In;

% IS
Sigma(2*para.n+1:3*para.n,1:para.n) = para.epsilon.*(1 - para.ha).*para.da.*In;
Sigma(2*para.n+1:3*para.n,2*para.n+1:3*para.n) = -para.gamma.*In;

% IPH
Sigma(3*para.n+1:4*para.n,1:para.n) = para.epsilon.*para.ha.*para.da.*In;
Sigma(3*para.n+1:4*para.n,3*para.n+1:4*para.n) = -para.zeta.*In;

% IH
Sigma(4*para.n+1:5*para.n,3*para.n+1:4*para.n) = para.zeta.*In;
Sigma(4*para.n+1:5*para.n,4*para.n+1:5*para.n) = -para.delta.*In;

% Next generation matrix
K = -T*inv(Sigma);
R = max(eig(K));
