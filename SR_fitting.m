%% Fitting the parameter curve

% define
syms x x0 y0
eta_x= 1.8928*exp(-2.4225*x)+1.2287*exp(-0.3899*x);
epsilon_x= 0.2174*exp(1.4357*x)+7.1217e-08*exp(13.9179*x);
xc_x=0.6178*x^1.6375+0.3952;

% diff
eta_df=diff(eta_x,x);
epsilon_df=diff(epsilon_x,x);
xc_df=diff(xc_x,x);

% equation
% (y0-y(x))*y'(x)+(x0-x)=0
eq_eta=(y0-eta_x)*eta_df+(x0-x);
eq_epsilon=(y0-epsilon_x)*epsilon_df+(x0-x);
eq_xc=(y0-xc_x)*xc_df+(x0-x);

% eq_eta=matlabFunction(eq_eta);
% eq_epsilon=matlabFunction(eq_epsilon);
% eq_xc=matlabFunction(eq_xc);
