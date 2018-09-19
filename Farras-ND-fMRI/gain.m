function f=gain(VAR,C1,C2,C3)
%nonlinear gain function for coupled ODE neural lump
f=C3*(1+tanh((VAR-C1)./C2));