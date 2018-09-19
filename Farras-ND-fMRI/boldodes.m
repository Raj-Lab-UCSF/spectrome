function fn=boldodes(t,x)
% Haemodynamic model embedding the Balloon-Windkessel model
% x(1) = s = vasodilatory signal
% x(2) = f = inflow
% x(3) = v = blood volume
% x(4) = q = deoxyhaemoglobin content
% z = neuronal activity
global kappa gamma tau alpha rho z
%if t<0.05, z=1; else z=0; end;
%t
%pause
tind = round(t*1000+1);
%if (t==0)
%    zt = 0; 
%else
    zt = z(tind);
%end;
fn(1) = zt-kappa*x(1) - gamma*(x(2)-1);
fn(2) = x(1);
fn(3) = (x(2) - (x(3)).^(1/alpha))/tau;
fn(4) = (x(2).*(1-(1-rho).^(1/x(2)))/rho - (x(3)).^((1-alpha)/alpha).*x(4))/tau;
fn=fn';