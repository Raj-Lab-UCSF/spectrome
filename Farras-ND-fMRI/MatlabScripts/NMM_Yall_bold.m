function  Ybold = NMM_Yall_bold(Yall)

Yall = Yall';
tfrom = 1;
tto = size(Yall,1);
N = size(Yall,2);

% Solve BOLD ODEs
% x(1) = s = vasodilatory signal
% x(2) = f = inflow
% x(3) = v = blood volume
% x(4) = q = deoxyhaemoglobin content
global kappa gamma tau alpha rho z

% parameters
kappa = 0.65;   % rate of signal decay, s^-1                0.65
gamma = 0.41;   % rate of flow-dependent elimination, s^-1  0.41
tau = 0.98;     % haemodynamic transit time, s              0.98
alpha = 0.32;   % Grubb's exponent                          0.32
rho = 0.34;     % resting oxygen extraction fraction        0.34
V0 = 0.02;      % resting blood volume fraction             0.02

% initialize bold signal
Ybold = zeros(tto,N);

% set tolerances
options = odeset('RelTol',1e-3,'AbsTol',1e-6);

% loop over time series
for n=1:N
    
    % get time series
    z = Yall(tfrom:tto,n);
    
    % get first 30 secs and append in front in order to get rid of initial transient
    ttrns = 30000;
    zt = Yall(tfrom:ttrns+1,n);
    
    % get abs(diff(glu))
    z = [abs(diff([zt;z])); 0];

    % ICs
    ics = [0 1 1 1];

    % solve
    tzero = 0;
    tend = (tto+ttrns-1)/1000;
    tsteps = 100;
    %[t,x] = rgk4('boldodes',tzero,tend,ics,(tend-tzero)*tsteps);
    [t,x] = ode45('boldodes',[tzero:0.001:tend],ics,options);

    % BOLD signal
    k1 = 7*rho;
    k2 = 2;
    k3 = 2*rho - 0.2;
    y = V0*(k1*(1-x(:,4))+k2*(1-x(:,4)./x(:,3))+k3*(1-x(:,3)));

    % remove transient
    Ybold(:,n) = y(ttrns+1:end);

end;
