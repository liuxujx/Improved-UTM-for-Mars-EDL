function outp=entry_tauf_bcs(Y0,Yf,tauf) % 
global R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 s0 Vf lambda_rf lambda_gammaf lambda_sf Hf ...
    qmax Qmax gmax ec eq eQ eg

% 求H(tauf)，但忽略下标tauf
r=Yf(1);
V=Yf(2);
gamma=Yf(3);
% s=Yf(4);
lambda_r=Yf(5-1);
lambda_V=Yf(6-1);
lambda_gamma=Yf(7-1);
% lambda_s=Yf(8);

h=(r-1)*R0;
rho=rho0*exp(-h/hs);
q=rho*(V*Vs)^2/2;
aq=q/qmax;
Q=kQ/sqrt(Rn)*sqrt(rho)*(V*Vs)^3; % 用3.15导致Qdot急剧增大
aQ=Q/Qmax;
g=R0*rho*V^2*Sr*sqrt(CL^2+CD^2)/(2*mass);
ag=g/gmax;

L=R0*rho*V^2*Sr*CL/(2*mass);
D=R0*rho*V^2*Sr*CD/(2*mass);

u=atan(-lambda_gamma*L*c1/V/ec);

rdot=V*sin(gamma);
Vdot=-D-sin(gamma)/r^2;
gammadot=L*(c0+c1*sin(u))/V + V*cos(gamma)/r - cos(gamma)/r^2/V;
sdot=V*cos(gamma)/r;

%H(tauf)表达式，注意要乘以 tauf_init*(
    % 有约束
H = tauf*( lambda_r*rdot + lambda_V*Vdot + lambda_gamma*gammadot + ...
    ec*cos(u) + eq*sec(pi/2*aq) + eQ*sec(pi/2*aQ) + eg*sec(pi/2*ag) ); % 
outp=[Y0(1)-r0;
    Y0(2)-V0;
    Y0(3)-gamma0;
    Yf(2)-Vf;
    Yf(4)-lambda_rf;
    Yf(6)-lambda_gammaf;
    H-Hf];

% H = tauf*( lambda_r*rdot + lambda_V*Vdot + lambda_gamma*gammadot + lambda_s*sdot + ...
%     ec*cos(u) + eq*sec(pi/2*aq) + eQ*sec(pi/2*aQ) + eg*sec(pi/2*ag) ); % 
% outp=[Y0(1)-r0;
%     Y0(2)-V0;
%     Y0(3)-gamma0;
%     Y0(4)-s0;
%     Yf(2)-Vf;
%     Yf(5)-lambda_rf;
%     Yf(7)-lambda_gammaf;
%     Yf(8)-lambda_sf;
%     H-Hf];
end