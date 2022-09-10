function outp=entry_tauf_odes(t,inp,tauf) % 
global R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 s0 Vf lambda_rf lambda_gammaf lambda_sf Hf ...
    qmax Qmax gmax ec eq eQ eg

%%
r=inp(1);
V=inp(2);
gamma=inp(3);
% s=inp(4);
lambda_r=inp(5-1);
lambda_V=inp(6-1);
lambda_gamma=inp(7-1);
% lambda_s=inp(8);

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

u=atan(-lambda_gamma*L*c1/V/ec); %arctan

rdot=V*sin(gamma);
Vdot=-D-sin(gamma)/r^2;
gammadot=L*(c0+c1*sin(u))/V + V*cos(gamma)/r - cos(gamma)/r^2/V;
% sdot=V*cos(gamma)/r;

    % 有约束
lambda_rdot =...
    lambda_gamma*((V*cos(gamma))/r^2 - (2*cos(gamma))/(V*r^3) + (CL*R0*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*((CL^3*R0^4*Sr^3*V^3*c1^4*lambda_gamma^3*rho0^3*exp(-(3*R0*(r - 1))/hs))/(8*ec^3*hs*mass^3*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) - (CL*R0^2*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*hs*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*mass) + (CL*R0^2*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*(c0 - (CL*R0*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*hs*mass)) - lambda_V*((2*sin(gamma))/r^3 + (CD*R0^2*Sr*V^2*rho0*exp(-(R0*(r - 1))/hs))/(2*hs*mass)) + (R0^2*V^2*eq*g0*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax)))/(4*hs*qmax*cos((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax))^2) - (CL^2*R0^3*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec*hs*mass^2*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) + (R0^2*Sr*V^2*eg*rho0*pi*sin((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*hs*mass*cos((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))^2) + (R0*V^3*eQ*kQ*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))*(R0*g0)^(3/2))/(4*Qmax*Rn^(1/2)*hs*cos((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))^2*(rho0*exp(-(R0*(r - 1))/hs))^(1/2));
lambda_Vdot =...
    (CD*R0*Sr*V*lambda_V*rho0*exp(-(R0*(r - 1))/hs))/mass - lambda_r*sin(gamma) - lambda_gamma*(cos(gamma)/r + cos(gamma)/(V^2*r^2) + (CL*R0*Sr*rho0*exp(-(R0*(r - 1))/hs)*(c0 - (CL*R0*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*mass) - (CL*R0*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*((CL*R0*Sr*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2)) - (CL^3*R0^3*Sr^3*V^2*c1^4*lambda_gamma^3*rho0^3*exp(-(3*R0*(r - 1))/hs))/(8*ec^3*mass^3*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2))))/(2*mass)) - (3*V^2*eQ*kQ*pi*sin((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)*cos((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))^2) - (R0*V*eq*g0*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax)))/(2*qmax*cos((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax))^2) + (CL^2*R0^2*Sr^2*V*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec*mass^2*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) - (R0*Sr*V*eg*rho0*pi*sin((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(2*gmax*mass*cos((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))^2);
lambda_gammadot =...
    lambda_gamma*((V*sin(gamma))/r - sin(gamma)/(V*r^2)) + (lambda_V*cos(gamma))/r^2 - V*lambda_r*cos(gamma);
outp=tauf*[rdot, Vdot, gammadot, lambda_rdot, lambda_Vdot, lambda_gammadot].'; %, sdot , lambda_sdot tauf_init*

% lambda_rdot =...
%     lambda_gamma*((V*cos(gamma))/r^2 - (2*cos(gamma))/(V*r^3) + (CL*R0*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*((CL^3*R0^4*Sr^3*V^3*c1^4*lambda_gamma^3*rho0^3*exp(-(3*R0*(r - 1))/hs))/(8*ec^3*hs*mass^3*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) - (CL*R0^2*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*hs*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*mass) + (CL*R0^2*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*(c0 - (CL*R0*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*hs*mass)) - lambda_V*((2*sin(gamma))/r^3 + (CD*R0^2*Sr*V^2*rho0*exp(-(R0*(r - 1))/hs))/(2*hs*mass)) + (V*lambda_s*sin(gamma))/r^2 + (R0^2*V^2*eq*g0*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax)))/(4*hs*qmax*cos((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax))^2) - (CL^2*R0^3*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec*hs*mass^2*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) + (R0^2*Sr*V^2*eg*rho0*pi*sin((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*hs*mass*cos((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))^2) + (R0*V^3*eQ*kQ*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))*(R0*g0)^(3/2))/(4*Qmax*Rn^(1/2)*hs*cos((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))^2*(rho0*exp(-(R0*(r - 1))/hs))^(1/2));
% lambda_Vdot =...
%     (CD*R0*Sr*V*lambda_V*rho0*exp(-(R0*(r - 1))/hs))/mass - lambda_r*sin(gamma) - (lambda_s*sin(gamma))/r - lambda_gamma*(cos(gamma)/r + cos(gamma)/(V^2*r^2) + (CL*R0*Sr*rho0*exp(-(R0*(r - 1))/hs)*(c0 - (CL*R0*Sr*V*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2))))/(2*mass) - (CL*R0*Sr*V*rho0*exp(-(R0*(r - 1))/hs)*((CL*R0*Sr*c1^2*lambda_gamma*rho0*exp(-(R0*(r - 1))/hs))/(2*ec*mass*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(1/2)) - (CL^3*R0^3*Sr^3*V^2*c1^4*lambda_gamma^3*rho0^3*exp(-(3*R0*(r - 1))/hs))/(8*ec^3*mass^3*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2))))/(2*mass)) - (3*V^2*eQ*kQ*pi*sin((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)*cos((V^3*kQ*pi*(R0*g0)^(3/2)*(rho0*exp(-(R0*(r - 1))/hs))^(1/2))/(2*Qmax*Rn^(1/2)))^2) - (R0*V*eq*g0*rho0*pi*exp(-(R0*(r - 1))/hs)*sin((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax)))/(2*qmax*cos((R0*V^2*g0*rho0*pi*exp(-(R0*(r - 1))/hs))/(4*qmax))^2) + (CL^2*R0^2*Sr^2*V*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec*mass^2*((CL^2*R0^2*Sr^2*V^2*c1^2*lambda_gamma^2*rho0^2*exp(-(2*R0*(r - 1))/hs))/(4*ec^2*mass^2) + 1)^(3/2)) - (R0*Sr*V*eg*rho0*pi*sin((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(2*gmax*mass*cos((R0*Sr*V^2*rho0*pi*exp(-(R0*(r - 1))/hs)*(CD^2 + CL^2)^(1/2))/(4*gmax*mass))^2);
% lambda_gammadot =...
%     lambda_gamma*((V*sin(gamma))/r - sin(gamma)/(V*r^2)) + (lambda_V*cos(gamma))/r^2 - V*lambda_r*cos(gamma) - (V*lambda_s*cos(gamma))/r;
% lambda_sdot = 0;
% outp=tauf*[rdot, Vdot, gammadot, lambda_rdot, sdot, lambda_Vdot, lambda_gammadot, lambda_sdot].'; % tauf_init*
end