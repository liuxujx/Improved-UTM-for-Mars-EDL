%% 复现utm，复现出来间接法+同伦+bvp4c
%%
% clc;clear all; close all;
global R0 g0 Vs Ts hs rho0 mass Sr CD CL Rn kQ ...
    u_min u_max c0 c1 h0 r0 V0 gamma0 s0 Vf lambda_rf lambda_gammaf lambda_sf Hf ...
    qmax Qmax gmax ec eq eQ eg

GM=4.284e13;% 2832e13;
R0=3397000; %m
g0=GM/R0^2; %m/s^2
Vs=sqrt(R0*g0);
Ts=sqrt(R0/g0);

hs=9354.5;
rho0=0.0158;

mass=3300;
Sr=15.9;
CD=1.45;
CL=0.348;
Rn=0.6;
kQ=1.9027e-4;

u_min=cosd(120);
u_max=cosd(30);
c0=(u_max+u_min)/2;
c1=(u_max-u_min)/2;

lambda_rf=-1;
lambda_gammaf=0;
lambda_sf=0;
Hf=0;

qmax=100e3;
Qmax=200e4;
gmax=50*9.80665/g0;

%% bvp4c 无约束
% ec=1;
% eq=0;
% eQ=0;
% eg=0;
% 
% h0=50e3;
% r0=1+h0/R0;
% V0=6e3/Vs;
% gamma0=-11.5*pi/180;
% s0=0;
% Vf=540/Vs;
% 
% Nt=500;
% tauf_init=200/Ts;
% tau=linspace(0,1,Nt);
% Y_init = [r0 V0 gamma0 0.1 0 0].';
% solinit = bvpinit(tau, Y_init, tauf_init);
% options = bvpset('Nmax',10*Nt,'RelTol',1e-4,'AbsTol',1e-6);% 默认Rel=1e-3 Abs=1e-6
% 
% tic
%     % 高度同伦
% h0h0=linspace(50,125,16)*1e3;
% for nh=1:length(h0h0)
%     nh
%     h0=h0h0(nh);
%     r0=1+h0/R0;
%     solinit.y = [[r0; solinit.y(2:6,1)],solinit.y(:,2:end)];% [[r0; Y(2:6,1)],Y(:,2:end)];
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit, options); % 
%     solinit = sol;
% end
%     % ec同伦
% ecec=[(9:-0.5:1)*1e-1, (9:-0.5:1)*1e-2, (9:-0.5:1)*1e-3, (9:-0.5:1)*1e-4, (9:-0.5:1)*1e-5, (9:-0.5:1)*1e-6, (9:-0.5:1)*1e-7];
% for nec=1:length(ecec)
%     nec
%     ec=ecec(nec);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% toc
% save sol1e-7.mat % 357.73s 
%% bvp4c 有约束
% ec=1;
% eq=1;
% eQ=1;
% eg=1;
% 
% h0=50e3;
% r0=1+h0/R0;
% V0=6e3/Vs;
% gamma0=-11.5*pi/180;
% s0=0;
% Vf=540/Vs;
% 
% Nt=500;
% tauf_init=200/Ts;
% tau=linspace(0,1,Nt);
% Y_init = [r0 V0 gamma0 0.1 0 0].';
% solinit = bvpinit(tau, Y_init, tauf_init);
% options = bvpset('Nmax',10*Nt,'RelTol',1e-4,'AbsTol',1e-6);% 默认Rel=1e-3 Abs=1e-6
% 
% tic
%     % 高度同伦
% h0h0=linspace(50,125,16)*1e3;
% for nh=1:length(h0h0)
%     nh
%     h0=h0h0(nh);
%     r0=1+h0/R0;
%     solinit.y = [[r0; solinit.y(2:6,1)],solinit.y(:,2:end)];% [[r0; Y(2:6,1)],Y(:,2:end)];
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit, options); % 
%     solinit = sol;
% end
% 
%     % Qmax同伦
% QQmax=linspace(200,70,131)*1e4;
% for nQ=1:length(QQmax)
%     nQ
%     Qmax=QQmax(nQ);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % gmax同伦
% ggmax=linspace(50,5,46)*9.80665/g0;
% for ng=1:length(ggmax)
%     ng
%     gmax=ggmax(ng);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % qmax同伦
% qqmax=linspace(100,10,10)*1e3;
% for nq=1:length(qqmax)
%     nq
%     qmax=qqmax(nq);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % eQ同伦
% eQeQ=[(9:-0.5:1)*1e-1, (9:-0.5:1)*1e-2, (9:-0.5:1)*1e-3, (9:-0.5:1)*1e-4, (9:-0.5:1)*1e-5, (9:-0.5:1)*1e-6, (9:-0.5:1)*1e-7];
% for neQ=1:length(eQeQ)
%     neQ
%     eQ=eQeQ(neQ);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % eg同伦
% egeg=[(9:-0.5:1)*1e-1, (9:-0.5:1)*1e-2, (9:-0.5:1)*1e-3, (9:-0.5:1)*1e-4, (9:-0.5:1)*1e-5, (9:-0.5:1)*1e-6, (9:-0.5:1)*1e-7];
% for neg=1:length(egeg)
%     neg
%     eg=egeg(neg);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % eq同伦
% eqeq=[(9:-0.5:1)*1e-1, (9:-0.5:1)*1e-2, (9:-0.5:1)*1e-3, (9:-0.5:1)*1e-4, (9:-0.5:1)*1e-5, (9:-0.5:1)*1e-6, (9:-0.5:1)*1e-7];
% for neq=1:length(eqeq)
%     neq
%     eq=eqeq(neq);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% 
%     % ec同伦 放到eQ, eg, eq后面，全部用-0.5+1e-3可以运行成功，结果可行，高度10.05km
% ecec=[(9:-0.2:1)*1e-1, (9:-0.2:1)*1e-2, (9:-0.2:1)*1e-3, (9:-0.2:1)*1e-4, (9:-0.2:1)*1e-5, ...
%     (9:-0.2:1)*1e-6, (9:-0.2:1)*1e-7, (9:-0.2:1)*1e-8, (9:-0.2:1)*1e-9, (9:-0.2:1)*1e-10];
% for nec=1:length(ecec)
%     nec
%     ec=ecec(nec);
%     sol = bvp4c(@entry_tauf_odes, @entry_tauf_bcs, solinit ,options); % 
%     solinit = sol;
% end
% toc
% save sol1e-10.mat % 1072.78s
%% plot
% load sol1e-7.mat % 无约束最终版
load sol1e-10.mat % 有约束最终版
tauf=sol.parameters;
tau=sol.x*tauf;
Y=sol.y;
    % 计算航程
y0=[r0,V0,gamma0,s0, Y(4:6,1).', 0];
[tout,yout]=ode45(@entry_lambda,tau,y0); % yout和Y还是有差距的
s = yout(:,4).';

Yf=[(Y(1,end)-1)*R0,Y(2,end)*Vs,Y(3,end)*180/pi,s(end)*R0];
r=Y(1,:);
V=Y(2,:);
gamma=Y(3,:);
lambda_r=Y(4,:);
lambda_V=Y(5,:);
lambda_gamma=Y(6,:);

h=(r-1)*R0;
rho=rho0*exp(-h/hs);
q=rho.*(V*Vs).^2/2;
aq=q/qmax;
Q=kQ/sqrt(Rn)*sqrt(rho).*(V*Vs).^3; % 用3.15导致Qdot急剧增大
aQ=Q/Qmax;
g=R0*rho.*V.^2*Sr*sqrt(CL^2+CD^2)/(2*mass);
ag=g/gmax;

L=R0*rho.*V.^2*Sr*CL/(2*mass);
D=R0*rho.*V.^2*Sr*CD/(2*mass);

u=atan(-lambda_gamma.*L.*c1./V/ec); %arctan
bank=acosd(c0+c1*sin(u));

figure(1);
subplot(211); hold on; box on;
plot(Y(2,:)*Vs/1e3, (Y(1,:)-1)*R0/1e3,'k:','linewidth',1.5);
subplot(212); hold on; box on;
plot(s*R0/1e3,(Y(1,:)-1)*R0/1e3,'k:','linewidth',1.5);

figure(2);
subplot(211); hold on; box on;
plot(tau*Ts,Y(3,:)/pi*180,'k:','linewidth',1.5);
subplot(212); hold on; box on;
plot(tau*Ts,bank,'k:','linewidth',1.5);
tUTM=tau*Ts; uUTM=bank;

figure(4);
subplot(311); plot(tau*Ts,lambda_r,'k:','linewidth',1.5);
subplot(312); plot(tau*Ts,lambda_V,'k:','linewidth',1.5);
subplot(313); plot(tau*Ts,lambda_gamma,'k:','linewidth',1.5);

figure(3);
subplot(311); plot(tau*Ts,q/1e3,'k:','linewidth',1.5);
subplot(312); plot(tau*Ts,Q/1e4,'k:','linewidth',1.5);
subplot(313); plot(tau*Ts,g*g0/9.80665,'k:','linewidth',1.5);


