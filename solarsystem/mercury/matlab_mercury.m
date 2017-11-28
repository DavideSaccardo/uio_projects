%% masses
MS = 2 * 10^30;
% in terms of muss of the sun
ME = 6 * 10^24/MS  ;
MJ = 1.9 * 10^27 /MS ;
MMa = 6.6 * 10^23/MS ;
MV = 4.9 * 10^24 /MS ;
MSa = 5.5 * 10^26 /MS ;
MMe = 3.3 * 10^23/MS ;
MU = 8.8 * 10^25 /MS ;
MN = 1.03 * 10^26/MS ; 
MP = 1.31 * 10^22/MS; 
% speed of light
c = 299792458; %m/s
AU = 149597870700; %km
y = 3600*24*365; 
c = c*y/AU;
%% animation
filename = 'positions_E102.dat';
A=importdata(filename);

%t=A(:,1);
%sun
x1_E = A(:,1);
y1_E = A(:,2);
z1_E = A(:,3);
%earth
x2_E = A(:,4);
y2_E = A(:,5);
z2_E = A(:,6);

%% animation
filename = 'positions_N102.dat';
A=importdata(filename);

%t=A(:,1);
%sun
x1_N = A(:,1);
y1_N = A(:,2);
z1_N = A(:,3);
%earth
x2_N = A(:,4);
y2_N = A(:,5);
z2_N = A(:,6);
%% position of mercury

%% orbits
hold on;
plot(x1_E,y1_E,'o')
plot(x2_N(1:10000),y2_N(1:10000))
plot(x2_E(1:10000),y2_E(1:10000))

%% oscillation
r_N = sqrt(x2_N.*x2_N + y2_N.*y2_N+ z2_N.*z2_N);
t = 1:length(r_N);
t = t.*1e-4;
r_E = sqrt(x2_E.*x2_E + y2_E.*y2_E+ z2_E.*z2_E);
%%
% looking for the angle
start = 100e4;
[Y_N,I_N] = min(r_N(start:100.16e4));
[Y_E,I_E] = min(r_E(start:100.16e4));
index_N =start+I_N;
index_E = start+I_E;
%%
hold on;
xlim([100 100.2]);
plot([t(index_N) t(index_N)],[0.3 0.48]);
plot(t,r_N,'.','MarkerSize',3);
plot(t,r_E,'.','MarkerSize',3);

title({'Distance r between Mercury and Sun',' as a function of time t'},'Interpreter','teX','Fontsize',18);
xlabel('$t$ $[years]$','interpreter','latex','fontsize',18);
ylabel('$r$ $[AU]$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','tex');
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times');
grid on;
ll=legend(gca,'show','Perihelion','location','northeast');
set(ll,'Interpreter','Latex');
%%
yp_N =  y2_N(index_N)-y1_N(index_N);
xp_N = x2_N(index_N)-x1_N(index_N);
%%
yp_E =  y2_E(index_E)-y1_E(index_E);
xp_E = x2_E(index_E)-x1_E(index_E);
%%
theta_N = atan2(yp_N,xp_N);
theta_E = atan2(yp_E,xp_E);
% diff = (theta_N - theta_E)
exp = 43*pi/(3600*180)
diff/exp
%% terra sole stabile
MS = 2 * 10^30;
% in terms of muss of the sun
MMe = 3.3 * 10^23/MS  ;
G = 6.16e-11;
L = 5.63377e-07;
E = -8.42e-5;
r = 0:1e-4:1;
V = ((((L.^2)./(2.*r.^2)) - 4*pi.^2./r)).*MMe;
[A,B]=min(V)
ylim([-2.5e-4 0]);
xlim([0 2])
hold on;
plot(r,V,'.');
plot([0 1.1],[E E])
plot([r(B) r(B)],[-2e-4 0])