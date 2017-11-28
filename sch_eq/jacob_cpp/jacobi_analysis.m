%% n and iterations
% number of points in the matrix are n-1
% the iterations --> how many times i should aplly the jacobi method to
% arrive with tolerance 1e-10
filename = 'prova.txt';
A=importdata(filename, ',');
%%
n = A.data(:,1);
it = A.data(:,2); %iterations
plot(n,it);
%% analysis
fil1 = 'rho_eig_500.txt';
B=importdata(fil1, ',');

fil2 = 'eigenkets_500.txt';
C=importdata(fil2, ',');
%
rho = B.data(:,1);
eigenval = B.data(:,2);
eigenf1 = C.data(:,1);
eigenf2 = C.data(:,2);
eigenf3 = C.data(:,3);
%parameters
beta_e = 1.44;%eV nm beta * e^2
h_bar = 6.582e-16; %eV s
m = 9.11e-31; %kg electron mass
alfa = h_bar^2 /(m*beta_e);
E =  (h_bar^2 /(m*(alfa)^2) ).* eigenval;
omega = 0.05;
% r_0 = (2*(omega)^2)^(-1/3);
%%
r = rho.*alfa;
hold on;
plot(r,(r).^2.*(eigenf1).^2,'-.');
plot(r,(r).^2.*(eigenf2).^2,'-.');
plot(r,(r).^2.*(eigenf3).^2,'-.');
%% 200
fil1 = 'rho_eig_v.txt';
A=importdata(fil1, ',');

fil2 = 'eigenkets.txt';
C=importdata(fil2, ',');

rho = A.data(:,1);
eigenval = A.data(:,2)
eigenf1 = C.data(:,1);
eigenf2 = C.data(:,2);
eigenf3 = C.data(:,3);

hold on;
plot(rho,(eigenf1).^2,'.');
plot(rho,(eigenf2).^2,'.');
plot(rho,(eigenf3).^2,'.');

