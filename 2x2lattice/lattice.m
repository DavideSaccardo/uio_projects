addpath 4cd
addpath 4b
%% 4b

file = 'dT1_L2.txt';
AA = importdata(file, ',');
file = 'dT1_L2_up.txt';
AAup = importdata(file, ',');
z = 12 + 2*exp(8) - 2*exp(-8);
t_e =2*8*exp(-8) + 2*(-8)*exp(8);
t_ee = 2*8*8*exp(-8) + 2*(-8)*(-8)*exp(8);
t_m = 2*4*exp(8)  + 2*4*2 ;
t_mm = 2*4*4*exp(+8) + 2*4*(2)*(2);
%% energy
hold on;
xlim([1.5 7.5])
ylim([-8.05 -7.75])
plot(log10(AA(:,1)),AA(:,2),'.');
plot(log10(AAup(:,1)),AAup(:,2),'.');
plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[t_e/z t_e/z])
plottalo({'Mean energy <E> as a function of';'the Monte Carlo cycles n';'T=1kT/J and L \times L = 2 \times 2'}, 20, ...
    '$log_{10}(n)$','$<E>$ unit J',18, 16);
ll=legend(gca,'show','Numerical data random initial condition',...
    'Numerical data all spins up','Theoretical value of < E > = -7.9839J','location','northeast');
%set(ll,'Interpreter','Latex');
%% |m|
hold on;
plot(log10(AA(:,1)),AA(:,4),'.');
plot(log10(AAup(:,1)),AAup(:,4),'.');
plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[t_m/z t_m/z]);
plottalo({'Mean absolute value of the magnetic momen <|M|>';' as a function of the Monte Carlo cycles n';'T=1kT/J and L \times L = 2 \times 2'}, 20, ...
    '$log_{10}(n)$','$<|M|>$',18, 16);
ll=legend(gca,'show','Numerical data random initial condition',...
    'Numerical data all spins up','Theoretical value of < |M| > =  3.9946','location','southeast');

%% cv
cv_e = AA(:,3);
cv_t = t_ee/z - t_e.*t_e./(z.*z);
hold on;
plot(log10(AA(:,1)),cv_e,'.');
plot(log10(AAup(:,1)),AAup(:,3),'.');
plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[cv_t cv_t]);
plottalo({'Heat capacity C_v as a function of';'the Monte Carlo cycles n';'T=1kT/J and L \times L = 2 \times 2'}, 20, ...
    '$log_{10}(n)$','$C_v$ unit $J^2/kT$',18, 16);
ll=legend(gca,'show','Numerical data random initial condition',...
    'Numerical data all spins up','Theoretical value of C_v =  0.1283 J^2/kT','location','northeast');

%% chi
chi_e = AA(:,5);
chi_t = t_mm/z ;
hold on;
plot(log10(AA(:,1)),chi_e,'.');
plot(log10(AAup(:,1)),AAup(:,5),'.');
plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[chi_t chi_t]);
plottalo({'Susceptability \chi as a function of';'the Monte Carlo cycles n';'T=1kT/J and L \times L = 2 \times 2'}, 20, ...
    '$log_{10}(n)$','$\chi$ unit $J^{-1}$',18, 16);
ll=legend(gca,'show','Numerical data random initial condition',...
    'Numerical data all spins up','Theoretical value of \chi =  -15.9572 J^{-1}','location','northeast');
%% m
hold on;
plot(log10(AA(:,1)),AA(:,6),'.');
plot(log10(AAup(:,1)),AAup(:,6),'.');
plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[0 0]);
plottalo({'Mean magnetic moment <M> as a function of';'the Monte Carlo cycles n';'T=1kT/J and L \times L = 2 \times 2'}, 20, ...
    '$log_{10}(n)$','$<M>$',18, 16);
ll=legend(gca,'show','Numerical data random initial condition',...
    'Numerical data all spins up','Theoretical value of < M > =  0','location','southeast');


%%
table(AA(1e5,2), t_e, cv(1e5), cv_t, AA(1e5,3), t_m, chi_e(1e5), chi_t);

%% 4c 
file = '950000_1.000000_20.txt';
AA1 = importdata(file, ',');
file = '950000_1.700000_20.txt';
AA2 = importdata(file, ',');
file = '950000_2.400000_20.txt';
AA3 = importdata(file, ',');
file = '950000_3.100000_20.txt';
AA4 = importdata(file, ',');
%% energy
hold on;
%xlim([1.5 7.5])
ylim([-805 -250])
plot(log10(AA1(:,1)),AA1(:,2),'.');
plot(log10(AA2(:,1)),AA2(:,2),'.');
plot(log10(AA3(:,1)),AA3(:,2),'.');
plot(log10(AA4(:,1)),AA4(:,2),'.');
%plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[t_e/z t_e/z])
plottalo({'Mean energy <E> as a function of';'the Monte Carlo cycles n'...
    ;'L \times L = 20 \times 20'}, 20, ...
    '$log_{10}(n)$','$<E>$ unit J',18, 16);
ll=legend(gca,'show','T=1kT/J',...
    'T=1.7kT/J','T=2.4kT/J','T=3.1kT/J','location','northeast');
%set(ll,'Interpreter','Latex');
%% |m|
hold on;
ylim([30 410])
plot(log10(AA1(:,1)),AA1(:,4),'.');
plot(log10(AA2(:,1)),AA2(:,4),'.');
plot(log10(AA3(:,1)),AA3(:,4),'.');
plot(log10(AA4(:,1)),AA4(:,4),'.');
%plot(log10(AAup(:,1)),AAup(:,4),'.');
%plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[t_m/z t_m/z]);
plottalo({'Mean absolute value of magnetization< <|M|> as a function of';'the Monte Carlo cycles n'
    ;'L \times L = 20 \times 20'}, 20, ...
    '$log_{10}(n)$','$<|M|>$',18, 16);
ll=legend(gca,'show','T=1kT/J',...
    'T=1.7kT/J','T=2.4kT/J','T=3.1kT/J','location','northeast');

%%
file = 'hist_1.000000.txt';
B1 = importdata(file, ',');
file = 'hist_1.700000.txt';
B2 = importdata(file, ',');
file = 'hist_2.400000.txt';
B3 = importdata(file, ',');
file = 'hist_3.100000.txt';
B4 = importdata(file, ',');
file = 'hist_3.800000.txt';
B5 = importdata(file, ',');

for ii=1:1:5
s = sprintf('m_%d=sum(B%d(:,1).*B%d(:,2));',ii,ii,ii);
eval(s);
s = sprintf('v_%d=sum((B%d(:,1)-m_%d).^2 .*B%d(:,2))',ii,ii,ii,ii);
eval(s);
end

%% accpeted moves
hold on;
%xlim([1.5 7.5])
%ylim([-805 -250])
plot((AA1(:,1)),(AA1(:,7)),'.');
plot((AA2(:,1)),(AA2(:,7)),'.');
plot((AA3(:,1)),(AA3(:,7)),'.');
plot((AA4(:,1)),(AA4(:,7)),'.');
%plot([log10(min(AA(:,1))) log10(max(AA(:,1)))],[t_e/z t_e/z])
plottalo({'Accepted configurations N_A as a function of';'the Monte Carlo cycles n'...
    ;'L \times L = 20 \times 20'}, 20, ...
    '$n$','$N_A$',18, 16);
ll=legend(gca,'show','T=1kT/J',...
    'T=1.7kT/J','T=2.4kT/J','T=3.1kT/J','location','northeast');
%set(ll,'Interpreter','Latex');







%% hist 

hold on
%xticks([-800:8:-700])
%xlim([-802 -782]);
%file = 'accepted_conf_T3.txt'
n_E = B2(:,2);
E = B2(:,1);

step=1;
color = 'c';
hold on
h_calibro =bar(E, n_E, step,'facecolor', color);
axs1 = gca; 
% Personalizzazione: trasparenza colonne
set(h_calibro, 'FaceAlpha', 0.3);

% Personalizzazione: colore bordi
set(h_calibro, 'EdgeColor', color);

% Accediamo agli assi
%gca ? il comando che da gli assi

% Scelta della colorazione degli istogrammi:
% true    usa i colori
% false   non usa i colori
hist_color_full = true;
%
%
plottalo({'Probability P(E) of the system';' to have an energy E at T=1kT/J'...
    ;' L \times L = 20 \times 20'}, 18, 'E unit J', 'P(E)', 16, 14);




















%% plot cv

figure();

for ii=1:1:4;
    hold on
    s1=sprintf('p=plot(A_%d(:,1), A_%d(:,5));',ii,ii);
    eval(s1);
    set(p,'Marker','.');
end

plottalo({'Susceptibility \chi';'as a function of temperature T'}, 20, ...
    'T unit kT/J','$\chi$ unit $J^{-1}$',18, 16);
ll=legend(gca,'show','L=40','L=60','L=80','L=100','location','southeast');

%% energy

for ii=1:1:4;
    hold on
    s1=sprintf('p=plot(A_%d(:,1), A_%d(:,2));',ii,ii);
    eval(s1);
    set(p,'Marker','.');
end

plottalo({'Mean energy <E> ';'as a function of temperature T'}, 20, ...
    'T unit kT/J','$<E>$',18, 16);
ll=legend(gca,'show','L=40','L=60','L=80','L=100','location','southeast');

%%
figure();
for ii=1:1:4;
    hold on
    s1=sprintf('p=plot(A_%d(:,1), A_%d(:,4));',ii,ii);
    eval(s1);
    set(p,'Marker','.');
end

plottalo({'Mean of the absolute value of magnetization <|M|> ';'as a function of temperature T'}, 20, ...
    'T unit kT/J','$<|M|>$',18, 16);
ll=legend(gca,'show','L=40','L=60','L=80','L=100','location','southeast');

%%

for ii=1:1:4;
    s1=sprintf('M_cv_%d = max(A_%d(:,3));',ii,ii,ii,ii);
    eval(s1);
    s2 = sprintf('M_chi_%d = max(A_%d(:,5));',ii,ii,ii,ii);
    eval(s2)
    s1=sprintf('index_cv_%d = find(abs(A_%d(:,3) - M_cv_%d)==0);',ii,ii,ii);
    eval(s1);
    s2 = sprintf('index_chi_%d = find(abs(A_%d(:,5)- M_chi_%d)==0);',ii,ii,ii);
    eval(s2);
    s1 = sprintf('Tc_cv_%d = A_%d(index_cv_%d,1);',ii,ii,ii);
    eval(s1);
    s2 = sprintf('Tc_chi_%d = A_%d(index_chi_%d,1);',ii,ii,ii);
    eval(s2); 
    
end
Tc_cv = [Tc_cv_1, Tc_cv_2,Tc_cv_3, Tc_cv_4];
Tc_chi = [Tc_chi_1, Tc_cv_2, Tc_cv_3, Tc_cv_4];
%%
L=40:20:100;
[L; Tc_cv; Tc_chi]

%%


for ii=1:1:4;
    s=sprintf('Tc_%d = mean([Tc_chi_%d, Tc_cv_%d]);',ii,ii,ii);
    eval(s)
end
diff_Tc_exp=zeros(size(40:20:180));

for ii=1:1:4;
    s=sprintf('diff_Tc_exp(%d) = Tc_%d - 2.269;',ii,ii);
    eval(s)
end


%% Tc function of L
L=40:20:100;
a=1;
Tc_t = a*L.^(-1);
hold on;
plot(L,Tc_t);
plot(L, Tc_cv - 2.269);

% 
% M_cv = max(A(start:end,3));
% M_chi = max(A(start:end,5));
% index_cv = find(A(start:end,3) == M_cv);
% index_chi = find(A(start:end,5) == M_chi);


%%
hold on;
%plot(A(6:701,1), A(6:701,3)./1e4,'.');
%plot([2.269 2.269],[0 0.01])
sta=30;
fin=181;
x_s = A120(sta:fin,1);
y_s =  ((A120(sta:fin,4)./1e4));
plot(x_s, y_s,'.')

%% derivative
d_y_s = zeros(size(y_s));
for tt=1:1:(length(y_s)-1)
    d_y_s(tt) = (y_s(tt+1)-y_s(tt))/(x_s(tt+1)-x_s(tt));
end
plot(x_s, d_y_s,'-.')
    
    
%hist(A(:,2))
%%
%xx=1.0:0.005:4;
f = (abs(+2.269-x_s)).^(1/8);
plot(x_s, -f)


%%
file = 'hist_1.700000.txt';
%%
file = 'hist_2.400000.txt';
%%
file = 'hist_3.100000.txt';
%%
file = 'hist_3.500000.txt';
%%
file = '4d_T3.txt';
%% 4b
file = 'dT1_L2.txt';
AA = importdata(file, ',');
z = 12 + 2*exp(8) - 2*exp(-8);
t_e =2*8*exp(-8) + 2*(-8)*exp(8);
t_ee = 2*8*8*exp(-8) + 2*(-8)*(-8)*exp(8);
t_m = 2*4*exp(8)  + 2*4*2 ;
t_mm = 2*4*4*exp(+8) + 2*4*(2)*(2);
%energy
figure();
hold on;
plot((AA(:,1)),AA(:,3));
plot([(min(AA(:,1))) (max(AA(:,1)))],[t_ee/z t_ee/z])

%%

file40 = 'l40_0_005.txt';
file60 = 'l60_0_005.txt';
file80 = 'l80_0_005.txt';
file100 = 'l100_0_005.txt';
file120 = 'l120_0_005.txt';
file140 = 'l140_0_005.txt';
file160 = 'l160_0_005.txt';
file180 = 'l180_0_005.txt';
file200 = 'l200_0_005.txt';

A40 = importdata(file40, ',');
A60 = importdata(file60, ',');
A80 = importdata(file80, ',');
A100 = importdata(file100, ',');
A120 = importdata(file120, ',');
A140 = importdata(file140, ',');
A160 = importdata(file160, ',');
A180 = importdata(file180, ',');
A200 = importdata(file200, ',');
%%

sta=30;
fin=181;

for ii=40:20:180;
    s1=sprintf('M_cv_%d = max(A%d(sta:fin,3));',ii,ii,ii,ii);
    eval(s1);
    s2 = sprintf('M_chi_%d = max(A%d(sta:fin,5));',ii,ii,ii,ii);
    eval(s2)
    s1=sprintf('index_cv_%d = find(abs(A%d(:,3) - M_cv_%d)<0.01);',ii,ii,ii);
    eval(s1);
    s2 = sprintf('index_chi_%d = find(abs(A%d(:,5)- M_chi_%d)<0.01);',ii,ii,ii);
    eval(s2);
    s1 = sprintf('Tc_cv_%d = A%d(index_cv_%d,1);',ii,ii,ii);
    eval(s1);
    s2 = sprintf('Tc_chi_%d = A%d(index_chi_%d,1);',ii,ii,ii);
    eval(s2); 
    
end



for ii=40:20:180;
    s=sprintf('Tc_%d = mean([Tc_chi_%d, Tc_cv_%d]);',ii,ii,ii);
    eval(s)
end
diff_Tc_exp=zeros(size(40:20:180));
for ii=40:20:180;
    s=sprintf('diff_Tc_exp(%d/20 -1) = Tc_%d - 2.269;',ii,ii);
    eval(s)
end


%% Tc 
hold on;
%plot(A160(sta:fin,1), A160(sta:fin,5)/(80*80));
%plot([Tc_chi_40 Tc_chi_40],[0 0.05]);
xx = A100(sta:fin,1);
yy = A100(sta:fin,3)/(100*100);
plot(xx, (yy),'.')
%%
rr = -(abs(xx-Tc_chi_40)).^(1/8);
plot(xx,rr)

%% Tc function of L
L=40:20:180;
a=1;
Tc_t = a*L.^(-1);
hold on;
plot(L,Tc_t);
plot(L, diff_Tc_exp);

% 
% M_cv = max(A(start:end,3));
% M_chi = max(A(start:end,5));
% index_cv = find(A(start:end,3) == M_cv);
% index_chi = find(A(start:end,5) == M_chi);


%%
hold on;
%plot(A(6:701,1), A(6:701,3)./1e4,'.');
%plot([2.269 2.269],[0 0.01])
sta=30;
fin=181;
x_s = A120(sta:fin,1);
y_s =  ((A120(sta:fin,4)./1e4));
plot(x_s, y_s,'.')

%% derivative
d_y_s = zeros(size(y_s));
for tt=1:1:(length(y_s)-1)
    d_y_s(tt) = (y_s(tt+1)-y_s(tt))/(x_s(tt+1)-x_s(tt));
end
plot(x_s, d_y_s,'-.')
    
    
%hist(A(:,2))
%%
%xx=1.0:0.005:4;
f = (abs(+2.269-x_s)).^(1/8);
plot(x_s, -f)



%%


% %%
% file40 = 'l40_0_005.txt';
% file60 = 'l60_0_005.txt';
% file80 = 'l80_0_005.txt';
% file100 = 'l100_0_005.txt';
% file120 = 'l120_0_005.txt';
% file140 = 'l140_0_005.txt';
% file160 = 'l160_0_005.txt';
% file180 = 'l180_0_005.txt';
% file200 = 'l200_0_005.txt';
% 
% A40 = importdata(file40, ',');
% A60 = importdata(file60, ',');
% A80 = importdata(file80, ',');
% A100 = importdata(file100, ',');
% A120 = importdata(file120, ',');
% A140 = importdata(file140, ',');
% A160 = importdata(file160, ',');
% A180 = importdata(file180, ',');
% A200 = importdata(file200, ',');
