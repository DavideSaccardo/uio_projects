%% omega analysis
%% omega=5
filename = 'rho_val_vec5.000000.txt';
A=importdata(filename, ',');

% eigenfunctin u(rho) = r * R(r)
rhoA = A.data(:,1);
valA = A.data(:,2);
f1A = A.data(:,3);
f2A = A.data(:,4);
f3A = A.data(:,5);
% probablity density function
pdf1A = (f1A).^2./trapz(rhoA,(f1A).^2);
pdf2A = (f2A).^2./trapz(rhoA,(f2A).^2);
pdf3A = (f3A).^2./trapz(rhoA,(f3A).^2);

% plot
hold on;
plot(rhoA,pdf1A,'.');
plot(rhoA,pdf2A,'.');
plot(rhoA,pdf3A,'.');

%% omega=1
filename = 'rho_val_vec1.000000.txt';
B=importdata(filename, ',');

%
rhoB = B.data(:,1);
valB = B.data(:,2);
f1B = B.data(:,3);
f2B = B.data(:,4);
f3B = B.data(:,5);
% probablity density function
pdf1B = (f1B).^2./trapz(rhoB,(f1B).^2);
pdf2B = (f2B).^2./trapz(rhoB,(f2B).^2);
pdf3B = (f3B).^2./trapz(rhoB,(f3B).^2);
% plot
hold on;
plot(rhoB,pdf1B,'.');
plot(rhoB,pdf2B,'.');
plot(rhoB,pdf3B,'.');

%% omega = 0.5
filename = 'rho_val_vec0.500000.txt';
C=importdata(filename, ',');

%
rhoC = C.data(:,1);
valC = C.data(:,2);
f1C = C.data(:,3);
f2C = C.data(:,4);
f3C = C.data(:,5);
% probablity density function
pdf1C = (f1C).^2./trapz(rhoC,(f1C).^2);
pdf2C = (f2C).^2./trapz(rhoC,(f2C).^2);
pdf3C = (f3C).^2./trapz(rhoC,(f3C).^2);
% plot
hold on;
plot(rhoC,pdf1C,'.');
plot(rhoC,pdf2C,'.');
plot(rhoC,pdf3C,'.');
%% omega = 0.01
filename = 'rho_val_vec0.010000.txt';
D=importdata(filename, ',');

%
rhoD = D.data(:,1);
valD = D.data(:,2);
f1D = D.data(:,3);
f2D = D.data(:,4);
f3D = D.data(:,5);
% probablity density function
pdf1D = (f1D).^2./trapz(rhoD,(f1D).^2);
pdf2D = (f2D).^2./trapz(rhoD,(f2D).^2);
pdf3D = (f3D).^2./trapz(rhoD,(f3D).^2);
% plot
hold on;
plot(rhoD,pdf1D,'.');
plot(rhoD,pdf2D,'.');
plot(rhoD,pdf3D,'.');
title({'Interacting solutions';'Probability density functions of the first three states';'with \omega_r=0.01'},'Interpreter','teX','Fontsize',18);
xlabel('$\rho$','interpreter','latex','fontsize',18);
ylabel('$|\psi(\rho)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times');
grid on;
xlim([0 50]);
plot([19.3 19.3],[0 0.08],'m');
ll=legend(gca,'show','ground state $\psi_0$','first state $\psi_1$','second state $\psi_2$','$\rho_{node}=19.3$ of $\psi_1$','location','northeast');
set(ll,'Interpreter','latex');

%% armonic oscillator 
filename = 'rho_val_vec_noninter.txt';
E=importdata(filename, ',');

%
rhoE = E.data(:,1);
valE = E.data(:,2);
f1E = E.data(:,3);
f2E = E.data(:,4);
f3E = E.data(:,5);
f4E = E.data(:,6);
% probablity density function
pdf1E = (f1E).^2./trapz(rhoE,(f1E).^2);
pdf2E = (f2E).^2./trapz(rhoE,(f2E).^2);
pdf3E = (f3E).^2./trapz(rhoE,(f3E).^2);
pdf4E = (f4E).^2./trapz(rhoE,(f4E).^2);

%%plot
hold on;
plot(rhoE,pdf1E,'r.');
plot(rhoE,pdf2E,'g.');
plot(rhoE,pdf3E,'k.');
plot(rhoE,pdf4E,'m.');

%% analytic solutions lambda = 2n + 3
hold on
xx = rhoE;
n= [1 2 3];
psi1 = xx.*exp(-xx.^2 ./2).* laguerreL(0,1/2, xx.^2);
pdf_psi1 = psi1.^2 ./trapz(xx,(psi1).^2);
psi2 = xx.*exp(-xx.^2 ./2).* laguerreL(1,1/2, xx.^2);
pdf_psi2 = psi2.^2 ./trapz(xx,(psi2).^2);
psi3 = xx.*exp(-xx.^2 ./2).* laguerreL(2,1/2, xx.^2);
pdf_psi3 = psi3.^2 ./trapz(xx,(psi3).^2);
psi4 = xx.*exp(-xx.^2 ./2).* laguerreL(3,1/2, xx.^2);
pdf_psi4 = psi4.^2 ./trapz(xx,(psi4).^2);
plot(xx, pdf_psi1);
plot(xx, pdf_psi2);
plot(xx, pdf_psi3);
plot(xx, pdf_psi4);
%% error analysis
hold on
eps1 = (abs((psi1-pdf1E)./psi1));
plot(xx, eps1);
%
eps2 = (abs((psi2-pdf2E)./psi2));
plot(xx, eps2);

eps3 = (abs((psi3-pdf3E)./psi3));
plot(xx, eps3);

eps4 = (abs((psi4-pdf4E)./psi4));
plot(xx, eps4);
%% three value fo omega
%% ground state
hold on
plot(rhoA,pdf1A,'.');
plot(rhoB,pdf1B,'.');
plot(rhoC,pdf1C,'.');
%plot(rhoD,pdf1D,'.');
title({'Interacting solutions';'Probability density functions of the ground state';'with three values of \omega_r'},'Interpreter','teX','Fontsize',18);
xlabel('$\rho$','interpreter','latex','fontsize',18);
ylabel('$|\psi(\rho)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','tex');
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times');
grid on;
xlim([0 4.5]);
ll=legend(gca,'show','$\omega_r$=5','$\omega_r$=1','$\omega_r$=0.5','location','northeast');
set(ll,'Interpreter','Latex');
%print('3valueomega', '-depsc');
%% first state
hold on
plot(rhoA,pdf2A,'.');
plot(rhoB,pdf2B,'.');
plot(rhoC,pdf2C,'.');
%plot(rhoD,pdf2D,'.');
title({'Interacting solutions';'Probability density functions of the first state';'with four values of \omega_r'},'Interpreter','teX','Fontsize',18);
xlabel('$\rho$','interpreter','latex','fontsize',18);
ylabel('$|u(\rho)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times');
grid on;
xlim([0 5.5]);
ll=legend(gca,'show','$\omega_r$=5','$\omega_r$=1','$\omega_r$=0.5','location','northeast');
set(ll,'Interpreter','Latex');
%print('3valueomega2', '-depsc');
%% interacting and non interacting case

%constants
beta_e = 1.44;%eV nm beta * e^2
h_bar = 6.582e-16; %eV s
m = 9.11e-31; %kg electron mass
%EE =  (h_bar^2 /(m*(alfa)^2) );

%omega quadro per m dell'oscillatore armonico in funzione di quello interagente
omega_r = [0.01 0.5 1 5];
alfa_i = h_bar^2 /(m*beta_e); %interacting
m_k = 4*h_bar^2 .*omega_r.^2 ./(alfa_i.^4);
alfa_o = (h_bar^2 ./(m_k)).^(1/4);

r_iA = rhoA.*alfa_i;
r_iB = rhoB.*alfa_i;
r_iC = rhoC.*alfa_i;
r_iD = rhoD.*alfa_i;
r_oA = (rhoE.*alfa_o(1,4));
r_oB = (rhoE.*alfa_o(1,3));
r_oC = (rhoE.*alfa_o(1,2));
r_oD = (rhoE.*alfa_o(1,1));

%% plots 5
plot(r_iA,pdf1A./trapz(r_iA,pdf1A),'.'); xlim([0 0.9]); ylim([0 6]); ...
    hold on; plot(r_oA,pdf1E./trapz(r_oA,pdf1E),'.'); 
title({'Probability density functions of the ground state';'with and without interactions'},'Interpreter','LateX','Fontsize',18);
xlabel('$r$ [nm]','interpreter','latex','fontsize',18);
ylabel('$|u(r)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontName', 'Times');
set(gca,'FontSize',18);
grid on;
xlim([0 0.5]);
ll=legend(gca,'show','with Coulomb potential $\omega_r = 5$','without Coulomb potential','location','northeast');
set(ll,'Interpreter','Latex');
%% plots 1
plot(r_iB,pdf1B./trapz(r_iB,pdf1B),'k.'); xlim([0 0.9]); ylim([0 6]); ...
    hold on; plot(r_oB,pdf1E./trapz(r_oB,pdf1E),'g.'); 
title({'Probability density functions of the ground state';'with and without interactions'},'Interpreter','LateX','Fontsize',18);
xlabel('$r$ [nm]','interpreter','latex','fontsize',18);
ylabel('$|u(r)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca, 'FontName', 'Times');
set(gca,'FontSize',18);
grid on;
xlim([0 1]); ylim([0 4]);
ll=legend(gca,'show','with Coulomb potential $\omega_r = 1$','without Coulomb potential','location','northeast');
set(ll,'Interpreter','Latex');
%% 0.5
plot(r_iC,pdf1C./trapz(r_iC,pdf1C),'c.'); xlim([0 0.9]); ylim([0 6]);  ...
    hold on; plot(r_oC,pdf1E./trapz(r_oC,pdf1E),'m.');
title({'Probability density functions of the ground state';'with and without interactions'},'Interpreter','LateX','Fontsize',18);
xlabel('$r$ [nm]','interpreter','latex','fontsize',18);
ylabel('$|u(r)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times');
grid on;
xlim([0 1]); ylim([0 3.2])
ll=legend(gca,'show','with Coulomb potential $\omega_r = 0.5$','without Coulomb potential','location','northeast');
set(ll,'Interpreter','Latex');
%% omega =0.01
plot(r_iD,pdf1D./trapz(r_iD,pdf1D),'r.'); ...
    hold on; plot(r_oD,pdf1E./trapz(r_oD,pdf1E),'b.'); 
title({'Probability density functions of the ground state';'with and without interactions'},'Interpreter','LateX','Fontsize',18);
xlabel('$r$ [nm]','interpreter','latex','fontsize',18);
ylabel('$|u(r)|^2$','interpreter','latex','fontsize',18);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',18);
set(gca, 'FontName', 'Times');
grid on;
xlim([0 12]); ylim([0 1.3])

ll=legend(gca,'show','with Coulomb potential $\omega_r = 0.01$','without Coulomb potential','location','northeast');
set(ll,'Interpreter','Latex');
%% A



subplot(2,2,1); plot(r_iA,pdf1A./trapz(r_iA,pdf1A),'r.'); xlim([0 0.5]); ylim([0 8.5]); ...
    hold on; plot(r_oA,pdf1E./trapz(r_oA,pdf1E),'b.'); 
title({'$\omega_r = 5$'},'Interpreter','LateX','Fontsize',14);

xlabel('$r$ [nm]','interpreter','latex','fontsize',12);
ylabel('$|\psi(r)|^2$','interpreter','latex','fontsize',12);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times');
grid on;


%ll=legend(gca,'show','with Coulomb potential $\omega_r = 5$','without Coulomb potential','location','northeast');
%set(ll,'Interpreter','Latex');
subplot(2,2,2); plot(r_iB,pdf1B./trapz(r_iB,pdf1B),'r.'); xlim([0 0.9]); ylim([0 4]); ...
    hold on; plot(r_oB,pdf1E./trapz(r_oB,pdf1E),'b.'); 
title({'$\omega_r = 1$'},'Interpreter','LateX','Fontsize',14);

xlabel('$r$ [nm]','interpreter','latex','fontsize',12);
ylabel('$|\psi(r)|^2$','interpreter','latex','fontsize',12);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times');
grid on;

subplot(2,2,3); plot(r_iC,pdf1C./trapz(r_iC,pdf1C),'r.'); xlim([0 1.5]); ylim([0 3.2]);  ...
    hold on; plot(r_oC,pdf1E./trapz(r_oC,pdf1E),'b.'); 
title({'$\omega_r = 0.5$'},'Interpreter','LateX','Fontsize',14);

xlabel('$r$ [nm]','interpreter','latex','fontsize',12);
ylabel('$|\psi(r)|^2$','interpreter','latex','fontsize',12);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times');
grid on;
subplot(2,2,4); plot(r_iD,pdf1D./trapz(r_iD,pdf1D),'r.'); ...
    hold on; plot(r_oD,pdf1E./trapz(r_oD,pdf1E),'b.'); ylim([0 0.4]);
xlim([0 11])
title({'$\omega_r = 0.01$'},'Interpreter','LateX','Fontsize',14);

xlabel('$r$ [nm]','interpreter','latex','fontsize',12);
ylabel('$|\psi(r)|^2$','interpreter','latex','fontsize',12);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
set(gca, 'FontName', 'Times');
grid on;
suptitle({'Probability density functions of the ground state','with and without Coulomb potential'})
%% uncertainty
m_A = trapz(r_iA,r_iA .*pdf1A./trapz(r_iA,pdf1A))
m_B = trapz(r_iB,r_iB .*pdf1B./trapz(r_iB,pdf1B))
m_C = trapz(r_iC,r_iC .*pdf1C./trapz(r_iC,pdf1C))
m_D = trapz(r_iD,r_iD .*pdf1D./trapz(r_iD,pdf1D))
m_oA = trapz(r_iA,r_oA.*pdf1E./trapz(r_oA,pdf1E))
m_oB = trapz(r_oB,r_oB.*pdf1E./trapz(r_oB,pdf1E))
m_oC = trapz(r_oC,r_oC.*pdf1E./trapz(r_oC,pdf1E))
m_oD = trapz(r_oD,r_oD.*pdf1E./trapz(r_oD,pdf1E))
table([m_A; m_oA],[m_B; m_oB],[m_C; m_oC],[m_D; m_oD])
%% varianza
var_A = trapz(r_iA,(r_iA-m_A).^2 .*pdf1A./trapz(r_iA,pdf1A))
var_B = trapz(r_iB,(r_iB-m_B).^2 .*pdf1B./trapz(r_iB,pdf1B))
var_C = trapz(r_iC,(r_iC-m_C).^2 .*pdf1C./trapz(r_iC,pdf1C))
var_D = trapz(r_iD,(r_iD-m_D).^2 .*pdf1D./trapz(r_iD,pdf1D))
table(var_A, var_B, var_C, var_D)
%%
dA = abs(max(pdf1A./trapz(r_iA,pdf1A)) - max(pdf1E./trapz(r_oA,pdf1E)));
dB = abs(max(pdf1B./trapz(r_iB,pdf1B)) - max(pdf1E./trapz(r_oB,pdf1E)));
dC = abs(max(pdf1C./trapz(r_iC,pdf1C)) - max(pdf1E./trapz(r_oC,pdf1E)));
dD = abs(max(pdf1D./trapz(r_iD,pdf1D)) - max(pdf1E./trapz(r_oD,pdf1E)));
maxx = [dD dC dB dA];
omega_r = [0.01 0.5 1 5];
plot(omega_r', maxx','-.')
table(dD, dC, dB, dA)
%% parameters interacting case
beta_e = 1.44;%eV nm beta * e^2
h_bar = 6.582e-16; %eV s
m = 9.11e-31; %kg electron mass
alfa = h_bar^2 /(m*beta_e);
E =  (h_bar^2 /(m*(alfa)^2) );%.* eigenval;
omega = 0.05;
% r_0 = (2*(omega)^2)^(-1/3);