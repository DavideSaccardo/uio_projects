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

%% 
filename = 'positions_E102.dat';
A=importdata(filename);

%%
filename = 'positions_N102.dat';
B=importdata(filename);


%% plot
hold on
plot(B(:,4),B(:,5),'-'); %newton
plot(A(:,4),A(:,5),'-'); %einstein
%plot(A(:,2),A(:,3),'.'); %einstein

%%
r_N = sqrt(B(:,5).*B(:,5) + B(:,6).*B(:,6) + B(:,7).*B(:,7));
r_E = sqrt(A(:,5).*A(:,5) + A(:,6).*A(:,6) + A(:,7).*A(:,7));


%% minimo

[Y_E,I_E] = min(r_E);
[Y_N,I_N] = min(r_N);
%

hold on;
%xlim([A(1,1) A(2000,1)])
plot(B(:,1),r_N,'r.');
plot(A(:,1),r_E,'b.');
plot([A(I_E,1) A(I_E,1)],[0.3 0.4])
plot([A(I_N,1) A(I_N,1)],[0.3 0.4])

ang_E = atan2d(A(I_E,6), A(I_E,5));
ang_N = atan2d(B(I_N,6), B(I_N,5));

a_e = atan2d(A(:,6), A(:,5));
a_n = atan2d(B(:,6), B(:,5));
d = a_e-a_n;
table(0.0119,(ang_E-ang_N))
% period of r is approx 0.3462

%%
x1 = A(:,4);
y1 = A(:,5);
z1 = A(:,6);


%%
dim = 1;

for ii=1:length(x1)
    
    
    plot3(x1(ii),y1(ii),z1(ii),'or','MarkerSize',1,'MarkerFaceColor','r')
    hold on;
    grid on;

    %plot3(x2(ii),y2(ii),z2(ii),'ob','MarkerSize',5,'MarkerFaceColor','b')
    %plot3(x3(ii),y3(ii),z3(ii),'ob','MarkerSize',5,'MarkerFaceColor','c')
    xlim([-dim dim]);
    ylim([-dim dim]);
    zlim([-dim dim]);
    pause(1e-11)
    %clf
end


%% interesting case
x = linspace(-6,6,1000);
y = sin(x);
plot(x,y)
axis manual
ax = gca;
h1 = hgtransform('Parent',ax);
hold on
plot(x(1),y(1),'o','Parent',h1);
hold off
t = text(x(1),y(1),num2str(y(1)),'Parent',h1,...
    'VerticalAlignment','top','FontSize',14);
for k = 2:length(x)
    m = makehgtform('translate',x(k)-x(1),y(k)-y(1),0);
    h1.Matrix = m;
    t.String = num2str(y(k));
    drawnow
end