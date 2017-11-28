%% three set of data x variable, s analytical solution, v numerical solution
% first set of n = 10
filename = 'A';
[A,delimiterOut]=importdata(filename);

xA = A.data(:,1);
sA = A.data(:,2);
vA = A.data(:,3);
% second set n=100
filename = 'B';
[B,delimiterOut]=importdata(filename);

xB = B.data(:,1);
sB = B.data(:,2);
vB = B.data(:,3);
% third set n=1000
filename = 'C';
[C,delimiterOut]=importdata(filename);

xC = C.data(:,1);
sC = C.data(:,2);
vC = C.data(:,3);
%% plot of the solution
hold on;
plot(xC,sC,'b.');
plot(xC,vC,'c.');
plot(xB,vB,'r.');
plot(xA,vA,'m.');
%% point d but the question isn't so clear
epsA = log10(abs( (sA - vA)./sA ) );
epsB = log10(abs( (sB - vB)./sB ) );
epsC = log10(abs( (sC - vC)./sC ) );
%log of h
lhA = log10(xA);
lhB = log10(xB);
lhC = log10(xC);

hold on;
plot(lhA,epsA);
plot(lhB,epsB);
plot(lhC,epsC);
%
