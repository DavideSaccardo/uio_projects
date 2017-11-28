file = 'pde.txt';
A = importdata(file, ',');
%%
%plot(0:0.01:1,A(:,100))
figure();

for ii = 1:1:3000000
    hold on
    axis([0 1 0 9])
    plot(0:0.001:1,A(ii,:));

    pause(0.0001)
    clf
end