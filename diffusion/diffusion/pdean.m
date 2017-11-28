file = 'pde.txt';
A = importdata(file, ',');
%%
%plot(0:0.01:1,A(:,100))
figure();
xx = 0:0.01:1;
for ii = 0:1:max(size(A));
    hold on
    q = 101*ii;
    qq = (ii+1)*101;
    %axis([0 1 0 9])
    surf(xx',xx',A((q+1):(qq),1:101));

    pause(0.0001)
    clf
end