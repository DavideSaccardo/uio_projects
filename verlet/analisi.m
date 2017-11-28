%% solution
filename = 'data.txt';
A=importdata(filename, ',');

t = A(:,1);
x = A(:,2);
y = A(:,3);

plot(A(:,3), A(:,2), '-.');
%%
h = animatedline;

pause on
for k = 1:length(x)
    
    addpoints(h,x,y);
    drawnow
    pause(1)
end
pause off
%%
h = animatedline;
axis([0 4*pi -1 1])
x = linspace(0,4*pi,10000);
y = sin(x);

for k = 1:length(x)
    addpoints(h,x(k),y(k));
    drawnow
end
%%
f = figure;
w = waitforbuttonpress;
if w == 0
    disp('Button click')
else
    disp('Key press')
end