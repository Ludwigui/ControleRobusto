clc
clear all
close all


x = [];
it = 1;
reta = [];

theta = 0:0.1:2*pi;
rad = 10;
xCenter=0;
yCenter=0;

xCoord = xCenter+rad*cos(theta);
yCoord = yCenter+rad*sin(theta);

figure
grid on
plot(xCoord, yCoord,'k')
hold on
xline(-3,'--r')
hold on
yline(0,'k')
hold on
xline(0,'k')
xlim([-11 0])
ylim([-11 11])
xlabel('Re')
ylabel('Im')
    