clc
clear all
close all

load initSig.txt
load p90.txt
load fid.txt
load p180.txt
load spinEcho.txt

load initSigTime.txt
load p90Time.txt
p90Time=p90Time+initSigTime(end);
load fidTime.txt
fidTime=fidTime+p90Time(end);
load p180Time.txt
p180Time=p180Time+fidTime(end);
load spinEchoTime.txt
spinEchoTime=spinEchoTime+p180Time(end);

gamma=1000;
B0=1;
omega=gamma*B0;

%for the first pulse: frame change
dt=p90Time(2)-p90Time(1);
t=dt*[0:(length(p90(:,1))-1)]';
V1=cos(omega*t);
V2=sin(omega*t);

%for the first pulse: frame change
dt2=p180Time(2)-p180Time(1);
t2=dt2*[0:(length(p180(:,1))-1)]';
V12=cos(omega*t2);
V22=sin(omega*t2);

MOM_MAG = [initSig; p90; fid; p180; spinEcho];
TIME = [initSigTime; p90Time; fidTime; p180Time; spinEchoTime];

%x
figure(1)
subplot(231)
plot(TIME,MOM_MAG(:,1))
subplot(232)
plot(initSigTime,initSig(:,1))
subplot(233)
plot(p90Time,p90(:,1))
subplot(234)
plot(fidTime,fid(:,1))
subplot(235)
plot(p180Time,p180(:,1))
subplot(236)
plot(spinEchoTime,spinEcho(:,1))

%y
figure(2)
subplot(231)
plot(TIME,MOM_MAG(:,2))
subplot(232)
plot(initSigTime,initSig(:,2))
subplot(233)
plot(p90Time,p90(:,2))
subplot(234)
plot(fidTime,fid(:,2))
subplot(235)
plot(p180Time,p180(:,2))
subplot(236)
plot(spinEchoTime,spinEcho(:,2))


%z
figure(3)
subplot(231)
plot(TIME,MOM_MAG(:,3))
subplot(232)
plot(initSigTime,initSig(:,3))
subplot(233)
plot(p90Time,p90(:,3))
subplot(234)
plot(fidTime,fid(:,3))
subplot(235)
plot(p180Time,p180(:,3))
subplot(236)
plot(spinEchoTime,spinEcho(:,3))



figure(4)
subplot(231)
plot(MOM_MAG(:,1),MOM_MAG(:,2))
subplot(232)
plot(initSig(:,1),initSig(:,2))
subplot(233)
plot(p90(:,1).*V1+p90(:,2).*V2,p90(:,2).*V2-p90(:,1).*V1)
subplot(234)
plot(fid(:,1),fid(:,2))
subplot(235)
plot(p180(:,1).*V12+p180(:,2).*V22,p180(:,2).*V22-p180(:,1).*V12)
subplot(236)
plot(spinEcho(:,1),spinEcho(:,2))


figure(5)
plot(spinEchoTime,spinEcho(:,1))
figure(6)
plot(spinEchoTime,spinEcho(:,2))
figure(7)
plot(spinEchoTime,spinEcho(:,3))



