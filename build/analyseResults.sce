close
clear

initSig= fscanfMat("initSig.txt");
p90= fscanfMat("p90.txt");
fid= fscanfMat("fid.txt");
p180= fscanfMat("p180.txt");
spinEcho= fscanfMat("spinEcho.txt");

initSigTime= fscanfMat("initSigTime.txt");
p90Time= fscanfMat("p90Time.txt");
p90Time=p90Time+initSigTime(length(initSigTime));
fidTime= fscanfMat("fidTime.txt");
fidTime=fidTime+p90Time(length(p90Time));
p180Time= fscanfMat("p180Time.txt");
p180Time=p180Time+fidTime(length(fidTime));
spinEchoTime= fscanfMat("spinEchoTime.txt");
spinEchoTime=spinEchoTime+p180Time(length(p180Time));

gamma=1000;
B0=1;
omega=gamma*B0;

//for the first pulse: frame change
dt=p90Time(2)-p90Time(1);
t=dt*[0:(length(p90(:,1))-1)]';
V1=cos(omega*t);
V2=sin(omega*t);

//for the first pulse: frame change
dt2=p180Time(2)-p180Time(1);
t2=dt2*[0:(length(p180(:,1))-1)]';
V12=cos(omega*t2);

MOM_MAG = [initSig; p90; fid; p180; spinEcho];
TIME = [initSigTime; p90Time; fidTime; p180Time; spinEchoTime];

//x
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

//y
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

//z
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




//Xtrue=fscanfMat("observedImage.png");
//Xtrue = Xtrue(3:(size(Xtrue,1)-2),3:(size(Xtrue,2)-2));

//proced image loading
//Y=fscanfMat("observedImage.png");
//X=fscanfMat("finalImage.png");

//NOISE = [50];//10 20 50 100];
//ITER = [50000];//10000 20000 50000
//LAMBDA = [100 200 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];// 20000 30000 40000];//

//for iter=1:length(ITER)
//    for noise=1:length(NOISE)
//        for l=1:length(LAMBDA)
//            filenameFinal = strcat(["noise" string(NOISE(noise)) "/" string(ITER(iter)) "/finalImage_N_" string(ITER(iter)) "_LAMBDA_" string(LAMBDA(l)) "_NOISE_" string(NOISE(noise)) ".png"]);
//            X=fscanfMat(filenameFinal);
//            X = X(3:(size(X,1)-2),3:(size(X,2)-2));
//            [ITER(iter) NOISE(noise) LAMBDA(l) sum((X(:)-Xtrue(:)).^2)/prod(size(Xtrue))]
//        end
//    end
//end
