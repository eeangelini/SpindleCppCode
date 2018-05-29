dataDir = '../data/histData/';
dataFile = 'rotationStabilizationWithTrans/2MTRotationStabilizationWithTrans';
run parameters.m
dataFile2 = 'rotationStabilizationWithTrans/20MTRotationStabilizationWithTrans';
dataFile3 = 'rotationStabilizationWithTrans/200MTRotationStabilizationWithTrans';
dataFile4 = 'rotationStabilizationWithTrans/2000MTRotationStabilizationWithTrans';
dataFileFull2 = [dataDir dataFile2 dataFileSuffix];
dataFileFull3 = [dataDir dataFile3 dataFileSuffix];
dataFileFull4 = [dataDir dataFile4 dataFileSuffix];

numRuns = 1000;
csvrange = [rowStart colStart rowStart+numRuns-1 colStart+2];DATA1 = csvread(dataFileFull,rowStart,colStart,csvrange);
DATA2 = csvread(dataFileFull2,rowStart,colStart,csvrange);
DATA3 = csvread(dataFileFull3,rowStart,colStart,csvrange);
DATA4 = csvread(dataFileFull4,rowStart,colStart,csvrange);
psi1 = abs((DATA1(:,3)-pi/2)*180/pi);
psi2 = abs((DATA2(:,3)-pi/2)*180/pi);
psi3 = abs((DATA3(:,3)-pi/2)*180/pi);
psi4 = abs((DATA4(:,3)-pi/2)*180/pi);

figure;
axes('Linewidth', 3.5);
nhist({psi1, psi2, psi3, psi4},'legend',{'2 MTs', '20 MTs', '200 MTs', '2000 MTs'},'box','location','EastOutside','fsize',16,'linewidth',3);
hold on;
StartingAng = plot(abs((startPsi-pi/2)*180/pi), 0, '*');
title('Rotation Stabilization with Translation')