%This code is written to simulate the K-alpha X-ray emission from Ar plasma

%Electron density maps provided by Dr Galata at LNL

%The Ar partial density is taken as constant, and reaction rate is 
%calculated using electron partial density. Distribution function for 
%electron partial density assumed as the sum of 3 separate Maxwell 
%distributions, one at extreme low temperature (kT<1), one at warmer 
%temperature (kT~1 keV) and the final at hot temperature (kT~10 keV) 

%Temperatures determined from Electron Energy Density maps provided by 
%Dr Galata at LNL

%Load all matrices

DensityMat=zeros(59,59,211,7);
EnerDensMat=zeros(59,59,211,7);

for i=1:7
    
    currentFile=sprintf('Dens%d_12_84_step3.mat',i);
    matdata=load(currentFile);
    DensityMat(:,:,:,i)=matdata.(sprintf('Densita%d',i));
    
end

for i=1:7
    
    currentFile=sprintf('Densita_energia%d_12_84_step3_coll.mat',i);
    matdata=load(currentFile);
    EnerDensMat(:,:,:,i)=matdata.(sprintf('Dist%d',i));
    
end
EnerDensMat(:,:,:,:)=EnerDensMat(:,:,:,:)./1000;


%Normalize to actual density in plasma

BackgroundDens=(10^5);     %in mm^-3
Scale=1;
ElectronDens=(10^6-BackgroundDens)*Scale+BackgroundDens;    %in mm^-3
ElTot=ElectronDens*(59*59*211);       %Expected total number of electrons in plasma

ElSimTot=sum(DensityMat,'all');      %Simulated total number of electrons in plasma 

Conv=ElTot/ElSimTot;
DensityMat(:,:,:,:)=DensityMat(:,:,:,:).*Conv;  %Normalised and converted to mm^-3


%Define Parameters of Maxwellian distributions

DensTot=sum(DensityMat,4);
ParDensTot(:,:,:)=zeros(59,59,211);
for i=2:6
    ParDensTot(:,:,:)=ParDensTot(:,:,:)+DensityMat(:,:,:,i);
end
NormCoeff=zeros(59,59,211,3);
NormCoeff(:,:,:,1)=DensityMat(:,:,:,1)./DensTot(:,:,:);
NormCoeff(:,:,:,2)=ParDensTot(:,:,:)./DensTot(:,:,:);
NormCoeff(:,:,:,3)=DensityMat(:,:,:,7)./DensTot(:,:,:);
NormCoeff(isnan(NormCoeff))=0;

Temperature=zeros(59,59,211,3);
Temperature(:,:,:,1)=sum(EnerDensMat.*DensityMat,4)./(DensTot(:,:,:).*1.5);
for i=2:5
    Temperature(:,:,:,2)=Temperature(:,:,:,2)+EnerDensMat(:,:,:,i).*DensityMat(:,:,:,i);
end
Temperature(:,:,:,2)=0.4.*Temperature(:,:,:,2)./((ParDensTot(:,:,:)-DensityMat(:,:,:,6)).*1.5);
Temperature(:,:,:,3)=1.2.*(EnerDensMat(:,:,:,7)./EnerDensMat(:,:,:,7)).*EnerDensMat(:,:,:,7)./1.5;
Temperature(isnan(Temperature))=2000;

%Define cross section matrix and calculate reaction rate
    
CrSection1=zeros(30,1);
CrSection2=fx.*(10.^-22);
CrSection=[CrSection1;CrSection2];

ReactionRate=0.25.*(DensTot(:,:,:).^2).*calarea_v5(NormCoeff,CrSection,Temperature);


%Plot images

figure;
for i=1:6
    subplot(2,3,i);
    imagesc(ReactionRate(:,:,25*(i+1)));
    title(sprintf("Z = %d",25*(i+1)));
end
                
TotalEmission=sum(ReactionRate,3);
        
figure;
imagesc(TotalEmission(:,:));
title('\bfFront view of Plasma Chamber')
        
