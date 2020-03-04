%This code is written to simulate the K-alpha X-ray emission from Ar plasma

%Electron density maps provided by Dr Galata at LNL

%The Ar partial density is taken as constant, and reaction rate is 
%calculated using electron partial density. Distribution function for 
%electron partial density assumed as a single Maxwell distribution 

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

Temperature=2.*sum(EnerDensMat.*DensityMat,4)./(DensTot(:,:,:).*1.5);
Temperature(isnan(Temperature))=2000;

%Define cross section matrix and calculate reaction rate
    
CrSection1=zeros(30,1);
CrSection2=fx.*(10.^-22);
CrSection=[CrSection1;CrSection2];

ReactionRate=0.25.*(DensTot(:,:,:).^2).*calarea_v3(CrSection,Temperature);


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
