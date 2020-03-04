%This code is written to simulate the K-alpha X-ray emission from Ar plasma

%Electron density maps provided by Dr Galata at LNL

%The Ar partial density is taken as constant, and reaction rate is 
%calculated using electron partial density. Distribution function for 
%electron partial density assumed as the sum of 2 separate Maxwell 
%distributions, one at low temperature (kT<1) and other at at higher 

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

%Calculate grand average temperature from 2-inf keV

AvTemp=zeros(59,59,211);
ParDensTot=zeros(59,59,211);

for i=2:7
    
    ParDensTot(:,:,:)=ParDensTot(:,:,:)+DensityMat(:,:,:,i);

end

for i=2:7
    
    AvTemp(:,:,:)=AvTemp(:,:,:)+(EnerDensMat(:,:,:,i).*DensityMat(:,:,:,i));
    
end

AvTemp(:,:,:)=AvTemp(:,:,:)./(ParDensTot(:,:,:).*1.5);                %Average temperature in terms of kT
AvTemp(isnan(AvTemp))=2000;

%Obtain coeffiecients of grand distribution functions

DensTot=sum(DensityMat,4);
LogicCoeff=zeros(59,59,211,2);
ActualCoeff=zeros(59,59,211,2);

LogicCoeff(:,:,:,1)=DensityMat(:,:,:,1)./DensityMat(:,:,:,1);
LogicCoeff(:,:,:,2)=ParDensTot(:,:,:)./ParDensTot(:,:,:);
LogicCoeff(isnan(LogicCoeff))=0;

ActualCoeff(:,:,:,1)=LogicCoeff(:,:,:,1).*DensityMat(:,:,:,1)./DensTot(:,:,:);
ActualCoeff(:,:,:,2)=LogicCoeff(:,:,:,2).*ParDensTot(:,:,:)./DensTot(:,:,:);
ActualCoeff(isnan(ActualCoeff))=0;


%Define cross section matrix and calculate reaction rate
    
CrSection1=zeros(30,1);
CrSection2=fx.*(10.^-22);
CrSection=[CrSection1;CrSection2];

ReactionRate=0.25.*(DensTot(:,:,:).^2).*calarea_v2_original(ActualCoeff,CrSection,AvTemp,EnerDensMat(:,:,:,1)./1.5);


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
        
