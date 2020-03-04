%This code is written to simulate the K-alpha X-Ray emission from Ar plasma

%Electron density maps calculated by Dr Galata at LNL

%Case 1: The Ar partial density is taken as constant, and reaction rate is 
%calculated using electron partial density
%Case 2: The Ar partial density is taken equal to electron partial density


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
EnerDensMat(:,:,:,:)=EnerDensMat(:,:,:,:).*1000;

%Normalize to actual density in plasma

BackgroundDens=(10^5);     %in mm^-3
Scale=1;
ElectronDens=(10^9-BackgroundDens)*Scale+BackgroundDens;    %in mm^-3
ElTot=sum(ElectronDens,'all');       %Expected total number of electrons in plasma

ElSimTot=sum(DensityMat,'all');      %Simulated total number of electrons in plasma 

Conv=ElTot/ElSimTot;
DensityMat(:,:,:,:)=DensityMat(:,:,:,:).*Conv;  %Normalised and converted to mm^-3


%Define K-alpha cross section matrix and corresponding energy

CrEnergy=[0 2 4 6 8 10 12 20];               
c=3*(10^11);
vel=c.*sqrt(2.*CrEnergy./511);
CrSection(1)=0;
CrSection(2)=0;
CrSection(3)=fx1(1)*(10^-22);
CrSection(4)=fx1(2)*(10^-22);
CrSection(5)=fx1(3)*(10^-22);
CrSection(6)=fx1(4)*(10^-22);
CrSection(7)=fx1(5)*(10^-22);
CrSection(8)=fx1(9)*(10^-22);


%Initiate main program sequence

choice=1;                           %Choice for type of approach to be taken

switch choice

    case 1
       
        %Calculate Ar ion density (taken as constant in this case)
        
        ArDensity=0.25.*sum(DensityMat,4);
        
        %Calculate reaction rate in each plasma cell (consider modification
        %of reaction rate using Trapezoidal area calculation)
        
        ReactionRate=zeros(59,59,211);
        
        for i=1:7
            ReactionRate(:,:,:)=ReactionRate(:,:,:)+0.5.*(ArDensity(:,:,:).*DensityMat(:,:,:,i).*(CrSection(i).*vel(i)+CrSection(i+1).*vel(i+1)));
        end
        
        figure;
        for i=1:6
            subplot(2,3,i);
            imagesc(ReactionRate(:,:,25*(i+1)));
            title(sprintf("Z = %d",25*(i+1)));
        end
        
        %Front view of X-ray emission in plasma chamber
        
        TotalEmission=zeros(59,59);
        
        for i=1:211
            TotalEmission(:,:)=TotalEmission+ReactionRate(:,:,i);
        end
        
        figure;
        imagesc(TotalEmission(:,:));
        title('\bfFront view of Plasma Chamber')
        
    case 2
        
        %Calculate reaction rate in each plasma cell
        
        ReactionRate=zeros(59,59,211);
        
        for i=1:7
            
            if (i<7)
                ReactionRate(:,:,:)=ReactionRate(:,:,:)+(0.25.*(DensityMat(:,:,:,i).^2).*CrSection(i).*CrEnergy(i).*2);
            else
                ReactionRate(:,:,:)=ReactionRate(:,:,:)+(0.25.*(DensityMat(:,:,:,i).^2).*CrSection(i).*CrEnergy(i).*8);
            end
        
        end
        
        figure;
        for i=1:6
            subplot(2,3,i);
            imagesc(ReactionRate(:,:,25*(i+1)));
            title(sprintf("Z = %d",25*(i+1)));
        end
        
        %Front view of X-ray emission in plasma chamber
        
        TotalEmission=zeros(59,59);
        
        for i=1:211
            TotalEmission(:,:)=TotalEmission+ReactionRate(:,:,i);
        end
        
        figure;
        imagesc(TotalEmission(:,:));
        title('\bfFront view of Plasma Chamber')
        
end


