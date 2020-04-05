%Code to study validity of 2-Maxwell population approximation
%Overlap correction added
%Druyvesteyn distribution used

%% Load all matrices

DensityMat=zeros(59,59,211,7);             %To hold densities of all cells in different energy intervals 
EnerDensMat=zeros(59,59,211,7);            %To hold average energies of all cells in different energy intervals

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
EnerDensMat(:,:,:,:)=EnerDensMat(:,:,:,:)./1000;     %Average energies in the files are in eV, need to convert into keV


%Normalize to actual density in plasma

BackgroundDens=(10^5);                                      %in mm^-3
Scale=1;
ElectronDens=(10^7-BackgroundDens)*Scale+BackgroundDens;    %in mm^-3
ElTot=ElectronDens*(59*59*211);                             %Expected total number of electrons in plasma

ElSimTot=sum(DensityMat,'all');                             %Simulated total number of electrons in plasma 

Conv=ElTot/ElSimTot;
DensityMat(:,:,:,:)=DensityMat(:,:,:,:).*Conv;              %Normalised and converted to mm^-3


%% Define Parameters of Maxwellian distributions

DensTot=sum(DensityMat,4);                                  %Density summed over all energy intervals for each cell
ParDensTot(:,:,:)=zeros(59,59,211);                         %Density summed over energy intervals from 2-inf for each cell
for i=2:7
    ParDensTot(:,:,:)=ParDensTot(:,:,:)+DensityMat(:,:,:,i);
end

Temperature=zeros(59,59,211,2);                                              %Define cold and warm electron temperatures
Temperature(:,:,:,1)=sum(EnerDensMat.*DensityMat,4)./(DensTot(:,:,:).*1.5);  %Cold electron temperature
for i=2:7
    Temperature(:,:,:,2)=Temperature(:,:,:,2)+EnerDensMat(:,:,:,i).*DensityMat(:,:,:,i);
end
Temperature(:,:,:,2)=0.85.*Temperature(:,:,:,2)./ParDensTot(:,:,:);         %Warm electron temperature
Temperature(isnan(Temperature))=2000;

NormCoeff=zeros(59,59,211,2);                               %Define normalisation coefficients
NormCoeff(:,:,:,1)=(DensityMat(:,:,:,1)-ParDensTot(:,:,:)./(1.513.*(Temperature(:,:,:,2).^0.816)-1))./DensTot(:,:,:); 
%Normalisation coefficient for cold electron population with overlap
NormCoeff(:,:,:,2)=(ParDensTot(:,:,:)./(1-0.661.*(Temperature(:,:,:,2).^-0.816)))./DensTot(:,:,:); 
%Normalisation coefficient for warm electron population with overlap
NormCoeff(isnan(NormCoeff))=0;


%% Calculate approximated density and average energy

DensApprox=zeros(59,59,211,7);                              %Numerically approximated density
EnerDensApprox=zeros(59,59,211,7);                          %Numerically approximated average energy
normint=zeros(59,59,211,7);                                 %Area under the distribution function

energy=0:0.001:35;                                            %Energy range to be considered (inf truncated at 35 keV)
j=1;

fnMax= @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^(3./2))).*exp(-y./x);        %Basic Maxwell distribution function (for 0-2 keV)
fnDr= @(x,y) 1.04.*(sqrt(y)./(x.^(3./2))).*exp(-0.55.*(y.^2)./(x.^2));  %Druyvesteyn distribution in 2-inf keV

for i=1:7
    
    if (i<7)
        while((j*0.001)<(2*i))                                      %calculate area under distribution function
            normint(:,:,:,i)=normint(:,:,:,i)+0.5.*0.001.*((NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j)))...
                                                        +(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j+1))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j+1))));
            j=j+1;
        end
    else
        while ((j*0.001)<35)
            normint(:,:,:,i)=normint(:,:,:,i)+0.5.*0.001.*((NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j)))...
                                                        +(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j+1))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j+1))));
            j=j+1;
        end
    end
    
end

j=1;
for i=1:7
              
    if (i<7)
        while((j*0.001)<(2*i))
            EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)+0.5.*0.001.*(energy(j).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j))...
                                                                                  +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j)))...
                                                                       +energy(j+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j+1))...
                                                                                  +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j+1))));
            j=j+1;
        end
    else
        while ((j*0.001)<35)
            EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)+0.5.*0.001.*(energy(j).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j))...
                                                                                  +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j)))...
                                                                       +energy(j+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),energy(j+1))...
                                                                                  +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),energy(j+1))));
            j=j+1;
        end
    end
        
     DensApprox(:,:,:,i)=DensTot(:,:,:).*normint(:,:,:,i);                 %Calculate densities
     EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)./normint(:,:,:,i);    %Calculate average energies
     EnerDensApprox(isnan(EnerDensApprox))=0;
    
end

%% Analyse goodness of fit

DensDev=abs(DensityMat-DensApprox)./DensityMat;
EnerDensDev=abs(EnerDensMat-EnerDensApprox)./EnerDensApprox;
DensDev(isnan(DensDev))=0;
EnerDensDev(isnan(EnerDensDev))=0;





