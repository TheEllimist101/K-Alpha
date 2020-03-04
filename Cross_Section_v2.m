%This code is written to simulate the K-alpha X-Ray emission from Ar plasma

%Electron density maps provided by Dr Galata at LNL

%The Ar partial density is taken as constant, and reaction rate is 
%calculated using electron partial density. Distribution function for 
%electron partial density assumed as the sum of 14 separate Maxwell 
%distributions, with temperature ranging from kT=2 keV to kT=28 keV, in
%steps of 2

%Obtain distribution functions

%Define temperature matrix 

kT=linspace(2,28,14);             %Maxwellian temperatures from 2 to 28 keV

%Obtain coefficient matrix for set of 14 equations

coeffmat=zeros(14,14);
for i=1:14
    for j=1:14
        coeffmat(i,j)=calarea(kT(j),i,1);
    end
end

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

%Define modified densities for constant matrix

DensTot=zeros(59,59,211);

DensTot(:,:,:)=sum(DensityMat,4);
DensityMat=DensityMat./DensTot;

%Define and assign values to cell arrays

Coeff_array=cell(59,59,211);
Const_array=cell(59,59,211);
Var_array=cell(59,59,211);

Coeff_array(:,:,:)={coeffmat};
Var_array(:,:,:)={zeros(14,1)};
Const_array(:,:,:)={zeros(14,1)};

for i=1:59
    for j=1:59
        for k=1:211
            Const_array{i,j,k}=[DensityMat(i,j,k,1);DensityMat(i,j,k,2);DensityMat(i,j,k,3);DensityMat(i,j,k,4);DensityMat(i,j,k,5);...
                                DensityMat(i,j,k,6);DensityMat(i,j,k,7);EnerDensMat(i,j,k,1);EnerDensMat(i,j,k,2);EnerDensMat(i,j,k,3);...
                                EnerDensMat(i,j,k,4);EnerDensMat(i,j,k,5);EnerDensMat(i,j,k,6);EnerDensMat(i,j,k,7)];
        end
    end
end

% Const_array{:,:,:}=[DensityMat(:,:,:,1);DensityMat(:,:,:,2);DensityMat(:,:,:,3);DensityMat(:,:,:,4);DensityMat(:,:,:,5);...
%                                 DensityMat(:,:,:,6);DensityMat(:,:,:,7);EnerDensMat(:,:,:,1);EnerDensMat(:,:,:,2);EnerDensMat(:,:,:,3);...
%                                  EnerDensMat(:,:,:,4);EnerDensMat(:,:,:,5);EnerDensMat(:,:,:,6);EnerDensMat(:,:,:,7)];

%Solve for every element of Var_array

for i=1:59
    for j=1:59
        for k=1:211
            Var_array{i,j,k}=linsolve(Coeff_array{i,j,k},Const_array{i,j,k});
        end
    end
end
Var_array{isnan{Var_array}}=0;    %Replace all NaN appearances with 0 to signify all Ai's are 0

%Create new density matrices (at intervals of 0.1 from 0 to 20)

GrandCell=cell(200,1);
TempDensity=zeros(59,59,211);

for i=1:200
    
    for m=1:59
        for n=1:59
            for p=1:211
                TempArray=Var_array{m,n,p};
                for q=1:length(TempArray)
                    TempDensity(m,n,p)=TempDensity(m,n,p)+DensTot(m,n,p)*TempArray(q)*calarea(kT(q),i*0.1,2); 
                end
            end
        end
    end
    GrandCell(i)={TempDensity};
    
end

%Calculate cross section matrix in new intervals and final reaction rate













