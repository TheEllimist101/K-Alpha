%This code will calculate the Local Geometrical Emission of the plasma
%cells by combining aspects of the Monte-Carlo and Perturbative codes

%The area over which the MC points are distributed will be restricted
%around the direction of emission sure to be detected as in Perturbation
%Method

%The code will be run over a large number of MC steps, each corresponding
%to a random direction of emission, but located within the restricted solid
%angle alone

%The LGE will be calculated with respect to the restricted space alone

%Code does not evaluate the photon counting algorithm or intrisic photon to
%charge cloud conversion process in the CCD

%% Define system parameters (all quantities expressed in SI units)

MeshTransparency=0.56;                                 %Mesh transparency at normal incidence (56%)
d0=1.191;                                              %Position of Al window + Pinhole (980 mm)
WindThick=3*(10^-6);                                   %Al window thickness (3 microns) 
WindDens=2699;                                         %Al density (2699 kg/m^3)
WindMassAbsCoeff=7.88*(10^1);                          %Al mass absorption coefficient for 2.96 keV X-ray photons (m^2/kg)
WindAbsCoeff=WindDens*WindMassAbsCoeff;                %Absorption coefficient in cm^-1
TungstenThick=(10^-3);                                 %W window thickness (1 mm) 
TungstenMassAbsCoeff=1.902*(10^2);                     %W mass absorption coefficient for 2.96 keV X-ray photons (m^2/kg)
TungstenDens=19300;                                    %W density (19300 kg/m^3)
TungstenAbsCoeff=TungstenDens*TungstenMassAbsCoeff;    %Absorption coefficient in cm^-1
LeadThick=0.2*(10^-3);                                 %Pb diaphragm thickness (0.2 mm) 
LeadDens=11350;                                        %Pb density (11.35 g/cm^3)
LeadMassAbsCoeff=1.965*(10^2);                         %Pb mass absorption coefficient for 2.96 keV X-ray photons (m^2/kg)
LeadAbsCoeff=LeadDens*LeadMassAbsCoeff;                %Absorption coefficient in cm^-1
PinholeDia=100*(10^-6);                                %Pinhole diameter (100 microns)
df=0.155+1.191;                                        %Z-Position of CCD camera
CCDxStart=-13.8*(10^-3);                               %Start X-Position of CCD camera
CCDxStop=13.8*(10^-3);                                 %Stop X-Position of CCD camera
CCDyStart=-3.45*(10^-3);                               %Start Y-Position of CCD camera
CCDyStop=3.45*(10^-3);                                 %Stop Y-Position of CCD camera
c=3*(10^8);                                            %Speed of light (m/s)

%% Declare required matrices and variables

PosMatX=zeros(59,59,211);        %Store X,Y and Z coordinates of the plasma cells
PosMatY=zeros(59,59,211);
PosMatZ=zeros(59,59,211);

AzimAngle_central=zeros(59,59,211);      %Store azimuthal angles (polar angles will be directly defined)

LGEMat=zeros(59,59,211);         %Store local geometrical efficieny (absorption taken into account)
ntcount=0;                       %Keep count of number of emissions for calculation of LGE
nt=1;                         %Define number of MC steps (effective number of photon emissions)
PolAngleLim=pi/6;                %Define parameters for restricted area (polar angles)
AzimAngleLim=pi/6;               %Define parameters for restricted area (azimuthal angles)

%% Define detectable directions

for i=1:59                       %Run loop to calculate initial positions and detectable azimuthal angles
    for j=1:59
        for k=1:211
            PosMatX(i,j,k)=(-29+(j-1))*(10^-3);          %Origin shifted to (30,30,1)
            PosMatY(i,j,k)=(29-(i-1))*(10^-3);
            PosMatZ(i,j,k)=k*(10^-3);
            if (i>=30&&j<=30)
                AzimAngle_central(i,j,k)=atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            elseif (i>=30&&j>30)
                AzimAngle_central(i,j,k)=pi-atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            elseif (i<30&&j>30)
                AzimAngle_central(i,j,k)=pi+atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            else
                AzimAngle_central(i,j,k)=(2*pi)-atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            end
        end
    end
end

AzimAngle_central(isnan(AzimAngle_central))=0;
PolAngle_central=atan(sqrt((PosMatX.^2)+(PosMatY.^2))./(d0-PosMatZ));   %Define detectable polar angles

%% Execute Monte-Carlo simulation

rng(0,'twister');      %Initialse random number generator
ntcount=0;

for i=1:nt                                      %Begin MC loop
    
    EmissionMat=ones(59,59,211);                %Setting all cells to emit one photon (direction will be randomized later)
    
    a=-1;
    b=1;
    
    %u=(b-a).*rand(59,59,211)+a;    %Define 3D matrix of random numbers between -1 and 1
    %v=rand(59,59,211);             %Define 3D matrix of random numbers between 0 and 1
    u=zeros(59,59,211);
    v=zeros(59,59,211);
    
    PolAngle=PolAngle_central+PolAngleLim.*(2.*v-1);       %Obtain uniformly distributed random polar angles for each plasma cell in restricted area
    AzimAngle=AzimAngle_central+AzimAngleLim.*u;           %Obtain uniformly distributed random azimuthal angles for each plasma cell in restricted area
    
    DirCosX=sin(PolAngle).*cos(AzimAngle);    %Define direction cosine matrices 
    DirCosY=sin(PolAngle).*sin(AzimAngle);
    DirCosZ=cos(PolAngle);
     
    Time=abs((d0-PosMatZ)./(c.*DirCosZ));      %Photon translation from position of emission to pinhole setup
    PosMatX=PosMatX+(c.*DirCosX.*Time);        %in X, Y and Z directions
    PosMatY=PosMatY+(c.*DirCosY.*Time);
    PosMatZ=PosMatZ+(c.*DirCosZ.*Time);
        
    RadMat=sqrt((PosMatX.^2)+(PosMatY.^2));
    EmissionMat(PosMatZ~=d0)=0;
    EmissionMat(RadMat>(PinholeDia/2))=EmissionMat(RadMat>(PinholeDia/2)).*exp(-TungstenAbsCoeff.*TungstenThick./DirCosZ(RadMat>(PinholeDia/2)));
    EmissionMat(RadMat>(PinholeDia/2))=EmissionMat(RadMat>(PinholeDia/2)).*exp(-LeadAbsCoeff.*LeadThick./DirCosZ(RadMat>(PinholeDia/2)));
        
    Time=abs((df-d0)./(c.*DirCosZ));
    PosMatX=PosMatX+(c.*DirCosX.*Time);
    PosMatY=PosMatY+(c.*DirCosY.*Time);
    PosMatZ=PosMatZ+(c.*DirCosZ.*Time);
    EmissionMat(PosMatZ~=df)=0;
    EmissionMat(PosMatX<CCDxStart)=0;
    EmissionMat(PosMatX>CCDxStop)=0;
    EmissionMat(PosMatY<CCDyStart)=0;
    EmissionMat(PosMatY>CCDyStop)=0;
        
    EmissionMat=EmissionMat.*exp(log(MeshTransparency)./DirCosZ);
    EmissionMat(isnan(EmissionMat))=0;
    EmissionMat=EmissionMat.*exp(-WindAbsCoeff.*WindThick./DirCosZ);
    EmissionMat(isnan(EmissionMat))=0;
        
    LGEMat=LGEMat+EmissionMat; 
    ntcount=ntcount+1; 
end

SolAngle=(cos(PolAngle_central-PolAngleLim)-cos(PolAngle_central+PolAngleLim)).*...
         ((AzimAngle_central+AzimAngleLim)-(AzimAngle_central-AzimAngleLim));
LGEMat=(LGEMat./ntcount).*(SolAngle./(4.*pi));

