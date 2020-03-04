%This code will calculate the Local Geometrical Efficiency of each cell in
%the plasma chamber using the same approach as the Restricted Monte Carlo
%method but the emission space (area and orientaation) will be defined as
%a function of the coordinates of the cells

%The code will be divided into 2 parts: the first will deal with
%calculation of the limits of the polar and azimuthal angles that will
%define the emission space, and the second will run multiple MC iterations
%to generate as many emissions as possible within the emission space

%This version differs from LGE_SelMC_v1 in the sense that the for-loops
%over MC will be vectorised to a 4D matrix

%The filtering algorithm will be optimised wtih respect to the perturbative
%code to account for fewer checks needed here

%The LGE will be calculated with respect to the selected emission space
%alone

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

%% Part 1: Selection of detectable emission space

PosMatX_init=zeros(59,59,211);           %Define 3D matrix to store X coordinate
PosMatY_init=zeros(59,59,211);           %Define 3D matrix to store Y coordinate
PosMatZ_init=zeros(59,59,211);           %Define 3D matrix to store Z coordinate

AzimAngleMin=zeros(59,59,211);      %Define 3D matrices to store the azimuthal angle
AzimAngleMax=zeros(59,59,211);      %limits of the emission space

for i=1:59                          %Run loop to recenter the origin to (30,30,1)
    for j=1:59
        for k=1:211
            PosMatX_init(i,j,k)=(-29+(j-1))*(10^-3);          
            PosMatY_init(i,j,k)=(29-(i-1))*(10^-3);
            PosMatZ_init(i,j,k)=k*(10^-3);
        end
    end
end

XdiffCorner1=PosMatX_init-CCDxStop;        %Matrices to store the distance along X axis from the corners of the CCD camera
XdiffCorner2=PosMatX_init-CCDxStart;
XdiffCorner3=XdiffCorner2;
XdiffCorner4=XdiffCorner1;

YdiffCorner1=PosMatY_init-CCDyStop;        %Matrices to store the distance along Y axis from the corners of the CCD camera    
YdiffCorner2=YdiffCorner1;
YdiffCorner3=PosMatY_init-CCDyStart;
YdiffCorner4=YdiffCorner3;

RadCorner1=sqrt(XdiffCorner1.^2+YdiffCorner1.^2);       %Matrices to store the complete distances from the corners of the CCD camera   
RadCorner2=sqrt(XdiffCorner2.^2+YdiffCorner2.^2);      
RadCorner3=sqrt(XdiffCorner3.^2+YdiffCorner3.^2);      
RadCorner4=sqrt(XdiffCorner4.^2+YdiffCorner4.^2);      

temp_arr=cat(4,RadCorner1,RadCorner2,RadCorner3,RadCorner4);    %Temporary array for purpose of comparison

PolAngleMax=atan(max(temp_arr,[],4)./(df-PosMatZ_init));      %Define and initialise 3D matrices to store polar angle
PolAngleMin=atan(min(temp_arr,[],4)./(df-PosMatZ_init));      %limits of the emission space

clear temp_arr    %clear memory

SlopeCorner1=abs(YdiffCorner1./XdiffCorner1);     %Matrices to store the slopes from the corners of the CCD camera
SlopeCorner2=abs(YdiffCorner2./XdiffCorner2);
SlopeCorner3=abs(YdiffCorner3./XdiffCorner3);
SlopeCorner4=abs(YdiffCorner4./XdiffCorner4);

%The X-Y plane of the plasma chamber is divided into 9 sections:
%C,D1,E1,D2,E2,D3,E3,D4 and E4, to facilitate understanding of emission
%angles

%Sector C
PolAngleMin(PosMatX_init>=CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<=CCDyStop)=0;      %The region parallel to the CCD has a 2pi azmuthal coverage  
AzimAngleMin(PosMatX_init>=CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<=CCDyStop)=0;     %and min polar angle is 0
AzimAngleMax(PosMatX_init>=CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<=CCDyStop)=2*pi;

%Sector D1
AzimAngleMax(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStop)=pi+atan(SlopeCorner4(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStop));
AzimAngleMin(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStop)=pi+atan(SlopeCorner2(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStop));

%Sector E1
AzimAngleMax(PosMatX_init>CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>CCDyStop)=2*pi-atan(SlopeCorner1(PosMatX_init>CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>CCDyStop));
AzimAngleMin(PosMatX_init>CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>CCDyStop)=pi+atan(SlopeCorner2(PosMatX_init>CCDxStart&PosMatX_init<=CCDxStop&PosMatY_init>CCDyStop));

%Sector D2
AzimAngleMax(PosMatX_init<=CCDxStart&PosMatY_init>CCDyStop)=2*pi-atan(SlopeCorner1(PosMatX_init<=CCDxStart&PosMatY_init>CCDyStop));
AzimAngleMin(PosMatX_init<=CCDxStart&PosMatY_init>CCDyStop)=2*pi-atan(SlopeCorner3(PosMatX_init<=CCDxStart&PosMatY_init>CCDyStop));

%Sector E2
AzimAngleMax(PosMatX_init<CCDxStart&PosMatY_init>CCDyStart&PosMatY_init<=CCDyStop)=2*pi+atan(SlopeCorner2(PosMatX_init<CCDxStart&PosMatY_init>CCDyStart&PosMatY_init<=CCDyStop));
AzimAngleMin(PosMatX_init<CCDxStart&PosMatY_init>CCDyStart&PosMatY_init<=CCDyStop)=2*pi-atan(SlopeCorner3(PosMatX_init<CCDxStart&PosMatY_init>CCDyStart&PosMatY_init<=CCDyStop));

%Sector D3
AzimAngleMax(PosMatX_init<CCDxStart&PosMatY_init<=CCDyStart)=2*pi+atan(SlopeCorner2(PosMatX_init<CCDxStart&PosMatY_init<=CCDyStart));
AzimAngleMin(PosMatX_init<CCDxStart&PosMatY_init<=CCDyStart)=2*pi+atan(SlopeCorner4(PosMatX_init<CCDxStart&PosMatY_init<=CCDyStart));

%Sector E3
AzimAngleMax(PosMatX_init>=CCDxStart&PosMatX_init<CCDxStop&PosMatY_init<CCDyStart)=3*pi-atan(SlopeCorner3(PosMatX_init>=CCDxStart&PosMatX_init<CCDxStop&PosMatY_init<CCDyStart));
AzimAngleMin(PosMatX_init>=CCDxStart&PosMatX_init<CCDxStop&PosMatY_init<CCDyStart)=2*pi+atan(SlopeCorner4(PosMatX_init>=CCDxStart&PosMatX_init<CCDxStop&PosMatY_init<CCDyStart));

%Sector D4
AzimAngleMax(PosMatX_init>=CCDxStop&PosMatY_init<CCDyStart)=3*pi-atan(SlopeCorner3(PosMatX_init>=CCDxStop&PosMatY_init<CCDyStart));
AzimAngleMin(PosMatX_init>=CCDxStop&PosMatY_init<CCDyStart)=3*pi-atan(SlopeCorner1(PosMatX_init>=CCDxStop&PosMatY_init<CCDyStart));

%Sector E4
AzimAngleMax(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<CCDyStop)=3*pi+atan(SlopeCorner4(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<CCDyStop));
AzimAngleMin(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<CCDyStop)=3*pi-atan(SlopeCorner1(PosMatX_init>CCDxStop&PosMatY_init>=CCDyStart&PosMatY_init<CCDyStop));

%% Part 2: Implementation of Monte Carlo

nt=100;                   %Define number of MC points on the emission surface
ntcount=0;                 %For averaging

EmissionMat=ones(59,59,211,nt);    %Setting all cells to emit one photon 

aPol=(0.5).*(1-cos(PolAngleMin));
bPol=(0.5).*(1-cos(PolAngleMax));
    
u=(bPol-aPol).*rand(59,59,211,nt)+aPol;     %Define 3D matrix of random numbers between cosine values of min and max Polar angle
v=rand(59,59,211,nt);                       %Define 3D matrix of random numbers between 0 and 1
     
PolAngle=acos(1-2.*u);                                    %Obtain uniformly distributed random polar angles for each plasma cell in restricted area
AzimAngle=(AzimAngleMax-AzimAngleMin).*v+AzimAngleMin;    %Obtain uniformly distributed random azimuthal angles for each plasma cell in restricted area
    
DirCosX=sin(PolAngle).*cos(AzimAngle);    %Define direction cosine matrices 
DirCosY=sin(PolAngle).*sin(AzimAngle);
DirCosZ=cos(PolAngle);    

PosMatX=zeros(59,59,211,nt);
PosMatY=zeros(59,59,211,nt);
PosMatZ=zeros(59,59,211,nt);

for i=1:nt
    PosMatX(:,:,:,i)=PosMatX_init;
    PosMatY(:,:,:,i)=PosMatY_init;
    PosMatZ(:,:,:,i)=PosMatZ_init; 
end

Time=abs((d0-PosMatZ)./(c.*DirCosZ));      %Photon translation from position of emission to pinhole setup
PosMatX=PosMatX+(c.*DirCosX.*Time);        %in X, Y and Z directions
PosMatY=PosMatY+(c.*DirCosY.*Time);
PosMatZ=PosMatZ+(c.*DirCosZ.*Time);
        
RadMat=sqrt((PosMatX.^2)+(PosMatY.^2));
EmissionMat(ismembertol(PosMatZ,d0)==0)=0;
EmissionMat(RadMat>(PinholeDia/2))=EmissionMat(RadMat>(PinholeDia/2)).*exp(-TungstenAbsCoeff.*TungstenThick./DirCosZ(RadMat>(PinholeDia/2)));
EmissionMat(RadMat>(PinholeDia/2))=EmissionMat(RadMat>(PinholeDia/2)).*exp(-LeadAbsCoeff.*LeadThick./DirCosZ(RadMat>(PinholeDia/2)));
           
Time=abs((df-d0)./(c.*DirCosZ));
PosMatX=PosMatX+(c.*DirCosX.*Time);
PosMatY=PosMatY+(c.*DirCosY.*Time);
PosMatZ=PosMatZ+(c.*DirCosZ.*Time);
EmissionMat(ismembertol(PosMatZ,df)==0)=0;
EmissionMat(PosMatX<CCDxStart)=0;
EmissionMat(PosMatX>CCDxStop)=0;
EmissionMat(PosMatY<CCDyStart)=0;
EmissionMat(PosMatY>CCDyStop)=0;
        
EmissionMat=EmissionMat.*exp(log(MeshTransparency)./DirCosZ);
EmissionMat(isnan(EmissionMat))=0;
EmissionMat=EmissionMat.*exp(-WindAbsCoeff.*WindThick./DirCosZ);
EmissionMat(isnan(EmissionMat))=0;

LGEMat=sum(EmissionMat,4);

SolAngle=(cos(PolAngleMin)-cos(PolAngleMax)).*(AzimAngleMax-AzimAngleMin);
LGEMat=(LGEMat./nt).*(SolAngle./(4.*pi));

