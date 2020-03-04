%This code will evaluate the effective geometrical efficiency (EGE) of the
%pinhole camera detection system

%A convergence loop will be run to calculate EGE and a Monte Carlo loop 
%nested inside the convergence loop to randomize the trajectory of each 
%emitted photon

%An average at the end of the convergence loop will determine the
%effective geometrical efficiency

%Parameters include pinhole diameter, relative position with respect to
%mesh and CCD, camera size, mesh transparency and window absorption. 

%Code does not evaluate the photon counting algorithm or intrisic photon to
%charge cloud conversion process in the CCD

%All quantites expressed in SI units unless present as a product that
%cancels dimensions

%% Define detector system geometry parameters

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

%% Define required matrices

EGE=zeros(59,59,211);

%% Define MC parameters

nt=2000;                                  %Number of convergence loops
ntcount=0;                                %For averaging

%% Initiate Monte Carlo Algorithm

rng(0,'twister');

for i=1:nt                                       %Begin convergence loop
    
    EmissionMat=ones(59,59,211);                %Setting all cells to emit one photon (direction will be randomized later)
    PosMatX=zeros(59,59,211);                   %Define initial X-coordinate array 
    PosMatY=zeros(59,59,211);                   %Define initial Y-coordinate array
    PosMatZ=zeros(59,59,211);                   %Define initial Z-coordinate array
    for m=1:59
        for j=1:59
            for k=1:211                         %Origin shifted to centre of X-Y plane
                PosMatX(m,j,k)=(-29+(m-1))*(10^-3);
                PosMatY(m,j,k)=(-29+(j-1))*(10^-3);
                PosMatZ(m,j,k)=k*(10^-3);
            end
        end
    end
    
    u=rand(59,59,211);    %Define 6000 random numbers between 0 and 1
    v=rand(59,59,211);    %Define 6000 random numbers between 0 and 1
    
    PolAngle=acos(2.*v-1);                                 %Obtain uniformly distributed random polar angles for each plasma cell
    AzimAngle=(2.*pi).*u;                                  %Obtain uniformly distributed random azimuthal angles for each plasma cell
    
    if (i==10)
        PolAngleTemp=PolAngle;
        AzimAngleTemp=AzimAngle;
    end
 
    Time=abs((MeshPos-PosMatZ)./(cos(PolAngle).*c));                 %Calculate time required to reach mesh
    PosMatX=PosMatX+(c.*sin(PolAngle).*cos(AzimAngle).*Time);        %Update X, Y and Z coordinates
    PosMatY=PosMatY+(c.*sin(PolAngle).*sin(AzimAngle).*Time);   
    PosMatZ=PosMatZ+(c.*cos(PolAngle).*Time);                           
    PosMatZ(PosMatZ<MeshPos)=NaN;
    PosMatX(PosMatX<=-0.029)=NaN;                                    %Apply filtering algorithm to remove undetected photons
    PosMatX(PosMatX>=0.029)=NaN;
    PosMatY(PosMatY<=-0.029)=NaN;
    PosMatY(PosMatY>=0.029)=NaN;
    PosMatZ(isnan(PosMatX))=NaN;
    PosMatZ(isnan(PosMatY))=NaN;
    EmissionMat(isnan(PosMatZ))=0;
    EmissionMat=EmissionMat.*exp(log(MeshTransparency)./cos(PolAngle));
    EmissionMat(isnan(EmissionMat))=0;
    
    Time=abs((d0-MeshPos)./(cos(PolAngle).*c));                      %Calculate time required to reach window
    PosMatX=PosMatX+(c.*sin(PolAngle).*cos(AzimAngle).*Time);        %Update X,Y and Z coordinates
    PosMatY=PosMatY+(c.*sin(PolAngle).*sin(AzimAngle).*Time);   
    PosMatZ=PosMatZ+(c.*cos(PolAngle).*Time);
    RadMat=sqrt((PosMatX.^2)+(PosMatY.^2));
    PosMatX(isnan(RadMat))=NaN;
    PosMatY(isnan(RadMat))=NaN;
    PosMatX(RadMat>=(PinholeDia/2))=NaN;
    PosMatY(RadMat>=(PinholeDia/2))=NaN;
    PosMatZ(isnan(PosMatX))=NaN;
    PosMatZ(isnan(PosMatY))=NaN;
    EmissionMat(isnan(PosMatZ))=0;                           %Apply filtering algorithm to remove undetected photons
    EmissionMat=EmissionMat.*exp(-WindAbsCoeff.*WindThick./cos(PolAngle));
    EmissionMat(isnan(EmissionMat))=0;
    
    Time=(df-d0)./(cos(PolAngle).*c);                           %Calculate time required to reach window
    PosMatY=PosMatY+(c.*sin(PolAngle).*sin(AzimAngle).*Time);   %Update X,Y and Z coordinates
    PosMatX=PosMatX+(c.*sin(PolAngle).*cos(AzimAngle).*Time);   
    PosMatZ=PosMatZ+(c.*cos(PolAngle).*Time);
    PosMatX(PosMatX<=CCDxStart)=NaN;                            %Apply filtering algorithm to remove undetected photons
    PosMatX(PosMatX>=CCDxStop)=NaN;
    PosMatY(PosMatY<=CCDyStart)=NaN;
    PosMatY(PosMatY>=CCDyStop)=NaN;
    PosMatZ(isnan(PosMatX))=NaN;
    PosMatZ(isnan(PosMatY))=NaN;
    EmissionMat(isnan(PosMatZ))=0;
    
    EGE=EGE+EmissionMat;
    ntcount=ntcount+1;
    
    
end

EGE=EGE/ntcount;
check=find(EGE);



