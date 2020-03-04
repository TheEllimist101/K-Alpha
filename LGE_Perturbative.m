%This code will evaluate the Local Geometrical Efficiency of the pinhole
%detection system

%This is based on the second approach wherein the polar and azimuthal
%angles are slightly perturbed from their known detection values to
%ascertain the detection limits and hence the solid angle into which the
%emission can take place for each plasma cell

%This is not a Monte Carlo randomisation approach since computational
%demands are too high for the former

%Code does not evaluate the photon counting algorithm or intrisic photon to
%charge cloud conversion process in the CCD

%All quantites expressed in SI units unless present as a product that
%cancels dimensions

%Define system parameters

MeshTransparency=0.56;                                 %Mesh transparency at normal incidence (56%)
d0=1.191;                                              %Position of Al window + Pinhole (980 mm)
WindThick=3*(10^-4);                                   %Al window thickness (3 microns) (expressed here in cm)
WindDens=2.699;                                        %Al density (2.699 g/cm^3)
WindMassAbsCoeff=7.88*(10^2);                          %Mass absorption coefficient for 2.96 keV X-ray photons (cm^2/g)
WindAbsCoeff=WindDens*WindMassAbsCoeff;                %Absorption coefficient in cm^-1
TungstenThick=0.1;                                     %Tungsten window thickness (1 mm) (expressed here in cm)
TungstenMassAbsCoeff=1.873*(10^3);                     %Mass absorption coefficient for 2.96 keV X-ray photons (cm^2/g)
TungstenDens=19.3;                                     %Tungsten density (19.3 g/cm^3)
TungstenAbsCoeff=TungstenDens*TungstenMassAbsCoeff;    %Absorption coefficient in cm^-1
LeadThick=0.02;                                        %Lead diaphragm thickness (0.2 mm) (expressed here in cm)
LeadDens=11.35;                                        %Lead density (11.35 g/cm^3)
LeadMassAbsCoeff=1.913*(10^3);                         %Mass absorption coefficient for 2.96 keV X-ray photons (cm^2/g)
LeadAbsCoeff=LeadDens*LeadMassAbsCoeff;                %Absorption coefficient in cm^-1
PinholeDia=100*(10^-6);                                %Pinhole diameter (100 microns)
df=0.155+1.191;                                        %Z-Position of CCD camera
CCDxStart=-13.8*(10^-3);                               %Start X-Position of CCD camera
CCDxStop=13.8*(10^-3);                                 %Stop X-Position of CCD camera
CCDyStart=-3.45*(10^-3);                               %Start Y-Position of CCD camera
CCDyStop=3.45*(10^-3);                                 %Stop Y-Position of CCD camera
c=3*(10^8);                                            %Speed of light (m/s)

%Define required matrices

PosMatX=zeros(59,59,211);        %Store X,Y and Z coordinates of the plasma cells
PosMatY=zeros(59,59,211);
PosMatZ=zeros(59,59,211);

AzimAngle=zeros(59,59,211);      %Store azimuthal angles

LGEMat=zeros(59,59,211);         %Store effective geometrical efficieny (absorption taken into account)
npertcount=0;                    %Keep count of number of emissions for calculation of LGE

for i=1:59                       %Define all initial positions and angles
    for j=1:59
        for k=1:211
            PosMatX(i,j,k)=(-29+(j-1))*(10^-3);
            PosMatY(i,j,k)=(29-(i-1))*(10^-3);
            PosMatZ(i,j,k)=k*(10^-3);
            if (i>=30&&j<=30)
                AzimAngle(i,j,k)=atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            elseif (i>=30&&j>30)
                AzimAngle(i,j,k)=pi-atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            elseif (i<30&&j>30)
                AzimAngle(i,j,k)=pi+atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            else
                AzimAngle(i,j,k)=(2*pi)-atan(abs(PosMatY(i,j,k))/abs(PosMatX(i,j,k)));
            end
        end
    end
end

PolAngle=atan(sqrt((PosMatX.^2)+(PosMatY.^2))./(d0-PosMatZ));   %Store polar angles

PolAngleLim1Act=zeros(59,59,211);                                  %Initiate angles
PolAngleLim2Act=zeros(59,59,211);
AzimAngleLim1Act=zeros(59,59,211);
AzimAngleLim2Act=zeros(59,59,211);

%Define perturbation angle vector

PertPol=0:1.5*(10^-12):1.5*(10^-11);
PertAzim=0:1.5*(10^-12):1.5*(10^-11);

%Initiate algorithm to evaluate scope of detection

for i=1:1%length(PertPol)
    for j=1:1%length(PertAzim)
        
        EmissionMatLim1=ones(59,59,211);
        EmissionMatLim2=ones(59,59,211);
        
        PolAngleLim1=PolAngle+PertPol(i);       %Define max and min polar and azimuthal angles 
        PolAngleLim2=PolAngle-PertPol(i);
        AzimAngleLim1=AzimAngle+PertAzim(j);
        AzimAngleLim2=AzimAngle-PertAzim(j);
        
        DirCosXLim1=sin(PolAngleLim1).*cos(AzimAngleLim1);    %Define max and min direction cosines 
        DirCosXLim2=sin(PolAngleLim2).*cos(AzimAngleLim2);
        DirCosYLim1=sin(PolAngleLim1).*sin(AzimAngleLim1);
        DirCosYLim2=sin(PolAngleLim2).*sin(AzimAngleLim2);
        DirCosZLim1=cos(PolAngleLim1);
        DirCosZLim2=cos(PolAngleLim2);
        
        TimeLim1=abs((d0-PosMatZ)./(c.*DirCosZLim1));
        TimeLim2=abs((d0-PosMatZ)./(c.*DirCosZLim2));
        PosMatXLim1=PosMatX+(c.*DirCosXLim1.*TimeLim1);
        PosMatXLim2=PosMatX+(c.*DirCosXLim2.*TimeLim2);
        PosMatYLim1=PosMatY+(c.*DirCosYLim1.*TimeLim1);
        PosMatYLim2=PosMatY+(c.*DirCosYLim2.*TimeLim2);
        PosMatZLim1=PosMatZ+(c.*DirCosZLim1.*TimeLim1);
        PosMatZLim2=PosMatZ+(c.*DirCosZLim2.*TimeLim2);
        
        RadMatLim1=sqrt((PosMatXLim1.^2)+(PosMatYLim1.^2));
        RadMatLim2=sqrt((PosMatXLim2.^2)+(PosMatYLim2.^2));
        EmissionMatLim1(PosMatZLim1~=d0)=0;
        EmissionMatLim2(PosMatZLim2~=d0)=0;
        EmissionMatLim1(RadMatLim1>(PinholeDia/2))=EmissionMatLim1(RadMatLim1>(PinholeDia/2)).*exp(-TungstenAbsCoeff.*TungstenThick./DirCosZLim1(RadMatLim1>(PinholeDia/2)));
        EmissionMatLim1(RadMatLim1>(PinholeDia/2))=EmissionMatLim1(RadMatLim1>(PinholeDia/2)).*exp(-LeadAbsCoeff.*LeadThick./DirCosZLim1(RadMatLim1>(PinholeDia/2)));
        EmissionMatLim2(RadMatLim2>(PinholeDia/2))=EmissionMatLim2(RadMatLim2>(PinholeDia/2)).*exp(-TungstenAbsCoeff.*TungstenThick./DirCosZLim2(RadMatLim2>(PinholeDia/2)));
        EmissionMatLim2(RadMatLim2>(PinholeDia/2))=EmissionMatLim2(RadMatLim2>(PinholeDia/2)).*exp(-LeadAbsCoeff.*LeadThick./DirCosZLim2(RadMatLim2>(PinholeDia/2)));
        
        TimeLim1=abs((df-d0)./(c.*DirCosZLim1));
        TimeLim2=abs((df-d0)./(c.*DirCosZLim2));
        PosMatXLim1=PosMatXLim1+(c.*DirCosXLim1.*TimeLim1);
        PosMatXLim2=PosMatXLim2+(c.*DirCosXLim2.*TimeLim2);
        PosMatYLim1=PosMatYLim1+(c.*DirCosYLim1.*TimeLim1);
        PosMatYLim2=PosMatYLim2+(c.*DirCosYLim2.*TimeLim2);
        PosMatZLim1=PosMatZLim1+(c.*DirCosZLim1.*TimeLim1);
        PosMatZLim2=PosMatZLim2+(c.*DirCosZLim2.*TimeLim2);
        EmissionMatLim1(PosMatZLim1~=df)=0;
        EmissionMatLim2(PosMatZLim2~=df)=0;
        EmissionMatLim1(PosMatXLim1<CCDxStart)=0;
        EmissionMatLim1(PosMatXLim1>CCDxStop)=0;
        EmissionMatLim2(PosMatXLim2<CCDxStart)=0;
        EmissionMatLim2(PosMatXLim2>CCDxStop)=0;
        EmissionMatLim1(PosMatYLim1<CCDyStart)=0;
        EmissionMatLim1(PosMatYLim1>CCDyStop)=0;
        EmissionMatLim2(PosMatYLim2<CCDyStart)=0;
        EmissionMatLim2(PosMatYLim2>CCDyStop)=0;
        PolAngleLim1Act(EmissionMatLim1~=0)=PolAngleLim1(EmissionMatLim1~=0);
        PolAngleLim2Act(EmissionMatLim2~=0)=PolAngleLim2(EmissionMatLim2~=0);
        AzimAngleLim1Act(EmissionMatLim1~=0)=AzimAngleLim1(EmissionMatLim1~=0);
        AzimAngleLim2Act(EmissionMatLim2~=0)=AzimAngleLim2(EmissionMatLim2~=0);
        
        EmissionMatLim1=EmissionMatLim1.*exp(log(MeshTransparency)./DirCosZLim1);
        EmissionMatLim2=EmissionMatLim2.*exp(log(MeshTransparency)./DirCosZLim2);
        EmissionMatLim1(isnan(EmissionMatLim1))=0;
        EmissionMatLim2(isnan(EmissionMatLim2))=0;
        EmissionMatLim1=EmissionMatLim1.*exp(-WindAbsCoeff.*WindThick./DirCosZLim1);
        EmissionMatLim2=EmissionMatLim2.*exp(-WindAbsCoeff.*WindThick./DirCosZLim2);
        EmissionMatLim1(isnan(EmissionMatLim1))=0;
        EmissionMatLim2(isnan(EmissionMatLim2))=0;
        
        if (i==1&&j==1)
            LGEMat(EmissionMatLim1~=0)=EmissionMatLim1(EmissionMatLim1~=0).*sin(PolAngleLim1Act(EmissionMatLim1~=0));
        else
            LGEMat(EmissionMatLim1~=0)=(1/4).*(LGEMat(EmissionMatLim1~=0)+EmissionMatLim1(EmissionMatLim1~=0).*sin(PolAngleLim1Act(EmissionMatLim1~=0))).*(10^-4);
            LGEMat(EmissionMatLim2~=0)=(1/4).*(LGEMat(EmissionMatLim2~=0)+EmissionMatLim2(EmissionMatLim2~=0).*sin(PolAngleLim2Act(EmissionMatLim2~=0))).*(10^-4);
        end
        
    end
end

LGEMat=LGEMat./(4.*pi);

%Plot Local Geometrical Efficiency

figure;
imagesc(LGEMat(:,:,105));
title("Z = 105 mm");
xlabel("X axis (mm)");
ylabel("Y axis (mm)")

figure;
imagesc(LGEMat(:,:,45));
title("Z = 45 mm");
xlabel("X axis (mm)");
ylabel("Y axis (mm)");

figure;
imagesc(LGEMat(:,:,200));
title("Z = 100 mm");
xlabel("X axis (mm)");
ylabel("Y axis (mm)");





