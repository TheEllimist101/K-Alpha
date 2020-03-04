function area = calarea_v4(NormCoeff,CrSection,Temperature,DensityMat)

area=zeros(59,59,211);
DensTot=sum(DensityMat,4);

%Calculate approximated density and average energy (for correction factor)

DensApprox=zeros(59,59,211,7);                              %Numerically approximated density
EnerDensApprox=zeros(59,59,211,7);                          %Numerically approximated average energy
normint=zeros(59,59,211,7);                                 %Area under the distribution function
correction=zeros(59,59,211,7);

delE=0.1;
enercount=0:delE:35;
c=3*(10^11);
vel=c.*sqrt(2.*enercount./511);
j=1;

fnMax= @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^1.5)).*exp(-y./x);         %Basic Maxwell distribution function
fnDr= @(x,y) 1.04.*(sqrt(y)./(x.^1.5)).*exp(-0.55.*(y.^2)./(x.^2));   %Basic Druyvesteyn distrbution function

for i=1:7
    
    if (i<7)
        while((j*0.1)<(2*i))                                      %calculate area under distribution function
            normint(:,:,:,i)=normint(:,:,:,i)+0.5.*0.1.*((NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(j)))...
                                                        +(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j+1))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(j+1))));
            j=j+1;
        end
    else
        while ((j*0.1)<35)
            normint(:,:,:,i)=normint(:,:,:,i)+0.5.*0.1.*((NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(j)))...
                                                        +(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j+1))...
                                                          +NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(j+1))));
            j=j+1;
        end
    end
    
end

j=1;
for i=1:7
              
    if (i<7)
        while((j*0.1)<(2*i))
            EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)+0.5.*0.1.*(enercount(j).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j))...
                                                                                  +NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(j)))...
                                                                       +enercount(j+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j+1))...
                                                                                  +NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(j+1))));
            j=j+1;
        end
        DensApprox(:,:,:,i)=DensTot(:,:,:).*normint(:,:,:,i);                                   %Calculate densities
        correction(:,:,:,i)=(DensApprox(:,:,:,i)-DensityMat(:,:,:,i))./(DensTot(:,:,:).*2);     %Calculate correction factor
    else
        while ((j*0.1)<35)
            EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)+0.5.*0.1.*(enercount(j).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j))...
                                                                                  +NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(j)))...
                                                                       +enercount(j+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(j+1))...
                                                                                  +NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(j+1))));
            j=j+1;
        end
        DensApprox(:,:,:,i)=DensTot(:,:,:).*normint(:,:,:,i);                                    %Calculate densities
        correction(:,:,:,i)=(DensApprox(:,:,:,i)-DensityMat(:,:,:,i))./(DensTot(:,:,:).*13);     %Calculate correction factor
    end
        
     EnerDensApprox(:,:,:,i)=EnerDensApprox(:,:,:,i)./normint(:,:,:,i);    %Calculate average energies
     EnerDensApprox(isnan(EnerDensApprox))=0;
    
end

j=1;
for i=1:(length(enercount)-1)
    
    area(:,:,:)=area(:,:,:)+0.5.*(CrSection(i).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(i))+...
                                                 NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(i))-...
                                                 correction(:,:,:,j)).*vel(i)+...
                                  CrSection(i+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(i+1))+...
                                                   NormCoeff(:,:,:,2).*fnMax(Temperature(:,:,:,2),enercount(i+1))-...
                                                   correction(:,:,:,j)).*vel(i+1)).*delE;
    if (enercount(i)>0&&(mod(enercount(i),2)==0)&&enercount(i)<=12)
        j=j+1;
    end
    
end
area(isnan(area))=0;

end