function area = calarea_v4(NormCoeff,CrSection,Temperature)

area=zeros(59,59,211);

delE=0.001;
enercount=0:delE:35;
c=3*(10^11);
vel=c.*sqrt(2.*enercount./511);
j=1;

fnMax= @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^1.5)).*exp(-y./x);         %Basic Maxwell distribution function
fnDr= @(x,y) 1.04.*(sqrt(y)./(x.^1.5)).*exp(-0.55.*(y.^2)./(x.^2));   %Basic Druyvesteyn distrbution function

for i=1:(length(enercount)-1)
    
    area(:,:,:)=area(:,:,:)+0.5.*(CrSection(i).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(i))+...
                                                 NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(i))).*vel(i)+...
                                  CrSection(i+1).*(NormCoeff(:,:,:,1).*fnMax(Temperature(:,:,:,1),enercount(i+1))+...
                                                   NormCoeff(:,:,:,2).*fnDr(Temperature(:,:,:,2),enercount(i+1))).*vel(i+1)).*delE;
    if (enercount(i)>0&&(mod(enercount(i),2)==0)&&enercount(i)<=12)
        j=j+1;
    end
    
end
area(isnan(area))=0;

end