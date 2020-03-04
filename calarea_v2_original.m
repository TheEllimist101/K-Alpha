function area = calarea_v2_original(ActualCoeff,CrSection,AvTemp,EnerDensMat1)

area=zeros(59,59,211);

delE=0.1;
enercount=0:delE:20;
c=3*(10^11);
vel=c.*sqrt(2.*enercount./511);

fn = @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^(3./2))).*exp(-y./x);

for i=1:(length(enercount)-1)
    
    area(:,:,:)=area(:,:,:)+0.5.*(CrSection(i).*(ActualCoeff(:,:,:,1).*fn(EnerDensMat1(:,:,:),enercount(i))+...
                                                 ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),enercount(i))).*vel(i)+...
                                  CrSection(i+1).*(ActualCoeff(:,:,:,1).*fn(EnerDensMat1(:,:,:),enercount(i+1))+...
                                                   ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),enercount(i+1))).*vel(i+1)).*delE;
    
end
area(isnan(area))=0;

end