function area = calarea_v5(NormCoeff,CrSection,Temperature)

area=zeros(59,59,211);

delE=0.1;
enercount=0:delE:35;
c=3*(10^11);
vel=c.*sqrt(2.*enercount./511);

fn = @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^(3./2))).*exp(-y./x);

for i=1:(length(enercount)-1)
    
    area(:,:,:)=area(:,:,:)+0.5.*(CrSection(i).*(NormCoeff(:,:,:,1).*fn(Temperature(:,:,:,1),enercount(i))+...
                                                 NormCoeff(:,:,:,2).*fn(Temperature(:,:,:,2),enercount(i))+...
                                                 NormCoeff(:,:,:,3).*fn(Temperature(:,:,:,3),enercount(i))).*vel(i)+...
                                  CrSection(i+1).*(NormCoeff(:,:,:,1).*fn(Temperature(:,:,:,1),enercount(i+1))+...
                                                   NormCoeff(:,:,:,2).*fn(Temperature(:,:,:,2),enercount(i+1))+...
                                                   NormCoeff(:,:,:,3).*fn(Temperature(:,:,:,3),enercount(i+1))).*vel(i+1)).*delE;
    
end
area(isnan(area))=0;

end