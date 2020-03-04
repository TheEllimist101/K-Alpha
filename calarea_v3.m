function area = calarea_v3(CrSection,Temperature)

area=zeros(59,59,211);

delE=0.1;
enercount=0:delE:20;
c=3*(10^11);
vel=c.*sqrt(2.*enercount./511);

fn = @(x,y) (2./sqrt(pi)).*(sqrt(y)./(x.^(3./2))).*exp(-y./x);

for i=1:(length(enercount)-1)
    
    area(:,:,:)=area(:,:,:)+0.5.*(CrSection(i).*fn(Temperature(:,:,:,1),enercount(i)).*vel(i)+...
                                  CrSection(i+1).*fn(Temperature(:,:,:,1),enercount(i+1)).*vel(i+1)).*delE;
    
end
area(isnan(area))=0;

end