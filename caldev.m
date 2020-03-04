function res = caldev(ActualCoeff,AvTemp,DensityMat,EnerDensMat)

fn= @(x,y) (2./sqrt(pi)).*((sqrt(y)./x.^(3./2))).*exp(-y./x);

densdev=zeros(59,59,211,7);
enerdensdev=zeros(59,59,211,7);

energy=0:0.1:20;
j=1;

DensTot=sum(DensityMat,4);

for i=1:7
    
    if (i<7)
        while((j*0.1)<(2*i))
            densdev(:,:,:,i)=densdev(:,:,:,i)+DensTot(:,:,:).*((ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+...
                             ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+(ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+...
                             ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))).*0.5.*0.1;
            enerdensdev(:,:,:,i)=enerdensdev(:,:,:,i)+(0.5.*0.1.*(energy(j).*(ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+...
                                 ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+energy(j+1).*...
                                (ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))))./...
                                (((ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+...
                                (ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))).*0.5.*0.1);
            j=j+1;
        end
    else
        while ((j*0.1)<20)
            densdev(:,:,:,i)=densdev(:,:,:,i)+DensTot(:,:,:).*((ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+...
                             ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+(ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+...
                             ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))).*0.5.*0.1;
            enerdensdev(:,:,:,i)=enerdensdev(:,:,:,i)+(0.5.*0.1.*(energy(j).*(ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+...
                                 ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+energy(j+1).*...
                                (ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))))./...
                                (((ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j)))+...
                                (ActualCoeff(:,:,:,1).*fn(EnerDensMat(:,:,:,1)./360,energy(j+1))+ActualCoeff(:,:,:,2).*fn(AvTemp(:,:,:),energy(j+1)))).*0.5.*0.1);
            j=j+1;

        end
    end
    
    densdev(:,:,:,i)=densdev(:,:,:,i)-DensityMat(:,:,:,i);
    enerdensdev(:,:,:,i)=enerdensdev(:,:,:,i)-EnerDensMat(:,:,:,i);
    
end

res=1;
    
figure;
for i=1:6
    subplot(2,3,i);
    imagesc(densdev(:,:,105,i));
end
title('\bfDensity Deviations');

figure;
imagesc(densdev(:,:,105,7));
title('\bfDensity Deviations');

figure;
for i=1:6
    subplot(2,3,i);
    imagesc(enerdensdev(:,:,105,i));
end
title('\bfAverage Energy Deviations');

figure;
imagesc(enerdensdev(:,:,105,7));
title('\bfAverage Energy Deviations');

    


