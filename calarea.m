function area = calarea(kT,index,operation)

area=0;

switch (operation)
    
    case 1
        
        Epts=linspace(0,12,7);

        if (index<=7)
            fn = @(x,y) (2/sqrt(pi))*(1/(x^(3/2)))*sqrt(y)*exp(-y/x);   %x stands for kT and y for E
            if (index<7)
                delE=0.1;
                Eint=linspace(Epts(index),Epts(index+1),(Epts(index+1)-Epts(index))/delE+1);
                for i=1:(length(Eint)-1)
                    area=area+0.5*(fn(kT,Eint(i))+fn(kT,Eint(i+1)))*delE;
                end
            else
                delE=0.5;
                Eint=linspace(Epts(index),20,(20-Epts(index))/delE+1);
                for i=1:(length(Eint)-1)
                    area=area+0.5*(fn(kT,Eint(i))+fn(kT,Eint(i+1)))*delE;
                end
            end        
        else
            fn = @(x,y) (2/sqrt(pi))*(1/(x^(3/2)))*(y^(3/2))*exp(-y/x);   %x stands for kT and y for E
            if ((index-7)<7)
               delE=0.1;
               Eint=linspace(Epts(index-7),Epts(index-6),(Epts(index-6)-Epts(index-7))/delE+1);
               for i=1:(length(Eint)-1)
                   area=area+0.5*(fn(kT,Eint(i))+fn(kT,Eint(i+1)))*delE;
               end
            else
               delE=0.5;
               Eint=linspace(Epts(index-7),20,(20-Epts(index-7))/delE+1);
               for i=1:(length(Eint)-1)
                   area=area+0.5*(fn(kT,Eint(i))+fn(kT,Eint(i+1)))*delE;
               end
            end  
        end
        
    case 2
        
        fn = @(x,y) (2/sqrt(pi))*(1/(x^(3/2)))*sqrt(y)*exp(-y/x);   %x stands for kT and y for E
        
        delE=0.02;
        Eint=linspace(index-0.1,index,0.1/delE+1);
        for i=1:(length(Eint)-1)
             area=area+0.5*(fn(kT,Eint(i))+fn(kT,Eint(i+1)))*delE;
        end
        
end

end