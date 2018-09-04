varst7=[];

for xx=1:mon_num
    varst2=zeros(1,cnum+mon_num+2+vercount+cnum);
    varst2(xx)=1;
    %varst2(mon_num+indj)=1; for weight
    varst7=[varst7;varst2];    
    
end

for xx=1:mon_num
    varst2=zeros(1,cnum+mon_num+2+vercount+cnum);
    varst2(xx)=2;
    varst7=[varst7;varst2];    
    
end

for ii=1:mon_num
    for jj=ii+1:mon_num
        varst2=zeros(1,cnum+mon_num+2+vercount+cnum);
        varst2(ii)=1;
        varst2(jj)=1;
        varst7=[varst7;varst2];
    end
end
