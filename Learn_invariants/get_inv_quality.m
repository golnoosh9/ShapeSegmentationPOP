maxp=0;
mine=0;
mon_num=15;
inv_num=2*mon_num+(mon_num-1)*mon_num/2;
while(mine-maxp<=0.1)
R=inv_coeffs;%rmvnrnd(exp,invs,1);
%R(mon_num+5)=100;
R2=reshape((R'*R),inv_num^2,1);
errors=A(:,1+pnum+nnum:pnum+nnum+inv_num^2)*R2;


perr=-errors(1:pnum);
nerr=errors(pnum+1:pnum+nnum);
%minp=min(perr);
maxp=mean(perr);
mine=mean(abs(nerr));
uuuu
end
tpnum=10000;
tnnum=10000;
count=1;
for i=1:tpnum
    temp=ones(1,mon_num+mon_num+(mon_num-1)*mon_num/2);
    rcoef=(rand(mon_num,1)-ones(mon_num,1)/2)*2;
     radi=rand(1)/2;
     %rcoef(5)=0;
     rcoef(6)=rcoef(4); 
     if(rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2<=0.1||rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2>=1)
         i=i-1;
         continue
     end
%      rcoef(6)=rcoef(4);
%      rcoef(1)=radi^2-0.25*(rcoef(2)^2)-0.25*(rcoef(3)^2);
     pcoeffs(:,i)=rcoef;
     temp(1:mon_num)=rcoef;
     pmat=rcoef*rcoef';
     temp(mon_num+1:mon_num+mon_num)=diag(pmat);
     row1=pmat(1,2:6);
    row2=pmat(2,3:6);
    row3=pmat(3,4:6);
    row4=pmat(4,5:6);
    row5=pmat(5,6);
    temp(1+2*mon_num:mon_num+mon_num+(mon_num-1)*mon_num/2)=[row1,row2,row3,row4,row5];
 terr_p(count)=temp*R';
     count=count+1;
    
end
max_test_p=mean((terr_p));

count=1;
for i=1+tpnum:tnnum+tpnum
   rcoef=(rand(mon_num,1)-ones(mon_num,1)/2)*2;
 % rcoef(5)=0.5;
 if abs(rcoef(6)==rcoef(4))<=0.1
     
     continue;
 end
 if(rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2<=0.99||rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2>=1)
         i=i-1;
         continue
     end
   temp(1:mon_num)=rcoef;
     pmat=rcoef*rcoef';
     temp(mon_num+1:mon_num+mon_num)=diag(pmat);
     row1=pmat(1,2:6);
    row2=pmat(2,3:6);
    row3=pmat(3,4:6);
    row4=pmat(4,5:6);
    row5=pmat(5,6);
    temp(2*mon_num+1:mon_num+mon_num+(mon_num-1)*mon_num/2)=[row1,row2,row3,row4,row5];
     terr_n(count)=temp*R';
     count=count+1;
end
min_test_n=mean(abs(terr_n));