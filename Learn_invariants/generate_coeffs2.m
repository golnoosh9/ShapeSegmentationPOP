function [A,b]=generate_coeffs(pnum,nnum,mon_num,inv_num,negs,poss)
hhhh
%%%pcoef:circle, ncoef:other conics
eps=0.01;

%mon_num=6;
pcoeffs=zeros(mon_num,pnum);
ncoefs=zeros(mon_num,nnum);
%inv_num=1+mon_num+mon_num+(mon_num-1)*mon_num/2;
A=zeros(pnum+nnum,inv_num^2+pnum+nnum);
b=zeros(pnum+nnum,1);
i=1;
while (i<=pnum)
   % i=i+1;
    temp=ones(1,mon_num+mon_num+(mon_num-1)*mon_num/2);
    rcoef=poss(:,i);
%      radi=rand(1)/2;
%      rcoef(5)=0;
%        rcoef(6)=rcoef(4);
%    %   rcoef(1)=radi^2-0.25*(rcoef(2)^2)-0.25*(rcoef(3)^2);
%      if(rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2<=0.99||rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2>=1)
%         
%          continue
%      end
    

     pcoeffs(:,i)=rcoef;
     temp(1:mon_num)=rcoef;
     pmat=rcoef*rcoef';
     temp(1+mon_num:mon_num+mon_num)=diag(pmat);
     row1=[];
     for ii=1:mon_num-1
         
             row1=[row1,pmat(ii,ii+1:mon_num)];
    end
%      row1=pmat(1,2:6);
%     row2=pmat(2,3:6);
%     row3=pmat(3,4:6);
%     row4=pmat(4,5:6);
%     row5=pmat(5,6);
    temp(2*mon_num+1:mon_num+mon_num+(mon_num-1)*mon_num/2)=[row1];
   temp_m=temp'*temp;
   A(i,1:pnum+nnum)=0;
   A(i,i)=0.000001;
     A(i,pnum+nnum+1:pnum+nnum+inv_num^2)=-reshape(temp_m,1,inv_num^2);
    
     b(i)=0;
     
    i=i+1;
end
i=pnum+1
% 
while(i<=nnum+pnum)
   % size(negs)
    rcoef=negs(:,i-pnum);;
    size(rcoef)
%     if(abs(rcoef(6)-rcoef(4))<=0.2 && abs(rcoef(5))<=0.2 )
% %        i=i-1;
%          continue;
%     end
% 
%      if(rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2<=0.99||rcoef(1)^2+rcoef(2)^2+rcoef(3)^2+rcoef(4)^2+rcoef(5)^2+rcoef(6)^2>=1)
%         
%          continue
%      end
   temp(1:mon_num)=rcoef;
     pmat=rcoef*rcoef';
     temp(1+mon_num:mon_num+mon_num)=diag(pmat);
     row1=[];
     for ii=1:mon_num-1
         
             row1=[row1,pmat(ii,ii+1:mon_num)];
        
    end
%      row1=pmat(1,2:6);
%     row2=pmat(2,3:6);
%     row3=pmat(3,4:6);
%     row4=pmat(4,5:6);
%     row5=pmat(5,6);
% size(row1)
% mon_num+mon_num+(mon_num-1)*mon_num/2
    temp(2*mon_num+1:mon_num+mon_num+(mon_num-1)*mon_num/2)=row1;
   temp_m=temp'*temp;
   A(i,1:pnum+nnum)=0;
   A(i,i)=-0.001;
   
   inv_num^2+pnum+nnum
     A(i,1+pnum+nnum:inv_num^2+pnum+nnum)=reshape(temp_m,1,inv_num^2);
     size(A)
    b(i)=1000;
    i=i+1;
end




end