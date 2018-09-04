function exp=learn_inv(pnum,nnum,poss,negs)
% pnum=4500;
% nnum=4500;
mon_num=28;
inv_num=2*mon_num+(mon_num-1)*mon_num/2;
[A,b]=generate_coeffs(pnum,nnum,mon_num,inv_num,negs,poss);
%%%learn second degree inv

c=zeros(inv_num^2+pnum+nnum,1);
c(1:pnum)=100;
c(pnum+1:pnum+nnum)=1;
K.l=pnum+nnum;
K.s=inv_num;
% size(A)
% size(b)
% size(c)
% K.l
% K.s
[x,y,info]=sedumi(A,b,c',K);
invs=reshape(x(pnum+nnum+1:pnum+nnum+inv_num^2),inv_num,inv_num);
elems=chol(invs);
exp=sum(elems,1)/(inv_num);

maxp=0;
mine=0;

while(mine-maxp<=0.1)
R=exp;%rmvnrnd(exp,invs,1);
%R(mon_num+5)=100;
R2=reshape((R'*R),inv_num^2,1);
errors=A(:,1+pnum+nnum:pnum+nnum+inv_num^2)*R2;


perr=-errors(1:pnum);
nerr=errors(pnum+1:pnum+nnum);
%minp=min(perr);
maxp=mean(perr);
mine=mean(abs(nerr));
save('errors.mat','perr','nerr');
break;
end

