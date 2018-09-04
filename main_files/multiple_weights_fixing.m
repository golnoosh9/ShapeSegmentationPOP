% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   clearvars -except invs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  % a1=invs(1);a2=invs(2);a3=invs(3);
% % % %5
%profile -memory on
%profile clear
%
tic
for i=1:1
ttt=0;
if (ttt==0)
    clearvars -except begin
mset clear;
mset clearmeas;  
run_inv=0;
global yyy 
global MMM
negatives_temp=[];
tic

mon_num=15; 
%file='rreal1.jpg';
file='aircraft_1test.png';
 sample_num=0; 
 cnum=79;
im = imread(file);;
sim=size(im);
sx=sim(1);
sy=sim(2);
max_dir=max(sx,sy);
maximum=max_dir/2;2
%im=rgb2gray(im);   

% imr=im(:,:,1);
% img=im(:,:,2);
% im=(im2double(imr)./im2double(imr+img));

BW=edge(im,'canny');
%imshow(BW);

varsum=0;
 bs = bwboundaries(BW,'noholes');
 
l0=concat(bs,ones(size(bs)));
segs=size(bs);
for j=1:segs
    temp=size(cell2mat(bs(j)));
    sarray(j)=temp(1);
end
dsize=sarray*ones(segs);
cnum=dsize;

csize=ceil(dsize/cnum);
[xp,yp]=cluster_neighbors(bs,cnum,csize);

im(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
end


meanx=sx/2;
meany=sy/2;
nn=ones(size(l0));
h=[meanx;meany];
nn(:,1)=nn(:,1).*meanx;
nn(:,2)=nn(:,2).*meany;

l0_new=l0-nn;
tl0=l0;
tl0x=tl0(:,1);
tl0y=tl0(:,2);
m=size(l0);%total points
m=m(1)
cluster_size=1;
% maximum=0;
% for a=1:size(l0_new)
%      norm_factor=max(abs(l0_new(a,1)),abs(l0_new(a,2)));;%normalize data
%       if(norm_factor>maximum)
%           maximum=norm_factor;
%       end
% end
 maximum=maximum/1;
    l0_new=l0_new/(maximum);
sl0=size(l0);
l0_s=size(l0_new); 

 xs=l0_new(1:m,1);
 ys=l0_new(1:m,2);

lcsize=floor(dsize/cnum);

syms x y
sumerr=0;
VEC = monomials([x; y],[0 1 2 3 4 5 6]);
%[coeffs,w,exp_mat,r1,mv,ll0,ll0_new,mu,ww,ee,cc]=
%[coeffs,e,miu,ws,prog]=dsos_test(im,file,0,0,csize,cnum,dsize);x
cnum=size(xp)
cnum=cnum(2)

[coeffs,ws,exp_mat,R1,mv,l00,l0_new,mu,w,ee,c,indxx1,indyy1,cnum,cluster_num]=write_Gams_levelset2_arplane(im,file,0,1,csize,dsize,0,begin);
toc

%gggg
% 

% for i=1:100
% %     xxx=full(basisSupports{i});
% %     mmm=size(xxx);
% %     bbb(i)=mmm(1);
% yyyy=full(varList(:,1));x
% % b(i-999,1:755)=xxx;     
% end
% % % 
% 
% 
% % % uuuuuuu
% % 
% % %%%%TODO: campare with wb multiply e
% % % % mon_num=6
% % % % segsize=cnum;
% % %  %mw=exp_mat(mon_num+1:mon_num+cnum);
% % % % meanw=0.5;
% % % % new_mw=mw-ones(size(mw))*meanw;
% % % % inds=find(new_mw<=0);
% % % % new_mw(inds)=0;
% % % % inds=find(new_mw>0);
% % % % new_mw(inds)=1;
% % % % im2=cluster_show(file,new_mw,lcsize,cnum,dsize);
% % % 
% % % % imshow(im2);
% % % 
% % % % 
% % % 
% % % % coeffs=c;
% % % % exp_mat=y;
% % % % mom_mat=double(mmat(mu));
% % % % cov_mat=mom_mat;
% % % % cov_mat(1,:)=[];
% % % % cov_mat(:,1)=[];
% % % % [V,D]=eig(cov_mat);
% % % % %recon=V*D*V';
% % % % Dsqr=sqrt(D);
% % % % U=V*Dsqr;
% % % % %exp_mat=U(:,72);
% % % % coeffs=exp_mat(1:6);
% % % % e=0.2
% % % % % mon_num=6
% % % % % load('affterr.mat');
% % % 
% % % % clear coeffs;
% % % %  coeffs=POP.xVectL(1:5);
% % % %  coeffs(6)=POP.xVect(4);
% % % %  coeffs(4)=POP.xVect(4);
% % e=0.01%POP.xVect(cnum+5+1)%POP.xVect(5+cnum+1)%exp_mat(mon_num+cnum+1)
% % % % % % % % % % % % %im2=cluster_show(file,new_mw,lcsize,cnum,dsize);
% % % % % % % % % % % imshow(im2);
% % % % % % % % % % % figure;
% % % l0=concat(bs,ones(size(bs)));
% % % im5=imread(file);;
% % % im5(:,:)=0;
% % % for(i=1:s2(1))
% % %     im5(l0(i,1),l0(i,2))=255;
% % % end
% % % clear im;
% % % for(i=1:s2(1))
% % %     im(l0(i,1),l0(i,2))=255;
% % % end
% % % im5=imread(file);
% % % 
% % % %coeffs=exp_mat(1:mon_num);
% % % sl0=size(l0);err=0;
% % % clear added;
% % % syms x y
% % % im3=im;
% % % tim=(im);
% % %  im(:,:,3)=0;
% % %  im(:,:,2)=0;
% % % VEC = [1;x;y;x^2;y^2];
% % % %VEC = monomials([x; y],[0 1 2]);count=1;
% % % for j=1:cnum
% % %     aerr=0;
% % %     xs=l0(j,1);
% % % 
% % %     ys=l0(j,2);
% % % 
% % %     csizes=size(xs);
% % %     xss=(xs-(ones(csizes)*meanx))/maximum;
% % %     yss=(ys-(ones(csizes)*meany))/maximum;
% % %     csizes=csizes(1);
% % %  
% % %    
% % % for i=1:csizes%cluster_size:size(xs)-cluster_size
% % % 
% % %    eval=(subs(VEC,[x,y],[xss(i),yss(i)]));
% % %     % eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
% % %      evald=double(eval);
% % %      err=power(evald'*coeffs,1);
% % % %     if (j==150)
% % % %        im() 
% % % %     end
% % %      if(err<=e && err>=-e)%381 381
% % %          j
% % %        im(xs(i),ys(i),2)=255;
% % % 
% % %      end
% % % end
% % % end
% % % % 
% % % 
% % % inds=find(im3==0);
% % % berr=0;
% % % [indx,indy]=ind2sub(size(tim),inds);
% % % indx1=(indx-ones(size(indx))*meanx)/(maximum);
% % % indy1=(indy-ones(size(indy))*meany)/(maximum);
% % % mmy=max(indy)
% % % my=min(indy)
% % % s=size(indx);
% % % for i=1:s(1)
% % %      eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
% % %      evald=double(eval);
% % %      err=power(evald'*coeffs,1);
% % %      if(err<=e &&err>=-e)
% % %          im(indx(i),indy(i),2)=255;
% % %         % disp('llll');
% % %        % added(count,:)=l0(t,:);
% % %        % count=count+1;
% % %      end
% % % end
% % % 
% % % 
% % % imshow(im);
% % % hhh
% % % % mom_mat=double(mmat(mu));
% % % % cov_mat=mom_mat;
% % % % cov_mat(1,:)=[];
% % % % cov_mat(:,1)=[];
% % % % [V,D]=eig(cov_mat);
% % % %recon=V*D*V';
% % % % Dsqr=sqrt(D);
% % % % U=V*Dsqr;
% % % % recon=U*U';
% % % % nnnnn
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%Sparse
%      load('middle.mat');
     %
   
     load('affterr.mat');
   
 

% mon_num=15;
% numvars=cnum+mon_num+vercount-8;
% wfixed=90;
% moments=SDPinfo.y;
% prob_w=moments(wfixed+mon_num);
% weights=POP.xVect(mon_num+1:mon_num+cnum);
% 
% msup=full(xIdxVec);
% smom=size(moments);
% msup(1,:)=[];
% degs=msup*ones(numvars,1);
% sec_terms=find(degs==2);
% first_terms=find(degs==1);
% % wfixed_column=msup(:,wfixed+mon_num-2);
% % w_nzinds=find(wfixed_column==1);
% % w_secterms=intersect(w_nzinds,sec_terms);
% % clique_num=1;
% % clear coeff_n;
% % for i=1:mon_num-2
% %     i
% %      wfixed_mon=msup(:,i);
% %      w_nzinds=find(wfixed_mon==1);
% %     temp=find(wfixed_column==1);
% %     w_mon1=intersect(w_nzinds,sec_terms);
% %     w_mon=intersect(w_mon1,temp);
% %     w_mon=w_mon(1);
% %     coeff_n(i)=moments(w_mon)/weights(wfixed);;
% % end
% % coeff_n(4)=1;
% % coeff_n(5)=1;
% % count=6;
% % coeff=coeff_n;
%  %weights=POP.xVect(5+cnum+1:5+cnum+bnum);
% % for i=1:mon_num
% %     for j=i+1:mon_num
% %         i 
% %         j
% %         wcol1=msup(:,i);
% %         wcol2=msup(:,j);
% %         terms=sec_terms;
% %         wtemp=find(wcol1==1);
% %         wtemp2=find(wcol2==1);
% %         winters=intersect(wtemp,wtemp2);
% %         if i>=4
% %             winters=wtemp2;
% %             terms=first_terms;
% %         end
% %         if j>=4
% %             terms=first_terms;
% %             winters=wtemp;
% %         end
% %         winters=intersect(winters,terms);
% %         if i>=4 && j>=4
% %             coeff_n(count)=1;
% %         else
% %             winters(1)
% %             coeff_n(count)=moments(winters(1));
% %         end
% %         count=count+1;
% %     end
% % end
% 
% %plot_implicit_2(coeff_n);
% 
% %coeffs=coeff_n';;
% weights_ds=POP.xVect(mon_num+cnum+vercount+3:mon_num+cnum+vercount+cnum+2);

   coeffs=POP.xVect(1:mon_num);
   get_img_distance;
   save('coeffs.mat','coeffs');
  %plot_implicit_six(coeffs);
  toc
return

 % pppp
if(run_inv==1)
adddd;
add_negative;
end
end

end
 save('coeffs.mat','coeffs');
 
 e=POP.xVect(mon_num+cnum+1);
%e=0.008;
i_epsil=e;
%i_epsil=i_epsil+0.003;

 e=POP.xVect(mon_num+cnum+1);
l0=concat(bs,ones(size(bs)));
im5=imread(file);;

image(im5)



%im5=rgb2gray(im5);
im5(:,:)=0;
s2=size(l0);
for(i=1:s2(1))
    im5(l0(i,1),l0(i,2))=255;
end
im=imread(file);;
im5=imread(file);;
im(:,:)=0;
for(i=1:s2(1))
    im(l0(i,1),l0(i,2))=255;
    
end


%coeffs=exp_mat(1:mon_num);
sl0=size(l0);err=0;
clear added;
syms x y
im3=im;
tim=(im);
 im(:,:,3)=0;
 im(:,:,2)=0;


VEC = [1;x;y;x^2;x*y;y^2];
%VEC = [1;x;y];
countc=1;
countn=1;
clear eee_chosen;
clear eee_nchosen
eee=zeros(cnum,3);
%get_higher_coeffs;
%e=e*3
VEC = monomials([x; y],[0 1 2 3 4]);count=1;
unit=1/(maximum*4);
meanx=sx/2;
meany=sy/2;
for j=1:cnum
    aerr=0;
    xs=l0(j,1);

    ys=l0(j,2);

    csizes=size(xs);
    xss=(xs-(ones(csizes)*meanx))/maximum;
    yss=(ys-(ones(csizes)*meany))/maximum;
    csizes=csizes(1);
 

for i=1:csizes%cluster_size:size(xs)-cluster_size
    
%     statem=VEC'*coeffs;
%     Derx=diff(statem,'x');
%     Dery=diff(statem,'y');
%     DerNorm=Derx^2+Dery^2;
%     Dereval=subs(DerNorm,[x,y],[xss(i),yss(i)]);

 for xjj=-1*unit:unit:1*unit
     for yjj=-1*unit:unit:1*unit
    xss(i)=xss(i)+xjj;
    yss(i)=yss(i)+yjj;
   eval=(subs(VEC,[x,y],[xss(i),yss(i)]));
    % eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
    evald=double(eval);
    
    if(l0(j,2)<=sy)%%%% point's y greater than 2
             neighbor_eval=(subs(VEC,[x,y],[(l0(j,1)-meanx)/maximum,(l0(j,2)+1-meany)/maximum]));
              neighbor_evald=double(neighbor_eval);
              
%       eee(j)=(neighbor_evald'*coeffs-evald'*coeffs)^2  ;     
    end
%     
%     create_coef;
%     evals_total=[evald(1)^2;evald(2)^2;evald(3)^2;evald(4)^2;evald(5)^2;evald(6)^2];
%     evals_total=[evals_total',coefs_sc];
    
 %err=coeffs_total*evals_total'
    %    err=coeffs'*evald
%     mcoef=zeros(15,1);
%       mcoef(1)=evald(1)^2;
%     mcoef(2)=evald(2)^2;
%     mcoef(3)=evald(3)^2;
%     mcoef(4)=evald(4)^2;
%     mcoef(5)=evald(5)^2;
%     coef_mat=evald*evald';
% coef_mat=coef_mat*2;
% row1=coef_mat(1,2:5);
% row2=coef_mat(2,3:5);
% row3=coef_mat(3,4:5);
% row4=coef_mat(4,5);
% mcoef(6:15)=[row1,row2,row3,row4];
% coeff_n*mcoef
%plot_implicit_six(coeffs);

     err=power(evald'*coeffs,1);
   
      if((err^2<=e^2 ))%381 381
          
         j;
        i=1;
       im(xs(i),ys(i),2)=255;
       im5(xs(i),ys(i),1:3)=0;
       
       break;
       % plot(ys(i),xs(i),'b.','MarkerSize',15)
     else
      
     end
     
     end
 end
end
end


inds=find(im==0);
berr=0;
[indx,indy]=ind2sub(size(im),inds);
indx1=(indx-ones(size(indx))*meanx)/(maximum);
indy1=(indy-ones(size(indy))*meany)/(maximum);
mmy=max(indy);
my=min(indy);
s=size(indx);
mcoef=zeros(15,1);

%i_epsil=i_epsil/2;
for i=1:s(1)
%      for xjj=-1*unit:unit:1*unit
%      for yjj=-1*unit:unit:1*unit
%     indx1(i)=indx1(i)+xjj;
%     indy1(i)=indy1(i)+yjj;
     eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
     evald=double(eval);
%      statem=VEC'*coeffs;
%     Derx=diff(statem,'x');
%     Dery=diff(statem,'y');
%     DerNorm=Derx^2+Dery^2;
%     Dereval=subs(DerNorm,[x,y],[indx1(i),indy1(i)]);
%     err=power(evald'*coeffs,1);
%      mcoef(1)=evald(1)^2;
%     mcoef(2)=evald(2)^2;
%     mcoef(3)=evald(3)^2;
%     mcoef(4)=evald(4)^2;
%     mcoef(5)=evald(5)^2;
%     coef_mat=evald*evald';
% coef_mat=coef_mat*2;
% row1=coef_mat(1,2:5);
% row2=coef_mat(2,3:5);
% row3=coef_mat(3,4:5);
% row4=coef_mat(4,5);
% mcoef(6:15)=[row1,row2,row3,row4];
%      

      
%     create_coef;
%     evals_total=[evald(1)^2;evald(2)^2;evald(3)^2;evald(4)^2;evald(5)^2;evald(6)^2];
%     evals_total=[evals_total',coefs_sc];
 
% err=coeffs_total*evals_total';
    err=coeffs'*evald;
      if(err^2<=e^2)
    %     eee_chosen(countc,1)=weights(i);
        eee_chosen(countc,2)=(err/(e^2));
        countc=countc+1;
         im(indx(i),indy(i),2)=255;
         im5(indx(i),indy(i),1:3)=0;
        % plot(indy(i),indx(i),'b.','MarkerSize',15)
       %  eee_chosen(countc,1)=weights(j);
       % eee_chosen(countc,2)=coeff_n*mcoef;
       % eee_chosen(countc,3)=(e^2-(weights(j)*err^2));
       % countc=countc+1;
        % disp('llll');
       % added(count,:)=l0(t,:);
       % count=count+1;
       
       
     else
         %coeff_n*mcoef
%         eee_nchosen(countn,1)=weights(i);
       eee_nchosen(countn,2)=(err/(e^2));
%        % eee_nchosen(countn,3)=(e^2-(weights(j)*err^2));
%         countn=countn+1;
      end
     end
%      end
%      i
% end


imshow(im);
cluster_num
 hhh
% hhhh






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

segsize=cnum;
mv_num = double(mv);
exp_mat=POP.xVect(1:mon_num);%%% just sample coefficients
get_higher_coeffs;
cov_mat=cpower2;
R1 = rmvnrnd(exp_mat',cov_mat,1);
plot_implicit(R1);
llll
[V,D]=eig(cov_mat);
     
mon_num=6;
 neg=1;count=0;

coeffs=zeros(6,1);
coeffs(1:mon_num)=exp_mat(1:mon_num);
sl0=size(l0);err=0;
syms x y
VEC = monomials([y; x],[0 1 2]);count=1;
 A1=eye(mon_num+cnum+7);
    A2=eye(mon_num+cnum+7).*(-1);
     B1=ones(mon_num+cnum+7,1);
    B2=zeros(mon_num+cnum+7,1);
    A3=eye(mon_num+cnum+2);
    A4=eye(mon_num+cnum+2).*(-1);
     B3=ones(mon_num+cnum+2,1);
    B4=zeros(mon_num+cnum+2,1);
    for i=1:mon_num
       A1(i,i)=0;
       A2(i,i)=0;
       B1(i,1)=0;
       B2(i,1)=0;
       A3(i,i)=0;
       A4(i,i)=0;
       B3(i,1)=0;
       B4(i,1)=0;
    end
    for i=mon_num+cnum+1:mon_num+cnum+1
       A1(i,i)=0;
       A2(i,i)=0;
       B1(i,1)=0;
       B2(i,1)=0;
       A3(i,i)=0;
       A4(i,i)=0;
       B3(i,1)=0;
       B4(i,1)=0;
    end
   
%      B1(mon_num+cnum+1)=0.01;%%for e
%      B2(mon_num+cnum+1)=0;
%      B1(mon_num+cnum+2)=1;%%for miu
%      B2(mon_num+cnum+2)=0.1;
    
Ws=zeros(cnum,1);    
    
    while(true)
R1 = rmvnrnd(exp_mat',cov_mat,1);
mw=R1(mon_num+1:mon_num+cnum);
meanw=mean(mw);
new_mw=mw-ones(size(mw))*meanw;
inds=find(new_mw<=0);
new_mw(inds)=0;
inds=find(new_mw>0);
new_mw(inds)=1;
im2=cluster_show(file,new_mw,lcsize,cnum,dsize);
imshow(im2);
coeffs(1:mon_num)=R1(1:mon_num);
Ws=R1(mon_num+1:mon_num+cnum);
    end
mmmm
% clear im;
% for(i=1:s2(1))
%     im(l0(i,1),l0(i,2))=255;
% end
% 
% im4=imread(file);
% im4(:,:)=0;
% %coeffs=R1(1:mon_num);
% sl0=size(l0);err=0;
% clear added;
% syms x y
% im3=im;
% tim=(im);
%  im(:,:,3)=0;
%  im(:,:,2)=0;
% VEC = monomials([x; y],[0 1 2]);count=1;
% for t=1:sl0(1)
%      eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
%      evald=double(eval);
%      err=power(evald'*coeffs,1);
%      if(err<=R1(cnum+mon_num+1) && err>=0)%381 381
%        im4(l0(t,1),l0(t,2))=255;
%        % added(count,:)=l0(t,:);
%         count=count+1;
%      end
% end
% 
% % inds=find(im3==0);
% % berr=0;
% % [indx,indy]=ind2sub(size(tim),inds);
% % indx1=(indx-ones(size(indx))*meanx)/sqrt(maximum);
% % indy1=(indy-ones(size(indy))*meany)/sqrt(maximum);
% % mmy=max(indy)
% % my=min(indy)
% % s=size(indx);
% % for i=1:s(1)
% %      eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
% %      evald=double(eval);
% %      err=power(evald'*coeffs,1);
% %      if(err<=0.01 && err>=-0.01)
% %          im(indx(i),indy(i),2)=255;
% %         % disp('llll');
% %        % added(count,:)=l0(t,:);
% %         count=count+1;
% %      end
% % end
% im4=rgb2gray(im4);
% imshow(im4);
% BW=edge(im4,'canny');
% bs = bwboundaries(BW);
% bs=size(bs);
% bs=bs(1)
% mw=R1(mon_num+1:mon_num+cnum);
% meanw=mean(mw);
% new_mw=mw-ones(size(mw))*meanw;
% inds=find(new_mw<=0);
% new_mw(inds)=0;
% inds=find(new_mw>0);
% new_mw(inds)=1;
% im2=cluster_show(file,new_mw,lcsize,cnum,dsize);
% if(bs==1)
%     break;
% end
% 
% end
% 
% 
% mw=R1(mon_num+1:mon_num+cnum);
% meanw=mean(mw);
% new_mw=mw-ones(size(mw))*meanw;
% inds=find(new_mw<=0);
% new_mw(inds)=0;
% inds=find(new_mw>0);
% new_mw(inds)=1;
% im2=cluster_show(file,new_mw,lcsize,cnum,dsize);
% imshow(im2);
% 
% 
% vvv
% 
% clear added; 
% while(true)
%     disp('hereee');
%     sumerr=0;
%     A1=ones(mon_num+cnum+52+4,1);
%     A2=ones(mon_num+cnum+52+4,1).*(-1);
%     for i=1:mon_num
%        A1(i,1)=0;
%        A2(i,1)=0;
%     end
%     for i=mon_num+cnum+52+1:mon_num+cnum+52+4
%        A1(i,1)=0;
%        A2(i,1)=0;
%     end
% R1 = mvnrnd(exp_mat',cov_mat,1,[A1;A2],[A1,zeros(mon_num+cnum+52+4,1)]);
% 
% coeffs=zeros(mon_num,1);
% coeffs(1:mon_num)=R1(1:mon_num);
% sr=size(R1);
% ebound=R1(sr(1));
% ws=R1(mon_num+1:mon_num+segsize(1));
% ws=(ws)-ones(size(ws))*0.5;
% inds=find(ws<=0);
% ws(inds)=0;
% inds=find(ws>0);
% ws(inds)=1;
% img=cluster_show(file,ws,lcsize,cnum,dsize);
% % if(img==0)
% %     continue;
% % end
% % sl0=size(img);err=0;
% % sl0=size(l0);err=0;
% % syms x y
% % VEC = monomials([y; x],[0 1 2]);count=1;
% % im(:,:,3)=0;
% % im(:,:,2)=0;
% % disp('uuuuu');
% % VEC = monomials([y; x],[0 1 2]);count=1;
% % for t=1:sl0(1)
% %      eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
% %      evald=double(eval);
% %      err=power(evald'*coeffs,1);
% %      if(err<=R1(mon_num+cnum+1)&& err>=-R1(mon_num+cnum+1) )
% %          im(l0(t,1),l0(t,2),2)=255;
% %          
% %         added(count,:)=l0(t,:);
% %         count=count+1;
% %      end
% % end
% % if(count==1)
% %     continue;
% % end
% % 
% % inds=find(im3==0);
% % berr=0;
% % [indx,indy]=ind2sub(size(tim),inds);
% % indx1=(indx-ones(size(indx))*meanx)/sqrt(maximum);
% % indy1=(indy-ones(size(indy))*meany)/sqrt(maximum);
% % mmy=max(indy)
% % my=min(indy)
% % s=size(indx);
% % for i=1:s(1)
% %      eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
% %      evald=double(eval);
% %      err=power(evald'*coeffs,1);
% %      if(err<=R1(mon_num+cnum+1)&& err>=-R1(mon_num+cnum+1) )
% %          im(indx(i),indy(i),2)=255;
% %         % disp('llll');
% %        % added(count,:)=l0(t,:);
% %         count=count+1;
% %      end
% % end
% % im(:,:)=0;
% % s2=size(added);
% % for(i=1:s2(1))
% %     im(added(i,1),added(i,2))=255;
% % end
% im2=img;
% imshow(im2);
% disp('uummaadd');
% BW=edge(im2,'canny');
% bs = bwboundaries(BW);
% bs=size(bs);
% bs=bs(1)
% % o=sum(ws);
% % sws=size(ws);
% % d=o*csize;
% % err=err/d;
% % im=cluster_show(file,ws,lcsize,cnum,dsize);
% % imshow(im);
% clear added;
% if(bs==1 )
%     break;
% end
% end
% mmmm
% im=cluster_show(file,ws,lcsize,cnum,dsize);
% imshow(im);
% pp
% % if(o==0)
% %     continue;
% % end
% % if(o==sws(2))
% %     continue;
% % end
% % innum=sws(2)-o;
% % present_dots=cluster_concat(l0_new,ws,lcsize,cnum,dsize);
% % k=convhull(present_dots(:,1),present_dots(:,2));
% % inds=find(ws==0);
% % sind=size(inds);
% % taccept=0;
% % for cc=1:cnum
% %     if(ws(cc)==1)
% %         continue;
% %     end
% %     
% % sample=cc;
% % if(sample~=cnum)
% % inpointsx=(l0_new((sample-1)*lcsize+1:sample*lcsize,1));
% % inpointsy=(l0_new((sample-1)*lcsize+1:sample*lcsize,2));
% % else
% %     sl0=size(l0_new);
% %     inpointsx=(l0_new((sample-1)*csize+1:sl0(1),1));
% %     inpointsy=(l0_new((sample-1)*csize+1:sl0(1),2));
% % end
% % sinp=size(inpointsx);
% % if(sinp(1)==0)
% %     continue;
% % end
% % ins=inpolygon(inpointsx,inpointsy,present_dots(k,1),present_dots(k,2));
% % sins=size(ins);
% % if(sins(1)==0)
% %     bbbbbb
% % end
% % accept=sum(ins)/sins(1);
% % taccept=accept+taccept;
% % end
% % 
% % xcluster=reshape(l0(1:cnum*lcsize,1),lcsize,cnum);%% have values l0 so we can iterate on pixels
% % ycluster=reshape(l0(1:cnum*lcsize,2),lcsize,cnum);
% % %xcluster(lcsize+1:sl0(1),cnum)=l0(cnum*lcsize+1:sl0(1),1);
% % sel_inds=find(ws==1);
% % xsel=xcluster(sel_inds);
% % ysel=ycluster(sel_inds);%%%% clusters with weights 1
% % minx=min(min(xsel));
% % maxx=max(max(xsel));
% % nopex=0;
% % for x=minx:maxx
% %    foundx=find(xsel==x);
% %    sfound=size(foundx);
% %    if(sfound(2)==2)
% %        nopex=nopex+1;
% %    end
% %     
% % end
% % 
% % miny=min(min(ysel));
% % maxy=max(max(ysel));
% % nopey=0;
% % for y=miny:maxy
% %    foundy=find(ysel==y);
% %    sfoundy=size(foundy);
% %    if(sfoundy(2)>=1)
% %        nopey=nopey+1;
% %    end
% %     
% % end
% % 
% % c=R1(1:mon_num);
% % smooth=0;count=1;maxe=-10000000;
% % 
% % 
% % 
% % taccept
% % d
% % 
% % if(d>=500  &&d<=900&&taccept>=0.004*innum)%%%R1(mon_num+cnum+2)=>eres:dist of rectangle
% %     break;%%cond::d>=500  &&d<=900&&ebound<=-100 &&nopex>=(maxx-minx-15)&&nopey>=(maxy-miny-25)&&taccept>=0.5*innum
% % end
% % end
% % % end
% % 
% % im=cluster_show(file,ws,lcsize,cnum,dsize);
% % imshow(im);
% % 
% % figure;
% % o=sum(ws);
% % sws=size(ws);
% % d=o*csize;
% % 
% % innum=sws(2)-o;
% % 
% % new_w=edit_weights_woinv(ws,l0,l0_new,cnum,csize,dsize,lcsize,file,10,innum);
% % im=cluster_show(file,new_w,lcsize,cnum,dsize);
% % imshow(im);
% % figure;
% % sum(new_w)
% 
% [innew_w,ds,tws,is]=edit_weights_BFS_prior(new_w,l0,l0_new,cnum,csize,dsize,lcsize,file,10,innum,invs);
%  im=cluster_show(file,innew_w,lcsize,cnum,dsize);
%  figure
%  imshow(im);
% bbbb
% 
% %%%%%%%%%%%
% %%%%%5second stage
% %%%%%%%%%%%%%
% file='H2.png';
% cnum=20;
%  im=imread(file);
% im=rgb2gray(im);
% BW=edge(im,'canny');
% bs = bwboundaries(BW);
% l0=concat(bs,ones(size(bs)));
% im(:,:)=0;
% s2=size(l0);
% for(i=1:s2(1))
%     im(l0(i,1),l0(i,2))=255;
% end
% meanx=mean(l0(:,1));
% meany=mean(l0(:,2));
% nn=ones(size(l0));
% h=[meanx;meany];
% nn(:,1)=nn(:,1).*meanx;
% nn(:,2)=nn(:,2).*meany;
% l0_new=l0-nn;
% m=size(l0);%total points
% m=m(1)
% cluster_size=1;
% maximum=0;
% for a=1:size(l0_new)
%      norm_factor=(power(l0_new(a,1),2)+power(l0_new(a,2),2));%normalize data
%       if(norm_factor>maximum)
%           maximum=norm_factor;
%       end
% end
%     l0_new=l0_new/sqrt(maximum);
% sl0=size(l0);
% l0_s=size(l0_new);
% r=randperm(sl0(1));
% r=r(1:m);
%  xs=l0_new(r,1);
%  ys=l0_new(r,2);
% segs=size(bs);
% for j=1:segs
%     temp=size(cell2mat(bs(j)));
%     sarray(j)=temp(1);
% end
% dsize=sarray*ones(segs);
% csize=ceil(dsize/cnum);
% lcsize=floor(dsize/cnum);
% ws=ones(cnum,1);
% [coeffs,ws,exp_ma,R1,mv]=cluster_weight_deg4(file,0,ws,csize,cnum,dsize)
% mon_num=28;
% mv_num = double(mv);
% ss=size(mv);
% ss=ss(1);
% exp_mat=mv_num(1+1:mon_num+4);%%% 3 for miu,2 for e: (w.o. e is 1
% coef=exp_mat(1:mon_num);
% 
% cov_elem=mv_num(mon_num+5:ss);
% mv_cov=mv(mon_num+5:ss);
% cov_mat=ones(mon_num);
% cov_cont=1;
% m=3;
% for i=1:mon_num+m
%         
%         mv_cov(cov_cont:cov_cont+mon_num+m-i)
%         cov_mat(i,i:mon_num+m)=cov_elem(cov_cont:cov_cont+mon_num+m-i)';
%         cov_mat(i:mon_num+m,i)=cov_elem(cov_cont:cov_cont+mon_num+m-i)
%     cov_cont=cov_cont+mon_num+m-i+1;
% 
% end
% ccov=cov_mat(1:mon_num,1:mon_num);
% [V,D]=eig(ccov);
%  neg=1;count=0;
% 
% R1 = mvnrnd(exp_mat',cov_mat,1);
% newc=R1(1:mon_num);
% 
% plot_implicit(coef);
% bbbb
