%get_image_dist
%clear all;
%load('rand_air_example');

im=imread('aircraft_1.png');

sim=size(im);
sx=sim(1);
sy=sim(2);
meanx=sx/2;
meany=sy/2;
max_dir=max(sx,sy);
maximum=max_dir/2;
BW=edge(im,'canny');
BW=imdilate(BW,[1,1,1,1;1,1,1,1]);
imshow(BW)

% imshow(im)
% hh
 bs = bwboundaries(BW,'noholes');
l0=concat(bs,ones(size(bs)));
x1=l0(:,1);
y1=l0(:,2);
cnum=size(x1,1)

syms x y
sumerr=0;
VEC = monomials([x; y],[0 1 2 3 4 5 6]);
 e=0.001;
unit=1/(maximum*4);
im(:,:)=0;
%neg_hs=[];
pos_hs=[];
i=1;
%for i=35:354
    i
    im(:,:)=0;
    x2=[];
    y2=[];
   % coeffs=neg_new(:,i);
   for j=1:cnum
    aerr=0;
    xs=l0(j,1);

    ys=l0(j,2);

    csizes=size(xs);
    xss=(xs-(ones(csizes)*meanx))/maximum;
    yss=(ys-(ones(csizes)*meany))/maximum;
    csizes=csizes(1);
 

for i=1:csizes%cluster_size:size(xs)-cluster_size

 for xjj=-1*unit:unit:1*unit
     for yjj=-1*unit:unit:1*unit
    xss(i)=xss(i)+xjj;
    yss(i)=yss(i)+yjj;
   eval=(subs(VEC,[x,y],[xss(i),yss(i)]));
    % eval=(subs(VEC,[x,y],[l0_new(t,1),l0_new(t,2)]));
    evald=double(eval);
  
     err=power(evald'*coeffs,1);
   
      if((err^2<=e^2 ))%381 381
          
         j;
        i=1;
        x2=[x2,xs(i)];
        y2=[y2,ys(i)];
       im(xs(i),ys(i))=255;
       
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

     eval=(subs(VEC,[x,y],[indx1(i),indy1(i)]));
     evald=double(eval);

    err=coeffs'*evald;
      if(err^2<=e^2)
    %     eee_chosen(countc,1)=weights(i);
%         eee_chosen(countc,2)=(err/(e^2));
%         countc=countc+1;
         im(indx(i),indy(i))=255;
          x2=[x2,indx(i)];
        y2=[y2,indy(i)];

      end
     end






%imshow(im);
% if size(x2,1)==0
%     continue
% end
[hd D] = HausdorffDist([x1,y1],[x2',y2']);
%neg_hs=[neg_hs,hd];
%end