% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   clearvars -except invs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  % a1=invs(1);a2=invs(2);a3=invs(3);
% % % %5
tic
diary
diary('mnew90000.txt');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SOSTOOLS-3.300');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/pics');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/sedumi-master');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/poly_helpers');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/glop/gloptipoly3/');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SOSTOOLS-3.300/internal');
addpath('/home/research/g.dehghanpoor/libraries/CSDP/Csdp-6.1.1/solver');
addpath('/cluster/cloud/mpich2/gnu-mpd/bin');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Examples');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/HSDSolver');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/HSDSolver/etc');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun/Oldfiles');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun/mexfun71');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun/pre7.5');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun_old');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun_old/Oldfiles');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Mexfun_old/mexfun71');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/Solver/Oldmfiles');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/dimacs');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SDPT3-4.0/sdplib');
addpath('/home/research/g.dehghanpoor/libraries/sdpa-c/mex')
addpath('/home/research/g.dehghanpoor/libraries/sdpa-c');
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mfiles')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mex')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V200')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V210SubPrograms')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V260SubPrograms')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mfiles/Solvers')


profile -memory on
profile clear

    clear all;clear classes;
 
run_inv=0;

negatives_temp=[];


mon_num=15; 
file='rrlimages-7.jpg';
 sample_num=0; 
 cnum=79;
im = imread(file);;
sim=size(im);
sx=sim(1);
sy=sim(2);
max_dir=max(sx,sy);
maximum=max_dir/2;
im=rgb2gray(im);
BW=edge(im,'canny');


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

 maximum=maximum/1;
    l0_new=l0_new/(maximum);
sl0=size(l0);
l0_s=size(l0_new); 

 xs=l0_new(1:m,1);
 ys=l0_new(1:m,2);

lcsize=floor(dsize/cnum);

cnum=size(xp)
cnum=cnum(2)

[coeffs,ws,exp_mat,R1,mv,l00,l0_new,mu,w,ee,c,indxx1,indyy1,cnum,cluster_num]=write_Gams_levelset2(im,file,0,1,csize,dsize,0);
toc
profile report
