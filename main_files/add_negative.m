function []=add_negative(begin,pos,POP,mon_num)
%neg_add=coeffs;
add=1;
%begin=0;
if begin==1
pos_new=[];
neg_new=[];
pos_count=1;
ex_count=1;

else
 %load('rand_four_example.mat');
 load('rand_air_example.mat');
end
  % pos=0;
%   neg_new(:,400:1218)=[];
  % pos_new(:,1:194)=[]];
 %  neg_new(:,100:183)=[];
%    pos_new(:,20:79)=[];
 %   neg_new(:,40:45)=[];
 %  neg_new(:,2)=[];
% %  ex_count=1;
% %  pos_count=1;
% %  neg_new=[];
% %  pos_new=[];5
% %  %load('invariant.mat')
if add==1
if pos==1
     pos_new(:,pos_count)=POP.xVect(1:mon_num);
     pos_count=pos_count+1;
else

            neg_new(:,ex_count)=POP.xVect(1:mon_num);
            ex_count=ex_count+1;
end
end
%   pos_new(:,163:166)=[];
% neg_new(:,513:517)=[];
% neg_new(:,455:459)=[];
% % pos_new(:,163:164)=[];
%    pos_new(:,1:12)=[];
%    pos_new(:,100:107)=[];
%    pos_new(:,150:168)=[];
% neg_new(:,1:188)=[];
% % % % % % % pos_new(:,163:166)=[];
% % % % % % % neg_new(:,513:517)=[];
% % % % % % % neg_new(:,455:459)=[];
% % % % % % % % pos_new(:,163:164)=[];
% % % % % % %    pos_new(:,1:12)=[];
% % % % % % %    pos_new(:,100:107)=[];
% % % % % % %    pos_new(:,150:168)=[];
% % % % % % % neg_new(:,1:188)=[];

%neg_new(:,8:15)=[];
% %pos_new(:,1)=[];    
%neg_new(:,70:84)=[];
%  pos_new(:,40:79)=[];
% % %        pos_new(:,13)=[];
% % % %    
%    neg_new(:,150:226)=[];
% %    neg_new(:,208)=[];
% %    neg_new(:,209)=[];
% %    neg_new(:,210)=[];
%    %count=size([206:210]);
%      %neg_new(:,20:226)=[];
  % save('rand_four_example.mat','pos_new','neg_new','ex_count','pos_count');
    save('rand_air_example.mat','pos_new','neg_new','ex_count','pos_count');
   hhhhh
% neg_new(:,60:70)=[];

pnum=size(pos_new);
nnum=size(neg_new)
inv_coeffs=learn_inv(pnum(2),nnum(2),pos_new,neg_new);
%save('four_random_inv','inv_coeffs');
save('airplane_random_inv','inv_coeffs');
