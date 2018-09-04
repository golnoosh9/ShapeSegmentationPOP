begin=0;
%begin=0;
for i=1:1
    multiple_weights_fixing;
    figure
    plot_implicit_six(coeffs)
     pos=0;
    if hd<5
        pos=1;
    end
    add_negative(begin,pos,POP,mon_num);
    begin=0;
%     
end

hhh

for j=1:1+50
    plot_implicit(neg_new(:,j));
    figure
    
    
end