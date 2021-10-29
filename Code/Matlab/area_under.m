function [result5, result10] = area_under(x,x2,f,f2,n,n2,m,path)

%% for 5 years

result5 = 0;
result10 = 0;
test1 = x(:) < n;
test2 = x2(:) < n;

if (ismember(0,test1(:)) == 0 || ismember(0,test2(:)) == 0)
    disp('Unable to calculate RMST for 5 years, no patients with such survival time for one of both of the groups')
    return
end

[~,indx] = (min(abs(x - n)));
after_tresholding = x(x <= n);
after_tresholding(indx) = n;
after_tresholding_f = f(f >= f(indx));
area1 = trapz(after_tresholding,after_tresholding_f);

[~,indx] = (min(abs(x2 - n)));
after_tresholding2 = x2(x2 <= n);
after_tresholding2(indx) = n;
after_tresholding_f2 = f2(f2 >= f2(indx));
area2 = trapz(after_tresholding2,after_tresholding_f2);

figure
subplot(2,1,1)
high5 = stairs(after_tresholding,after_tresholding_f);
hold on
low5 = stairs(after_tresholding2,after_tresholding_f2);
title(sprintf('%s for 5 years', m))
xlabel('Lifetime (days)')
ylabel('Survival probability')
legend([high5,low5],{sprintf('Above treshhold for %d patients', length(after_tresholding)), ...
    sprintf('Below treshhold for %d patients', length(after_tresholding2))}, 'Location', 'best')

result5 = area1 - area2;

%% for 10 years

test3 = x(:) < n2;
test4 = x2(:) < n2;

if (ismember(0,test3(:)) == 0 || ismember(0,test4(:)) == 0)
    disp('Unable to calculate RMST for 10 years, no patients with such time')    
    return 
end


[~,indx] = (min(abs(x - n2)));
after_tresholding = x(x <= n2);
after_tresholding(indx) = n2;
after_tresholding_f = f(f >= f(indx));
area1 = trapz(after_tresholding,after_tresholding_f);
[~,indx] = (min(abs(x2 - n2)));
after_tresholding2 = x2(x2 <= n2);
after_tresholding2(indx) = n2;
after_tresholding_f2 = f2(f2 >= f2(indx));
area2 = trapz(after_tresholding2,after_tresholding_f2);

subplot(2,1,2)
high5 = stairs(after_tresholding,after_tresholding_f);
hold on
low5 = stairs(after_tresholding2,after_tresholding_f2);
title(sprintf('%s for 10 years', m))
xlabel('Lifetime (days)')
ylabel('Survival probability')
legend([high5,low5],{sprintf('Above treshhold for %d patients', length(after_tresholding)), ...
    sprintf('Below treshhold for %d patients', length(after_tresholding2))}, 'Location', 'best')
savefig(append(path,m,'510'))

result10 = area1 - area2;

end

