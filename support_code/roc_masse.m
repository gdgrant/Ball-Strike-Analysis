function AUC = roc_masse(x1, x2)

AUC = 0;
unique_vals = unique(x1);

for i = 1:length(unique_vals)
    
    p1 = mean(x1 == unique_vals(i));
    p2 = mean(x2 > unique_vals(i)) + 0.5*mean(x2 == unique_vals(i));
    AUC = AUC + p1*p2;
    
end

