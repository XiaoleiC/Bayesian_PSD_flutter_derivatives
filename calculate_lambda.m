function lambda0 = calculate_lambda(y1,y2)
n = length(y1);
a1 = (n*sum(y1.*y2) - sum(y1)*sum(y2))^2;
a2 = (n*sum(y1.^2) - (sum(y1))^2)*(n*sum(y2.^2) - (sum(y2))^2);
lambda0 = a1/a2;
end