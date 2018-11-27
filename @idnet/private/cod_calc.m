function codk = cod_calc(y, A, w, k)
% calculate the k-step-ahead coefficient of determination, COD_k
% k == 1, using R_D^2, instead of R^2, adjusted COD_1 for time series

if k == 1
    yhat = A*w;
    num = mean((y(2:end) - yhat(2:end)).^2);
    % dom = mean((y - mean(y)).^2); % only works for stationary ts
    dy = y(2:end) - y(1:end-1);
    dom = mean((dy - mean(dy)).^2);
else
    error('COD_k (k>1) not yet available.')
end

codk = 1 - num/dom;

end
