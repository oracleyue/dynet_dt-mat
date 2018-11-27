function rmse = rmse_calc(y, A, w)
% calculate the root-mean-square-error of the model

yhat = A*w;
mse = mean((y - yhat).^2);
rmse = sqrt(mse);

end