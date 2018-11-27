figure;
idxExp = randi(num_exp);
time = 1:size(inputs{idxExp},1);
subplot(4,1,1)
idxChannel = randi(num_inputs);
plot(time, inputs{idxExp}(:,idxChannel), '-');
ylabel(sprintf('input(%u)',idxChannel))
title('Examples of input and output signals')
subplot(4,1,2)
plot(time, outputs{idxExp}(:,1), '-');
ylabel('output(1)')
subplot(4,1,3)
idxChannel = randi(num_nodes);
plot(time, outputs{idxExp}(:,idxChannel), '-');
ylabel(['output(' num2str(idxChannel) ')'])
subplot(4,1,4)
plot(time, outputs{idxExp}(:,end), '-');
ylabel(sprintf('output(%u)', num_nodes))
xlabel(sprintf('time (Ts = %g)', Ts))