function [idx_train, idx_test] = get_test_train_trials(td,n)

if 0
    itest       = ishuf( 1:bins_p_fold );
    itrain      = ishuf(~ismember(ishuf,itest));
else % get 'n' random trial(s) for each target direction
    u = unique([td.target_direction]);
    idx_train = 1:length(td);
    
    idx_test = zeros(length(u),n);
    for i = 1:length(u)
        idx_test(i,:) = getTDidx(td,'target_direction',u(i),'rand',n);
    end
    idx_test = reshape(idx_test,1,numel(idx_test));
    
    % now remove the test ones
    idx_train = idx_train(~ismember(idx_train,idx_test));
end
