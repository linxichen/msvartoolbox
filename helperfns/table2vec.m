function [y_vec,X_vec] = table2vec(y_table,X_table,model)
% convert data from the usual table form (think Stata) into a stacked
% vector form required by my math notes
% input: y_table the T-by-N depdent variable table
%        X_table the T-by-K indie variable table, K the number all regressors
%        model = contains the info about the model, mostly sizes
% assumption: each depdent variable shares the same set of regressor!
%             if not, follow the note in how to transform different set of
%             regressors into one big mommy set in a sparese matrix

%% basic error checking
if isequal(size(y_table),[model.T,model.N])
else
	error('dependent vector is not expected size');
end

if isequal(size(X_table),[model.T,model.K])
else
	error('regressor matrix is not expected size');
end
%% actuall processing
yprime = y_table';
y_vec = yprime(:);

X_vec = kron(X_table,ones(model.N,1)); % repeat each row N times because same RHS

end