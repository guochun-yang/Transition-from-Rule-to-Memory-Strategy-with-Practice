function [eta_p] = geteffectfromfitlme(lme,myvarname)
%Usage: [eta_p] = geteffectfromfitlme(lme,myvarname)

% get index for myvarname
varnames = lme.Coefficients.Name;
idx_var = find(ismember(varnames,myvarname));
residual = lme.residuals;

% get the table
tbl = lme.Variables;
tbvarnames = tbl.Properties.VariableNames;
if contains(myvarname,":") % interaction effect
    intvarnames = split(myvarname,':'); %'a:b' will be {'a','b'}
    nint = length(intvarnames);% how many way interaction
    for iint = 1:nint
        idx_tbl(iint) = find(ismember(tbvarnames,intvarnames{iint}));
        predictors(:,iint) = table2array(tbl(:,idx_tbl(iint)));
    end
    predictor = prod(predictors,2);
else
    idx_tbl = find(ismember(tbvarnames,myvarname));
    predictor = table2array(tbl(:,idx_tbl));
end

y_fit_var = predictor*lme.Coefficients.Estimate(idx_var);

% get partial eta
eta_p = var(y_fit_var,'omitnan')/(var(y_fit_var,'omitnan')+var(residual,'omitnan'));

end

