function val = get_frot_pdim_multilinear_form(myAssembly, axisR, mu0, ...
        parameterMultiIndex, inputx)
%GET_FROT_PDIM_MULTILINEAR_FORM Multilinear Frot coefficient for PDIM.
%
% inputx is a cell containing the state directions. Its length determines
% the state degree. parameterMultiIndex selects the coefficient after the
% expansion mu = mu0 + delta, with mu = [omega; dotOmega].

if ~iscell(inputx) || isempty(inputx)
    error('RotFramePDIM:InvalidStateInputs', ...
        'inputx must be a nonempty cell array of state directions.');
end

stateOrder = numel(inputx);
F_handle = @(z) myAssembly.genSysFrotPdimPart( ...
    z, axisR, mu0, stateOrder, parameterMultiIndex);

if stateOrder == 1
    val = F_handle(inputx{1});
else
    val = get_multilinear_form(F_handle, inputx{:});
end
end
