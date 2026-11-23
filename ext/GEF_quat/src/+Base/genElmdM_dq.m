function dM_dq = genElmdM_dq(N, Jrho, G, numDof)
    dG_de = Base.genElmdG_de();
    
    JG  = Jrho * G;
    GTJ = G.' * Jrho;

    K_stack = zeros(numDof, numDof, 4);
    for alpha = 1:4
        dGa = dG_de(:,:,alpha);

        core_term = dGa.' * JG + GTJ * dGa;
        K_stack(:,:,alpha) = 4 * (N.' * core_term * N);
    end

    K_flat = reshape(K_stack, [], 4); 

    dM_dq_flat = K_flat * N;
    dM_dq = reshape(dM_dq_flat, numDof, numDof, numDof);
end