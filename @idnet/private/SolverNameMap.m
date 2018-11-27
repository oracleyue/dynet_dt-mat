function name = SolverNameMap(solver)

    switch solver
      case 'gl1'
        name = 'Group LASSO';
        
      case 'gil1'
        name = 'Group Iterative Reweighted l1 Method';
        
      case 'gil2'
        name = 'Group Iterative Reweighted l2 Method';
        
      case 'gsbl'
        name = '(group) Sparse Bayesian Learning';
        
      case 'gspmc'
        name = 'Group sparsity using MCMC methods';
        
      case 'l1'
        name = 'LASSO';
        
      case 'il1'
        name = 'Iterative Reweighted l1 Method';
        
      case 'il2'
        name = 'Iterative Reweighted l2 Method';
    end

end