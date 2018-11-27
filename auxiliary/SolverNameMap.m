function names = SolverNameMap(solver)

    switch solver
      case 'gl1'
        names.fullName = 'Group LASSO';
        names.abbrName = 'GL';

      case 'gil1'
        names.fullName = 'Group Iterative Reweighted l1 Method';
        names.abbrName = 'GIRL1';

      case 'gil2'
        names.fullName = 'Group Iterative Reweighted l2 Method';
        names.abbrName = 'GIRL2';

      case 'gsbl'
        names.fullName = '(Group) Sparse Bayesian Learning';
        names.abbrName = 'GSBL';

      case 'l1'
        names.fullName = 'LASSO';
        names.abbrName = 'LASSO';

      case 'il1'
        names.fullName = 'Iterative Reweighted l1 Method';
        names.abbrName = 'IRL1';

      case 'il2'
        names.fullName = 'Iterative Reweighted l2 Method';
        names.abbrName = 'IRL2';
    end

end