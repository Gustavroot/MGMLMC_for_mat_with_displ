function [mgh] = two_level_mg(mgh,level_nr,from_pgmres)

  % post-smoothing only
  if level_nr==1
    smooth_iters = 5;
  else
    smooth_iters = 3;
  end

  if level_nr>1
    if from_pgmres==1
      w_cycle_iters = 1;
    else
      w_cycle_iters = 3;
    end
  else
    w_cycle_iters = 1;
  end

  % initial guess zero
  mgh.xs{level_nr}(:) = 0;

  for i=1:w_cycle_iters

    % coarse grid correction
    mgh.rs{level_nr} = mgh.bs{level_nr} - mgh.D{level_nr}*mgh.xs{level_nr};
    mgh.bs{level_nr+1} = mgh.R{level_nr}*mgh.rs{level_nr};
    if level_nr==length(mgh.D)-1
      mgh.xs{level_nr+1} = mgh.invD{level_nr+1} * mgh.bs{level_nr+1};
    else
      [mgh] = two_level_mg(mgh,level_nr+1,0);
    end
    mgh.xs{level_nr} = mgh.xs{level_nr} + mgh.P{level_nr}*mgh.xs{level_nr+1};

    % compute the residual
    mgh.rs{level_nr} = mgh.bs{level_nr} - mgh.D{level_nr}*mgh.xs{level_nr};
    % apply smoother
    mgh.es{level_nr} = smooth(mgh.D{level_nr},mgh.rs{level_nr},smooth_iters,mgh.omega{level_nr});
    % update solution
    mgh.xs{level_nr} = mgh.xs{level_nr} + mgh.es{level_nr};

  end

end
