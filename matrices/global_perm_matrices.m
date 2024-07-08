function [GPM] = global_perm_matrices()

  % returns the matrix which is a reordering into blocks, a-la
  % DDalphaAMG

  % this is how Dreord was built:
  % GPM = global_perm_matrices();
  % save('GPM.mat','GPM','-v7.3');
  % load GPM.mat
  % load D.mat
  % Dreord = GPM*(D*GPM');

  % these parameters are for a 3-level method for D16
  aggr_dim = [4,2];
  latt_dim = [16,4,2];
  nr_levels = 3;
  dofs = 12;

  % set of matrices
  GP = {{}};

  % the overall matrix
  GPM = speye(12 * ( latt_dim(1)^4 ));

  eye_entries = dofs;

  for level=1:(nr_levels-1)

    GP{level} = sparse(latt_dim(level)^4,latt_dim(level)^4);

    ctr = 1;

    % go over the coordinates, in block-order
    for x1=1:latt_dim(level)/aggr_dim(level)
      for x2=1:latt_dim(level)/aggr_dim(level)
        for x3=1:latt_dim(level)/aggr_dim(level)
          for x4=1:latt_dim(level)/aggr_dim(level)

            offset_x1 = (x1-1)*aggr_dim(level);
            offset_x2 = (x2-1)*aggr_dim(level);
            offset_x3 = (x3-1)*aggr_dim(level);
            offset_x4 = (x4-1)*aggr_dim(level);

            % each iteration at this level corresponds to one block
          
            for i1=1:aggr_dim(level)
              for i2=1:aggr_dim(level)
                for i3=1:aggr_dim(level)
                  for i4=1:aggr_dim(level)

                    i = i4+offset_x4 + (i3-1+offset_x3)*latt_dim(level) + ...
                        (i2-1+offset_x2)*(latt_dim(level)*latt_dim(level)) + ...
                        (i1-1+offset_x1)*(latt_dim(level)*latt_dim(level)*latt_dim(level));
                    GP{level}(ctr,i) = 1;
                    ctr = ctr + 1;

                  end
                end
              end
            end

          end
        end
      end
    end

    GP{level} = kron(GP{level},speye(eye_entries));
    eye_entries = eye_entries * (aggr_dim(level)^4);
    GPM = GP{level}*GPM;

  end

end