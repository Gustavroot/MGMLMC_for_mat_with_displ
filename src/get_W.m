function [W] = get_W(mgh)

  % Test the construction of derivative based operators.
  % Set up the lattice geometry,
  L = 16;
  T = 16;
  %Nc = 3;
  %Ns = 4;
  %N = L*L*L*Nc*Ns;

  load('matrices/D16/Us.mat');

  % Fix a value of time in Matlab index convention: t = 1,2,...
  t = 1;

  % Select a spatial direction.
  %k = 3;

  % Create the 3D gauge covariant derivative matrix. It has spatial and color indices.
  %D1 = deriv3d(k,t,Us,L);
  % Visualize.
  %figure();
  %spy(D1);

  % Create the full 3D gauge covariant derivative matrix. It has spatial, color
  % indices and spin indices.
  %D1f = derivFull(k,t,Us,L);
  % Visualize.
  %figure();
  %spy(D1f);

  % %Testing: If we move along the spin diagonal the matrices should remain the
  % %same.
  % hop = makeHop(L);
  % X = 45;
  % Y = hop(45,k);
  % disp('All non-zero values should be equal:');
  % err = 0;
  % for ix = 1:L*L*L
  %     for jx = 1:L*L*L
  %         if jx == hop(ix,k) || jx == hop(ix,k+3)
  %             for mu = 1:4
  %                 for a = 1:3
  %                     for b = 1:3
  %                         term1 = D1f(a + 3*(mu-1) + 12*(ix-1),b + 3*(mu-1) + 12*(jx-1));
  %                         term2 = D1(a + 3*(ix-1),b + 3*(jx-1));
  %                         err = err + abs(term1 - term2);
  %                     end
  %                 end
  %             end
  %         end
  %     end
  % end
  % disp(err);

  U1 = Us;
  for n = 1:5
      U2 = smearAPE(U1,L,T,0.35);
      U1 = U2;
      %disp(n);
  end

  %n = size(mgh.D{1},1);

  dim4D = nthroot( size(mgh.D{1},1)/12,4 );
  size3D = size(mgh.D{1},1)/dim4D;
  n = size3D;

  W = sparse(n,n);
  for k=1:3
    W = W + mgh.gs{k+1}*derivFull(k,t,U1,L);
  end

  %W = derivFull(k,t,U1,L);

end

% ---------------------------------------------------------------
% AUXILIARY FUNCTIONS
% --------------------------------------------------------------------------

% Calculate lexicographic index for 3D coordinates
function newMat = projectSU3(mat)
    H = sqrtm(mat'*mat);
    U = mat/H;
    d = det(U);
    newMat = U./(d^(1/3));
end

function index = lexic(x,y,z,L)
    index = z + L*y + L*L*x;
end

function hop = makeHop(L)
    %Array of 1D extension.
    pos = 0:L-1;
    %Positively shifted array.
    posP = circshift(pos,-1);
    %Negatively shifted array.
    posM = circshift(pos,1);
    %Build the hopping matrix. 
    hop = zeros(L*L*L,6);
    for x = 0:L-1
        for y = 0:L-1
            for z = 0:L-1
                ind = z + L*y + L*L*x + 1;
                %Calculate lexic. index of neighbors in all 6 possible 3D
                %directions.             
                %Positive x:
                xNew = posP(x+1);
                hop(ind,1) = lexic(xNew,y,z,L) + 1;
                %Positive y:
                yNew = posP(y+1);
                hop(ind,2) = lexic(x,yNew,z,L) + 1;
                %Positive z:
                zNew = posP(z+1);
                hop(ind,3) = lexic(x,y,zNew,L) + 1;
                %Negative x:
                xNew = posM(x+1);
                hop(ind,4) = lexic(xNew,y,z,L) + 1;
                %Negative y:
                yNew = posM(y+1);
                hop(ind,5) = lexic(x,yNew,z,L) + 1;
                %Negative z:
                zNew = posM(z+1);
                hop(ind,6) = lexic(x,y,zNew,L) + 1;
            end
        end
    end
end

function Unew = smearAPE(U,L,T,a)
    hop = makeHop(L);
    %Copy the gauge field.
    Unew = U;
    %Loop over all values of time.
    for t = 1:T
        %Loop over 3D volume.
        for ix = 1:L*L*L
            %Loop over 4 directions.
            for mu = 2:4
                %Extract link variable.
                Link = U{4*(L*L*L*(t-1) + ix-1) + mu};
                %Loop over the perpendicular spatial directions.
                staple = zeros(3);
                for nu = 2:4
                    if nu ~= mu
                        %First staple.
                        %U_{nu}(x)
                        U1 = U{4*(L*L*L*(t-1) + ix-1) + nu};
                        %U_{mu}(x + nu).
                        xPnu = hop(ix,nu-1);
                        U2 = U{4*(L*L*L*(t-1) + xPnu-1) + mu};
                        %U_{nu}{x+mu}^{+}.
                        xPmu = hop(ix,mu-1);
                        U3 = U{4*(L*L*L*(t-1) + xPmu-1) + nu}';
                        termA = U1*U2*U3;
                        %Second staple.
                        %U_{nu}(x-nu)^{+}.
                        xMnu = hop(ix,nu-1 + 3);
                        U1 = U{4*(L*L*L*(t-1) + xMnu-1) + nu}';
                        %U_{mu}(x-nu)
                        U2 = U{4*(L*L*L*(t-1) + xMnu-1) + mu};
                        %U_{nu}(x + mu - nu)
                        xPmuMnu = hop(xPmu, nu-1);
                        U3 = U{4*(L*L*L*(t-1) + xPmuMnu-1) + nu};
                        termB = U1*U2*U3;
                        staple = staple + termA + termB;
                    end
                end
                %Sum of old link and weighted staples.
                newLink = Link + a.*staple;
                %Project new link back to SU(3).
                Unew{4*(L*L*L*(t-1) + ix-1) + mu} = projectSU3(newLink);
            end
        end
    end
    
end

function Dk = deriv3d(k,t,G,L)
    %We create the 3D gauge-covariant derivative in direction k at time t
    %as a sparse matrix.
    %Create hopping matrix.
    hop = makeHop(L);
    %Array to store the non-zero entries.
    TermVals = zeros(2*9*L*L*L,1);
    %Array to store the row indices of the non-zero values.
    TermRows = zeros(2*9*L*L*L,1);
    %Array to store the column indices of the non-zero values.
    TermCols = zeros(2*9*L*L*L,1);
    %Index to enumerate the non-zero entries.
    idx = 1;
    %Loop over 3D volume rows.
    for ix = 1:L*L*L
        %Index of positive neighbour in direction k.
        xPk = hop(ix,k);
        %Index of negative neighbour in direction k.
        xMk = hop(ix,k+3);
        %Loop over 3D volume columns.
        for jx = 1:L*L*L
            %Check if we are in the positive 3D off-diagonals.
            if xPk == jx
                %Loop over color matrix.
                for a = 1:3
                    for b = 1:3
                        %Fix row index.
                        rowInd = a + 3*(ix-1);
                        %Fix column index.
                        colInd = b + 3*(jx-1);
                        %Extract U_k(x)_{a,b}.
                        Ukx = G{4*(L*L*L*(t-1) + ix-1) + k + 1};
                        %Store the non-zero entry,
                        TermVals(idx) = Ukx(a,b)/2;
                        %Store the row.
                        TermRows(idx) = rowInd;
                        TermCols(idx) = colInd;
                        %Advance the counter.
                        idx = idx + 1;
                    end
                end
            end
            %Check if we are in the negative 3D off-diagonals.
            if xMk == jx
                %Loop over color matrix.
                for a = 1:3
                    for b = 1:3
                        %Fix row index.
                        rowInd = a + 3*(ix-1);
                        %Fix column index.
                        colInd = b + 3*(jx-1);
                        %Extract U_k(x - k)^{+}_{ab}
                        UkxmkD = G{4*(L*L*L*(t-1) + xMk - 1) + k + 1}';
                        %Store the non-zero entry.
                        TermVals(idx) = -UkxmkD(a,b)/2;
                        %Store row and column of the non-zero entry.
                        TermRows(idx) = rowInd;
                        TermCols(idx) = colInd;
                        %Advance the counter.
                        idx = idx + 1;
                    end
                end
            end
        end
    end
    Dk = sparse(TermRows,TermCols,TermVals,3*L*L*L,3*L*L*L);
end

function Dk = derivFull(k,t,G,L)
    %We create the 3D gauge-covariant derivative in direction k at time t
    %with spin indices as a sparse matrix.
    %Create hopping matrix.
    hop = makeHop(L);
    %Array to store the non-zero entries.
    TermVals = zeros(4*2*9*L*L*L,1);
    %Array to store the row indices of the non-zero values.
    TermRows = zeros(4*2*9*L*L*L,1);
    %Array to store the column indices of the non-zero values.
    TermCols = zeros(4*2*9*L*L*L,1);
    %Index to enumerate the non-zero entries.
    idx = 1;
    %Loop over 3D volume rows.
    for ix = 1:L*L*L
        %Index of positive neighbour in direction k.
        xPk = hop(ix,k);
        %Index of negative neighbour in direction k.
        xMk = hop(ix,k+3);
        %Loop over 3D volume columns.
        for jx = 1:L*L*L
            %Loop over spin diagonal index.
            for mu = 1:4
                %Check if we are in the positive 3D off-diagonals.
                if xPk == jx
                    %Loop over color matrix.
                    for a = 1:3
                        for b = 1:3
                            %Fix row index.
                            rowInd = a + 3*(mu-1) + 12*(ix-1);
                            %Fix column index.
                            colInd = b + 3*(mu-1) + 12*(jx-1);
                            %Extract U_k(x)_{a,b}.
                            Ukx = G{4*(L*L*L*(t-1) + ix-1) + k + 1};
                            %Store the non-zero entry,
                            TermVals(idx) = Ukx(a,b)/2;
                            %Store the row.
                            TermRows(idx) = rowInd;
                            TermCols(idx) = colInd;
                            %Advance the counter.
                            idx = idx + 1;
                        end
                    end
                end
                %Check if we are in the negative 3D off-diagonals.
                if xMk == jx
                    %Loop over color matrix.
                    for a = 1:3
                        for b = 1:3
                            %Fix row index.
                            rowInd = a + 3*(mu-1) + 12*(ix-1);
                            %Fix column index.
                            colInd = b + 3*(mu-1) + 12*(jx-1);
                            %Extract U_k(x - k)^{+}_{ab}
                            UkxmkD = G{4*(L*L*L*(t-1) + xMk - 1) + k + 1}';
                            %Store the non-zero entry.
                            TermVals(idx) = -UkxmkD(a,b)/2;
                            %Store row and column of the non-zero entry.
                            TermRows(idx) = rowInd;
                            TermCols(idx) = colInd;
                            %Advance the counter.
                            idx = idx + 1;
                        end
                    end
                end
            end
        end
    end
    Dk = sparse(TermRows,TermCols,TermVals,4*3*L*L*L,4*3*L*L*L);
end