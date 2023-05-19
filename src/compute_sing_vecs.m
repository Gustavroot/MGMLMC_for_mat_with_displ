function [] = compute_sing_vecs(mgh)

  nr_sing_vecs = 28;
  Vsing = rand(size(mgh.D{1},1),nr_sing_vecs);
  C = ones(size(Vsing));
  for i=1:5
    for j=1:nr_sing_vecs
      C(:,j) = pgmres(mgh.g5{1}*Vsing(:,j),mgh,1,1.0e-12);
      fprintf(".");
    end
    [C,~] = qr(C,0);
    Vsing(:,:) = C(:,:);
    fprintf("\n");
  end

  save('Vsing.mat','Vsing');

end
