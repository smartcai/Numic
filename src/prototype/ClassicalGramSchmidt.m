function [Q R] = ClassicalGramSchmidt(A)
  [m n] = size(A);
  Q = zeros(m,n);
  R = zeros(n,n);

  V = A;
  for i = 1:n
    for j = 1:i-1
      R(j,i) = Q(:,j)' * A(:,i);
      V(:,i) = V(:,i) - R(j,i) * Q(:,j);
    end
    R(i,i) = norm(V(:,i));
    Q(:,i) = V(:,i) / R(i,i);
  end
end
