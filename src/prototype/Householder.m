function [Q R] = Householder(A)
  [m n] = size(A);
  Q = eye(m,m);
  R = A;

  m1 = max([m n]);
  n1 = min([m n]);
  m = m1;
  n = n1;

  for k = 1:n
    x = R(k:m,k);
    e1 = zeros(m-k+1,1);
    e1(1) = 1;
    vk = sign(x(1)) * norm(x) * e1 + x;
    vk = vk / norm(vk);
    Qk = eye(m);
    Qk(k:m,k:m) = eye(m-k+1) - 2 * vk * vk';
    Q = Q * Qk;
    R = Qk * R;
  end
end
