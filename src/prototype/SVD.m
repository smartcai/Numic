clear; close all; clc;

function [B1, G] = oneiter_right(B, k)
  B1 = B;
  [c, s] = givens(B1(k, k), B1(k, k+1));
  G = [c s; -s c];
  disp([c, s]);
  B1(:, k:k+1) = B1(:, k:k+1) * G';
end

function [B1, G] = oneiter_left(B, k)
  B1 = B;
  [c, s] = givens(B1(k, k), B1(k+1, k));
  G = [c s; -s c];
  B1(k:k+1, :) = G * B1(k:k+1, :);
end

B0 = [
      0.538516  3.145988 0.000000 0.000000 -0.000000 ;
      0.000000  0.700351 1.822189 0.000000 -0.000000 ;
     -0.000000 -0.000000 0.806397 0.040211  0.000000 ;
];

B1 = [
       0.538516  3.145988 0.000000 ;
       0.000000  0.700351 1.822189 ;
      -0.000000 -0.000000 0.806397 ;
];

Bout = B0;

m = 3;
n = 4;
r = min(m, n);
iter = 30;
for i = 1:iter
  for k = 1:r - 1
    [Bout, G] = oneiter_right(Bout, k);
    [Bout, G] = oneiter_left(Bout, k);
  end
  if n > m
    [Bout, G] = oneiter_right(Bout, r);
  end
end
