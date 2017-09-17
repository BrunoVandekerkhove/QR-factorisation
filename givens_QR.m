function [Q,R] = givens_QR(A)
% Calculates the QR factorisation of a given matrix using 
% Givens rotations.

[n,m] = size(A); % Dimensies van A
I = eye(n); Q = I; % Eenheidsmatrix (wordt orthogonale matrix)
R = A; % Wordt rechterbovendriehoeksmatrix

for col = 1:m
    for row = n:-1:col+1
        if R(row,col) ~= 0 % Geen rotatie nodig voor dit element
            G = planerot([R(row-1,col) R(row,col)]'); % Givens-rotatie
            R([row-1 row],:) = G*R([row-1 row],:);
            % R(row,col) = 0;
            Q(:,[row-1 row]) = Q(:,[row-1 row]) * G';
        end
    end
end

end