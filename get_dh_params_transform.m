% Breno Niehues dos Santos - 2022/11/05

function [dh, T] = get_transform(a, alpha, d, theta)
%Symbolic homogenous transform matrix from the DH params.
% 𝑎_𝑖: a distancia entre os eixos 𝑍_(𝑖−1) e 𝑍_𝑖 ao longo do eixo 𝑋_(𝑖−1)
% 𝛼_𝑖: o angulo entre os eixos 𝑍_(𝑖−1) e 𝑍_𝑖 ao longo do eixo 𝑋_(𝑖−1)
% 𝑑_𝑖: a distancia entre os eixos 𝑋_(𝑖−1) e 𝑋_𝑖 ao longo do eixo 𝑍_𝑖
% 𝜃_𝑖: o angulo entre os eixos 𝑋_(𝑖−1) e 𝑋_𝑖 ao longo do eixo 𝑍_𝑖

    dh = [transpose(a), transpose(alpha), transpose(d), transpose(theta)];

    n = length(a);
    T = cell(n+1, 4);
    T{1, 1} = sym(eye(4,4));
    for k = 1:n
        row1 = [cos(theta(k)), -sin(theta(k)), 0, a(k)];
        row2 = [sin(theta(k))*cos(alpha(k)), cos(theta(k))*cos(alpha(k)), -sin(alpha(k)), -sin(alpha(k))*d(k)];
        row3 = [sin(theta(k))*sin(alpha(k)), cos(theta(k))*sin(alpha(k)), cos(alpha(k)), cos(alpha(k))*d(k)];
        row4 = [0, 0, 0, 1];

        T{k, k+1} = [row1; row2; row3; row4];
        T{k+1, k+1} = sym(eye(4, 4));
    end

    for k = 1:n+1
        for m = k+2:n+1
            [k m];
            T{k, m} = simplify( T{k, m-1} * T{m-1, m});
        end
    end
end

