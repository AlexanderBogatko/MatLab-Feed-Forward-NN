% функция для вычисления расстояния между точками (X, Y, Z) и (Xp, Yp, Zp)
function [result] = Rasst(X, Y, Z, Xp, Yp, Zp)
    result = sqrt((X - Xp)^2 + (Y - Yp)^2 + (Z - Zp)^2);
end

