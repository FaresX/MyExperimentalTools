let
    f(z) = 2(1 + 2z^2) * atanh(2z * √((1 + z^2) / (1 + 6z^2 + 4z^4))) * (z * √((1 + z^2) * (1 + 6z^2 + 4z^4)))^-1 - 4 / 3
    D(z) = 1 / (1 + z^2)
    g(D) = D == 1 ? π : π * (1 - D)^(1 // 4) / √((1 + √(1 - D)) * (1 + √(1 - D) - D))
    global function solveD(Ic, Iexc)
        # z = find_zero(z -> Iexc * π * Ic0(D(z)) / 2 - Ic * f(z), 0.5)
        z = find_zero(z -> Iexc / Ic - f(z) / g(D(z)), 0.5)
        # z = find_zero(z -> Ic / Iexc - π * Ic0(D(z)) / (2f(z)), 0.5)
        return D(z), Ic / g(D(z))
    end
end