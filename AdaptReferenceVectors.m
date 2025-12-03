function V_new = AdaptReferenceVectors(Population, V0)

    Objs = Population.objs;
    Zmin = min(Objs, [], 1);
    Zmax = max(Objs, [], 1);
    Range = Zmax - Zmin;
    Range(Range < 1e-6) = 1e-6;
    V_new = V0 .* Range;
    NormV = vecnorm(V_new, 2, 2);
    NormV(NormV < 1e-6) = 1e-6;
    V_new = V_new ./ NormV;
end