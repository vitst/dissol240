//  Cubic perm model //

{
    Kinv = (nu/K0)*pow(P, -3);   // Perm model
    Kinv.correctBoundaryConditions();
}
