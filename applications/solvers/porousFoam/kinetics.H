//  Calculate R and dR for a species; include mineral conditions

{
    #include "reacModel.H"              // Get R & dR for all reactions

    // Constraints from individual mineral volumes

    for (int f=0; f < nF; f++){         // Calculate dF
        dFList[f]().internalField() = scalar(0);
        for (int r=0; r < nR; r++)
            dFList[f]() += RList[r]()*stoiTable[r][f]*deltaT;
    }

    if (!FRAC2D){                       // Limit mineral dissolution
        double Fn, dF, fac;
        for (int f=0; f < nF; f++){
            forAll(FList[f](), n){
                Fn =   FList[f]()[n];
                dF = -dFList[f]()[n];   // - delta F
                if (dF < SMALL)         // No limit
                    fac = 1;
                else if (Fn < SMALL)
                    fac = 0;
                else if (Fn < dF)
                    fac = Fn/dF; 
                else
                    fac = 1;
                for (int r=0; r < nR; r++)
                if (stoiTable[r][f] != 0){  // Scale rates in f
                    RList[r]()[n] *= fac;
                   dRList[r]()[n] *= fac;
                }
            }
            dFList[f]().internalField() = scalar(0);
            for (int r=0; r < nR; r++)
                dFList[f]() += RList[r]()*stoiTable[r][f]*deltaT;
        }
    }

    // Constraints from porosity

    for (int f=1; f < nF; f++)          // Calculate total dF
        dFList[0]() += dFList[f]();

    double Pn, dP, fac;
    forAll(dFList[0](), n){
        if (!FRAC2D){                   // Porous
            Pn = P[n] - Pmin;           // Pn > 0
            dP = dFList[0]()[n];        // - delta P
        }
        else{                           // Fracture
            Pn = Pmax - P[n];           // Pn > 0
            dP =-dFList[0]()[n];        // + delta P
        }
        if (dP < SMALL)                 // No limit
            fac = 1;
        else if (Pn < SMALL)            // Limit porosity
            fac = 0;
        else if (Pn < dP)
            fac = Pn/dP; 
        else
            fac = 1.0;
        for (int f=0; f < nF; f++)     // Scale all rates
        for (int r=0; r < nR; r++)
        if (stoiTable[r][f] != 0){     // Scale rates in f
            RList[r]()[n] *= fac;
           dRList[r]()[n] *= fac;
        }
    }

    //  Calculate reaction rate for Case

     R.internalField() = scalar(0);
    dR.internalField() = scalar(0);
    for (int r=0; r < nR; r++){
        R +=  RList[r]()*stoiTable[r][Case];
       dR += dRList[r]()*stoiTable[r][Case];
    }
     R.correctBoundaryConditions();
    dR.correctBoundaryConditions();
}
