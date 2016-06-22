//  Linear reaction kinetics - transport corrected //


    double stoiTable[1][2] = {{-1.0,-1.0}};  // Stoichiometric coeffs

    double G=0.187;

    RList[0]()  == k0*CList[0]()/(1+G*P);
    dRList[0]() == k0/(1+G*P);

