//  Nonlinear reaction kinetics fracture; threshold = Ct  //

    double stoiTable[1][2] = {{-1.0,-1.0}};     // Stoichiometric coeffs

    {
        double Sh = 8.0;                        // Parameters
        double p = 4.0, Ct = 0.2, len = 0.02;
        double w, w1, w2, x, xp, r, dr, drr, g;
        volScalarField& C = CList[0]();

        R  = 2*k0;
        dR = 2*k0;
        forAll (C, i)
        {
            w  = std::tanh((C[i]-Ct)/len);
            w1 = 0.5*(1+w);
            w2 = 0.5*(1-w);
            xp = std::pow(C[i]/Ct,p-2)/Ct;
            x  = xp*C[i];
            r  = (w1 + w2*x)*C[i];
            dr = (w1 + p*w2*x) + (1-w*w)*(1-x)*C[i]/(2*len);
            drr= p*(p-1)*w2*xp + (1-w*w)*(1-p*x)/len - w*(1-w*w)*(1-x)*C[i]/(len*len);
            g  = 2*k0.value()*P[i]/(D0.value()*Sh);
            R[i]  *= r/(1+g*dr);
            dR[i] *= dr/(1+g*dr) -  g*r*drr/((1+g*dr)*(1+g*dr));
        }
        RList[0]()  == R;
        dRList[0]() == dR;
    }
