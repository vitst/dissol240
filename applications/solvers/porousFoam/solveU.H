//  NS-Brinkman flow solver: SIMPLE with Darcy initial condition

    U   == U*Uscale;                // Undo scaling
    phi == phi*Uscale;              // Improves solver convergence

    if (DARCY)                      //  Darcy-flow solver
    {
        forAll(Kinv, i)
        {
            if (Kinv[i] < 1.0/Kmax)  Kinv[i] = 1.0/Kmax;
        }
        Kinv.correctBoundaryConditions();

        Info << "Start Darcy solver" << endl;

        int counter=0;
        while(counter < nCycle)
        {
            counter++;
            double residual = solve
            (
                fvm::laplacian(1./Kinv,p)
            ).initialResidual();

            if (DEBUG)
                Info << "Iteration "    << counter 
                     << ": residual = " << residual << endl;

            if (residual < tolerance)
            {
                Info << "Darcy solver: Converged in "
                     << counter << " steps" << endl;
                break;                              // Done
            }                
            else if (counter == nCycle)
            {
                Info << "Darcy: FAILED to converge in "
                     << nCycle << " steps" << nl << endl;
                     exit(1);                       // No convergence 
            }
        }

        phi == -fvc::interpolate(1/Kinv)*fvc::snGrad(p)*mag(mesh.Sf());
        U   == -(1/Kinv)*fvc::grad(p);  //  Note: U is not conservative
    }
    else
    {
        Info << "Start Brinkman solver" << endl;

        while (simple.loop())
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }
    }
    if (DEBUG)
    {
        #include "continuityErrs.H"
    }

    #include "avgU.H"           // Average flow velocity (Ux)

    U   == U/Uscale;            // Normalize U
    phi == phi/Uscale;          // Normalize surface volume flux

    Info << "Flow solver done: Uscale = " << Uscale.value() << nl
         << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
    	 << "    ClockTime = " << runTime.elapsedClockTime() << " s"
    	 << endl;
