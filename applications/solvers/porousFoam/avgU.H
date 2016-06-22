//  Calculate average flow velocity <Ux>

{
    dimensionedScalar intV 
    (
        "intV",
        dimVol,
        gSum(mesh.V())
    );

    dimensioned<vector> intU
    (
        "intU",
        U.dimensions()*dimLength*dimLength*dimLength,
        pTraits<vector>::zero
     );
    intU = fvc::domainIntegrate(U);
    dimensionedScalar avgU = intU[0]/intV;


//  Calculate velocity scaling factor (Uscale)

    if (Uscale.value() == 0)
        Uscale.value() =  avgU.value();     // Initial value of <U>

    if (QCON)
        Uscale.value() = avgU.value();
    else if (avgU.value() > Qmax*Uscale.value())
        Uscale.value() = avgU.value()/Qmax;
}
