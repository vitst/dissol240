/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
  surfRoughGen

Description
  Preprocessing utility to modify the fracture surface represented by
  two parallel plates. It overwrites 0 directory with modified surface.
  
  Limitations now are the next:
    - a geometry should include two parallel plates
    - the number of faces on each surface should be 
      a power of 2 (for Fourier transform)

Usage
  - surfRoughGen

Needs dictionary
  system/surfRoughGenDict
    
\*---------------------------------------------------------------------------*/

#include <complex>
#include <fftw3.h>

#include "meshRelax.H"

/*
 #######################################################################################
 *    Main program body
 #######################################################################################
*/

class RoughnessGenerator
{
protected:
  // Protected data
  int seed;
  int M;
  int N;
  double rgh;
  word wayToApply;
  double fractalParam; // Set Hausdorff dimension H = 3 - D_f
  double smoothing;
  
public:
  
  // Constructor
  RoughnessGenerator
  (
    int seed_,
    int M_,
    int N_,
    double rgh_,
    word wayToApply_,
    double fractalParam_,
    double smoothing_
  )
  :
    seed(seed_),
    M(M_),
    N(N_),
    rgh(rgh_),
    wayToApply(wayToApply_),
    fractalParam(fractalParam_),
    smoothing(smoothing_)
  {
  }

  void getSurfaceDisplacement(dynamicFvMesh& mesh,scalarField& wd,label& wallID)
  {
    scalarField sFn(M*N, 0.0);
    fftDisp(sFn);
    scalarField sFp(M*N, 0.0);
    if(wayToApply=="asymmetric")
    {
      // shift the seed to get two different numbers
      seed += 125522;
      fftDisp(sFp);
    }
    else
    {
      sFp = sFn;
    }

    pointField pointFace = mesh.boundaryMesh()[wallID].faceCentres();
    scalar maxX = max( pointFace.component(vector::X) );
    scalar maxZ = max( pointFace.component(vector::Z) );
    scalar minX = min( pointFace.component(vector::X) );
    scalar minZ = min( pointFace.component(vector::Z) );
    double Lx = maxX - minX;
    double Lz = maxZ - minZ;

    forAll(pointFace, i)
    {
      scalar curX = pointFace[i].x() - minX;
      scalar curZ = pointFace[i].z() - minZ;

      scalar sign = pointFace[i].y() / mag(pointFace[i].y());

      int curm = std::floor(curZ / Lz * (M-1));
      int curn = std::floor(curX / Lx * (N-1));
      int ind = curn + N * curm;

      if(wayToApply=="symmetric")
        wd[i] = sFn[ind];

      if(wayToApply=="synchronous")
        wd[i] = sign * sFn[ind];

      if(wayToApply=="asymmetric")
      {
        if(sign<0)
          wd[i] = sFn[ind];
        else
          wd[i] = sFp[ind];
      }
    }
  }

private:
  // converts indexes
  label index(label m, label n){ return n + N * m; }
  
  scalar pspec(int u)
  {
    scalar p = Foam::pow(u, -0.5 * (fractalParam+1) );
    p *= Foam::exp(-smoothing*u);
    return p;
  }

  void fftDisp(scalarField& disp)
  {
    unsigned int MN = M*N;

    if ( MN & (MN - 1) )
    {
        FatalErrorIn
        (
             "getSurfaceDisplacement  "
        )   << "number of elements is not a power of 2" << endl
            << "    Number of elements = " << MN
            << abort(FatalError);
    }

    std::vector<std::complex<double> > f, F;
    f.resize(MN);
    F.resize(MN);

    Random rnd( seed );
    scalar TwoPi = constant::mathematical::twoPi;

    Info<< "Displacement calc starts...."<<nl;
    /*
     *   --- ---
     *  | 1 | 2 |
     *   --- ---
     *  | 3 | 4 |
     *   --- ---
     */
    // calculating 1 and 4
    for(int m=0; m<M/2+1; ++m)
    {
      for(int n=0; n<N/2+1; ++n)
      {
        scalar p = TwoPi * rnd.scalar01();
        int u = N * m*m / static_cast<double>(M) + M * n*n / static_cast<double>(N);

        scalar rad = 0.0;
        if(u == 0)
          rad = 0.0;
        else
          rad = pspec(u) * rnd.GaussNormal();

        f[ index(m,n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(((M-m)%M),(N-n)%N) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    f[ index(M/2,0)   ].imag() = 0.0;
    f[ index(0,  N/2) ].imag() = 0.0;
    f[ index(M/2,N/2) ].imag() = 0.0;

    for(int m=1; m<M/2; ++m)
    {
      for(int n=1; n<N/2; ++n)
      {
        scalar p = TwoPi * rnd.scalar01();
        int u = N * m*m / static_cast<double>(M) + M * n*n / static_cast<double>(N);
        scalar rad = pspec(u) * rnd.GaussNormal();

        f[ index(  m, N-n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(M-m,   n) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }

    fftw_plan plan;
    plan = fftw_plan_dft_2d(M, N,
                             reinterpret_cast<fftw_complex*>(&f[0]),
                             reinterpret_cast<fftw_complex*>(&F[0]),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);


    scalarField sF(MN);
    forAll(sF, ii)
    {
      sF[ii] = F[ii].real();
    }

    scalarField sF2 = sqr(sF);
    scalar avSF     = average(sF);
    scalar avSF2    = average(sF2);

    scalar factor = rgh / Foam::sqrt( mag(avSF2 - sqr(avSF)) );

    sF *= factor;

    disp = sF;
  }
};


int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  // reading dictionary surfRoughGenDict
  IOdictionary surfRoughGenDict
  (
    IOobject
    (
      "surfRoughGenDict",
      runTime.system(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );
  
  int seed;
  if( !surfRoughGenDict.readIfPresent<int>("seed", seed) ){
    SeriousErrorIn("main")
            <<"There is no `seed` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }

  int M;
  if( !surfRoughGenDict.readIfPresent<int>("sizeZ", M) ){
    SeriousErrorIn("main")
            <<"There is no `sizeZ` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  int N;
  if( !surfRoughGenDict.readIfPresent<int>("sizeX", N) ){
    SeriousErrorIn("main")
            <<"There is no `sizeX` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  
  scalar rgh;
  if( !surfRoughGenDict.readIfPresent<scalar>("roughness", rgh) ){
    SeriousErrorIn("main")
            <<"There is no `roughness` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  // symmetric, synchronous, asymmetric
  word wayToApply;
  if( !surfRoughGenDict.readIfPresent<word>("apply", wayToApply) ){
    SeriousErrorIn("main")
            <<"There is no `synchronous` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  
  double fractalParam;
  if( !surfRoughGenDict.readIfPresent<double>("fractalParam", fractalParam) ){
    SeriousErrorIn("main")
            <<"There is no `fractalParam` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  
  double smoothing;
  if( !surfRoughGenDict.readIfPresent<double>("smoothing", smoothing) ){
    SeriousErrorIn("main")
            <<"There is no `smoothing` parameter in surfRoughGenDict dictionary"
            <<exit(FatalError);
  }
  
  Info<< "Seed:       " << seed << endl;
  Info<< "M:          " << M << endl;
  Info<< "N:          " << N << endl;
  Info<< "roughness:  " << rgh << endl;
  Info<< "apply:      " << wayToApply << endl;
  Info<< "fractalParam:  " << fractalParam << endl;
  Info<< "smoothing:     " << smoothing << endl;

  Info<< "Setup RoughnessGenerator class" << endl;
  RoughnessGenerator rg(seed, M, N, rgh, wayToApply, fractalParam, smoothing);
  
  Info<< "Setup mesh relaxation class" << endl;
  
  double cpuTime = runTime.elapsedCpuTime();
  meshRelax mesh_rlx(mesh, args);
  
  // Get patch ID for boundaries we want to move ("walls" "inlet")
  label wallID  = mesh.boundaryMesh().findPatchID("walls");
  label inletID = mesh.boundaryMesh().findPatchID("inlet");
  label outletID = mesh.boundaryMesh().findPatchID("outlet");
  
  coupledPatchInterpolation patchInterpolator( mesh.boundaryMesh()[wallID], mesh );
  
  const pointField& boundaryPoints = mesh.boundaryMesh()[wallID].localPoints();
  vectorField pointDispWall(boundaryPoints.size(), vector::zero);
  
  vectorField pointNface = mesh.boundaryMesh()[wallID].faceNormals();
  vectorField motionN = patchInterpolator.faceToPointInterpolate(pointNface);
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  scalarField faceDisp(pointNface.size(), 0.0);
  rg.getSurfaceDisplacement(mesh, faceDisp, wallID);
  scalarField pointDisp = patchInterpolator.faceToPointInterpolate(faceDisp);
  
  Info<<nl<< "Maximum and minimum face displacement in Y direction:"<<nl;
  Info<<nl<< "Max: " << max(faceDisp) << "  min: "<<min(faceDisp)<<nl<<nl;
  Info<<nl<< "Maximum and minimum point displacement in Y direction:"<<nl;
  Info<<     "Max: " << max(pointDisp) << "  min: "<<min(pointDisp)<<nl<<nl;
  
  forAll( pointDispWall, i )
    pointDispWall[i] = pointDisp[i] * motionN[i];
  
  pointDispWall /= runTime.deltaTValue();
  
  bool auxSw = mesh_rlx.get_fixInletWallEdgeDispl();
  mesh_rlx.set_fixInletWallEdgeDispl(false);
  mesh_rlx.meshUpdate(pointDispWall, runTime);
  mesh_rlx.set_fixInletWallEdgeDispl(auxSw);
  
  cpuTime = runTime.elapsedCpuTime() - cpuTime;
  
  Info<<nl<<"Time statistics:"<<nl;
  
  int wlNP = mesh.boundaryMesh()[wallID].nPoints();
  int inNP = mesh.boundaryMesh()[inletID].nPoints();
  int ouNP = mesh.boundaryMesh()[outletID].nPoints();
  Info<<"Total number of points:                   "<<mesh.nPoints()<<nl;
  Info<<"Total (approx) number of boundary points: "<<wlNP+inNP+ouNP<<nl;
  Info<<"Number of points on the walls:            "<<wlNP<<nl;
  Info<<"Number of points on the inlet:            "<<inNP<<nl;
  Info<<"Number of points on the outlet:           "<<ouNP<<nl;
  Info<<"Running time:                             "<<cpuTime<<nl<<endl;

  Info<<"Overwriting points in current time directory."<<nl;
  runTime.writeNow();
  
  Info<<"End"<<nl;
  return 0;
}

// **************************** End of the solver ******************************** //
