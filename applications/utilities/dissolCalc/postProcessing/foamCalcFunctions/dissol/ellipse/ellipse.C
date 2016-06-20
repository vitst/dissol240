/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.


Description
    contact info: vitaliy.starchenko@gmail.com
    
SourceFiles
    ellipse.C

\*---------------------------------------------------------------------------*/

#define NUMBER_OF_COLUMNS 10

#include "ellipse.H"

#include "addToRunTimeSelectionTable.H"
#include "memInfo.H"

#include "OFstreamMod.H"

#include "interpolation.H"

#include "triSurfaceSearch.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(ellipse, 0);
        addToRunTimeSelectionTable(calcType, ellipse, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::ellipse::ellipse()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::ellipse::~ellipse()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::ellipse::init()
{
  argList::validArgs.append("ellipse");
  argList::validArgs.append("processingType");
  argList::validArgs.append("N");
  argList::validArgs.append("M");
  argList::validArgs.append("K");
}


void Foam::calcTypes::ellipse::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
  #ifdef FOAM_DEV
    const word& processingType = args.additionalArgs()[1];
    const word& Nword = args.additionalArgs()[2];
    const word& Mword = args.additionalArgs()[3];
    const word& Kword = args.additionalArgs()[4];

    N_ = std::atoi( Nword.c_str() );
    M_ = std::atoi( Mword.c_str() );
    K_ = std::atoi( Mword.c_str() );
    Info << "FOAM_DEV true: "<< FOAM_DEV <<nl;
  #else
    const word processingType = args[2];
    N_ = std::atoi( args[3].c_str() );
    M_ = std::atoi( args[4].c_str() );
    K_ = std::atoi( args[5].c_str() );
    Info << "FOAM_DEV false: " <<nl;
  #endif
  processingType_ = const_cast<word&>(processingType);
  
  N1_ = N_+1;
  M1_ = M_+1;
  K1_ = K_+1;

  Info << "Calculating the max and min coordinates"<<nl;
  const pointField &allPoints = mesh.points();
  
  maxPosX = max( allPoints.component(vector::X) );
  maxPosZ = max( allPoints.component(vector::Z) );
  
  minPosX = min( allPoints.component(vector::X) );
  minPosZ = min( allPoints.component(vector::Z) );
  
  // z becomes x; x becomes y
  dx = (maxPosZ - minPosZ) / static_cast<scalar>(N_);
  dy = (maxPosX - minPosX) / static_cast<scalar>(M_);

  pointsXY_2D_X.setSize( N1_ * M1_ );
  
  pointsXY_2D_Z.setSize( N1_ );

  scalar p0Z = minPosZ;
  
  // N+1 and M+1 because we need points at x=minPosX and x=maxPosX
  for(int i=0; i<N1_; i++){
    scalar zz = p0Z + i*dx;
    pointsXY_2D_Z[i] = zz;
  }
  
  pointsXYonsurface.setSize( 2.0 * N1_ * M1_ );
}

void Foam::calcTypes::ellipse::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
  Info << "Find the points on the surface"<<nl;
  //build_surface_points( mesh );

  fileName current_dissolpostproc_dir;
  current_dissolpostproc_dir = "dissolpostproc" / runTime.timeName();
  if ( !isDir(current_dissolpostproc_dir) ) mkDir(current_dissolpostproc_dir);

  if(processingType_ == "all"){
    write_q(mesh, runTime);
    write_Ccup(mesh, runTime);
  }
  else if(processingType_ == "aspect"){
    write_aspect(mesh, runTime);
  }
  else if(processingType_ == "U"){
    write_q(mesh, runTime);
  }
  else if(processingType_ == "p"){
    //write_p(mesh, runTime);
    FatalError<<"p processing is not implemented yet "<<nl<<nl<<exit(FatalError);
  }
  else if(processingType_ == "C"){
    write_Ccup(mesh, runTime);
  }
  else if(processingType_ == "temp"){
    write_temp(mesh, runTime);
  }
  else{
    FatalError<<"Unable to process "<<processingType_<<nl<<nl<<exit(FatalError);
  }
      
}

scalar Foam::calcTypes::ellipse::primitive_simpson_integration
(
  scalarField& x,
  scalarField& y  
)
{
  scalar result = 0.0;
  
  // TODO check len y == len x
  int step = 2;
  
  scalarField dlt_x( x.size()-1 );
  
  forAll(dlt_x, i){
    dlt_x = x[i+1]-x[i];
  }
  
  int len_lab_list = static_cast<int>( std::floor( dlt_x.size() / 2.0 ) );
  
  labelList slice0( len_lab_list );
  labelList slice1( len_lab_list );
  labelList slice2( len_lab_list );
  
  int count = 0;
  forAll(slice0, i){
    slice0[i] = count;
    slice1[i] = count+1;
    slice2[i] = count+2;
    count += step;
  }
  
  scalarField dx0(len_lab_list);
  scalarField dx1(len_lab_list);
  forAll(dx0, j){
    dx0[j] = dlt_x[slice0[j]];
    dx1[j] = dlt_x[slice1[j]];
  }

  scalarField hsum(len_lab_list); //= h0 + h1
  scalarField hprod(len_lab_list); //= h0 * h1
  scalarField h0divh1(len_lab_list); // = h0 / h1

  
  forAll(dx0, i){
    hsum[i] = dx0[i] + dx1[i];
    hprod[i] = dx0[i] * dx1[i];
    h0divh1[i] = dx0[i] / dx1[i];
  }
  
  forAll(hsum, i){
    result += hsum[i]/6.0 *
                          (
                            y[slice0[i]] * (2-1.0/h0divh1[i]) +
                            y[slice1[i]] * hsum[i]*hsum[i]/hprod[i] +
                            y[slice2[i]] * (2-h0divh1[i])
                          );
  }
  
  return result;
}



void Foam::calcTypes::ellipse::build_surface_points
(
  const fvMesh& mesh
)
{
  
  oppositeWallUniform.clear();
  
  // triangulation
  // @TODO in case of parallel calculations see surfaceMeshTriangulate.C
  // @TODO add the error handling for "walls"
  label wallID = mesh.boundaryMesh().findPatchID("walls");
  labelHashSet includePatches(1);
  includePatches.insert(wallID);
  
  triSurface wallTriSurface
  (
    triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
  );
  
  // intersection
  //const vectorField& normals = wallTriSurface.faceNormals();
  const triSurfaceSearch querySurf(wallTriSurface);
  const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
  
  forAll(pointsXY_2D_X, i){
    point searchStart(0.0, 0.0, pointsXY_2D_Z[i]);
    while(searchStart.x() >= maxPosX ) searchStart.x() = searchStart.x() - 0.001;
    while(searchStart.x() <= minPosX ) searchStart.x() = searchStart.x() + 0.001;
    while(searchStart.z() >= maxPosZ ) searchStart.z() = searchStart.z() - 0.001;
    while(searchStart.z() <= minPosZ ) searchStart.z() = searchStart.z() + 0.001;
    point searchEnd = searchStart;
    searchStart.y() = 500.0;
    
    point hitPoint(0.0, 0.0, 0.0);
    
    label posWallLab = 2*i;
    label negWallLab = 2*i+1;
    
    // @TODO description
    bool pm = false, mp = false;
    
    pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        pm = true;
    }
    else{
    }
    
    pointsXYonsurface[posWallLab] = hitPoint;
    
    // search symmetric point
    searchStart.y() = -500.0;
    pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        mp = true;
    }
    else{
    }

    pointsXYonsurface[negWallLab] = hitPoint;
    
    if(pm && !mp){
      pointsXYonsurface[negWallLab] = pointsXYonsurface[posWallLab];
      pointsXYonsurface[negWallLab].y() *= -1.0;
    }
    else if(!pm && mp){
      pointsXYonsurface[posWallLab] = pointsXYonsurface[negWallLab];
      pointsXYonsurface[posWallLab].y() *= -1.0;
    }
    
    oppositeWallUniform[posWallLab] = negWallLab;
  }
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// an empty function for something extremely not general
void Foam::calcTypes::ellipse::write_temp
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  Info << "temporary function..." << endl;
  
  fileName current_file_path;
  current_file_path = "dissolpostproc" / runTime.timeName() / "csurf";
  OFstream aFile( current_file_path );
  
  aFile << N_ << "   " << M_ << "   " << runTime.value() 
          << "   " << dx << "   " << dy << endl;
  
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  fieldC field_c(headerC, mesh);
  autoPtr<interpolation<scalar> >interpolatorC(interpolation<scalar>::New("cellPointFace",field_c));
  
  int count = 0;
  meshSearch searchEngine(mesh);
  for (std::map<int,int>::const_iterator it  = oppositeWallUniform.begin(); 
                                   it != oppositeWallUniform.end(); ++it ){
    point ap = pointsXYonsurface[it->first];
    label cellI = searchEngine.findNearestCell( ap );
    label faceI = searchEngine.findNearestFace( ap );
    scalar ac = interpolatorC->interpolate(ap, cellI, faceI);
    
    aFile << ac << "  ";

    count++;
    if(count>=NUMBER_OF_COLUMNS){ 
      aFile <<"\n";
      count=0;
    }
  }  
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Foam::calcTypes::ellipse::write_aspect
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  //Info << "Processing the  aperture of the fracture..." << endl;
  
  // triangulation
  // @TODO in case of parallel calculations see surfaceMeshTriangulate.C
  // @TODO add the error handling for "walls"
  label wallID = mesh.boundaryMesh().findPatchID("walls");
  labelHashSet includePatches(1);
  includePatches.insert(wallID);
  
  triSurface wallTriSurface
  (
    triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches )
  );
  
  // intersection
  //const vectorField& normals = wallTriSurface.faceNormals();
  const triSurfaceSearch querySurf(wallTriSurface);
  const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

  
  fileName current_file_path;
  current_file_path = "dissolpostproc" / runTime.timeName() / "aspect";
  OFstream apertureFile( current_file_path );
  
  forAll(pointsXY_2D_Z, i){
    point x1, x2, y1, y2;
    
    point searchStart(0.0, 0.0, pointsXY_2D_Z[i]);
    while(searchStart.z() >= maxPosZ ) searchStart.z() = searchStart.z() - 0.001;
    while(searchStart.z() <= minPosZ ) searchStart.z() = searchStart.z() + 0.001;
    point searchEnd = searchStart;
    searchStart.x() = 100.0;
    
    point hitPoint(0.0, 0.0, 0.0);
    
    // @TODO description
    bool pm = false, mp = false;
    
    pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        pm = true;
    }
    else{
    }
    x1 = hitPoint;
    
    // search symmetric point
    searchStart.x() = -100.0;
    pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        mp = true;
    }
    else{
    }
    x2 = hitPoint;

    
    if(pm && !mp){
      x2 = x1;
      x2.x() *= -1.0;
    }
    else if(!pm && mp){
      x1 = x2;
      x1.x() *= -1.0;
    }
    // ----------------------------------------
    searchStart = point(0.0, 0.0, pointsXY_2D_Z[i]);
    while(searchStart.z() >= maxPosZ ) searchStart.z() = searchStart.z() - 0.001;
    while(searchStart.z() <= minPosZ ) searchStart.z() = searchStart.z() + 0.001;
    searchEnd = searchStart;
    searchStart.y() = 100.0;
    
    hitPoint = point(0.0, 0.0, 0.0);
    
    // @TODO description
    pm = false, mp = false;
    
    pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        pm = true;
    }
    else{
    }
    y1 = hitPoint;
    
    // search symmetric point
    searchStart.y() = -100.0;
    pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
        hitPoint = pHit.hitPoint();
        mp = true;
    }
    else{
    }
    y2 = hitPoint;

    
    if(pm && !mp){
      y2 = y1;
      y2.y() *= -1.0;
    }
    else if(!pm && mp){
      y1 = y2;
      y1.y() *= -1.0;
    }
    
    scalar xdist = mag( x2 - x1 );
    scalar ydist = mag( y2 - y1 );
    scalar asp = 0.0;
    if(ydist==0){
      Info<<"ydist=0"<<nl;
    }else{
      asp = xdist/ydist;
    }
    apertureFile<<pointsXY_2D_Z[i]<<"  "<<asp<<"  "<<xdist<<"  "<<ydist<<"\n";
  }
  
}

void Foam::calcTypes::ellipse::write_q
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  Info << "Processing the flux qx and qy..." << endl;
  
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;

  IOobject header
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  if (header.headerOk())
  {
    fieldU field(header, mesh);

    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx = "dissolpostproc" / runTime.timeName() / "qx";
    current_file_path_qy = "dissolpostproc" / runTime.timeName() / "qy";

    autoPtr<interpolation<vector> >interpolator(interpolation<vector>::New("cellPointFace",field));

    OFstream mapXYqx( current_file_path_qx );
    OFstream mapXYqy( current_file_path_qy );

    mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
            << "   " << dx << "   " << dy << endl;
    
    mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
            << "   " << dx << "   " << dy << endl;
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    for (std::map<int,int>::iterator it  = oppositeWallUniform.begin(); 
                                     it != oppositeWallUniform.end(); ++it ){
      vector dirc =  pointsXYonsurface[it->first] - pointsXYonsurface[it->second];
      
      scalarField variable(K1_);
      scalarField Ux(K1_);
      scalarField Uy(K1_);

      scalar integratedQx = 0.0;
      scalar integratedQy = 0.0;
      if( mag(dirc) > 0.0000000001 ){
        vector dr = dirc / static_cast<scalar>(K_);

        pointField samp_points(K1_);
        
        forAll(samp_points, i){
          samp_points[i] = pointsXYonsurface[it->second] + i*dr;
          label cellI = searchEngine.findNearestCell( samp_points[i] );
          label faceI = searchEngine.findNearestFace( samp_points[i] );
          vector interp_field = interpolator->interpolate(samp_points[i], cellI, faceI);
          
          Ux[i] = interp_field.z(); // z becomes x; x becomes y
          Uy[i] = interp_field.x();
          variable[i] = samp_points[i].y();
        }
        
        integratedQx = primitive_simpson_integration(variable, Ux);
        integratedQy = primitive_simpson_integration(variable, Uy);
      }
      int ind = std::distance( oppositeWallUniform.begin(), it);
      Info<<"\r"<< 100.0 * ind / static_cast<double>(N1_*M1_) << "%          ";
      
      mapXYqx << integratedQx << "  ";
      mapXYqy << integratedQy << "  ";
      
      count++;
      if(count>=NUMBER_OF_COLUMNS){ 
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        count=0;
      }
    }  
    Info<<endl;
  }
  else{
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
}

void Foam::calcTypes::ellipse::write_p
(
  const fvMesh& mesh,
  const Time& runTime
)
{

}

void Foam::calcTypes::ellipse::write_Ccup
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  Info << "Processing the concentration field..." << endl;
  
  typedef GeometricField<vector, fvPatchField, volMesh> fieldU;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldC;

  IOobject headerU
  (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  IOobject headerC
  (
      "C",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ
  );
  
  if (headerU.headerOk() && headerC.headerOk())
  {
    fieldU field_u(headerU, mesh);
    fieldC field_c(headerC, mesh);

    fileName current_file_path_cup;
    current_file_path_cup = "dissolpostproc" / runTime.timeName() / "c";

    autoPtr<interpolation<vector> >interpolatorU(interpolation<vector>::New("cellPointFace",field_u));
    autoPtr<interpolation<scalar> >interpolatorC(interpolation<scalar>::New("cellPointFace",field_c));

    OFstream mapXYccup( current_file_path_cup );

    mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
    
    //meshSearch searchEngine(mesh, true);
    meshSearch searchEngine(mesh);
        
    int count = 0;
    for (std::map<int,int>::iterator it  = oppositeWallUniform.begin(); 
                                     it != oppositeWallUniform.end(); ++it ){
      vector dirc =  pointsXYonsurface[it->first] - pointsXYonsurface[it->second];
      
      scalarField variable(K1_);
      scalarField Ux(K1_);
      scalarField Uy(K1_);
      
      scalarField U_C3d(K1_);

      scalar integratedQx = 0.0;
      scalar integratedQy = 0.0;
      scalar integratedU_C3d = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > 0.0000000001 ){
        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        pointField samp_points(K1_);
        forAll(samp_points, i){
          samp_points[i] = pointsXYonsurface[it->second] + i*dr;
          label cellI = searchEngine.findNearestCell( samp_points[i] );
          label faceI = searchEngine.findNearestFace( samp_points[i] );
          vector interp_fieldU = interpolatorU->interpolate(samp_points[i], cellI, faceI);
          Ux[i] = interp_fieldU.z();
          Uy[i] = interp_fieldU.x();
          variable[i] = samp_points[i].y();
          
          scalar interp_fieldC = interpolatorC->interpolate(samp_points[i], cellI, faceI);
          U_C3d[i] = mag( interp_fieldU ) * interp_fieldC;
        }
        
        integratedQx = primitive_simpson_integration(variable, Ux);
        integratedQy = primitive_simpson_integration(variable, Uy);
        integratedU_C3d = primitive_simpson_integration(variable, U_C3d);
        scalar q_mag = integratedQx*integratedQx+integratedQy*integratedQy;
        if( q_mag != 0 ){
          Ccup = integratedU_C3d / std::sqrt(q_mag);
        }
      }
      int ind = std::distance(oppositeWallUniform.begin(), it);
      Info<<"\r"<< 100.0 * ind / static_cast<double>(N1_*M1_) << "%          ";
      
      mapXYccup << Ccup << "  ";
      
      count++;
      if(count>=NUMBER_OF_COLUMNS){ 
        mapXYccup <<"\n";
        count=0;
      }
    }  
    Info<<endl;
  }
  else if( !headerU.headerOk() && headerC.headerOk() ){
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
  else if( headerU.headerOk() && !headerC.headerOk() ){
    FatalError<<"There is no C field"<<nl<<nl<<exit(FatalError);
  }
  else{
    FatalError<<"There is no U and C field"<<nl<<nl<<exit(FatalError);
  }
}

// ************************************************************************* //

