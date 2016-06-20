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
    fieldMap2D.C

\*---------------------------------------------------------------------------*/

#define NUMBER_OF_COLUMNS 10

#include "fieldMap2D.H"

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
        defineTypeNameAndDebug(fieldMap2D, 0);
        addToRunTimeSelectionTable(calcType, fieldMap2D, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::fieldMap2D::fieldMap2D()
:
    calcType()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::fieldMap2D::~fieldMap2D()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::fieldMap2D::init()
{
  argList::validArgs.append("fieldMap2D");
  argList::validArgs.append("processingType");
  argList::validArgs.append("N");
  argList::validArgs.append("M");
  argList::validArgs.append("K");
}

void Foam::calcTypes::fieldMap2D::preCalc
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
  
  N1M1 = N1_*M1_;

  // need temporary arguments in order to add a 0 time directory
  argList argsTmp = args;
  argsTmp.setOption("time", "0");
  Foam::Time timeTmp(Foam::Time::controlDictName, argsTmp);
  Foam::instantList timeDirs = Foam::timeSelector::select0(timeTmp, argsTmp);
  timeTmp.setTime(timeDirs[0], 0);
  
  Foam::fvMesh meshTmp
  (
      Foam::IOobject
      (
          Foam::fvMesh::defaultRegion,
          timeTmp.timeName(),
          timeTmp,
          Foam::IOobject::MUST_READ
      )
  );
  if( timeTmp.timeName()!="0" ){
    SeriousErrorIn("fieldOperations::getInletAreaT0")
            <<"There is no 0 time directory. Check your decomposition as well!"<<nl
            <<"TimeTmp: "<< timeTmp.timeName()<<nl
            <<exit(FatalError);
  }
  
  Info << "Calculating the max and min coordinates"<<nl;
  const pointField &allPoints = meshTmp.points();
  
  maxPosX = max( allPoints.component(vector::X) );
  maxPosY = max( allPoints.component(vector::Y) );
  maxPosZ = max( allPoints.component(vector::Z) );
  
  minPosX = min( allPoints.component(vector::X) );
  minPosY = min( allPoints.component(vector::Y) );
  minPosZ = min( allPoints.component(vector::Z) );
  
  
  // z becomes x; x becomes y
  dx = (maxPosZ - minPosZ) / static_cast<scalar>(N_);
  dy = (maxPosX - minPosX) / static_cast<scalar>(M_);
  
  // calculate the size of the current pattern being not bigger then 10^7
  thisTimeSize = min(N1_*M1_, maxNumProcPoints);
  curNum = 0;
  curEnd = thisTimeSize;
  
  totNumLoop = (N1_*M1_-1) / maxNumProcPoints + 1;
}


// TODO make it more general
void Foam::calcTypes::fieldMap2D::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
  Info << "Find the points on the surface"<<nl;
  // coordinates of points on the surface of the fracture walls
  //memInfo mf;
  
  // print what is calculated
  if(processingType_ == "all")
    Info << "Processing all fields..." << endl;
  else if(processingType_ == "h")
    Info << "Processing the  aperture of the fracture..." << endl;
  else if(processingType_ == "U")
    Info << "Processing the flux qx and qy..." << endl;
  else if(processingType_ == "p")
    FatalError<<"p processing is not implemented yet "<<nl<<nl<<exit(FatalError);
  else if(processingType_ == "C")
    Info << "Processing the concentration field..." << endl;
  else if(processingType_ == "temp")
    Info << "Running temporary function..." << endl;
  else
    FatalError<<"Unable to process "<<processingType_<<nl<<nl<<exit(FatalError);
  
  for(int cI=0; cI<totNumLoop; cI++){
    
    curNum = cI;
    curBlock = thisTimeSize * curNum;
    
    sizeAA = thisTimeSize;
    if(cI==totNumLoop-1){
      sizeAA = N1M1 - (totNumLoop-1)*thisTimeSize;
    }

    pointsXYonsurface.clear();
    pointsXYonsurface.setSize( 2.0 * sizeAA );
    build_surface_points( mesh );

    fileName current_dissolpostproc_dir;
    current_dissolpostproc_dir = "dissolpostproc" / runTime.timeName();
    if ( !isDir(current_dissolpostproc_dir) ) mkDir(current_dissolpostproc_dir);

    if(processingType_ == "all"){
      write_all(mesh, runTime);
    }
    else if(processingType_ == "h"){
      write_h(mesh, runTime);
    }
    else if(processingType_ == "U"){
      write_q(mesh, runTime);
    }
    else if(processingType_ == "p"){
      //write_p(mesh, runTime);
      //FatalError<<"p processing is not implemented yet "<<nl<<nl<<exit(FatalError);
    }
    else if(processingType_ == "C"){
      write_Ccup(mesh, runTime);
    }
    else if(processingType_ == "temp"){
      write_temp(mesh, runTime);
    }
    else{
      //FatalError<<"Unable to process "<<processingType_<<nl<<nl<<exit(FatalError);
    }
  }
}

scalar Foam::calcTypes::fieldMap2D::primitive_simpson_integration
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

  scalarField hsum(len_lab_list);     //= h0 + h1
  scalarField hprod(len_lab_list);    //= h0 * h1
  scalarField h0divh1(len_lab_list);  //= h0 / h1

  
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



void Foam::calcTypes::fieldMap2D::build_surface_points
(
  const fvMesh& mesh
)
{
  const pointField &allPoints = mesh.points();
  maxPosY = max( allPoints.component(vector::Y) );
  minPosY = min( allPoints.component(vector::Y) );
  
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
  const triSurfaceSearch querySurf(wallTriSurface);
  const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();
  
  scalar curX = 0.0;
  scalar curZ = 0.0;
  
  // N+1 and M+1 because we need points at x=minPosX and x=maxPosX
  for(int ijk=0; ijk<sizeAA; ijk++)
  {
    int totIJK = ijk + curNum * thisTimeSize;
    label i = totIJK / M1_;
    label j = totIJK % M1_;
    curZ = minPosZ + i*dx;
    curX = minPosX + j*dy;
    label ind = ijk;
      
    point searchStart(curX, 0.0, curZ);
    while(searchStart.x() >= maxPosX ) searchStart.x() = searchStart.x() - 0.001;
    while(searchStart.x() <= minPosX ) searchStart.x() = searchStart.x() + 0.001;
    while(searchStart.z() >= maxPosZ ) searchStart.z() = searchStart.z() - 0.001;
    while(searchStart.z() <= minPosZ ) searchStart.z() = searchStart.z() + 0.001;
    point searchEnd = searchStart;
    searchStart.y() = maxPosY+1.0;

    point hitPoint(0.0, 0.0, 0.0);

    label posWallLab = 2*ind;
    label negWallLab = 2*ind+1;

    // @TODO description
    bool pm = false, mp = false;

    pointIndexHit pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
      hitPoint = pHit.hitPoint();
      pm = true;
    }

    pointsXYonsurface[posWallLab] = hitPoint;

    // search symmetric point
    searchStart.y() = minPosY-1.0;
    pHit = tree.findLine(searchStart, searchEnd);
    if ( pHit.hit() )
    {
      hitPoint = pHit.hitPoint();
      mp = true;
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
void Foam::calcTypes::fieldMap2D::write_temp
(
  const fvMesh& mesh,
  const Time& runTime
)
{
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

void Foam::calcTypes::fieldMap2D::write_h
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  fileName current_file_path;
  current_file_path = "dissolpostproc" / runTime.timeName() / "h";
  
  ios_base::openmode mode = (curNum==0) ? ios_base::out|ios_base::trunc : ios_base::out|ios_base::app;  
  OFstreamMod apertureMapFile( current_file_path, mode);
  
  if( curNum==0 ){
    apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
  }

  int count = 0;
  for (std::map<int,int>::const_iterator it  = oppositeWallUniform.begin(); 
                                   it != oppositeWallUniform.end(); ++it ){
    scalar dist = mag( pointsXYonsurface[it->first] - pointsXYonsurface[it->second] );
    apertureMapFile << dist << "  ";

    count++;
    if(count>=NUMBER_OF_COLUMNS){ 
      apertureMapFile <<"\n";
      count=0;
    }
  }
  
}

void Foam::calcTypes::fieldMap2D::write_q
(
  const fvMesh& mesh,
  const Time& runTime
)
{
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
    
    ios_base::openmode mode = (curNum==0) ? ios_base::out|ios_base::trunc : ios_base::out|ios_base::app;  
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );

    if( curNum==0 ){
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;

      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
    }
    
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
  }
  else{
    FatalError<<"There is no U field"<<nl<<nl<<exit(FatalError);
  }
}

void Foam::calcTypes::fieldMap2D::write_p
(
  const fvMesh& mesh,
  const Time& runTime
)
{
  Info<< "Not implemented yet"<<nl;
}

void Foam::calcTypes::fieldMap2D::write_Ccup
(
  const fvMesh& mesh,
  const Time& runTime
)
{
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

    ios_base::openmode mode = (curNum==0) ? ios_base::out|ios_base::trunc : ios_base::out|ios_base::app;  
    OFstreamMod mapXYccup( current_file_path_cup, mode );

    if(curNum==0){
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
    }
    
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



void Foam::calcTypes::fieldMap2D::write_all
(
  const fvMesh& mesh,
  const Time& runTime
)
{
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
    fileName current_file_path_qx, current_file_path_qy;
    current_file_path_qx = "dissolpostproc" / runTime.timeName() / "qx";
    current_file_path_qy = "dissolpostproc" / runTime.timeName() / "qy";
    fileName current_file_path_h;
    current_file_path_h = "dissolpostproc" / runTime.timeName() / "h";
  
    autoPtr<interpolation<vector> >interpolatorU(interpolation<vector>::New("cellPoint",field_u));
    autoPtr<interpolation<scalar> >interpolatorC(interpolation<scalar>::New("cellPoint",field_c));

    ios_base::openmode mode = (curNum==0) ? ios_base::out|ios_base::trunc : ios_base::out|ios_base::app;  
    OFstreamMod mapXYccup( current_file_path_cup, mode );
    OFstreamMod mapXYqx( current_file_path_qx, mode );
    OFstreamMod mapXYqy( current_file_path_qy, mode );
    OFstreamMod apertureMapFile( current_file_path_h, mode);
    
    if(curNum==0){
      mapXYccup << N_ << "   " << M_ << "   " << runTime.value() 
            << "   " << dx << "   " << dy << endl;
      mapXYqx << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;

      mapXYqy << N_ << "   " << M_ << "   " << runTime.value()
              << "   " << dx << "   " << dy << endl;
      
      apertureMapFile << N_ << "   " << M_ << "   " << runTime.value() 
              << "   " << dx << "   " << dy << endl;
    }
    
    meshSearch searchEngine(mesh);
        
    int count = 0;
    for(
            std::map<int,int>::iterator it  = oppositeWallUniform.begin(); 
                                        it != oppositeWallUniform.end();
                                      ++it 
       )
    {
      vector dirc =  pointsXYonsurface[it->first] - pointsXYonsurface[it->second];
      scalar dist = mag( dirc );
      
      scalar integratedQx = 0.0;
      scalar integratedQy = 0.0;
      scalar integratedU_C3d = 0.0;
      scalar Ccup = 0.0;
      if( mag(dirc) > SMALL ){
        scalarField variable(K1_);
        scalarField Ux(K1_);
        scalarField Uy(K1_);

        scalarField U_C3d(K1_);

        vector dr = 1.0 / static_cast<scalar>(K_) * dirc;

        for(int i=0; i<K1_; i++){
          point samp_point = pointsXYonsurface[it->second] + i*dr;
          
          label cellI = searchEngine.findCell( samp_point );
          
          if (cellI==-1) cellI = searchEngine.findNearestCell( samp_point );
          label faceI = -1;
          vector interp_fieldU = interpolatorU->interpolate(samp_point, cellI, faceI);
          
          if(i==0 || i==K1_-1){
            interp_fieldU = vector::zero;
          }
          
          Ux[i] = interp_fieldU.z();
          Uy[i] = interp_fieldU.x();
          variable[i] = samp_point.y();
          
          scalar interp_fieldC = interpolatorC->interpolate(samp_point, cellI, faceI);
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
      int ind1 = ind + curBlock;
      int prcnt  = 100 * (ind1    ) / N1M1;
      int prcnt1 = 100 * (ind1 - 1) / N1M1;
      if( prcnt%1==0 && prcnt!=prcnt1 ){
        Info<<"\r"<< prcnt << "%  "<<flush;
      }
      
      mapXYccup << Ccup << "  ";
      mapXYqx << integratedQx << "  ";
      mapXYqy << integratedQy << "  ";
      apertureMapFile << dist << "  ";

      count++;
      if(count>=NUMBER_OF_COLUMNS){ 
        mapXYccup <<"\n";
        mapXYqx <<"\n";
        mapXYqy <<"\n";
        apertureMapFile <<"\n";
        count=0;
      }
    }  
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

