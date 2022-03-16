#include <iostream>
#include <mpi.h>
#include "Mesh.h"

using namespace std;

int main(int argc, char **argv){

  if ( argc < 4 ){
    cerr << "Informar 3 valores: " << endl;
    cerr << argv[0] << " meshSize tipMeshSize meshFieldRadius" << endl;
    return -1;
  }

  // MPI setup
  MPI_Status status;
  MPI_Init(&argc, &argv);

  double meshSize = atof(argv[1]);
  double tipMeshSize = atof(argv[2]);
  double meshFieldRadius = atof(argv[3]);

  // Inicia um novo modelo na malha GMSH
  gmshr::Model model(meshSize, tipMeshSize, meshFieldRadius); 

  MPI_Finalize();

}
