#include <iostream>
#include <mpi.h>
#include "Mesh.h"

using namespace std;

gmshr::Model::Model(double meshSize, double tipMeshSize, double meshFieldRadius){

  _meshSize = meshSize;
  _tipMeshSize = tipMeshSize;
  _meshFieldRadius = meshFieldRadius;

  // Process rank
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  _init_cracks();

  // Apenas 1 processo edita a malha
  if ( rank == 0 )
  {
    gmsh::initialize(0, NULL);

    // Level of information printed on the terminal
    // 0: silent except for fatal errors
    // 1: +errors
    // 2: +warnings
    // 3: +direct
    // 4: +information
    // 5: +status
    // 99: +debug
    gmsh::option::setNumber("General.Verbosity", 5);

    _model = "block";
    
    gmsh::model::add(_model);

    // Configurações gerais da malha
    //if ( CFG._rj["mesh"].HasMember("AngleToleranceFacetOverlap") )
    //  gmsh::option::setNumber("Mesh.AngleToleranceFacetOverlap", CFG.dbl("mesh", "AngleToleranceFacetOverlap"));

    // Zerar tag de marcação de physical groups. O 
    // espaço de tags é o mesmo para todas as dimensões.
    _phy_group_id = 0;
    _phy_frac_face = 0;
    _phy_frac_open = 0;

    if ( _model == "block" ){
      block();

      // Cria fraturas no primeiro (único) volume
      crack();
      
      // Marca as entidades na malha
      label();

      // Campos de refino
      refine();

      // Gera malha e roda crack plugin
      render();
    }

    if ( _model == "poco" ) {
      block();
      poco();
    }

    // Salvar em disco
    save();

    gmsh::finalize();
  }

  // Todos os processos devem aguardar a malha
	MPI_Barrier(MPI_COMM_WORLD);

}

void gmshr::Model::_init_cracks() {

	gmshr::Crack crack;

	// id é um sufixo único para identificar as entidades
	// correspondentes à fratura: id_point, id_line, id_face.
	crack.id = "crack";
	// length define o comprimento da fratura
	crack.length = 50;
	// depth define a profundidade da fratura
	crack.depth = 15;
	// Altura dos elementos planos acima e abaixo da fratura
	crack.faceLayer = 0.5;
	// faceMeshSize define o tamanho aproximado dos elementos
	// na face da fratura
	crack.faceMeshSize = 1.25;
	// tipMeshSize define o tamanho aproximado dos elementos
	// na ponta da fratura
	crack.tipMeshSize = _tipMeshSize;
	// Raio de refino
	crack.meshFieldRadius = _meshFieldRadius;

  cracks.push_back(crack);

}

void gmshr::Model::block() {
  
  _x_min = 0;
  _y_min = 0;
  _z_min = 0;
  _x_max = 100;
  _y_max = 100;
  _z_max = 50;

  // Pontos do bloco
  vector<int> ptags;

  using namespace gmsh::model::occ;
  ptags.push_back(addPoint(_x_min, _y_min, _z_min));
  ptags.push_back(addPoint(_x_min, _y_min, _z_max));
  ptags.push_back(addPoint(_x_min, _y_max, _z_max));
  ptags.push_back(addPoint(_x_min, _y_max, _z_min));

  ptags.push_back(addPoint(_x_max, _y_min, _z_min));
  ptags.push_back(addPoint(_x_max, _y_min, _z_max));
  ptags.push_back(addPoint(_x_max, _y_max, _z_max));
  ptags.push_back(addPoint(_x_max, _y_max, _z_min));

  // Linhas
  vector<int> ltags;
  ltags.push_back(addLine(ptags[0], ptags[1]));
  ltags.push_back(addLine(ptags[1], ptags[2]));
  ltags.push_back(addLine(ptags[2], ptags[3]));
  ltags.push_back(addLine(ptags[3], ptags[0]));

  ltags.push_back(addLine(ptags[0], ptags[4]));
  ltags.push_back(addLine(ptags[1], ptags[5]));
  ltags.push_back(addLine(ptags[2], ptags[6]));
  ltags.push_back(addLine(ptags[3], ptags[7]));

  ltags.push_back(addLine(ptags[4], ptags[5]));
  ltags.push_back(addLine(ptags[5], ptags[6]));
  ltags.push_back(addLine(ptags[6], ptags[7]));
  ltags.push_back(addLine(ptags[7], ptags[4]));

  // Criar superfícies a partir das linhas
  vector<int> loops;
  loops.push_back(addCurveLoop({ltags[0], ltags[1], ltags[2], ltags[3]}));
  loops.push_back(addCurveLoop({ltags[8], ltags[9], ltags[10], ltags[11]}));
  loops.push_back(addCurveLoop({ltags[0], ltags[5], -ltags[8], -ltags[4]}));
  loops.push_back(addCurveLoop({-ltags[2], ltags[6], ltags[10], -ltags[7]}));
  loops.push_back(addCurveLoop({ltags[3], ltags[4], -ltags[11], -ltags[7]}));
  loops.push_back(addCurveLoop({-ltags[1], ltags[5], ltags[9], -ltags[6]}));

  vector<int> surfaces;
  for (int l = 0; l < 6; l++)
    surfaces.push_back(addSurfaceFilling(loops[l]));

  _vol_tag = addVolume({gmsh::model::occ::addSurfaceLoop(surfaces)});

  // Sincronizar occ e representação interna
  synchronize();

}

double gmshr::Model::distance(const vector<double> p1, const vector<double> p2)
{
  double s = 0;
  for (uint i = 0; i < p1.size(); i++)
    s += pow(p2[i] - p1[i], 2.);
  return sqrt(s);
}

vector<double> gmshr::Model::moveScale(const vector<double> pos, const vector<double> v, const double scalar)
{
  vector<double> r(v.size());
  for (uint i = 0; i < v.size(); i++)
    r[i] = pos[i] + scalar * v[i];
  return r;
}

void gmshr::Model::projection(vector<double> p, int surf, vector<double>& proj, vector<double>& par)
{
  proj.clear();
  par.clear();
  gmsh::model::getClosestPoint(2, surf, p, proj, par);
  gmsh::model::getParametrization(2, surf, proj, par);
}

void gmshr::Model::crack()
{

  gmsh::model::setCurrent(_model);

  for ( auto& crack : cracks )
  {

    // Tag da face fraturada
    //int cracked = _faces[crack.face].tag;
    int cracked = 1;

    // Placeholder para coordenadas paramétricas ao projetar ponto no plano
    vector<double> par;

    // Pontos da fenda de entrada da fratura
    vector<vector<double>> slit({
      {0, (_y_min + _y_max)/2., (_z_min + _z_max)/2. - crack.length/2.},
      {0, (_y_min + _y_max)/2., (_z_min + _z_max)/2. + crack.length/2.}
    });

    // Verificar se fratura corta todo o volume
    if ( crack.length >= (_z_max-_z_min) ){
      // Garantir extremos na fenda na fronteira
      slit[0][2] = _z_min;
      slit[1][2] = _z_max;
    }

    // Projeção do primeiro ponto da fenda no plano fraturado
    vector<double> slit_p1;
    projection(slit[0], cracked, slit_p1, par);

    // Projeção do segundo ponto da fenda no plano fraturado
    vector<double> slit_p2;
    projection(slit[1], cracked, slit_p2, par);

    // Vetor normal da face fraturada
    vector<double> cracked_normal;
    gmsh::model::getNormal(cracked, par, cracked_normal);

    // Pontos do plano de entrada da fratura
    vector<vector<double>> points({
      slit_p1,
      slit_p2,
      moveScale(slit_p2, cracked_normal, crack.depth),
      moveScale(slit_p1, cracked_normal, crack.depth)
    });

    // Coletar tags dos pontos inseridos
    vector<int> ptags;

    // Inserção dos pontos
    for ( int i = 0; i < 4; i++ ) {
      ptags.push_back(gmsh::model::occ::addPoint(points[i][0],points[i][1],points[i][2]));
    }

    // Coletar tags dos pontos da camada acima
    vector<int> ptagsUp;
    for ( int i = 0; i < 4; i++ ) {
      ptagsUp.push_back(gmsh::model::occ::addPoint(points[i][0],points[i][1]+crack.faceLayer,points[i][2]));
    }

    // Coletar tags dos pontos da camada abaixo
    vector<int> ptagsDown;
    for ( int i = 0; i < 4; i++ ) {
      ptagsDown.push_back(gmsh::model::occ::addPoint(points[i][0],points[i][1]-crack.faceLayer,points[i][2]));
    }

    vector<int> ltags;

    // Inserção de linhas
    for ( int i = 0; i < 4; i++ ) {
      int p1 = i;
      int p2 = (i+1) % 4;

      // Inserir linha na malha
      ltags.push_back(gmsh::model::occ::addLine(ptags[p1], ptags[p2]));

    }

    vector<int> ltagsUp;
    for ( int i = 0; i < 4; i++ ) {
      int p1 = i;
      int p2 = (i+1) % 4;
      ltagsUp.push_back(gmsh::model::occ::addLine(ptagsUp[p1], ptagsUp[p2]));
    }

    vector<int> ltagsDown;
    for ( int i = 0; i < 4; i++ ) {
      int p1 = i;
      int p2 = (i+1) % 4;
      ltagsDown.push_back(gmsh::model::occ::addLine(ptagsDown[p1], ptagsDown[p2]));
    }

    vector<int> loops;
    loops.push_back(gmsh::model::occ::addCurveLoop({ltags[0], ltags[1], ltags[2], ltags[3]}));
    loops.push_back(gmsh::model::occ::addCurveLoop({ltagsUp[0], ltagsUp[1], ltagsUp[2], ltagsUp[3]}));
    loops.push_back(gmsh::model::occ::addCurveLoop({ltagsDown[0], ltagsDown[1], ltagsDown[2], ltagsDown[3]}));

    vector<int> surfaces;
    surfaces.push_back(gmsh::model::occ::addSurfaceFilling(loops[0]));
    surfaces.push_back(gmsh::model::occ::addSurfaceFilling(loops[1]));
    surfaces.push_back(gmsh::model::occ::addSurfaceFilling(loops[2]));

    // ov contains all the generated entities of the same dimension as the input entities
    vector<pair<int, int>> ov;
    // ovv contains the parent-child relationships for all the input entities:
    vector<vector<pair<int, int>>> ovv;
    gmsh::model::occ::fragment({{2, surfaces[2]}}, {{3, _vol_tag}}, ov, ovv);
    gmsh::model::occ::fragment({{2, surfaces[1]}}, {{3, _vol_tag}}, ov, ovv);
    gmsh::model::occ::fragment({{2, surfaces[0]}}, {{3, _vol_tag}}, ov, ovv);

    gmsh::model::occ::synchronize();

    // Campo de refino na superfície da fratura
    //_fields.push_back(refine_layer(crack.meshFieldRadius, crack.faceMeshSize, _meshSize, surfaces));

    // Campo de refino nas pontas internas
    _fields.push_back(refine_ball(crack.meshFieldRadius, crack.tipMeshSize, _meshSize, points[2]));
    _fields.push_back(refine_ball(crack.meshFieldRadius, crack.tipMeshSize, _meshSize, points[3]));

    // Campo de refino cilíndrico no corte interno
    vector<double> axis({(points[3][0]-points[2][0])/2., (points[3][1]-points[2][1])/2., (points[3][2]-points[2][2])/2.});
    vector<double> center({points[2][0]+axis[0], points[2][1]+axis[1], points[2][2]+axis[2]});
    _fields.push_back(refine_cylinder(crack.meshFieldRadius, crack.tipMeshSize, _meshSize, center, axis));

    // Por hora, apenas uma fratura
    break;
  }

  gmsh::model::occ::synchronize();
}

void gmshr::Model::label()
{

  gmsh::model::setCurrent(_model);

  // Vetor de tags a serem marcadas
  vector<int> tags;

  // Pontos do bloco
  tags = {14,16,15,13,20,19,18,17};
  for ( int i = 0; i < tags.size(); i++ ){
    gmsh::model::addPhysicalGroup(0, {tags[i]}, ++_phy_group_id);
    gmsh::model::setPhysicalName(0, _phy_group_id, _model + "_point_" + to_string(i+1));
  }

  // Linhas do bloco
  tags = {19,20,18,17,26,28,24,22,21,27,25,23};
  for ( int i = 0; i < tags.size(); i++ ){
    gmsh::model::addPhysicalGroup(1, {tags[i]}, ++_phy_group_id);
    gmsh::model::setPhysicalName(1, _phy_group_id, _model + "_line_" + to_string(i+1));
  }

  // Faces do bloco
  tags = {8,13,9,12,11,10};
  for ( int i = 0; i < tags.size(); i++ ){
    gmsh::model::addPhysicalGroup(2, {tags[i]}, ++_phy_group_id);
    gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_" + to_string(i+1));
  }

  // Fratura
  for ( auto& crack : cracks )
  {
    // Pontos da fratura
    tags = {9,10,11,12};
    for ( int i = 0; i < tags.size(); i++ ){
      gmsh::model::addPhysicalGroup(0, {tags[i]}, ++_phy_group_id);
      gmsh::model::setPhysicalName(0, _phy_group_id, crack.id + "_point_" + to_string(i));
    }

    // Linhas da fratura

    // Verificar se fratura corta todo o volume
    if ( crack.length >= (_z_max-_z_min) ){
      gmsh::model::addPhysicalGroup(1, {13, 14, 16}, ++_phy_group_id);
      gmsh::model::setPhysicalName(1, _phy_group_id, "frac_open");
      //tags = {15};
    }
    else {
      gmsh::model::addPhysicalGroup(1, {13}, ++_phy_group_id);
      gmsh::model::setPhysicalName(1, _phy_group_id, "frac_open");
      //tags = {14,15,16};
    }
    _phy_frac_open = _phy_group_id;

    //for ( int i = 0; i < tags.size(); i++ ){
    //  gmsh::model::addPhysicalGroup(1, {tags[i]}, ++_phy_group_id);
    //  gmsh::model::setPhysicalName(1, _phy_group_id, crack.id + "_line_" + to_string(i));
    //}

    // Face da fratura
    gmsh::model::addPhysicalGroup(2, {7}, ++_phy_group_id);
    gmsh::model::setPhysicalName(2, _phy_group_id, crack.id + "_face");
    _phy_frac_face = _phy_group_id;

    // Por hora, apenas uma fratura
    break;
  }

  // Marcação do volume total
  gmsh::model::addPhysicalGroup(3, {_vol_tag}, ++_phy_group_id);
  gmsh::model::setPhysicalName(3, _phy_group_id, _model);

}

int gmshr::Model::refine_ball(double radius, double VIn, double VOut, vector<double> center)
{
  int field = gmsh::model::mesh::field::add("Ball");

  gmsh::model::mesh::field::setNumber(field, "Radius", radius);
  gmsh::model::mesh::field::setNumber(field, "Thickness", 0);
  gmsh::model::mesh::field::setNumber(field, "VIn", VIn);
  gmsh::model::mesh::field::setNumber(field, "VOut", VOut);

  gmsh::model::mesh::field::setNumber(field, "XCenter", center[0]);
  gmsh::model::mesh::field::setNumber(field, "YCenter", center[1]);
  gmsh::model::mesh::field::setNumber(field, "ZCenter", center[2]);

  return field;
}

int gmshr::Model::refine_cylinder(double radius, double VIn, double VOut, vector<double> center, vector<double> axis)
{
  int field = gmsh::model::mesh::field::add("Cylinder");

  gmsh::model::mesh::field::setNumber(field, "Radius", radius);
  gmsh::model::mesh::field::setNumber(field, "VIn", VIn);
  gmsh::model::mesh::field::setNumber(field, "VOut", VOut);

  gmsh::model::mesh::field::setNumber(field, "XCenter", center[0]);
  gmsh::model::mesh::field::setNumber(field, "YCenter", center[1]);
  gmsh::model::mesh::field::setNumber(field, "ZCenter", center[2]);

  gmsh::model::mesh::field::setNumber(field, "XAxis", axis[0]);
  gmsh::model::mesh::field::setNumber(field, "YAxis", axis[1]);
  gmsh::model::mesh::field::setNumber(field, "ZAxis", axis[2]);

  return field;
}

int gmshr::Model::refine_layer(double d, double VIn, double VOut, vector<int> tags)
{

  int distance_field = gmsh::model::mesh::field::add("Distance");

  vector<double> dtags;
  for ( auto i : tags ) dtags.push_back((double)i);

  //gmsh::model::mesh::field::setNumbers(distance_field, "PointsList", dtags);
  //gmsh::model::mesh::field::setNumbers(distance_field, "CurvesList", dtags);
  gmsh::model::mesh::field::setNumbers(distance_field, "SurfacesList", dtags);
  gmsh::model::mesh::field::setNumber(distance_field, "Sampling", 100);

  int field = gmsh::model::mesh::field::add("Threshold");
  gmsh::model::mesh::field::setNumber(field, "InField", distance_field);
  gmsh::model::mesh::field::setNumber(field, "SizeMin", VIn);
  gmsh::model::mesh::field::setNumber(field, "SizeMax", VOut);
  gmsh::model::mesh::field::setNumber(field, "DistMin", d/2.0);
  gmsh::model::mesh::field::setNumber(field, "DistMax", d);

  return field;
}

void gmshr::Model::refine()
{

  gmsh::model::setCurrent(_model);

  // Definir tamanho padrão dos elementos na malha
  gmsh::option::setNumber("Mesh.MeshSizeMax", _meshSize);
  //gmsh::option::setNumber("Mesh.MeshSizeMin", 0.01);

  // Refino nos pontos de fratura
  vector<pair<int, int>> entities;
  gmsh::model::getEntities(entities);

  double crackFaceMeshSize;
  for ( auto& c : cracks ){
    crackFaceMeshSize = c.faceMeshSize;
    break;
  }
  for(auto e : entities)
    if ( e.first == 0 && e.second < 13 )
      gmsh::model::mesh::setSize({{0 , e.second}}, crackFaceMeshSize);

  // Field refinement
  int minField = gmsh::model::mesh::field::add("Min");
  gmsh::model::mesh::field::setNumbers(minField, "FieldsList", _fields);
  gmsh::model::mesh::field::setAsBackgroundMesh(minField);

  // Definir métodos padrão de refino
  gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 1);
  gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 1);
  gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);
}

void gmshr::Model::render()
{
  gmsh::model::setCurrent(_model);

  // Gerar a malha tridimensional
  gmsh::model::mesh::generate(3);

  for ( auto& crack : cracks ) {

    // Crack plugin
    gmsh::plugin::setNumber("Crack", "Dimension", 2);

    if ( _phy_frac_face != 0 && _phy_frac_open != 0 )
    {
      gmsh::plugin::setNumber("Crack", "PhysicalGroup", _phy_frac_face);
      gmsh::plugin::setNumber("Crack", "OpenBoundaryPhysicalGroup", _phy_frac_open);
      gmsh::plugin::run("Crack");

      // Agora o grupo físico da face da fratura está vinculado
      // a duas faces: a original e a recém criada pelo crack plugin
      vector<int> crack_face_tags;
      gmsh::model::getEntitiesForPhysicalGroup(2, _phy_frac_face, crack_face_tags);

      // Remover grupo físico com as duas faces
      vector<pair<int, int>> old_phy_group;
      old_phy_group.push_back({2, _phy_frac_face});
      gmsh::model::removePhysicalGroups(old_phy_group);

      // Incluir grupo físico exclusivo para cada face
      gmsh::model::addPhysicalGroup(2, {crack_face_tags[0]}, ++_phy_group_id);
      gmsh::model::setPhysicalName(2, _phy_group_id, "frac_p");
      gmsh::model::addPhysicalGroup(2, {crack_face_tags[1]}, ++_phy_group_id);
      gmsh::model::setPhysicalName(2, _phy_group_id, "frac_m");

      if ( crack.length >= (_z_max-_z_min) ){
        // Após o corte fora-a-fora, a face de entrada está dividida
        // em duas e precisa ser agregada num único physical group

        // Remover Physical groups marcados incorretamente
        vector<pair<int, int>> old_phy_faces;
        old_phy_faces.push_back({2, 21});
        old_phy_faces.push_back({2, 24});
        old_phy_faces.push_back({2, 25});
        old_phy_faces.push_back({2, 26});
        gmsh::model::removePhysicalGroups(old_phy_faces);

        // Surface 8,10 -> metades desmarcadas -> block_face_1
        gmsh::model::addPhysicalGroup(2, {8,10}, ++_phy_group_id);
        gmsh::model::setPhysicalName(2, _phy_group_id, "block_face_1");

        gmsh::model::addPhysicalGroup(2, {11}, ++_phy_group_id);
        gmsh::model::setPhysicalName(2, _phy_group_id, "block_face_4");
        gmsh::model::addPhysicalGroup(2, {14}, ++_phy_group_id);
        gmsh::model::setPhysicalName(2, _phy_group_id, "block_face_5");
        gmsh::model::addPhysicalGroup(2, {12}, ++_phy_group_id);
        gmsh::model::setPhysicalName(2, _phy_group_id, "block_face_6");

        // O crack plugin o grupo _phy_frac_open com três curvas,
        // mas o injection point deve ser somente a curva que
        // cortou a _phy_frac_face (tag 13)
        vector<pair<int, int>> phy_fo;
        phy_fo.push_back({1, 31});
        gmsh::model::removePhysicalGroups(phy_fo);
        gmsh::model::addPhysicalGroup(1, {13}, 31);
        gmsh::model::setPhysicalName(1, 31, "frac_open");

      }

      // Marcar demais linhas da fratura
      vector<int> tags = {14,15,16};
      for ( int i = 0; i < tags.size(); i++ ){
        gmsh::model::addPhysicalGroup(1, {tags[i]}, ++_phy_group_id);
        gmsh::model::setPhysicalName(1, _phy_group_id, crack.id + "_line_" + to_string(i));
      }

    }
    
    // Apenas uma fratura por hora
    break;
  }
}

void gmshr::Model::save()
{
  gmsh::model::setCurrent(_model);

  double tetrahedra;
  gmsh::option::getNumber("Mesh.NbTetrahedra", tetrahedra);

  gmsh::write("mesh_"+to_string((int)tetrahedra)+".msh");
}

void gmshr::Model::poco() {

  gmsh::model::setCurrent(_model);

  _x_min = 0;
  _y_min = 0;
  _z_min = 0;
  _x_max = 50;
  _y_max = 50;
  _z_max = 5;

  const double res0 = 0;
  const double res1 = 0;
  const double wradius = 0;
  const double sizeOut = _meshSize;
  const double sizeIn = 0;

  // Tags dos pontos
  vector<int> ptags;

  // Início do reservatório
  ptags.push_back(gmsh::model::occ::addPoint(_x_min, _y_min, res0));
  ptags.push_back(gmsh::model::occ::addPoint(_x_max, _y_min, res0));
  ptags.push_back(gmsh::model::occ::addPoint(_x_max, _y_max, res0));
  ptags.push_back(gmsh::model::occ::addPoint(_x_min, _y_max, res0));

  // Final do reservatório
  ptags.push_back(gmsh::model::occ::addPoint(_x_min, _y_min, res1));
  ptags.push_back(gmsh::model::occ::addPoint(_x_max, _y_min, res1));
  ptags.push_back(gmsh::model::occ::addPoint(_x_max, _y_max, res1));
  ptags.push_back(gmsh::model::occ::addPoint(_x_min, _y_max, res1));

  // Linhas
  vector<int> ltags;

  ltags.push_back(gmsh::model::occ::addLine(ptags[0], ptags[1]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[1], ptags[2]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[2], ptags[3]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[3], ptags[0]));

  ltags.push_back(gmsh::model::occ::addLine(ptags[4], ptags[5]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[5], ptags[6]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[6], ptags[7]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[7], ptags[4]));

  // Criar superfícies a partir das linhas
  vector<int> loops;
  loops.push_back(gmsh::model::occ::addCurveLoop({ltags[0], ltags[1], ltags[2], ltags[3]}));
  loops.push_back(gmsh::model::occ::addCurveLoop({ltags[4], ltags[5], ltags[6], ltags[7]}));

  vector<int> surfaces;
  surfaces.push_back(gmsh::model::occ::addSurfaceFilling(loops[0]));
  surfaces.push_back(gmsh::model::occ::addSurfaceFilling(loops[1]));

  // ov contains all the generated entities of the same dimension as the input entities
  vector<pair<int, int>> ov;
  // ovv contains the parent-child relationships for all the input entities:
  vector<vector<pair<int, int>>> ovv;
  gmsh::model::occ::fragment({{2, surfaces[0]}}, {{3, 1}}, ov, ovv);
  gmsh::model::occ::fragment({{2, surfaces[1]}}, {{3, 2}}, ov, ovv);

  gmsh::model::occ::synchronize();

  // Definição da linha do poço

  // Quantidade de pontos ao redor do ponto de entrada
  int PTS = (int)(wradius / sizeIn);

  // Ângulo áureo
  double PHI = (1.0 + (1.0-sqrt(5.0))/2.0) * 2.0 * M_PI;

  // Posição angular do ponto
  double A;

  // Passo de crescimento do raio
  double RHO = wradius / (double)PTS;

  // Raio crescente
  double R;

  // Coordenada dos pontos em z
  double z;

  ptags.clear();
  for ( int i = 0; i < 4*PTS; ++i ){
    if ( i % PTS == 0 ){
      A = PHI;
      R = RHO;
    }
    if ( i / PTS == 0 ) z = _z_min;
    if ( i / PTS == 1 ) z = res0;
    if ( i / PTS == 2 ) z = res1;
    if ( i / PTS == 3 ) z = _z_max;
    ptags.push_back(gmsh::model::occ::addPoint((_x_max-_x_min)/2.0 + R*cos(A), (_y_max-_y_min)/2.0 + R*sin(A), z));
    A = A + PHI;
    R = R + RHO;
  }

  gmsh::model::occ::synchronize();

  for ( int i = 0; i < 4*PTS; ++i ){
    if ( i / PTS == 0 ) gmsh::model::mesh::embed(0, {ptags[i]}, 2, 10);
    if ( i / PTS == 1 ) gmsh::model::mesh::embed(0, {ptags[i]}, 2, 7);
    if ( i / PTS == 2 ) gmsh::model::mesh::embed(0, {ptags[i]}, 2, 8);
    if ( i / PTS == 3 ) gmsh::model::mesh::embed(0, {ptags[i]}, 2, 20);
  }

  // Pontos para as linhas do poço
  ptags.clear();
  ptags.push_back(gmsh::model::occ::addPoint((_x_max-_x_min)/2.0, (_y_max-_y_min)/2.0, _z_min));
  ptags.push_back(gmsh::model::occ::addPoint((_x_max-_x_min)/2.0, (_y_max-_y_min)/2.0, res0));
  ptags.push_back(gmsh::model::occ::addPoint((_x_max-_x_min)/2.0, (_y_max-_y_min)/2.0, res1));
  ptags.push_back(gmsh::model::occ::addPoint((_x_max-_x_min)/2.0, (_y_max-_y_min)/2.0, _z_max));

  gmsh::model::occ::synchronize();
  gmsh::model::mesh::embed(0, {ptags[0]}, 2, 10);
  gmsh::model::mesh::embed(0, {ptags[1]}, 2, 7);
  gmsh::model::mesh::embed(0, {ptags[2]}, 2, 8);
  gmsh::model::mesh::embed(0, {ptags[3]}, 2, 20);

  // Fundir linhas do poço aos volumes
  ltags.clear();
  ltags.push_back(gmsh::model::occ::addLine(ptags[0], ptags[1]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[1], ptags[2]));
  ltags.push_back(gmsh::model::occ::addLine(ptags[2], ptags[3]));

  gmsh::model::occ::synchronize();

  gmsh::model::mesh::embed(1, {ltags[0]}, 3, 1);
  gmsh::model::mesh::embed(1, {ltags[1]}, 3, 2);
  gmsh::model::mesh::embed(1, {ltags[2]}, 3, 3);

  // Pontos nas linhas de poço modificam refino
  gmsh::model::mesh::setTransfiniteCurve(ltags[0], (res0  - _z_min) / sizeIn);
  gmsh::model::mesh::setTransfiniteCurve(ltags[1], (res1  - res0 ) / sizeIn);
  gmsh::model::mesh::setTransfiniteCurve(ltags[2], (_z_max - res1 ) / sizeIn);

  // Marcar nomes de regiões e entidades

  // Faces do bloco
  gmsh::model::addPhysicalGroup(2, {9,14,18}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_1");

  gmsh::model::addPhysicalGroup(2, {13,17,22}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_2");

  gmsh::model::addPhysicalGroup(2, {11,15,19}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_3");

  gmsh::model::addPhysicalGroup(2, {12,16,21}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_4");

  gmsh::model::addPhysicalGroup(2, {10}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_5");

  gmsh::model::addPhysicalGroup(2, {20}, ++_phy_group_id);
  gmsh::model::setPhysicalName(2, _phy_group_id, _model + "_face_6");

  gmsh::model::addPhysicalGroup(0, {ptags[3]}, ++_phy_group_id);
  gmsh::model::setPhysicalName(0, _phy_group_id, "injection_point");

  gmsh::model::addPhysicalGroup(1, {ltags[0]}, ++_phy_group_id);
  gmsh::model::setPhysicalName(1, _phy_group_id, "well0");
  gmsh::model::addPhysicalGroup(1, {ltags[1]}, ++_phy_group_id);
  gmsh::model::setPhysicalName(1, _phy_group_id, "well1");
  gmsh::model::addPhysicalGroup(1, {ltags[2]}, ++_phy_group_id);
  gmsh::model::setPhysicalName(1, _phy_group_id, "well2");

  gmsh::model::addPhysicalGroup(3, {1}, ++_phy_group_id);
  gmsh::model::setPhysicalName(3, _phy_group_id, "res0");
  gmsh::model::addPhysicalGroup(3, {2}, ++_phy_group_id);
  gmsh::model::setPhysicalName(3, _phy_group_id, "res1");
  gmsh::model::addPhysicalGroup(3, {3}, ++_phy_group_id);
  gmsh::model::setPhysicalName(3, _phy_group_id, "res2");

  // Campo de refino
  int field = gmsh::model::mesh::field::add("Cylinder");

  gmsh::model::mesh::field::setNumber(field, "Radius", wradius);
  gmsh::model::mesh::field::setNumber(field, "VIn", sizeIn);
  gmsh::model::mesh::field::setNumber(field, "VOut", sizeOut);

  gmsh::model::mesh::field::setNumber(field, "XCenter", (_x_max-_x_min)/2.0);
  gmsh::model::mesh::field::setNumber(field, "YCenter", (_y_max-_y_min)/2.0);
  gmsh::model::mesh::field::setNumber(field, "ZCenter", (_z_max-_z_min)/2.0);

  gmsh::model::mesh::field::setNumber(field, "XAxis", 0);
  gmsh::model::mesh::field::setNumber(field, "YAxis", 0);
  gmsh::model::mesh::field::setNumber(field, "ZAxis", (_z_max-_z_min)/2.0);

  //gmsh::model::mesh::field::setAsBackgroundMesh(field);

  // Gerar a malha tridimensional
  gmsh::model::mesh::generate(3);
}

