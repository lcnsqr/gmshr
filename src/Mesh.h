#ifndef gmshr_h
#define gmshr_h

#include <vector>
#include <set>
#include <iostream>
#include <map>

#include <gmsh.h>

namespace gmshr {

  // Propriedades de uma fratura
  typedef struct {
    // id é um sufixo único para identificar as entidades
    // correspondentes à fratura: id_point, id_line, id_face.
    std::string id;
    // length define o comprimento da fratura
    double length;
    // depth define a profundidade da fratura
    double depth;
    // Altura dos elementos planos acima e abaixo da fratura
    double faceLayer;
    // faceMeshSize define o tamanho aproximado dos elementos
    // na face da fratura
    double faceMeshSize;
    // tipMeshSize define o tamanho aproximado dos elementos
    // na ponta da fratura
    double tipMeshSize;
    // Raio de refino
    double meshFieldRadius;
  } Crack;

  /**
   * Classe de modelo dentro da malha.
   */
  class Model {

    public:

    /**
     * Construtor
     */
    Model(double meshSize, double tipMeshSize, double meshFieldRadius);

    /**
     * Criar bloco simples
     */
    void block();

    /**
     * Poço 1D
     */
    void poco();

    /**
     * Salva a malha. \p filename é o arquivo onde a malha será gravada.
     */
    void save();

    /**
     * Distância euclidiana entre o ponto \p p1 e \p p2
     */
    double distance(const std::vector<double> p1, const std::vector<double> p2);

    /**
     * Deslocar e escalonar
     */
    std::vector<double> moveScale(const std::vector<double> pos, const std::vector<double> v, const double scalar);

    /**
     * Gerar malha tridimensional
     */
    void render();

    /**
     * Cria fraturas no volume definido na geometria.
     */
    void crack();

    /**
     * Identifica as entidades na malha.
     */
    void label();

    /**
     * Refino da malha. 
     */
    void refine();

    /**
     * Projeção de ponto no plano
     * \p p ponto, \p surf tag do plano, \p proj projeção, \p par coordenadas paramétricas
     */
    void projection(std::vector<double> p, int surf, std::vector<double>& proj, std::vector<double>& par);

    // Fraturas
    std::vector<gmshr::Crack> cracks;

    private:

    // Tamanho padrão de elemento 3d
    double _meshSize;

    // tipMeshSize define o tamanho aproximado dos elementos
    // na ponta da fratura
    double _tipMeshSize;
    // Raio de refino
    double _meshFieldRadius;

    /**
     * Identifica as entidades na malha com uma tolerância de distância;
     */
    void label(std::vector<std::pair<int, int>>& entities, double tol);

    /**
     * Entidades geométricas em campo de refino
     */
    int refine_layer(double d, double VIn, double VOut, std::vector<int> tags);
    int refine_ball(double radius, double VIn, double VOut, std::vector<double> center);
    int refine_cylinder(double radius, double VIn, double VOut, std::vector<double> center, std::vector<double> axis);

		// Nome do modelo
		std::string _model;

    // Tags de campos de refino
    std::vector<double> _fields;

    // Id acumulada de physical groups
    int _phy_group_id;

    // Physical group da face da fratura
    int _phy_frac_face;

    // Physical group da fenda da fratura
    int _phy_frac_open;

    // Tag do volume final
    int _vol_tag;

		// Coordenadas do bloco
		double _x_min;
		double _x_max;
		double _y_min;
		double _y_max;
		double _z_min;
		double _z_max;

    void _init_cracks();
  };

}

#endif
