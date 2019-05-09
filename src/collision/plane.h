#ifndef COLLISIONOBJECT_PLANE_H
#define COLLISIONOBJECT_PLANE_H

#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "collisionObject.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "../spring.h"

using namespace nanogui;
using namespace CGL;
using namespace std;

enum e_orientation2 { HORIZONTAL2 = 0, VERTICAL2 = 1 };

//struct Plane : public CollisionObject {
//public:
//  Plane(const Vector3D &point, const Vector3D &normal, double friction)
//      : point(point), normal(normal.unit()), friction(friction) {}
//
//  void render(GLShader &shader);
//  void collide(PointMass &pm);
//
//  Vector3D point;
//  Vector3D normal;
//
//  double friction;
//};

struct PlaneParameters {
	PlaneParameters() {}
	PlaneParameters(bool enable_structural_constraints,
		bool enable_shearing_constraints,
		bool enable_bending_constraints, double damping,
		double density, double ks)
		: enable_structural_constraints(enable_structural_constraints),
		enable_shearing_constraints(enable_shearing_constraints),
		enable_bending_constraints(enable_bending_constraints),
		damping(damping), density(density), ks(ks) {}
	~PlaneParameters() {}

	// Global simulation parameters

	bool enable_structural_constraints;
	bool enable_shearing_constraints;
	bool enable_bending_constraints;

	double damping;

	// Mass-spring parameters
	double density;
	double ks;
};

struct Plane {
	Plane() {}
	Plane(double width, double height, int num_width_points,
		int num_height_points, float thickness);
	~Plane();

	void buildGrid();



	double distance3D(double x1, double y1, double z1, double x2, double y2, double z2);


	void simulate(double frames_per_sec, double simulation_steps, PlaneParameters *cp,
		vector<Vector3D> external_accelerations,
		vector<CollisionObject *> *collision_objects);

	void reset();
	void buildClothMesh();

	void build_spatial_map();
	void self_collide(PointMass &pm, double simulation_steps);
	float hash_position(Vector3D pos);

	// Cloth properties
	double width;
	double height;
	int num_width_points;
	int num_height_points;
	double thickness;
	int windSpeed = 1000000;
	e_orientation2 orientation;

	// Cloth components
	vector<PointMass> point_masses;
	vector<PointMass> wind_masses;
	vector<vector<int>> pinned;
	vector<Spring> springs;
	ClothMesh *clothMesh;

	// Spatial hashing
	unordered_map<float, vector<PointMass *> *> map;
};

#endif /* COLLISIONOBJECT_PLANE_H */


