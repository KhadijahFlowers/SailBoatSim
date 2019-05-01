#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
	// TODO (Part 3): Handle collisions with spheres.
	//whoevers direction i want to point in goes first
	double dist1 = sqrt(pow((origin.x - pm.position.x), 2) +
		pow((origin.y - pm.position.y), 2) + pow((origin.z - pm.position.z), 2));

	bool inDist1 = dist1 < radius;

	if (inDist1) {
		Vector3D pos = pm.position - origin;
		Vector3D tangent = pos.unit() * radius + origin;

		Vector3D correctionVector = tangent - pm.last_position;
		Vector3D newPos = (pm.last_position + (1 - friction) * correctionVector);
		pm.position = newPos;
	}
	

	/*"/scene/sphere.json"*/
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
