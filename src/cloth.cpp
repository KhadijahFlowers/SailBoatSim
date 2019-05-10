#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
	int num_height_points, float thickness) {
	this->width = width;
	this->height = height;
	this->num_width_points = num_width_points;
	this->num_height_points = num_height_points;
	this->thickness = thickness;

	buildGrid();
	buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

//PART 1
//c++ map insert pait
void Cloth::buildGrid() {
	// TODO (Part 1): Build a grid of masses and springs.
	  /*Start by creating an evenly spaced grid of masses.
	  There should be num_width_points by num_height_points
	  total masses spanning width and height lengths,
	  respectively.
	  */


	for (int x = 0; x < num_width_points; x++) {
		for (int y = 0; y < num_height_points; y++) {
			
			float pn = rand();
			float pnq;
			if (pn <= RAND_MAX / 2) {
				pnq = 1.f;
			}
			else {
				pnq = -1.f;
			}
			float z = pnq * rand() / RAND_MAX / 1000.f;
			

			Vector3D pos = Vector3D(x * width / num_width_points, y * height / num_height_points + 1, z);
			PointMass newGuy = PointMass(pos, false);
			newGuy.LeftMost = false;
			newGuy.RightMost = false;
			point_masses.push_back(newGuy);
			
		}
	}


	//within pinned? Always pinned for now

	for (int x = 0; x < num_width_points; x++) {
		for (int y = 0; y < num_height_points; y++) {
			if (x == 0 && y == 0) {
				//bottom left
				point_masses.at(num_width_points * y + x).pinned = true;
				point_masses.at(num_width_points * y + x).LeftMost = true;
			}
			else if (x == (num_width_points - 1) && y == 0) {
				//top left
				point_masses.at(num_width_points * y + x).pinned = true;
				point_masses.at(num_width_points * y + x).LeftMost = true;
			}
			else if (x == (num_width_points - 1) && y == (num_height_points - 1)) {
				//top right
				point_masses.at(num_width_points * y + x).pinned = true;
				point_masses.at(num_width_points * y + x).RightMost = true;
			}
			else if (x == 0 && y == (num_height_points - 1)) {
				//bottom right
				point_masses.at(num_width_points * y + x).pinned = true;
				point_masses.at(num_width_points * y + x).RightMost = true;
			}
		}
	}


	//structural loop
	for (int pm = 0; pm < point_masses.size(); pm++) {
		int leftRow, leftCol, topRow, topCol, row, col;
		row = pm / num_width_points;
		col = pm % num_width_points;
		leftRow = (pm - 1) / num_width_points;
		leftCol = (pm - 1) % num_width_points;
		topRow = (pm - num_width_points) / num_width_points;
		topCol = (pm - num_width_points) % num_width_points;

		if (!(leftRow < 0) && !(leftRow != row) && !(abs(col - leftCol) != 1) && (num_width_points * leftRow + leftCol < point_masses.size())) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * leftRow + leftCol), STRUCTURAL));
		}
		if (!(topRow < 0) && !(abs(row - topRow) != 1) && !(abs(col - topCol) != 0) && (num_width_points * topRow + topCol) < point_masses.size()) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * topRow + topCol), STRUCTURAL));
		}
	}
	
	//shear loop
	for (int pm = 0; pm < point_masses.size(); pm++) {
		int topLeftRow, topLeftCol, topRightRow, topRightCol, row, col;
		row = pm / num_width_points;
		col = pm % num_width_points;
		topLeftRow = (pm - num_width_points - 1) / num_width_points;
		topLeftCol = (pm - num_width_points - 1) % num_width_points;
		topRightRow = (pm - num_width_points + 1) / num_width_points;
		topRightCol = (pm - num_width_points + 1) % num_width_points;

		if (!(topLeftRow < 0) && !(abs(topLeftRow - row) != 1) && !(abs(col - topLeftCol) != 1) && (num_width_points * topLeftRow + topLeftCol < point_masses.size())) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * topLeftRow + topLeftCol), SHEARING));
		}
		if (!(topRightRow < 0) && !(abs(row - topRightRow) != 1) && !(abs(col - topRightCol) != 1) && (num_width_points * topRightRow + topRightCol < point_masses.size())) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * topRightRow + topRightCol), SHEARING));
		}
	}

	//bending
	for (int pm = 0; pm < point_masses.size(); pm++) {
		int leftRow, leftCol, topRow, topCol, row, col;
		row = pm / num_width_points;
		col = pm % num_width_points;
		leftRow = (pm - 2) / num_width_points;
		leftCol = (pm - 2) % num_width_points;
		topRow = (pm - 2 * num_width_points) / num_width_points;
		topCol = (pm - 2 * num_width_points) % num_width_points;

		if (!(leftRow < 0) && !(leftRow != row) && !(abs(col - leftCol) != 1) && (num_width_points * leftRow + leftCol < point_masses.size())) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * leftRow + leftCol), BENDING));
		}
		if (!(topRow < 0) && !(abs(row - topRow) != 2) && !(abs(col - topCol) != 0) && (num_width_points * leftRow + leftCol < point_masses.size())) {
			springs.emplace_back(Spring(&point_masses.at(pm), &point_masses.at(num_width_points * topRow + topCol), BENDING));
		}
	}

}


double Cloth::distance3D(double x1, double y1, double z1, double x2, double y2, double z2) {
	return (double) sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
}



void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
	vector<Vector3D> external_accelerations,
	vector<CollisionObject *> *collision_objects) {

	double mass = width * height * cp->density / num_width_points / num_height_points;
	double delta_t = 1.0f / frames_per_sec / simulation_steps;

	for (int pm = 0; pm < point_masses.size(); pm++) {
		point_masses.at(pm).forces = Vector3D(0, 0, 0);
	}

	// TODO (Part 2): Compute total force acting on each point mass.

	//apply external force to every point mass
	vector<Vector3D> forcesVector = vector<Vector3D>();
	for (int acc = 0; acc < external_accelerations.size(); acc++) {
		Vector3D f = external_accelerations.at(acc) * mass;
		forcesVector.push_back(f);
	}

	for (int pm = 0; pm < point_masses.size(); pm++) {
		for (int i = 0; i < forcesVector.size(); i++) {
			if (!point_masses.at(pm).pinned) {
				point_masses.at(pm).forces += forcesVector.at(i);
			}
		}
	}

	//apply spring correction forces
	for (int sp = 0; sp < springs.size(); sp++) {
		Vector3D vect = (springs.at(sp).pm_a->position - springs.at(sp).pm_b->position);

		double norm = vect.norm();
		double l = springs.at(sp).rest_length;
		double ks = cp->ks;
		double fs = ks * (norm - l);

		if (cp->enable_structural_constraints) {
			Vector3D f1 = (springs.at(sp).pm_a)->forces - fs * vect.unit();
			Vector3D f2 = (springs.at(sp).pm_b)->forces + fs * vect.unit();
			(springs.at(sp).pm_a)->forces = f1;
			(springs.at(sp).pm_b)->forces = f2;
		}
		else if (cp->enable_shearing_constraints) {
			Vector3D f1 = (springs.at(sp).pm_a)->forces - fs * vect.unit();
			Vector3D f2 = (springs.at(sp).pm_b)->forces + fs * vect.unit();
			(springs.at(sp).pm_a)->forces = f1;
			(springs.at(sp).pm_b)->forces = f2;
		}
		else if (cp->enable_bending_constraints) {
			Vector3D f1 = (springs.at(sp).pm_a)->forces - fs * vect.unit(); //0.2
			Vector3D f2 = (springs.at(sp).pm_b)->forces + fs * vect.unit();
			(springs.at(sp).pm_a)->forces = f1;
			(springs.at(sp).pm_b)->forces = f2;
		}


	}
	//Ideas: check accelerations
	//

	// TODO (Part 2): Use Verlet integration to compute new point mass positions
	// TODO (Part 2): Use Verlet integration to compute new point mass positions
	for (int pm = 0; pm < point_masses.size(); pm++) {


		if (point_masses.at(pm).pinned) {
			if (rotationChanges < 6000) {

				Vector3D gravity = external_accelerations.at(0);
				//WIND SPEED 10000000
				rotationChanges++;
				Vector3D n, s, e, w;
				n = Vector3D(0, 0, -0.0001);
				e = Vector3D(0.0001, 0, 0);
				s = Vector3D(0, 0, 0.0001);
				w = Vector3D(-0.0001, 0, 0);

				if (dot(lastDir, gravity) != 0 || lastDir.y == -9.8) {

					point_masses.at(pm).position = point_masses.at(pm).position + external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm();
					continue;
				}
				else {

					//NORTH
					if (lastDir.z == -9.8) {
						//W
						if (gravity.x == -9.8) {
							if (point_masses.at(pm).RightMost) {
								
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + w + 5 * n;
								continue;
							}
							else if (point_masses.at(pm).LeftMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * e;
								continue;
							}
						}
						//E
						else if (gravity.x == 9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * w;
								continue;
							}
							else
								if (point_masses.at(pm).LeftMost) {
									point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + e + 5 * n;
									continue;
								}
						}
					}
					//EAST
					else if (lastDir.x == 9.8) {
						//N
						if (gravity.z == -9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * e + n;
								continue;
							}
							else if (point_masses.at(pm).LeftMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * s;
								continue;
							}
						}
						//S
						else if (gravity.z == 9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * w + s;
								continue;
							}
							else
								if (point_masses.at(pm).LeftMost) {
									point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * s;
									continue;
								}
						}
					}
					//SOUTH
					else if (lastDir.z == 9.8) {
						//E
						if (gravity.x == 9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + e + 5 * s;
								continue;
							}
							else if (point_masses.at(pm).LeftMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * w;
								continue;
							}
						}
						//W
						else if (gravity.x == -9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * e;
								continue;
							}
							else
								if (point_masses.at(pm).LeftMost) {
									point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + w + 5 * s;
									continue;
								}
						}
					}
					//WEST
					else if (lastDir.x == -9.8) {
						//N
						if (gravity.z == -9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * e;
								continue;
							}
							else
								if (point_masses.at(pm).LeftMost) {
									point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * n;
									continue;
								}

							//S
						}
						else if (gravity.z == 9.8) {
							if (point_masses.at(pm).RightMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + s + 5 * w;
								continue;
							}
							else if (point_masses.at(pm).LeftMost) {
								point_masses.at(pm).position = point_masses.at(pm).position + (external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm()) + 5 * e;
								continue;
							}
						}
					}
					
				}
			}
			else {
				point_masses.at(pm).position = point_masses.at(pm).position + external_accelerations.at(0) / windSpeed * external_accelerations.at(0).norm();
				continue;
			}
		}
		else {

			Vector3D newPos = point_masses[pm].position + (1.0 - (cp->damping / 100.0)) *
				(point_masses[pm].position - point_masses[pm].last_position)
				+ (external_accelerations.at(0) + (point_masses[pm].forces / mass)) * pow(delta_t, 2);


			for (int i = 1; i < external_accelerations.size(); i++) {
				newPos += point_masses[pm].position + (1.0 - (cp->damping / 100.0)) *
					(point_masses[pm].position - point_masses[pm].last_position)
					+ (external_accelerations.at(i) + (point_masses[pm].forces / mass)) * pow(delta_t, 2);
			}


			Vector3D oldPos = point_masses[pm].position;
			point_masses.at(pm).last_position = oldPos;

			point_masses.at(pm).position = newPos;


		}
	}



	// TODO (Part 3): Handle collisions with other primitives.
	/*for (int pm = 0; pm < point_masses.size(); pm++) {
		for (int c = 0; c < collision_objects->size(); c++) {
			collision_objects->at(c)->collide(point_masses.at(pm));
		}
	}*/



	// TODO (Part 2): Constrain the changes to be such that the spring does not change
	// in length more than 10% per timestep [Provot 1995].
	/*for (int sp = 0; sp < springs.size(); sp++) {

		Vector3D direction = (*springs.at(sp).pm_a).position - (*springs[sp].pm_b).position;
		direction = direction.unit();
		double atMost = springs.at(sp).rest_length * 1.1;

		double at = distance3D((*springs.at(sp).pm_a).position.x, (*springs.at(sp).pm_a).position.y, (*springs.at(sp).pm_a).position.z,
			(*springs.at(sp).pm_b).position.x, (*springs.at(sp).pm_b).position.y, (*springs.at(sp).pm_b).position.z);
		double lastAt = 0;
		
			double diff = at - atMost;
		
			if ((*springs.at(sp).pm_a).pinned && !(*springs.at(sp).pm_b).pinned) {
				(*springs.at(sp).pm_b).position += diff * direction;
			}
			else if (!(*springs.at(sp).pm_a).pinned && (*springs.at(sp).pm_b).pinned) {
				(*springs.at(sp).pm_a).position -= diff * direction;
			}
			else {
				(*springs.at(sp).pm_b).position += (diff / 2.0) * direction;
				(*springs.at(sp).pm_a).position -= (diff / 2.0) * direction;
			}
			at = distance3D((*springs.at(sp).pm_a).position.x, (*springs.at(sp).pm_a).position.y, (*springs.at(sp).pm_a).position.z,
				(*springs.at(sp).pm_b).position.x, (*springs.at(sp).pm_b).position.y, (*springs.at(sp).pm_b).position.z);

	}*/

}















//DO NOT NEED!!!!!!!!!!!!!!!!!!!!!!
void Cloth::build_spatial_map() {
	for (const auto &entry : map) {
		delete(entry.second);
	}
	map.clear();

	
	//map[float you found] = vector
	//map[f].emplace_back(pm)
	// TODO (Part 4): Build a spatial map out of all of the point masses.
	for (int pm = 0; pm < point_masses.size(); pm++) {

		float box = hash_position(point_masses[pm].position);
		
		if (map[box] == NULL) {
			map[box] = new vector<PointMass*>();
			//issue
			map[box]->push_back(&point_masses[pm]);
		} else{
			map[box]->push_back(&point_masses[pm]);
		}
	}
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
	// TODO (Part 4): Handle self-collision for a given point mass.
	Vector3D finalCorrection = Vector3D(0, 0, 0);
	int howManyCorrections = 0;

	double ref = hash_position(pm.position);
	vector<PointMass*>* guy = map[ref];

	for (PointMass* neighbor : *guy) {

		//PointMass *neighbor = (*guy).at(i);

		//distance is bad, not going in
		if (neighbor != &pm) {

			Vector3D direction = (pm.position - (*neighbor).position);
			float dist = direction.norm();
			if (dist < 2 * thickness) {
				finalCorrection = finalCorrection + (2 * thickness - dist) * direction.unit();
				howManyCorrections += 1;
			}
		}
	
	}
	
	if (howManyCorrections != 0) {
		pm.position = pm.position + (finalCorrection / (double)howManyCorrections) / simulation_steps;
	}
	
}

float Cloth::hash_position(Vector3D pos) {
	// TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
	float w = num_width_points;
	float h = num_height_points;
	float t = max(w, h);
	float dimensions = w * h * t;

	float newXBox = fmod(pos.x * 3, w * width);
	float newYBox = fmod(pos.y * 3, h * height);
	float newZBox;

	if (h > w) {
		newZBox = fmod(pos.z * 3, w * width);
	}
	else {
		newZBox = fmod(pos.z * 3, h * height);
	}
	newZBox = fmod(newZBox, 10);
	string s1 = std::to_string((int)newXBox);
	string s2 = std::to_string((int)newYBox);
	string s3 = std::to_string((int)newZBox);

	string ultimate = s1 + s2 + s3;
	float ret = std::stof(s1 + s2 + s3, NULL);
	//turn to string
	//concatenate
	//this is number for box
	return ret;
}

//DO NOT NEED!!!!!!!!!!!!!!!!!!!!!!!!










///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}

