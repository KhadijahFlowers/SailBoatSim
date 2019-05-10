#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
  void BezierCurve::evaluateStep()
  {
    // TODO Part 1.
    // Perform one step of the Bezier curve's evaluation at t using de Casteljau's algorithm for subdivision.
    // Store all of the intermediate control points into the 2D vector evaluatedLevels.
	  std::vector<Vector2D> result;
	  std::vector<Vector2D>& lastLevel = evaluatedLevels[evaluatedLevels.size() - 1];
	  //t;
	  for (int i = 0; i < lastLevel.size() - 1; ++i) {
		  Vector2D p_i = lastLevel[i];
		  Vector2D p_i_1 = lastLevel[i+1];
		  result.push_back( (1 - t)*p_i + t * p_i_1 );		
	  }
	  evaluatedLevels.push_back(result);
    return;
  }


  Vector3D BezierPatch::evaluate(double u, double v) const
  {
    // TODO Part 2.
    // Evaluate the Bezier surface at parameters (u, v) through 2D de Casteljau subdivision.
    // (i.e. Unlike Part 1 where we performed one subdivision level per call to evaluateStep, this function
    // should apply de Casteljau's algorithm until it computes the final, evaluated point on the surface)

	  std::vector<Vector3D> q;
	  for (int row = 0; row < controlPoints.size(); ++row) {
		  q.push_back(evaluate1D(controlPoints[row], u));
	  }
	  
	  return evaluate1D(q, v);
  }

  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> points, double t) const
  {
    // TODO Part 2.
    // Optional helper function that you might find useful to implement as an abstraction when implementing BezierPatch::evaluate.
    // Given an array of 4 points that lie on a single curve, evaluates the Bezier curve at parameter t using 1D de Casteljau subdivision.
	  std::vector<Vector3D> first = points;
	  std::vector<Vector3D> result = points;
	  // Run until only 1 point left
	  while (result.size() > 1) {
		  result.clear(); //this is like the evaluateStep() func p1
		  for (int i = 0; i < first.size() - 1; ++i) {
			  Vector3D p_i = first[i];
			  Vector3D p_i_1 = first[i + 1];
			  result.push_back((1 - t)*p_i + t * p_i_1);
		  }
		  // End of func p1. One less point was added. Recursion point to repeat
		  first = result;
	  }	  
	  
	  return result[result.size() -1];
 }



  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // TODO Returns an approximate unit normal at this vertex, computed by
    // TODO taking the area-weighted average of the normals of neighboring
    // TODO triangles, then normalizing.
	  HalfedgeCIter h = halfedge();
	  HalfedgeCIter h_orig = h;
	  //area weighted norm n
	  Vector3D n = Vector3D(0,0,0);
	  do {
		  n += cross(h->vertex()->position, h->twin()->next()->twin()->vertex()->position);
		  h = h->next();
		  h = h->twin();
	  } while (h != h_orig);
    return n.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // TODO This method should flip the given edge and return an iterator to the flipped edge.
	  // Check that it is not a boundary edge
	  if (e0->halfedge()->face()->isBoundary() || e0->halfedge()->twin()->face()->isBoundary())
		  return e0;

	  //The current
	  HalfedgeIter bc = e0->halfedge();
	  HalfedgeIter ca = bc->next();
	  HalfedgeIter ab = ca->next();
	  HalfedgeIter cb = bc->twin();
	  HalfedgeIter bd = cb->next();
	  HalfedgeIter dc = bd->next();
	  //reverse - unnecessary since works in any direction
	  HalfedgeIter ac = ca->twin();
	  HalfedgeIter ba = ab->twin();
	  HalfedgeIter cd = dc->twin();
	  HalfedgeIter db = bd->twin();

	  FaceIter face_abc = bc->face(); //becomes the top one
	  FaceIter face_bcd = bd->face(); //becomes the bottom one
	  VertexIter a = ab->vertex();
	  VertexIter b = bc->vertex();
	  VertexIter c = cb->vertex();
	  VertexIter d = dc->vertex();
	  
	  //The flip
	  // Starting at ca
	  ca->next() = bc;
	  // Doing ad new edge
	  HalfedgeIter ad = ca->next();
	  ad->edge() = bc->edge();
	  ad->vertex() = a;
	  ad->face() = face_abc;
	  ad->next() = dc;
	  //missing twin. will be after;
	  // Doing dc
	  dc->next() = ca;
	  dc->face() = face_abc;
	  // Done with top triangle
	  
	  // Starting at bd
	  bd->next() = cb;
	  // Doing da new edge
	  HalfedgeIter da = bd->next();
	  da->edge() = cb->edge();
	  da->vertex() = d;
	  da->face() = face_bcd;
	  da->next() = ab;
	  da->twin() = ad;
	  ad->twin() = da;
	  // Doing ab
	  ab->next() = bd;
	  ab->face() = face_bcd;
	  // Done with bottom triangle

	  //// Top triangle in reverse starting with cd
	  //cd->next() = da;
	  ////already did da
	  //ac->next() = cd;

	  //// Bottom triangle in reverse starting with db
	  //db->next() = ba;
	  //ba->next() = ad;
	  ////already did ad

	  // Change vertices
	  b->halfedge() = bd;
	  c->halfedge() = ca;

	  // And faces
	  face_abc->halfedge() = ca;
	  face_bcd->halfedge() = bd;
	  
    return ad->edge();
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
    // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.
	  // Check that it is not a boundary edge
	  if (e0->halfedge()->face()->isBoundary() || e0->halfedge()->twin()->face()->isBoundary())
		  return e0->halfedge()->vertex();

	  // The current
	  HalfedgeIter bc = e0->halfedge();
	  HalfedgeIter ca = bc->next();
	  HalfedgeIter ab = ca->next();
	  HalfedgeIter cb = bc->twin();
	  HalfedgeIter bd = cb->next();
	  HalfedgeIter dc = bd->next();

	  FaceIter face_abc = bc->face(); //becomes the top left one
	  FaceIter face_bcd = bd->face(); //becomes the top right one
	  VertexIter a = ab->vertex();
	  VertexIter b = bc->vertex();
	  VertexIter c = cb->vertex();
	  VertexIter d = dc->vertex();

	  // Create new vertex
	  VertexIter m = newVertex();
	  m->position = (b->position + c->position)/2; //average

	  // Start upper left with mc
	  HalfedgeIter mc = newHalfedge();
	  mc->edge() = newEdge();
	  mc->edge()->halfedge() = mc;
	  mc->vertex() = m;
	  m->halfedge() = mc;
	  mc->face() = face_abc;
	  face_abc->halfedge() = ca;
	  //mc->isBoundary() = false;
	  mc->next() = ca;
	  //did not do twin
	  HalfedgeIter am = newHalfedge();
	  ca->next() = am;	  
	  am->edge() = newEdge();
	  am->edge()->halfedge() = am;
	  am->vertex() = a;
	  am->face() = face_abc;
	  am->next() = mc;
	  //again, no twin

	  // Next bottom left with ab
	  ab->face() = newFace();
	  FaceIter face_abm = ab->face();
	  face_abm->halfedge() = ab;
	  ab->next() = bc;
	  c->halfedge() = ca;
	  //bc is now the same thing as bm
	  HalfedgeIter bm = bc;
	  HalfedgeIter ma = newHalfedge();
	  bm->face() = face_abm;
	  bm->next() = ma;
	  bm->vertex() = b;
	  bm->edge() = bc->edge();
	  bc = bm;
	  bc->next() = bm->next();
	  ab->next() = bm;
	  ma->vertex() = m;
	  ma->edge() = am->edge();
	  ma->face() = face_abm;
	  ma->next() = ab;
	  ma->twin() = am;
	  am->twin() = ma;

	  // Next bottom right with bd
	  HalfedgeIter dm = newHalfedge();
	  bd->next() = dm;
	  bd->face() = newFace();
	  FaceIter face_bdm = bd->face();
	  face_bdm->halfedge() = bd;
	  dm->face() = face_bdm;
	  dm->vertex() = d;
	  dm->edge() = newEdge();
	  dm->edge()->halfedge() = dm;
	  dm->next() = cb;
	  cb->next() = bd;
	  HalfedgeIter mb = dm->next();
	  mb->twin() = bm;
	  bm->twin() = mb;
	  mb->next() = bd;
	  mb->face() = face_bdm;
	  mb->edge() = bm->edge();
	  cb = mb;
	  mb->vertex() = m;

	  // Last triangle - top right with dc
	  HalfedgeIter cm = newHalfedge();
	  cm->twin() = mc;
	  mc->twin() = cm;
	  dc->next() = cm;
	  cm->vertex() = c;
	  cm->edge() = mc->edge();
	  HalfedgeIter md = newHalfedge();
	  cm->face() = face_bcd;
	  face_bcd->halfedge() = dc;
	  cm->next() = md;
	  md->face() = face_bcd;
	  md->vertex() = m;
	  md->edge() = dm->edge();
	  md->next() = dc;
	  dm->twin() = md;
	  md->twin() = dm;

	  b->halfedge() = bd;
	  
    return m;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
    // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse. We will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.


    // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // TODO a vertex of the original mesh.
	  float n, u;
	  Vector3D original_position, neighbor_position_sum;
	  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
		  v->isNew = false;
		  original_position = v->position;
		  neighbor_position_sum = Vector3D(0,0,0);
		  n = 0;
		  
		  // need to look at the neighbors
		  HalfedgeCIter h = v->halfedge();
		  do {
			  neighbor_position_sum += h->twin()->vertex()->position;
			  n += 1; 
			  h = h->twin()->next();
		  } while (h != v->halfedge());

		  if (n == 3)
			  u = 3.0f/16.0f;
		  else
			  u = 3.0f/(8.0f * n);
		  
		  v->newPosition = (1.0f - n * u) * original_position + (u * neighbor_position_sum);
		  //std::cout << n << "\s";
		  //std::cout << u << "\n";
	  }



    // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
	  EdgeIter end = mesh.edgesEnd();
	  Vector3D A, B, C, D;
	  for (EdgeIter e = mesh.edgesBegin(); e != end; e++) {
		  A = e->halfedge()->vertex()->position;
		  B = e->halfedge()->twin()->vertex()->position;

		  C = e->halfedge()->next()->vertex()->position;
		  D = e->halfedge()->twin()->next()->vertex()->position;
		  e->newPosition = 3.0f/8.0 * (A + B) + 1.0f/8.0 * (C + D);
		  e->isNew = false;
		  //std::cout << e->newPosition << "\s";
	  }


    // TODO Next, we're going to split every edge in the mesh, in any order.  For future
    // TODO reference, we're also going to store some information about which subdivided
    // TODO edges come from splitting an edge in the original mesh, and which edges are new,
    // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
    // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
    // TODO just split (and the loop will never end!)

	  end = mesh.edgesEnd();
	  VertexIter new_v;
	  std::list<EdgeIter> edges;
	  // e;
	  for (EdgeIter e = mesh.edgesBegin(); e != end; e++) {
		  edges.push_back(e);
	  }

	  for (std::list<EdgeIter>::iterator e = edges.begin(); e != edges.end(); e++) {
		  new_v = mesh.splitEdge(*e);
		  new_v->newPosition = (*e)->newPosition;
		  new_v->isNew = true;
		  new_v->halfedge()->next()->next()->edge()->isNew = true;
		  new_v->halfedge()->twin()->next()->edge()->isNew = true;
	  }
	  //should only be splitting 3?


    // TODO Now flip any new edge that connects an old and new vertex.
	  end = mesh.edgesEnd();
	  bool ab, ba;
	  for (EdgeIter e = mesh.edgesBegin(); e != end; e++) {
		  if (e->isNew) {
			  ab = e->halfedge()->vertex()->isNew;
			  ba = e->halfedge()->twin()->vertex()->isNew;
			  if ((ab&&!ba) || (!ab&&ba))
				  mesh.flipEdge(e);
		  }
			
	  }
	  //should only be flipping 1?

    // TODO Finally, copy the new vertex positions into final Vertex::position.
	  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
		  v->position = v->newPosition;
	  }

    return;
  }
}
