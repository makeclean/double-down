
// MOAB
#include "moab/EntityType.hpp"
#include "moab/CartVect.hpp"
#include "moab/GeomUtil.hpp"

// Local
#include "double-down/primitives.hpp"

void intersectionFilter(void* ptr, RTCDRayHit &rayhit)
{
  switch(rayhit.ray.rf_type)
    {
    case RayFireType::RF: //if this is a typical ray_fire, check the dot_product
      if (0 > rayhit.dot_prod()) { rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; }
      break;
    case RayFireType::PIV: //if this is a point_in_vol fire, accept all hits
      break;
    }
}

void DblTriBounds(const RTCBoundsFunctionArguments* args)
{
  // get the array of DblTri's stored on the geometry
  void* tris_i = args->geometryUserPtr;
  // referencde to the returned triangle bounds
  RTCBounds& bounds_o = *args->bounds_o;

  // convert the void user data pointer to the DblTri pointer
  const DblTri* tris = (const DblTri*) tris_i;
  // select the DblTri in question based on the primitive ID
  const DblTri& this_tri = tris[args->primID];

  // get the MOAB direct access manager from the triangle (it's present on each one)
  MBDirectAccess* mdam = (MBDirectAccess*) this_tri.mdam;

  // compute the bounding box using the direct access manager and the triangle's handle
  bounds_o = DblTriBounds(mdam, this_tri.handle);
}

void DblTriIntersectFunc(RTCIntersectFunctionNArguments* args) {
  // get the array of DblTri's stored on the geometry
  void* tris_i = args->geometryUserPtr;
  // convert the void user data pointer to the DblTri pointer
  const DblTri* tris = (const DblTri*) tris_i;
  // select the DblTri in question based on the primitive ID
  const DblTri& this_tri = tris[args->primID];

  // get the rayhit object from the args
  RTCDRayHit* rayhit = (RTCDRayHit*)args->rayhit;
  RTCDRay& ray = rayhit->ray;
  RTCDHit& hit = rayhit->hit;

  // get the direct access manager from the triangle (it's present on each one)
  MBDirectAccess* mdam = (MBDirectAccess*) this_tri.mdam;

  // get the triangle coordinates from the direct access manager
  std::array<moab::CartVect, 3> coords = mdam->get_mb_coords(this_tri.handle);


  double dist; // local variable for the distance to the triangle intersection
  // convert the ray origin and directory to MOAB CartVect
  moab::CartVect ray_org(ray.dorg.x, ray.dorg.y, ray.dorg.z);
  moab::CartVect ray_dir(ray.ddir.x, ray.ddir.y, ray.ddir.z);

  double nonneg_ray_len = 1e17;
  const double* ptr = &nonneg_ray_len;
  bool hit_tri = moab::GeomUtil::plucker_ray_tri_intersect(&(coords[0]), ray_org, ray_dir, dist, ptr);

  if ( hit_tri ) {
    // only accept hits that are closer than the current one
    if (dist > ray.dtfar) {
      hit.geomID = -1;
      return;
    }
    // set both the single and double ray distances
    ray.set_len(dist);
    // zero out the barycentric coordinates of the hit
    hit.u = 0.0f;
    hit.v = 0.0f;
    // set the hit
    hit.geomID = this_tri.geomID;
    hit.primID = args->primID;

    // compute the triangle normal
    moab::CartVect normal = (coords[1] - coords[0]) * (coords[2] - coords[0]);
    normal.normalize();

    // flip the triangle normal if the sense is reversed
    if( -1 == this_tri.sense ) normal *= -1;

    // set the double precision normal of the hit
    hit.dNg[0] = normal[0];
    hit.dNg[1] = normal[1];
    hit.dNg[2] = normal[2];
  } else {
    // otherwise set the hit's geomID as invalid to ensure we don't consider this triangle as hit
    hit.geomID = RTC_INVALID_GEOMETRY_ID;
  }
}

void DblTriOccludedFunc(RTCOccludedFunctionNArguments* args) {
  // get the array of DblTri's stored on the geometry
  void* tris_i = args->geometryUserPtr;
  // convert the void user data pointer to the DblTri pointer
  const DblTri* tris = (const DblTri*) tris_i;
  // select the DblTri in question based on the primitive ID
  const DblTri& this_tri = tris[args->primID];

  // get the double precision ray from the args
  RTCDRay* ray = (RTCDRay*)&(args->ray);

  MBDirectAccess* mdam = (MBDirectAccess*) mdam;

  std::array<moab::CartVect, 3> coords = mdam->get_mb_coords(this_tri.handle);

  double dist;
  double nonneg_ray_len = 1e37;
  double* ptr = &nonneg_ray_len;
  Vec3da& dorg = ray->dorg;
  moab::CartVect ray_org(dorg[0], dorg[1], dorg[2]);
  Vec3da& ddir = ray->ddir;
  moab::CartVect ray_dir(ddir[0], ddir[1], ddir[2]);

  bool hit_tri = moab::GeomUtil::plucker_ray_tri_intersect(&(coords[0]), ray_org, ray_dir, dist, ptr);
  // if we hit the triangle, we hit the triangle
  if ( hit_tri ) {
    ray->set_len(neg_inf);
  }
}

bool DblTriPointQueryFunc(RTCPointQueryFunctionArguments* args) {
  RTCDPointQuery* pq = (RTCDPointQuery*) args->query;

  double pnt[3];
  pnt[0] = pq->dx;
  pnt[1] = pq->dy;
  pnt[2] = pq->dz;

  RTCGeometry g = rtcGetGeometry(*(RTCScene*)args->userPtr, args->geomID);

  // get the array of DblTri's stored on the geometry
  void* tris_i = rtcGetGeometryUserData(g);
  // convert the void user data pointer to the DblTri pointer
  const DblTri* tris = (const DblTri*) tris_i;
  // select the DblTri in question based on the primitive ID
  const DblTri& this_tri = tris[args->primID];
  // compute the distance to the nearest point on the triangle
  double dist = DblTriClosestFunc(this_tri, pnt);
  // update the point query values if this distance is less than the current one
  if (dist < pq->dradius) {
    pq->radius = dist;
    pq->dradius = dist;
    pq->primID = args->primID;
    pq->geomID = args->geomID;
    return true;
  } else {
    return false;
  }
}

double DblTriClosestFunc(const DblTri& tri, const double loc[3]) {
  // get the direct access manager
  MBDirectAccess* mdam = (MBDirectAccess*) tri.mdam;
  // get the triangle coordinates
  std::array<moab::CartVect, 3> coords = mdam->get_mb_coords(tri.handle);
  // convert query location to a MOAB CartVect
  moab::CartVect location(loc[0], loc[1], loc[2]);
  moab::CartVect result;
  // compute the nearest location on the triangle
  moab::GeomUtil::closest_location_on_tri(location, coords.data(), result);
  // return the distance to the triangle location
  return (result - location).length();
}
