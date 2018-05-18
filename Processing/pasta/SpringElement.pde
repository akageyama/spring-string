
class SpringElement
{
  //  A spring, with two ends called alpha and beta.
  //           alpha         beta
  //              O=========O

  int alpha;
  int beta;
  float springConst;

  SpringElement(float springConst, int alpha, int beta)
  {
    this.springConst = springConst;
    this.alpha = alpha;
    this.beta  = beta;
  }


  Vec3 getPullForce(int particleId, Vec3[] pos)
  {
    Vec3 force = new Vec3(0.0);

    float ax = pos[alpha].x;
    float ay = pos[alpha].y;
    float az = pos[alpha].z;
    float bx = pos[beta ].x;
    float by = pos[beta ].y;
    float bz = pos[beta ].z;

    float distance = dist(ax, ay, az, bx, by, bz);

    if ( distance >= EDGE_ELONGATION_CUT_LIMIT )
      springConst = 0.0;

    if ( distance <= EDGE_CONTRACTION_CUT_LIMIT )
      springConst = 0.0;

    float pullForceAmplitude = springConst * (distance - EDGE_LENGTH);

    if (particleId==alpha) {
      force = new Vec3((bx-ax)/distance,
                       (by-ay)/distance,
                       (bz-az)/distance);  // unit vec from alpha to beta.
      force.multiply(pullForceAmplitude);
    }
    else if (particleId==beta) {
      force = new Vec3((ax-bx)/distance,
                       (ay-by)/distance,
                       (az-bz)/distance);  // unit vec from beta to alpha.
      force.multiply(pullForceAmplitude);
    }
    else {
      println("??? SpringElement/getForce");
      exit();
    }
    return force;
  }


  void display(Vec3[] pos)
  {
    if ( springConst > 0.0) {
      float ax = mapx(pos[alpha].x);
      float ay = mapy(pos[alpha].y);
      float az = mapz(pos[alpha].z);
      float bx = mapx(pos[beta ].x);
      float by = mapy(pos[beta ].y);
      float bz = mapz(pos[beta ].z);
      line(ax,ay,az,bx,by,bz);
    }
  }


  float energy(Vec3[] pos)
  {
    float ax = pos[alpha].x;
    float ay = pos[alpha].y;
    float az = pos[alpha].z;
    float bx = pos[beta ].x;
    float by = pos[beta ].y;
    float bz = pos[beta ].z;

    float distance = dist(ax, ay, az, bx, by, bz);
    float elongation = distance - EDGE_LENGTH;

    return 0.5*springConst*pow(elongation,2);
  }

}
