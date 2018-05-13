
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


  Vec3 getPullForce(int particleId,
                    float[] posx,
                    float[] posy,
                    float[] posz)
  {
    Vec3 force = new Vec3();

    float ax = posx[alpha];
    float ay = posy[alpha];
    float az = posz[alpha];
    float bx = posx[beta ];
    float by = posy[beta ];
    float bz = posz[beta ];

    float distance = dist(ax, ay, az, bx, by, bz);

    if ( distance >= SPRING_CUT_LIMIT_LENGTH ) springConst = 0.0; // cut off.

    float pullForceAmplitude = springConst * (distance - EDGE_LENGTH);

//  if (nonlinearSpringFlag) {
//    float factor = springConst*EDGE_LENGTH;
//    pullForceAmplitude += factor * pow(distance/EDGE_LENGTH - 1, 3);
//  }

    if (particleId==alpha) {
      force = new Vec3((bx-ax)/distance,
                       (by-ay)/distance,
                       (bz-az)/distance);  // unit vector from alpha to beta.
      force.multiply(pullForceAmplitude);
    }
    else if (particleId==beta) {
      force = new Vec3((ax-bx)/distance,
                       (ay-by)/distance,
                       (az-bz)/distance);  // unit vector from beta to alpha.
      force.multiply(pullForceAmplitude);
    }
    else {
      println("??? SpringElement/getForce");
      exit();
    }
    return force;
  }


  void display(float[] posx,
               float[] posy,
               float[] posz)
  {

    if ( springConst > 0.0 ) {
      float ax = mapx(posx[alpha]);
      float ay = mapy(posy[alpha]);
      float az = mapz(posz[alpha]);
      float bx = mapx(posx[beta]);
      float by = mapy(posy[beta]);
      float bz = mapz(posz[beta]);
      line(ax,ay,az,bx,by,bz);
    }
  }


  float energy(float[] posx,
               float[] posy,
               float[] posz)
  {
    float ax = posx[alpha];
    float ay = posy[alpha];
    float az = posz[alpha];
    float bx = posx[beta ];
    float by = posy[beta ];
    float bz = posz[beta ];

    float distance = dist(ax, ay, az, bx, by, bz);
    float elongation = distance - EDGE_LENGTH;

    return 0.5*springConst*pow(elongation,2);
  }

}
