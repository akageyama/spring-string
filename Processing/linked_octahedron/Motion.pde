
class Motion
{

  void rungeKuttaIncrement(int num,
                           float[] posx,
                           float[] posy,
                           float[] posz,
                           float[] velx,
                           float[] vely,
                           float[] velz,
                           float[] posx1,
                           float[] posy1,
                           float[] posz1,
                           float[] velx1,
                           float[] vely1,
                           float[] velz1,
                           float[] dposx,
                           float[] dposy,
                           float[] dposz,
                           float[] dvelx,
                           float[] dvely,
                           float[] dvelz,
                           float factor)
  {
    for (int p=0; p<num; p++) {
      posx[p] = posx1[p] + factor*dposx[p];
      posy[p] = posy1[p] + factor*dposy[p];
      posz[p] = posz1[p] + factor*dposz[p];
      velx[p] = velx1[p] + factor*dvelx[p];
      vely[p] = vely1[p] + factor*dvely[p];
      velz[p] = velz1[p] + factor*dvelz[p];
    }
  }

  void zeroset(int n, Vec3[] v)
  {
    for (int i=0; i<n; i++) {
      v[i] = new Vec3(0.0, 0.0, 0.0); 
    }
  }

  void getForceSpring(float[] posx,
                      float[] posy,
                      float[] posz,
                      Vec3[] force)
  {
    for (int l=0; l<N_TRIANGLES; l++) {
      for (int j=0; j<3; j++) {
        int p = particles.id(l,j);
        int[] splist = new int[6];
        splist = particles.getConnectedSpingListForThisParticle(p);

        Vec3 forceSum = new Vec3(0.0, 0.0, 0.0);

        int numCounterpart
            = particles.numberOfConnectedSpringsToThisParticle(p);
        for (int s=0; s<numCounterpart; s++) {
          SpringElement aSpring = springs.element[splist[s]];
          Vec3 pullForceFromTheSpring = aSpring.getPullForce(p,posx,
                                                               posy,
                                                               posz);
          forceSum.add(pullForceFromTheSpring);
        }
        force[p] = forceSum;
      }
    }
  }

  void getForceFriction(float[] velx,
                        float[] vely,
                        float[] velz,
                        Vec3[] force)
  {
    for (int p=0; p<N_PARTICLES; p++) {
      force[p].x = -FRICTION_COEFF*velx[p];
      force[p].y = -FRICTION_COEFF*vely[p];
      force[p].z = -FRICTION_COEFF*velz[p];
    }
  }


  void getForceGravity(float[] posx,
                       float[] posy,
                       float[] posz,
                       Vec3[] force)
  {
    for (int p=0; p<N_PARTICLES; p++) {
     force[p].x = 0.0;
     force[p].y = -GRAVITY_ACCELERATION*PARTICLE_MASS;
     force[p].z = 0.0;
    }
  }

  void getForceTension(float[] posx,
                       float[] posy,
                       float[] posz,
                       Vec3[] force)
  {
    for (int k=0; k<2; k++) {
      int l = (N_TRIANGLES-1)*k;
      for (int j=0; j<3; j++) {
        int p = particles.id(l,j);
        if (k==0) {
          force[p].x = -0.01*SPRING_CONST*(posx[p] + ROPE_LENGTH/2);
        }
        else {
          force[p].x = -0.01*SPRING_CONST*(posx[p] - ROPE_LENGTH/2);
        }
      }
    }
  }


  void equationOfMotion(float posx[],
                        float posy[],
                        float posz[],
                        float velx[],
                        float vely[],
                        float velz[],
                        float dposx[],
                        float dposy[],
                        float dposz[],
                        float dvelx[],
                        float dvely[],
                        float dvelz[],
                        float dt)
  {
    float dtm = dt / PARTICLE_MASS;

    Vec3[] forceSpring   = new Vec3[N_PARTICLES];
    Vec3[] forceFriction = new Vec3[N_PARTICLES];
    Vec3[] forceTension  = new Vec3[N_PARTICLES];
    Vec3[] forceGravity  = new Vec3[N_PARTICLES];

    zeroset(N_PARTICLES,forceSpring);
    zeroset(N_PARTICLES,forceTension);
    zeroset(N_PARTICLES,forceFriction);
    zeroset(N_PARTICLES,forceGravity);

    getForceSpring(posx, posy, posz, forceSpring);
    if ( frictionFlag ) {
      getForceFriction(velx, vely, velz, forceFriction);
    }
    getForceTension(posx, posy, posz, forceTension);
    getForceGravity(posx, posy, posz, forceGravity);

    for (int l=1; l<N_TRIANGLES-1; l++) {
      for (int j=0; j<3; j++) {
        int p = particles.id(l,j); // particle id
        Vec3 force = new Vec3(0.0, 0.0, 0.0);

        force.add(forceSpring[p]);
        force.add(forceFriction[p]);
        force.add(forceTension[p]);
        force.add(forceGravity[p]);

        dposx[p] = velx[p] * dt;
        dposy[p] = vely[p] * dt;
        dposz[p] = velz[p] * dt;
        dvelx[p] = force.x * dtm;
        dvely[p] = force.y * dtm;
        dvelz[p] = force.z * dtm;
      }
    }

  }


  void boundaryCondition(float angle,
                         float[] posx,
                         float[] posy,
                         float[] posz,
                         float[] velx,
                         float[] vely,
                         float[] velz)
  {
    Vec3[] verts = new Vec3[3];

    for (int l=0; l<2; l++) {
      int tl = l*(N_TRIANGLES-1);
      for (int j=0; j<3; j++) {
        int p = particles.id(tl,j);
        verts[j] = new Vec3(posx[p], posy[p], posz[p]);
      }

      if (l==0)
        twist(angle/2, verts);
      else
        twist(-angle/2, verts);

      for (int j=0; j<3; j++) {
        int p = particles.id(tl,j);
        Vec3 vert = verts[j];
        posx[p] = vert.x;
        posy[p] = vert.y;
        posz[p] = vert.z;
      }
    }
  }


  void rungeKutta()
  {
    final float ONE_SIXTH = 1.0/6.0;
    final float ONE_THIRD = 1.0/3.0;
    final int NN = N_PARTICLES;

    float[] posxprev = new float[NN];
    float[] posxwork = new float[NN];
    float[]   dposx1 = new float[NN];
    float[]   dposx2 = new float[NN];
    float[]   dposx3 = new float[NN];
    float[]   dposx4 = new float[NN];
    float[] posyprev = new float[NN];
    float[] posywork = new float[NN];
    float[]   dposy1 = new float[NN];
    float[]   dposy2 = new float[NN];
    float[]   dposy3 = new float[NN];
    float[]   dposy4 = new float[NN];
    float[] poszprev = new float[NN];
    float[] poszwork = new float[NN];
    float[]   dposz1 = new float[NN];
    float[]   dposz2 = new float[NN];
    float[]   dposz3 = new float[NN];
    float[]   dposz4 = new float[NN];
    float[] velxprev = new float[NN];
    float[] velxwork = new float[NN];
    float[]   dvelx1 = new float[NN];
    float[]   dvelx2 = new float[NN];
    float[]   dvelx3 = new float[NN];
    float[]   dvelx4 = new float[NN];
    float[] velyprev = new float[NN];
    float[] velywork = new float[NN];
    float[]   dvely1 = new float[NN];
    float[]   dvely2 = new float[NN];
    float[]   dvely3 = new float[NN];
    float[]   dvely4 = new float[NN];
    float[] velzprev = new float[NN];
    float[] velzwork = new float[NN];
    float[]   dvelz1 = new float[NN];
    float[]   dvelz2 = new float[NN];
    float[]   dvelz3 = new float[NN];
    float[]   dvelz4 = new float[NN];

    arrayCopy(particles.posx, posxprev);
    arrayCopy(particles.posy, posyprev);
    arrayCopy(particles.posz, poszprev);
    arrayCopy(particles.velx, velxprev);
    arrayCopy(particles.vely, velyprev);
    arrayCopy(particles.velz, velzprev);

    //step 1
    equationOfMotion(posxprev,
                     posyprev,
                     poszprev,
                     velxprev,
                     velyprev,
                     velzprev,
                     dposx1,
                     dposy1,
                     dposz1,
                     dvelx1,
                     dvely1,
                     dvelz1,
                     dt);
    rungeKuttaIncrement(NN,
                        posxwork,
                        posywork,
                        poszwork,
                        velxwork,
                        velywork,
                        velzwork,
                        posxprev,
                        posyprev,
                        poszprev,
                        velxprev,
                        velyprev,
                        velzprev,
                        dposx1,
                        dposy1,
                        dposz1,
                        dvelx1,
                        dvely1,
                        dvelz1,
                        0.5);

    //step 2
    time += 0.5*dt;
    float angle = EDGE_TWIST_RATE_OMEGA * dt * 0.5;
    boundaryCondition(angle,
                      posxwork,
                      posywork,
                      poszwork,
                      velxwork,
                      velywork,
                      velzwork);
    equationOfMotion(posxwork,
                     posywork,
                     poszwork,
                     velxwork,
                     velywork,
                     velzwork,
                     dposx2,
                     dposy2,
                     dposz2,
                     dvelx2,
                     dvely2,
                     dvelz2,
                     dt);
    rungeKuttaIncrement(NN,
                        posxwork,
                        posywork,
                        poszwork,
                        velxwork,
                        velywork,
                        velzwork,
                        posxprev,
                        posyprev,
                        poszprev,
                        velxprev,
                        velyprev,
                        velzprev,
                        dposx2,
                        dposy2,
                        dposz2,
                        dvelx2,
                        dvely2,
                        dvelz2,
                        0.5);

    //step 3
    boundaryCondition(angle,
                      posxwork,
                      posywork,
                      poszwork,
                      velxwork,
                      velywork,
                      velzwork);
    equationOfMotion(posxwork,
                     posywork,
                     poszwork,
                     velxwork,
                     velywork,
                     velzwork,
                     dposx3,
                     dposy3,
                     dposz3,
                     dvelx3,
                     dvely3,
                     dvelz3,
                     dt);
    rungeKuttaIncrement(NN,
                        posxwork,
                        posywork,
                        poszwork,
                        velxwork,
                        velywork,
                        velzwork,
                        posxprev,
                        posyprev,
                        poszprev,
                        velxprev,
                        velyprev,
                        velzprev,
                        dposx3,
                        dposy3,
                        dposz3,
                        dvelx3,
                        dvely3,
                        dvelz3,
                        1.0);

    //step 4
    time += 0.5*dt;
    angle = EDGE_TWIST_RATE_OMEGA * dt * 0.5;
    boundaryCondition(angle,
                      posxwork,
                      posywork,
                      poszwork,
                      velxwork,
                      velywork,
                      velzwork);
    equationOfMotion(posxwork,
                     posywork,
                     poszwork,
                     velxwork,
                     velywork,
                     velzwork,
                     dposx4,
                     dposy4,
                     dposz4,
                     dvelx4,
                     dvely4,
                     dvelz4,
                     dt);
    // weighted sum
    for (int tl=0; tl<N_TRIANGLES; tl++) {
      for (int j=0; j<3; j++) { // three verteces in a triangle.
        int p = particles.id(tl,j);
        posxwork[p] =            posxprev[p] + (
                         ONE_SIXTH*dposx1[p]
                       + ONE_THIRD*dposx2[p]
                       + ONE_THIRD*dposx3[p]
                       + ONE_SIXTH*dposx4[p]
                       );
        posywork[p] =            posyprev[p] + (
                         ONE_SIXTH*dposy1[p]
                       + ONE_THIRD*dposy2[p]
                       + ONE_THIRD*dposy3[p]
                       + ONE_SIXTH*dposy4[p]
                       );
        poszwork[p] =            poszprev[p] + (
                         ONE_SIXTH*dposz1[p]
                       + ONE_THIRD*dposz2[p]
                       + ONE_THIRD*dposz3[p]
                       + ONE_SIXTH*dposz4[p]
                       );
        velxwork[p] =            velxprev[p] + (
                         ONE_SIXTH*dvelx1[p]
                       + ONE_THIRD*dvelx2[p]
                       + ONE_THIRD*dvelx3[p]
                       + ONE_SIXTH*dvelx4[p]
                       );
        velywork[p] =            velyprev[p] + (
                         ONE_SIXTH*dvely1[p]
                       + ONE_THIRD*dvely2[p]
                       + ONE_THIRD*dvely3[p]
                       + ONE_SIXTH*dvely4[p]
                       );
        velzwork[p] =            velzprev[p] + (
                         ONE_SIXTH*dvelz1[p]
                       + ONE_THIRD*dvelz2[p]
                       + ONE_THIRD*dvelz3[p]
                       + ONE_SIXTH*dvelz4[p]
                       );
      }
    }
    boundaryCondition(angle,
                      posxwork,
                      posywork,
                      poszwork,
                      velxwork,
                      velywork,
                      velzwork);

    arrayCopy(posxwork, particles.posx);
    arrayCopy(posywork, particles.posy);
    arrayCopy(poszwork, particles.posz);
    arrayCopy(velxwork, particles.velx);
    arrayCopy(velywork, particles.vely);
    arrayCopy(velzwork, particles.velz);
  }


  void display()
  {
    particles.display();
    springs.display(particles.posx,
                    particles.posy,
                    particles.posz);
  }


  Vec3 calcCenter(Vec3[] verts)
  {
    Vec3 center = new Vec3(0.0, 0.0, 0.0);
    for (int j=0; j<3; j++) {
      center.add(verts[j]);
    }
    center.divide(3.0);

    return center;
  }

  Vec3 vecSubtract(Vec3 a, Vec3 b) // a-b
  {
    Vec3 ans = new Vec3(0.0, 0.0, 0.0);

    ans.x = a.x - b.x;
    ans.y = a.y - b.y;
    ans.z = a.z - b.z;

    return ans;
  }


  void calcUnitVectors(Vec3 center,
                       Vec3[] verts,
                       Vec3[] unitVec)
  {
    Vec3 vecC0 = vecSubtract(verts[0], center);
    Vec3 vecC1 = vecSubtract(verts[1], center);

    Vec3 v = vecC0;
    v.normalize();
    unitVec[0] = new Vec3(v);

    v = vecC0.crossProduct(vecC1); // along the axis
    v.normalize();

    unitVec[1] = v.crossProduct(unitVec[0]);
  }


  void twist(float angle, Vec3[] verts)
  {
    Vec3 center;
    Vec3 vecC0, vecC1;
    Vec3[] unitVec = new Vec3[2];

    center = calcCenter(verts);
    calcUnitVectors(center,
                    verts,
                    unitVec);

    float deltaPhi = (PI*2) / 3;
    float r = TUBE_RADIUS;

    for (int j=0; j<3; j++) {
      float phi = angle + j*deltaPhi;
      float x, y, z;
      x = center.x + r*unitVec[0].x*cos(phi)
                   + r*unitVec[1].x*sin(phi);
      y = center.y + r*unitVec[0].y*cos(phi)
                   + r*unitVec[1].y*sin(phi);
      z = center.z + r*unitVec[0].z*cos(phi)
                   + r*unitVec[1].z*sin(phi);
      verts[j] = new Vec3(x,y,z);
    }
  }


  float totalEnergy()
  {
    float kinetic = particles.energy();
    float potentialSpring  = springs.energy(particles.posx,
                                            particles.posy,
                                            particles.posz);

    float ysum=0.0;
    for (int p=0; p<N_PARTICLES; p++)
      ysum += particles.posy[p];

    float potentialGravity = PARTICLE_MASS*GRAVITY_ACCELERATION*ysum;
    return kinetic + potentialSpring + potentialGravity;
  }

}
