
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

    for (int l=1; l<N_PAIRS; l++) { // ladder
      for (int j=0; j<2; j++) {
        int p = particles.id(l,j); // particle id
        int[] splist = new int[3];
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

        // gravity force
        Vec3 gForce = new Vec3(0.0,
                               0.0,
                               -GRAVITY_ACCELERATION*PARTICLE_MASS);
        forceSum.add(gForce);

        // viscous force
        if ( frictionFlag ) {
          float vForceX = -FRICTION_COEFF*velx[p];
          float vForceY = -FRICTION_COEFF*vely[p];
          float vForceZ = -FRICTION_COEFF*velz[p];
          forceSum.add(vForceX, vForceY, vForceZ);
        }

        dposx[p] = velx[p] * dt;
        dposy[p] = vely[p] * dt;
        dposz[p] = velz[p] * dt;
        dvelx[p] = forceSum.x * dtm;
        dvely[p] = forceSum.y * dtm;
        dvelz[p] = forceSum.z * dtm;
      }
    }
  }


  void boundaryCondition(float t,
                         float[] posx,
                         float[] posy,
                         float[] posz)
  {
    Vec3[] verts = new Vec3[3];

    particles.upperBoundaryConfiguration(t, verts);

    for (int j=0; j<2; j++) { // three vertices at the bottom.
      int p = particles.id(0,j);
      posx[p] = verts[j].x;
      posy[p] = verts[j].y;
      posz[p] = verts[j].z;
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
    boundaryCondition(time,
                      posxwork,
                      posywork,
                      poszwork);
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
    boundaryCondition(time,
                      posxwork,
                      posywork,
                      poszwork);
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
    boundaryCondition(time,
                      posxwork,
                      posywork,
                      poszwork);
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
    for (int l=1; l<N_PAIRS-1; l++) {
      for (int j=0; j<2; j++) { // three verteces in a triangle.
        int p = particles.id(l,j);
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
    boundaryCondition(time,
                      posxwork,
                      posywork,
                      poszwork);

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
