
class Motion
{

  Vec3 calcCenter(int size, Vec3[] verts)
  {
    float sumx = 0.0;
    float sumy = 0.0;
    float sumz = 0.0;
    for (int j=0; j<size; j++) {
      sumx += verts[j].x;
      sumy += verts[j].y;
      sumz += verts[j].z;
    }

    Vec3 ans = new Vec3(sumx/size,
                        sumy/size,
                        sumz/size);
    return ans;
  }


  Vec3[] calcRelativeVecFromCenter(Vec3[] triangleVerts)
  {
    Vec3[] ans = new Vec3[3];

    Vec3 centerOfEndTriangle = new Vec3(calcCenter(3, triangleVerts));
    for (int j=0; j<3; j++) {
      ans[j] = triangleVerts[j].ssubtract(centerOfEndTriangle);
    }
    return ans;
  }

  void boundaryCondition(float tim, Vec3[] pos)
  {
    Vec3[] triangleVerts = new Vec3[3];
    Vec3[] relativeVecFromCenter = new Vec3[3];

    float endPointPositionX;

    for (int l=0; l<2; l++) {
      if (l==0)
        endPointPositionX = - STICK_LENGTH/2
                            + STICK_END_POINT_MOVE_SPEED*tim;
      else
        endPointPositionX =   STICK_LENGTH/2
                            - STICK_END_POINT_MOVE_SPEED*tim;
      int t = (N_TRIANGLES-1)*l;
      for (int j=0; j<3; j++) {
        int p = particles.id(t,j);
        triangleVerts[j] = pos[p];
      }
      relativeVecFromCenter = calcRelativeVecFromCenter(triangleVerts);
      for (int j=0; j<3; j++) {
        int p = particles.id(t,l);
        pos[p].x = endPointPositionX + relativeVecFromCenter[j].x;
      }
    }
  }

  void rungeKuttaIncrement(int num,
                           Vec3[] pos,
                           Vec3[] vel,
                           Vec3[] pos1,
                           Vec3[] vel1,
                           Vec3[] dpos,
                           Vec3[] dvel,
                           float factor)
  {
    for (int p=0; p<num; p++) {
      pos[p].x = pos1[p].x + factor*dpos[p].x;
      pos[p].y = pos1[p].y + factor*dpos[p].y;
      pos[p].z = pos1[p].z + factor*dpos[p].z;
      vel[p].x = vel1[p].x + factor*dvel[p].x;
      vel[p].y = vel1[p].y + factor*dvel[p].y;
      vel[p].z = vel1[p].z + factor*dvel[p].z;
    }
  }


  void equationOfMotion(Vec3[] pos,
                        Vec3[] vel,
                        Vec3[] dpos,
                        Vec3[] dvel,
                        float dt)
  {
    float dtm = dt / PARTICLE_MASS;

      //  u2=2
      //   o
      //        .     upper triangle
      //     .      .
      //               o
      //       .     .
      //           .
      //         o u1=0          vertex index
      //                           in a triangle

      //            o m2=2           o 2
      //         .                  . .
      //      .      .             .   .
      //    o                     .     .
      //  me=0  .      .         o . . . o
      //            .           0         1
      //                o m1=1
      //  l2=2
      //   o
      //        .
      //     .      .
      //               o
      //       .     .
      //           .    lower triangle
      //         o
      //        l1=0


    for (int t=1; t<N_TRIANGLES-1; t++) { // skip end triangles.
      for (int j=0; j<3; j++) {
        int p = particles.id(t,j);
        int[] splist = new int[6];
        splist = particles.getConnectedSpingListForThisParticle(p);

        Vec3 forceSum = new Vec3(0.0, 0.0, 0.0);
        int numCounterpart
            = particles.numberOfConnectedSpringsToThisParticle(p);
        for (int s=0; s<numCounterpart; s++) {
          SpringElement aSpring = springs.element[splist[s]];
          Vec3 pullForceFromTheSpring = aSpring.getPullForce(p,pos);
          forceSum.add(pullForceFromTheSpring);
        }

        // friction force
        if ( frictionFlag ) {
          Vec3 frictionForce = vel[p].mmultiply(-FRICTION_COEFF);
          forceSum.add(frictionForce);
        }

        dpos[p].x = (vel[p].x) * dt;
        dpos[p].y = (vel[p].y) * dt;
        dpos[p].z = (vel[p].z) * dt;

        dvel[p].x = forceSum.x * dtm;
        dvel[p].y = forceSum.y * dtm;
        dvel[p].z = forceSum.z * dtm;
      }
    }
  }

  void copyIt(Vec3[] from, Vec3[] to)
  {
    for (int p=0; p<N_PARTICLES; p++) {
      to[p].x = from[p].x;
      to[p].y = from[p].y;
      to[p].z = from[p].z;
    }
  }

  void rungeKutta()
  {
    final float ONE_SIXTH = 1.0/6.0;
    final float ONE_THIRD = 1.0/3.0;
    final int NN = N_PARTICLES;

    Vec3[] posprev = new Vec3[NN];
    Vec3[] poswork = new Vec3[NN];
    Vec3[]   dpos1 = new Vec3[NN];
    Vec3[]   dpos2 = new Vec3[NN];
    Vec3[]   dpos3 = new Vec3[NN];
    Vec3[]   dpos4 = new Vec3[NN];
    Vec3[] velprev = new Vec3[NN];
    Vec3[] velwork = new Vec3[NN];
    Vec3[]   dvel1 = new Vec3[NN];
    Vec3[]   dvel2 = new Vec3[NN];
    Vec3[]   dvel3 = new Vec3[NN];
    Vec3[]   dvel4 = new Vec3[NN];

    copyIt(particles.pos, posprev);
    copyIt(particles.vel, velprev);

    //step 1
    equationOfMotion(posprev,
                     velprev,
                     dpos1,
                     dvel1,
                     dt);
    rungeKuttaIncrement(NN,
                        poswork,
                        velwork,
                        posprev,
                        velprev,
                        dpos1,
                        dvel1,
                        0.5);

    //step 2
    time += 0.5*dt;
    boundaryCondition(time,
                      poswork);
    equationOfMotion(poswork,
                     velwork,
                     dpos2,
                     dvel2,
                     dt);
    rungeKuttaIncrement(NN,
                        poswork,
                        velwork,
                        posprev,
                        velprev,
                        dpos2,
                        dvel2,
                        0.5);

    //step 3
    boundaryCondition(time,
                      poswork);
    equationOfMotion(poswork,
                     velwork,
                     dpos3,
                     dvel3,
                     dt);
    rungeKuttaIncrement(NN,
                        poswork,
                        velwork,
                        posprev,
                        velprev,
                        dpos3,
                        dvel3,
                        1.0);

    //step 4
    time += 0.5*dt;
    boundaryCondition(time,
                      poswork);
    equationOfMotion(poswork,
                     velwork,
                     dpos4,
                     dvel4,
                     dt);
    // weighted sum
    for (int t=1; t<N_TRIANGLES-1; t++) {
      for (int j=0; j<3; j++) { // three verteces in a triangle.
        int p = particles.id(t,j);
        poswork[p].x =           posprev[p].x + (
                         ONE_SIXTH*dpos1[p].x
                       + ONE_THIRD*dpos2[p].x
                       + ONE_THIRD*dpos3[p].x
                       + ONE_SIXTH*dpos4[p].x
                       );
        poswork[p].y =           posprev[p].y + (
                         ONE_SIXTH*dpos1[p].y
                       + ONE_THIRD*dpos2[p].y
                       + ONE_THIRD*dpos3[p].y
                       + ONE_SIXTH*dpos4[p].y
                       );
        poswork[p].z =           posprev[p].z + (
                         ONE_SIXTH*dpos1[p].z
                       + ONE_THIRD*dpos2[p].z
                       + ONE_THIRD*dpos3[p].z
                       + ONE_SIXTH*dpos4[p].z
                       );
        velwork[p].x =           velprev[p].x + (
                         ONE_SIXTH*dvel1[p].x
                       + ONE_THIRD*dvel2[p].x
                       + ONE_THIRD*dvel3[p].x
                       + ONE_SIXTH*dvel4[p].x
                       );
        velwork[p].y =           velprev[p].y + (
                         ONE_SIXTH*dvel1[p].y
                       + ONE_THIRD*dvel2[p].y
                       + ONE_THIRD*dvel3[p].y
                       + ONE_SIXTH*dvel4[p].y
                       );
        velwork[p].z =           velprev[p].z + (
                         ONE_SIXTH*dvel1[p].z
                       + ONE_THIRD*dvel2[p].z
                       + ONE_THIRD*dvel3[p].z
                       + ONE_SIXTH*dvel4[p].z
                       );
      }
    }
    boundaryCondition(time,
                      poswork);

    copyIt(poswork, particles.pos);
    copyIt(velwork, particles.vel);
  }


  void display()
  {
    particles.display();
    springs.display(particles.pos);
  }



  float totalEnergy()
  {
    float kinetic = particles.energy();
    float potential = springs.energy(particles.pos);
    return kinetic + potential;
  }

}
