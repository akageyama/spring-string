
class Particles
{
  float[] posx = new float[N_PARTICLES];
  float[] posy = new float[N_PARTICLES];
  float[] posz = new float[N_PARTICLES];
  float[] velx = new float[N_PARTICLES];
  float[] vely = new float[N_PARTICLES];
  float[] velz = new float[N_PARTICLES];

  int[][] sixSpringList = new int[N_PARTICLES][6];
            // each particle are connected with 6 springs.
            //
            //    upper layer
            //  (3) (2)
            //   |  /    (1)
            //   | /   o    same
            //   |/ o          layer
            //  (X) o  o  o  (0)
            //    \ \
            //     \  \
            //      \   \
            //     (4)   (5)
            //       lower layer


  final int NULL_MARK = -1; // used to count already set elements.

        //
        //                     o  p=8
        //                  .    .
        //               .         .
        //            .              .
        //         o .   .  .  .  .  . o
        //      p=6                     p=7
        //
        //     p=5                      p=4
        //        o x x x x x x x x x o
        //          x              x
        //            x         x
        //              x    x
        //                 o
        //              p=3
        //                     o  p=2
        //                  .    .
        //               .         .
        //            .              .
        //         o .   .  .  .  .  . o
        //      p=0                     p=1
        //


  Vec3 getBasicTriangleVertCoord(int nvert)
  {
    //              (0,a/sqrt(3))
    //                     .
    //                     .
    //                     o 2
    //                    / \
    //                   /   \
    //                  /     \
    //               0 o - - - o 1
    //                .         .
    //               .           .
    //  (-a/2,-a/(2*sqrt(3))     (a/2,-a/(2*sqrt(3))
    //
    final float C0 = EDGE_LENGTH / 2;
    final float C1 = EDGE_LENGTH / (2*sqrt(3.0));
    final float C2 = EDGE_LENGTH / sqrt(3);
    final float C3 = EDGE_LENGTH * sqrt(2.0/3.0);

    final float V0x =  0;
    final float V0y = -C0;
    final float V0z = -C1;
    final float V1x =  0;
    final float V1y =  C0;
    final float V1z = -C1;
    final float V2x =  0;
    final float V2y =  0;
    final float V2z =  C2;
    //
    //  (-a/2,a/(2*sqrt(3))     (a/2,a/(2*sqrt(3))
    //               .           .
    //                .         .
    //               5 o - - - o 4
    //                  \     /
    //                   \   /
    //                    \ /
    //                     o 3
    //                     .
    //                     .
    //                (0,-a/sqrt(3))
    //
    final float V3x =  C3;
    final float V3y =   0;
    final float V3z = -C2;
    final float V4x =  C3;
    final float V4y =  C0;
    final float V4z =  C1;
    final float V5x =  C3;
    final float V5y = -C0;
    final float V5z =  C1;

    assert nvert>=0 && nvert<6;

    Vec3 ans = new Vec3(0.0, 0.0, 0.0);

    switch(nvert)
    {
      case 0:
        ans = new Vec3(V0x, V0y, V0z);
        break;
      case 1:
        ans = new Vec3(V1x, V1y, V1z);
        break;
      case 2:
        ans = new Vec3(V2x, V2y, V2z);
        break;
      case 3:
        ans = new Vec3(V3x, V3y, V3z);
        break;
      case 4:
        ans = new Vec3(V4x, V4y, V4z);
        break;
      case 5:
        ans = new Vec3(V5x, V5y, V5z);
        break;
    }

    return ans;
  }


  void initialConfiguration()
  {
    Vec3 vert0, vert1, vert2, vert3, vert4, vert5;

    vert0 = getBasicTriangleVertCoord(0);
    vert1 = getBasicTriangleVertCoord(1);
    vert2 = getBasicTriangleVertCoord(2);
    vert3 = getBasicTriangleVertCoord(3);
    vert4 = getBasicTriangleVertCoord(4);
    vert5 = getBasicTriangleVertCoord(5);

    float xHeight = vert3.x;

    for (int n=0; n<N_TRIANGLES; n++) {
      if ( n%2==0 ) {
        posx[id(n,0)] = vert0.x;
        posy[id(n,0)] = vert0.y;
        posz[id(n,0)] = vert0.z;
        posx[id(n,1)] = vert1.x;
        posy[id(n,1)] = vert1.y;
        posz[id(n,1)] = vert1.z;
        posx[id(n,2)] = vert2.x;
        posy[id(n,2)] = vert2.y;
        posz[id(n,2)] = vert2.z;
      }
      else {
        posx[id(n,0)] = vert3.x;
        posy[id(n,0)] = vert3.y;
        posz[id(n,0)] = vert3.z;
        posx[id(n,1)] = vert4.x;
        posy[id(n,1)] = vert4.y;
        posz[id(n,1)] = vert4.z;
        posx[id(n,2)] = vert5.x;
        posy[id(n,2)] = vert5.y;
        posz[id(n,2)] = vert5.z;
      }
      for (int i=0; i<3; i++) { // shift in x-direction.
        posx[id(n,i)] += 2*xHeight*(n/2);
      }
    }

    for (int i=0; i<N_PARTICLES; i++) {
      velx[i] = 0.0;
      vely[i] = 0.0;
      velz[i] = 0.0;
    }
  }

  void shiftCenterOfGravityToOrigin()
  {
    float cogx = 0.0;  // center of gravity
    float cogy = 0.0;
    float cogz = 0.0;

    for (int p=0; p<N_PARTICLES; p++) {
      cogx += posx[p];
      cogy += posy[p];
      cogz += posz[p];
    }
    cogx /= float(N_PARTICLES);
    cogy /= float(N_PARTICLES);
    cogz /= float(N_PARTICLES);

    for (int p=0; p<N_PARTICLES; p++) {
      posx[p] -= cogx;
      posy[p] -= cogy;
      posz[p] -= cogz;
    }
  }

  Vec3[] leftBoundaryConfiguration(float t, float xOrg)
  {
    float twistFactor;

    if ( t<SPRING_CHAR_PERIOD*30 )
      twistFactor = 0.0;
    else
      twistFactor = 0.1;

    float angle = (PI*2 / SPRING_CHAR_PERIOD) * twistFactor * t;

    Vec3 vert0, vert1, vert2;

    vert0 = getBasicTriangleVertCoord(0);
    vert1 = getBasicTriangleVertCoord(1);
    vert2 = getBasicTriangleVertCoord(2);

    Vec3[] rotatedVerts = new Vec3[3];

    float y, z;

    y =  cos(angle)*vert0.y + sin(angle)*vert0.z;
    z = -sin(angle)*vert0.y + cos(angle)*vert0.z;
    rotatedVerts[0] = new Vec3(xOrg,y,z);

    y =  cos(angle)*vert1.y + sin(angle)*vert1.z;
    z = -sin(angle)*vert1.y + cos(angle)*vert1.z;
    rotatedVerts[1] = new Vec3(xOrg,y,z);

    y =   cos(angle)*vert2.y + sin(angle)*vert2.z;
    z = - sin(angle)*vert2.z + cos(angle)*vert2.y;
    rotatedVerts[2] = new Vec3(xOrg,y,z);

    return rotatedVerts;
  }


  Particles()
  {
    initialConfiguration();

    shiftCenterOfGravityToOrigin();

    for (int nt=0; nt<N_TRIANGLES; nt++) {
      for (int j=0; j<3; j++) {
        int pid = id(nt,j); // particle id
        for (int s=0; s<6; s++) {  // six springs.
          // each particles is connectd by 6 springs.
          sixSpringList[pid][s] = NULL_MARK;
        }
      }
    }
  }

  int id(int layerId, int vertexId)
  {
    //  each triangle's              each particle's
    //    vertex id                        id
    // - - - - - - - - - - - - - - - - - - - - - - - - - -
    //
    //        2                             8
    //      .   .       triangle          .   .
    //    0 . . . 1      layerId=2      6 . . . 7
    //
    //
    //   2 . . . 1                     5 . . . 4
    //    .   .         triangle        .   .
    //     0              layerId=1      3
    //
    //
    //        2                             2
    //      .   .       triangle          .   .
    //    0 . . . 1       layerId=0     0 . . . 1
    //
    return 3*layerId + vertexId;
  }


  int[] getSixSpringListForThisParticle(int particleId)
  {
    int[] list = new int[6];

    for (int i=0; i<6; i++) {
      list[i] = sixSpringList[particleId][i];
    }

    return list;
  }


  int numberOfAlreadyRegisteredSpring(int particleId)
  {
    int ans;

    for (int i=0; i<6; i++) {
      int val = sixSpringList[particleId][i];
      if ( val==NULL_MARK ) {
        ans = i;
        return ans;
      }
    }
    // alerady set all the six elements.
    ans = 6;
    return ans;
  }


  void sixSpringsAppend(int particleId, int springid)
  {
    int num = numberOfAlreadyRegisteredSpring(particleId);

    assert num >=0 && num<6;

    sixSpringList[particleId][num] = springid;
  }


  void display() {
    noStroke();
    fill(100,0,130);
    for (int p=0; p<N_PARTICLES; p++) {
      pushMatrix();
        float x = posx[p];
        float y = posy[p];
        float z = posz[p];
        translate(mapx(x), mapy(y), mapz(z));
        sphere(3);
      popMatrix();
    }
  }

}
