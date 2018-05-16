
class Particles
{
  float twistAngle=0.0;

  float[] posx = new float[N_PARTICLES];
  float[] posy = new float[N_PARTICLES];
  float[] posz = new float[N_PARTICLES];
  float[] velx = new float[N_PARTICLES];
  float[] vely = new float[N_PARTICLES];
  float[] velz = new float[N_PARTICLES];

  int[][] connectedSpringList = new int[N_PARTICLES][6];
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


  void boundaryConfiguration(float t, Vec3[] verts)
  {
    float r = STICK_LENGTH;
    float leftEdgeX = -STICK_LENGTH/2;
    float leftEdgeY = 0.0;
    float x = leftEdgeX + r*cos(EDGE_ROTATION_OMEGA*t);
    float y = leftEdgeY + r*cos(EDGE_ROTATION_OMEGA*t);
    for (int j=0; j<3; j++) {


    }
  }

  void initialConfiguration()
  {
    float deltaX = STICK_LENGTH / (N_TRIANGLES-1);
    float angle;

    for (int l=0; l<N_TRIANGLES; l++) {
      for (int j=0; j<3; j++) {
        int p = id(l,j);
        if ( l%2==0 )
          angle = (TWO_PI/3)*j;
        else
          angle = (TWO_PI/3)*j + TWO_PI/6;
        posx[p] = -STICK_LENGTH/2 + deltaX*l
        posy[p] = STICK_RADIUS*cos(angle);
        posz[p] = STICK_RADIUS*sin(angle);

        posx[p] += EDGE_LENGTH*(random(0.02)-0.01);
        posy[p] += EDGE_LENGTH*(random(0.02)-0.01);
        posz[p] += EDGE_LENGTH*(random(0.02)-0.01);
      }
    }

    for (int l=0; l<N_TRIANGLES; l++) {
      for (int j=0; j<3; j++) {
        int p = id(l,j);
        velx[p] = 0.0;
        vely[p] = 0.0;
        velz[p] = 0.0;
      }
    }
  }


  Particles()
  {
    initialConfiguration();

    for (int nt=0; nt<N_TRIANGLES; nt++) {
      for (int j=0; j<3; j++) {
        int pid = id(nt,j); // particle id
        for (int s=0; s<6; s++) {  // six springs.
          // each particles is connectd by 6 springs.
          connectedSpringList[pid][s] = NULL_MARK;
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


  int[] getConnectedSpingListForThisParticle(int particleId)
  {
    int[] list = new int[6];

    for (int i=0; i<6; i++) {
      list[i] = connectedSpringList[particleId][i];
    }

    return list;
  }


  int numberOfConnectedSpringsToThisParticle(int particleId)
  {
    int ans;

    for (int i=0; i<6; i++) {
      int val = connectedSpringList[particleId][i];
      if ( val==NULL_MARK ) {
        ans = i;
        return ans;
      }
    }
    ans = 6;
    return ans;
  }


  void connectedSpringsAppend(int particleId, int springid)
  {
    int num = numberOfConnectedSpringsToThisParticle(particleId);

    assert num >=0 && num<6;

    connectedSpringList[particleId][num] = springid;
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
        sphere(1);
      popMatrix();
    }
  }


  float energy()
  {
    float sum = 0.0;
    for (int p=0; p<N_PARTICLES; p++) {
      float vxsq = pow(velx[p],2);
      float vysq = pow(vely[p],2);
      float vzsq = pow(velz[p],2);
      float vsq = vxsq + vysq + vzsq;
      sum += 0.5*PARTICLE_MASS*vsq;
    }
    return sum;
  }

}
