
class Particles
{
  float twistAngle=0.0;

  final float END_POINTS_SEPARATION = LADDER_LENGTH*1.0;

  float[] posx = new float[N_PARTICLES];
  float[] posy = new float[N_PARTICLES];
  float[] posz = new float[N_PARTICLES];
  float[] velx = new float[N_PARTICLES];
  float[] vely = new float[N_PARTICLES];
  float[] velz = new float[N_PARTICLES];

  int[][] connectedSpringList = new int[N_PARTICLES][3];
            // each particle are connected with 6 springs.
            //
            //
            //   o---------o  upper layer
            //   |         |
            //   3         |
            //   |         |
            //   X----1----o  self layer
            //   |         |
            //   2         |
            //   |         |
            //   o---------o  lower layer

  final int NULL_MARK = -1; // used to count already set elements.

  void initialConfiguration()
  {
    for (int j=0; j<2; j++) {
      int l = 0;  // uppermost ladder
      int p = id(l,j);
      posx[p] = -SPRING_NATURAL_LENGTH/2 + SPRING_NATURAL_LENGTH*j;
      posy[p] = 0.0;
      posz[p] = 0.0;
    }

    float deltaZ = END_POINTS_SEPARATION / (N_PAIRS-1);

    for (int l=1; l<N_PAIRS; l++) { // ladder
      for (int j=0; j<2; j++) {
        int p = id(l,j);
        posx[p] = -SPRING_NATURAL_LENGTH/2 + SPRING_NATURAL_LENGTH*j;
        posy[p] = 0.0;
        posz[p] = - deltaZ * l;

//      posx[p] += SPRING_NATURAL_LENGTH*random(0.01);
//      posy[p] += SPRING_NATURAL_LENGTH*random(0.01);
//      posz[p] += SPRING_NATURAL_LENGTH*random(0.01);
      }
    }

    float omega = LADDER_TWIST_RATE_OMEGA;
    float vy0 = omega*SPRING_NATURAL_LENGTH/2;
    for (int l=0; l<N_PAIRS; l++) {
      for (int j=0; j<2; j++) {
        int p = id(l,j);
        velx[p] = 0.0;
  //    vely[p] = 0.0;
        velz[p] = 0.0;
        if (j%2==0)
          vely[p] = vy0;
        else
          vely[p] = -vy0;
      }
    }
  }

  void upperBoundaryConfiguration(float time, Vec3[] ans)
  {
    float angle = LADDER_TWIST_RATE_OMEGA * time;

    for (int j=0; j<2; j++) {
      float x0 = -SPRING_NATURAL_LENGTH/2 + SPRING_NATURAL_LENGTH*j;
      float y0 = 0.0;
      float z0 = 0.0;
      float x =  cos(angle)*x0 + sin(angle)*y0;
      float y = -sin(angle)*x0 + cos(angle)*y0;
      ans[j] = new Vec3(x,y,z0);
    }
  }


  Particles()
  {
    initialConfiguration();

    for (int l=0; l<N_PAIRS; l++) {
      for (int j=0; j<2; j++) {
        int pid = id(l,j); // particle id
        for (int s=0; s<3; s++) {
          // each particle is connected with 3 springs.
          connectedSpringList[pid][s] = NULL_MARK;
        }
      }
    }
  }

  int id(int layerId, int vertexId)
  {
    return 2*layerId + vertexId;
  }


  int[] getConnectedSpingListForThisParticle(int particleId)
  {
    int[] list = new int[3];

    for (int i=0; i<3; i++) {
      list[i] = connectedSpringList[particleId][i];
    }

    return list;
  }


  int numberOfConnectedSpringsToThisParticle(int particleId)
  {
    int ans;

    for (int i=0; i<3; i++) {
      int val = connectedSpringList[particleId][i];
      if ( val==NULL_MARK ) {
        ans = i;
        return ans;
      }
    }
    // alerady set all the six elements.
    ans = 3;
    return ans;
  }


  void connectedSpringsAppend(int particleId, int springid)
  {
    int num = numberOfConnectedSpringsToThisParticle(particleId);

    assert num >=0 && num<3;

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
