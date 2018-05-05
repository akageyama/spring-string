

class Particles 
{
  private float[] posx = new float[N_PARTICLES];
  private float[] posy = new float[N_PARTICLES];
  private float[] posz = new float[N_PARTICLES];
  private float[] velx = new float[N_PARTICLES];
  private float[] vely = new float[N_PARTICLES];
  private float[] velz = new float[N_PARTICLES];

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
  final float C0 = EDGE_LENGTH/2;
  final float C1 = EDGE_LENGTH/(2*sqrt(3.0));
  final float C2 = EDGE_LENGTH/sqrt(3);
  final float C3 = EDGE_LENGTH * sqrt(2.0/3.0);
  
  final float V0x = -C0;
  final float V0y = -C1;
  final float V0z =  0;
  final float V1x =  C0;
  final float V1y = -C1;
  final float V1z =  0;
  final float V2x =  0;
  final float V2y =  C2;
  final float V2z =  0;
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
  final float V3x =   0;
  final float V3y = -C2;
  final float V3z =  C3;
  final float V4x =  C0;
  final float V4y =  C1;
  final float V4z =  C3;
  final float V5x = -C0;
  final float V5y =  C1;
  final float V5z =  C3;  
    
  private int[][] sixSpringList = new int[N_PARTICLES][6];
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


  private final int NULL_MARK = -1; // used to count already set elements.
  
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

  private void initialConfiguration()
  {

    
    for (int n=0; n<N_TRIANGLES; n++) {
      if ( n%2==0 ) {        
        pos[id(n,0)] = new Vec3(V0x,V0y,V0z); 
        pos[id(n,1)] = new Vec3(V1x,V1y,V1z); 
        pos[id(n,2)] = new Vec3(V2x,V2y,V2z); 
      }
      else {
        pos[id(n,0)] = new Vec3(V3x,V3y,V3z); 
        pos[id(n,1)] = new Vec3(V4x,V4y,V4z); 
        pos[id(n,2)] = new Vec3(V5x,V5y,V5z); 
      }
      for (int i=0; i<3; i++) { // shift in z-direction.
        pos[id(n,i)].z += 2*C3*(n/2);
      }
    }
    
    for (int i=0; i<N_PARTICLES; i++) {
      vel[i] = new Vec3(0.0, 0.0, 0.0);
    }
  }
  
  
  void shiftCenterOfGravityToOrigin()
  {
    Vec3 cog = new Vec3(0.0, 0.0, 0.0); // center of gravity
    
    for (int p=0; p<N_PARTICLES; p++) {
      cog.add(pos[p]);
    }
    cog.divide(float(N_PARTICLES));
    
    for (int p=0; p<N_PARTICLES; p++) {
      pos[p].subtract(cog);
    }
  }
  
  Vec3[] lowerBoundaryConfiguration(float t, float z)
  {
    float factor = 0.1;
    float angle = (PI*2 / SPRING_CHAR_PERIOD) * factor * t;
    
    Vec3[] verts = new Vec3[3];
    
    float x, y;
    
    x = cos(angle)*V0x - sin(angle)*V0y;
    y = sin(angle)*V0x + cos(angle)*V0y;
    verts[0] = new Vec3(x,y,z);
    
    x = cos(angle)*V1x - sin(angle)*V1y;
    y = sin(angle)*V1x + cos(angle)*V1y;
    verts[1] = new Vec3(x,y,z);
    
    x = cos(angle)*V2x - sin(angle)*V2y;
    y = sin(angle)*V2x + cos(angle)*V2y;
    verts[2] = new Vec3(x,y,z);

    return verts;
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
  
  
  float getPosX(int p)
  {
    return posx[p];
  }
  
  float getPosY(int p)
  {
    return posy[p];
  }
  
  float getPosZ(int p)
  {
    return posz[p];
  }
  
  void setPosX(int p, float x)
  {
    posx[p] = x;
  }  
  
  void setPosX(float[] x)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      posx[i] = x[i];
    }
  }
  
  void setPosY(int p, float y)
  {
    posy[p] = y;
  }
  
  void setPosY(float[] y)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      posy[i] = y[i];
    }
  }
  
  void setPosZ(int p, float z)
  {
    posz[p] = z;
  }

  void setPosZ(float[] z)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      posz[i] = z[i];
    }
  }
    
  float getVelX(int p)
  {
    return velx[p];
  }
  
  float getVelY(int p)
  {
    return vely[p];
  }
  
  float getVelZ(int p)
  {
    return velz[p];
  }
  
  void setVelX(int p, float vx)
  {
    velx[p] = vx;
  }
  
  void setVelX(float[] vx)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      velx[i] = vx[i];
    }      
  }
  
  void setVelY(int p, float vy)
  {
    vely[p] = vy;
  }
  
  void setVelY(float[] vy)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      vely[i] = vy[i];
    }      
  }
  
  void setVelZ(int p, float vz)
  {
    velz[p] = vz;
  }
  
  void setVelZ(float[] vz)
  {
    for (int i=0; i<N_PARTICLES; i++) {
      velz[i] = vz[i];
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
  
  
  private int numberOfAlreadyRegisteredSpring(int particleId)
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
  
  
  void draw() {
    noStroke();
    fill(100,0,130);
    for (int i=0; i<N_PARTICLES; i++) {
      pushMatrix();
        float x = posx[i];
        float y = posy[i];
        float z = posz[i];
        translate(mapx(x), mapy(y), mapz(z));      
        sphere(3);
      popMatrix();
    }
  }

}
