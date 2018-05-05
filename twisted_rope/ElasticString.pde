
class ElasticString
{


  
  private void rungeKuttaIncrement(int num, 
                                   float[] p, 
                                   float[] p1, 
                                   float[] dp, 
                                   float factor)
  {
    for (int j=0; j<num; j++) {
      p[j] = p1[j] + factor*dp[j];
    }  
  }


  
  //private void equationOfMotion(float q[], float dq[], float dt)
  //{
     
  //}
  
  private void equationOfMotion(float q[], float dq[], float dt) 
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


    for (int tl=1; tl<N_TRIANGLES-1; tl++) { // triangle layer. skip the ends.   
      for (int me=0; me<3; me++) {
        int pid = particles.id(tl,me); // particle id
        int[] splist = new int[6];
        splist = particles.getSixSpringListForThisParticle(pid);

        Vec3 forceSum = new Vec3(0.0, 0.0, 0.0);
        for (int s=0; s<6; s++) {
          SpringElement aSpring = springs.getElement(splist[s]);
          Vec3 pullForceFromTheSpring = aSpring.getPullForce(pid,q);
          forceSum.add(pullForceFromTheSpring);          
        }
        
        dq[6*pid+0] = ( q[6*pid+3] ) * dt; // dx = vx * dt
        dq[6*pid+1] = ( q[6*pid+4] ) * dt; // dy = vy * dt
        dq[6*pid+2] = ( q[6*pid+5] ) * dt; // dz = vz * dt
        dq[6*pid+3] = ( forceSum.x ) * dtm; // dvx = (fx/m)*dt 
        dq[6*pid+4] = ( forceSum.y ) * dtm; // dvy = (fy/m)*dt 
        dq[6*pid+5] = ( forceSum.z ) * dtm; // dvz = (fz/m)*dt 
      }
    }
  }


  void boundaryCondition(float t, float[] generalCoord)
  {
    Vec3[] verts = new Vec3[3];
    
    //for (int j=0; j<3; j++) { // three vertices at the bottom.
    //  verts[j].x = generalCoord[6*j+0];
    //  verts[j].y = generalCoord[6*j+1];
    //  verts[j].z = generalCoord[6*j+2];
    //}
    
    float z0 = generalCoord[0*6+2];
    
    verts = particles.lowerBoundaryConfiguration(t,z0);
    
    for (int j=0; j<3; j++) { // three vertices at the bottom.
      generalCoord[6*j+0] = verts[j].x;
      generalCoord[6*j+1] = verts[j].y;
      generalCoord[6*j+2] = verts[j].z;
    }
    
  }

  
  
  void rungeKutta()
  {
    final float ONE_SIXTH = 1.0/6.0;
    final float ONE_THIRD = 1.0/3.0;
    final int NN = 6*N_PARTICLES;  // (pos.x,y,z) and (vel.x,y,z)
  
    float[] qprev = new float[NN];
    float[] qwork = new float[NN];
    float[] dq1 = new float[NN];
    float[] dq2 = new float[NN];
    float[] dq3 = new float[NN];
    float[] dq4 = new float[NN];
      
    for (int n=0; n<N_PARTICLES; n++) {
      Vec3 pos = particles.getPos(n);
      Vec3 vel = particles.getVel(n);
      qprev[6*n+0] = pos.x;
      qprev[6*n+1] = pos.y;
      qprev[6*n+2] = pos.z;
      qprev[6*n+3] = vel.x;
      qprev[6*n+4] = vel.y;
      qprev[6*n+5] = vel.z;
    }
  
    //step 1
    equationOfMotion(qprev, dq1, dt);
    rungeKuttaIncrement(NN, qwork, qprev, dq1, 0.5);

    time += 0.5*dt;
    boundaryCondition(time, qwork);
  
    //step 2
    equationOfMotion(qwork, dq2, dt);
    rungeKuttaIncrement(NN, qwork, qprev, dq2, 0.5);

    boundaryCondition(time, qwork);
  
    //step 3
    equationOfMotion(qwork, dq3, dt);
    rungeKuttaIncrement(NN, qwork, qprev, dq3, 1.0);

    time += 0.5*dt;
    boundaryCondition(time, qwork);
  
    //step 4
    equationOfMotion(qwork, dq4, dt);

    //the result
    for (int tl=1; tl<N_TRIANGLES-1; tl++) { 
      // See boundaryCondition() for end points.
      for (int j=0; j<3; j++) { // three verteces in a triangle.
        int pid = particles.id(tl,j);        
        for (int i=0; i<6; i++) { // x,y,z,vx,vy,vz
          qwork[6*pid+i] =         qprev[6*pid+i] + (
                           ONE_SIXTH*dq1[6*pid+i]
                         + ONE_THIRD*dq2[6*pid+i]
                         + ONE_THIRD*dq3[6*pid+i]
                         + ONE_SIXTH*dq4[6*pid+i]
                         );
        } 
      }
    }
    boundaryCondition(time, qwork);
  
    for (int p=0; p<N_PARTICLES; p++) { 
      for (int i=0; i<6; i++) { // x,y,z,vx,vy,vz
        float work = qwork[6*p+i];
        if (i==0)
          particles.setPosX(p,work);
        else if (i==1)
          particles.setPosY(p,work);
        else if (i==2)
          particles.setPosZ(p,work);
        else if (i==3)
          particles.setVelX(p,work);
        else if (i==4)
          particles.setVelY(p,work);
        else if (i==5)
          particles.setVelZ(p,work);
      } 
    }
    
  }
  


  void draw()
  {
    particles.draw();
    springs.draw(particles.getPos());
  }
  
  //float totalEnergy() 
  //{
  //  float sum_kinetic = 0.0;
  //  float sum_potential = 0.0;
  //  float springc = spring.getConst();
  //  for (int i=0; i<N_PARTICLES; i++) {
  //    float posx = pos[i].x;
  //    float posy = pos[i].y;
  //    float posz = pos[i].z;
  //    float velx = vel[i].x;
  //    float vely = vel[i].y;
  //    float velz = vel[i].z;
  //    sum_potential += 0.5*springc*util.squaredSum(
  //    sum_kinetic += 0.5*PARTICLE_MASS*(velx*velx
  //                            +vely*vely
  //                            +velz*velz);
  //  }
  //  float kinetic = 0.5*PARTICLE_MASS*(
  //}
}
