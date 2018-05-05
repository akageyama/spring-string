
class ElasticString
{


  
  private void rungeKuttaIncrement(int num, 
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
    for (int j=0; j<num; j++) {
      posx[j] = posx1[j] + factor*dposx[j];
      posy[j] = posy1[j] + factor*dposy[j];
      posz[j] = posz1[j] + factor*dposz[j];
      velx[j] = velx1[j] + factor*dvelx[j];
      vely[j] = vely1[j] + factor*dvely[j];
      velz[j] = velz1[j] + factor*dvelz[j];
    }  
  }


  
  private void equationOfMotion(float posx[],
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
          Vec3 pullForceFromTheSpring = aSpring.getPullForce(pid,posx,
                                                                 posy,
                                                                 posz);
          forceSum.add(pullForceFromTheSpring);          
        }
        
        dposx[pid] = velx[pid] * dt;
        dposy[pid] = vely[pid] * dt;
        dposz[pid] = velz[pid] * dt;
        dvelx[pid] = forceSum.x * dtm;
        dvely[pid] = forceSum.y * dtm;
        dvelz[pid] = forceSum.z * dtm;
      }
    }
  }


  void boundaryCondition(float t, 
                         float[] posx,
                         float[] posy,
                         float[] posz)
  {
    Vec3[] verts = new Vec3[3];
        
    float z0 = posz[0];
    
    verts = particles.lowerBoundaryConfiguration(t,z0);
    
    for (int j=0; j<3; j++) { // three vertices at the bottom.
      posx[j] = verts[j].x;
      posy[j] = verts[j].y;
      posz[j] = verts[j].z;
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
    
      
    for (int n=0; n<NN; n++) {
      posxprev[n] = particles.getPosX(n);
      posyprev[n] = particles.getPosY(n);
      poszprev[n] = particles.getPosZ(n);
      velxprev[n] = particles.getVelX(n);
      velyprev[n] = particles.getVelY(n);
      velyprev[n] = particles.getVelZ(n);
    }
  
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

    time += 0.5*dt;
    boundaryCondition(time, 
                      posxwork,
                      posywork,
                      poszwork);
  
    //step 2
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
                        dposx2,
                        dposy2,
                        dposz2,
                        dvelx2,
                        dvely2,
                        dvelz2,
                        0.5);

    boundaryCondition(time, 
                      posxwork,
                      posywork,
                      poszwork);
  
    //step 3
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

    boundaryCondition(time, 
                      posxwork,
                      posywork,
                      poszwork);
    

    time += 0.5*dt;
    boundaryCondition(time, 
                      posxwork,
                      posywork,
                      poszwork);
  
    //step 4
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
    
    //the result
    for (int tl=1; tl<N_TRIANGLES-1; tl++) { 
      for (int j=0; j<3; j++) { // three verteces in a triangle.
        int pid = particles.id(tl,j);        
        posxwork[pid] =            posxprev[pid] + (
                           ONE_SIXTH*dposx1[pid]
                         + ONE_THIRD*dposx2[pid]
                         + ONE_THIRD*dposx3[pid]
                         + ONE_SIXTH*dposx4[pid]
                         );
        posywork[pid] =            posyprev[pid] + (
                           ONE_SIXTH*dposy1[pid]
                         + ONE_THIRD*dposy2[pid]
                         + ONE_THIRD*dposy3[pid]
                         + ONE_SIXTH*dposy4[pid]
                         );
        poszwork[pid] =            poszprev[pid] + (
                           ONE_SIXTH*dposz1[pid]
                         + ONE_THIRD*dposz2[pid]
                         + ONE_THIRD*dposz3[pid]
                         + ONE_SIXTH*dposz4[pid]
                         );
        velxwork[pid] =            velxprev[pid] + (
                           ONE_SIXTH*dvelx1[pid]
                         + ONE_THIRD*dvelx2[pid]
                         + ONE_THIRD*dvelx3[pid]
                         + ONE_SIXTH*dvelx4[pid]
                         );
        velywork[pid] =            velyprev[pid] + (
                           ONE_SIXTH*dvely1[pid]
                         + ONE_THIRD*dvely2[pid]
                         + ONE_THIRD*dvely3[pid]
                         + ONE_SIXTH*dvely4[pid]
                         );
        velzwork[pid] =            velzprev[pid] + (
                           ONE_SIXTH*dvelz1[pid]
                         + ONE_THIRD*dvelz2[pid]
                         + ONE_THIRD*dvelz3[pid]
                         + ONE_SIXTH*dvelz4[pid]
                         );
      }
    }
    boundaryCondition(time, 
                      posxwork,
                      posywork,
                      poszwork);
  
    particles.setPosX(posxwork);
    particles.setPosY(posywork);
    particles.setPosZ(poszwork);
    particles.setVelX(velxwork);
    particles.setVelY(velywork);
    particles.setVelZ(velzwork);   
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
