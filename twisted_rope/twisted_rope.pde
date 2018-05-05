/* 

  twisted_rope.pde
 
 
 Developed 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.05
 
 
 */


//MouseCamera mouseCamera;

float time = 0.0;
int step = 0;

final int N_TRIANGLES = 6;
final int N_PARTICLES = N_TRIANGLES*3;
final float EDGE_LENGTH = 0.3;
final float PARTICLE_MASS = 0.1;
final float SPRING_CHAR_PERIOD = 0.1; // second

float dt = SPRING_CHAR_PERIOD*0.01;

float x_coord_min = -3.0;
float x_coord_max =  3.0;
float y_coord_min = x_coord_min;
float y_coord_max = x_coord_max;
float z_coord_min = x_coord_min;
float z_coord_max = x_coord_max;


class Rotor {
  float rotx = 0;
  float roty = 0;
  float rotz = 0;
  float delta = PI/100;
  boolean toggle_keep_rotate_x = false;
  boolean toggle_keep_rotate_y = false;
  boolean toggle_keep_rotate_z = false;

  Rotor(float rotx, float roty, float rotz) {
    this.rotx = rotx;
    this.roty = roty;
    this.rotz = rotz;
  }
  
  
  float fit(float angle) {
    float ans = angle;
    if ( ans > 2*PI ) ans -= 2*PI;
    if ( ans < 0    ) ans += 2*PI;
    return ans;
  }
  
  void update() {
    if ( toggle_keep_rotate_x ) rotx = fit(rotx + delta);
    if ( toggle_keep_rotate_y ) roty = fit(roty + delta);
    if ( toggle_keep_rotate_z ) rotz = fit(rotz + delta);
  }
  
  void toggle(char key) {
    switch (key) {
      case 'x':
        toggle_keep_rotate_x = !toggle_keep_rotate_x;
        break;
      case 'y':
        toggle_keep_rotate_y = !toggle_keep_rotate_y;
        break;
      case 'z':
        toggle_keep_rotate_z = !toggle_keep_rotate_z;
        break;
    }           
  }
}

Rotor rotor = new Rotor(0,0,0); 



class Vec3 {
  float x, y, z;
  
  Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  
  Vec3(Vec3 rhs) {
    this.x = rhs.x;
    this.y = rhs.y;
    this.z = rhs.z;
  }
  
  void add(Vec3 v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
  }
  
  void multiply(float a) {
    this.x *= a;
    this.y *= a;
    this.z *= a;
  }
}

class Util
{
  Util() 
  {
  }
  
  
  float squaredSum(float a, float b, float c)
  {
    return a*a+b*b+c*c;
  }
}


Util util = new Util();



class Particles 
{
  private Vec3[] pos = new Vec3[N_PARTICLES];
  private Vec3[] vel = new Vec3[N_PARTICLES];
    
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
  
  
  Particles()  
  {
    initialConfiguration();
    
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
  
  
  Vec3 getPos(int particleId)
  {
    return pos[particleId];
  }
  
  void setPosX(int particleId, float x)
  {
    pos[particleId].x = x;
  }  
  
  void setPosY(int particleId, float y)
  {
    pos[particleId].y = y;
  }
  
  void setPosZ(int particleId, float z)
  {
    pos[particleId].z = z;
  }
  
  Vec3 getVel(int particleId)
  {
    return vel[particleId];
  }
  
  void setVelX(int particleId, float vx)
  {
    vel[particleId].x = vx;
  }
  
  void setVelY(int particleId, float vy)
  {
    vel[particleId].y = vy;
  }
  
  void setVelZ(int particleId, float vz)
  {
    vel[particleId].z = vz;
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
        float x = pos[i].x;
        float y = pos[i].y;
        float z = pos[i].z;
        translate(mapx(x), mapy(y), mapz(z));      
        sphere(3);
      popMatrix();
    }
  }

}

Particles particles = new Particles();



class SpringElement
{
  //  A spring, with two ends called alpha and beta.
  //           alpha         beta
  //              O=========O 
  
  private int endParticleIdAlpha;
  private int endParticleIdBeta;
  private float springConst;
  
  SpringElement()
  {
  }
  
  SpringElement(float springConst, int alpha, int beta)
  {
    this.springConst = springConst;
    endParticleIdAlpha = alpha;
    endParticleIdBeta  = beta;
  }
  

  int getAlpha()
  {
    return endParticleIdAlpha;
  }
  
  int getBeta()
  {
    return endParticleIdBeta;
  }
  
  Vec3 getPullForce(int particleId, float[] generalCoord) 
  {
    Vec3 force = new Vec3(0.0, 0.0, 0.0);
    
    float ax = generalCoord[6*endParticleIdAlpha+0];
    float ay = generalCoord[6*endParticleIdAlpha+1];
    float az = generalCoord[6*endParticleIdAlpha+2];
    float bx = generalCoord[6*endParticleIdBeta +0];
    float by = generalCoord[6*endParticleIdBeta +1];
    float bz = generalCoord[6*endParticleIdBeta +2];
    
    float distance = dist(ax, ay, az, bx, by, bz);
    
    float pullForceAmplitude = springConst * (distance - EDGE_LENGTH);
    
    if (particleId==endParticleIdAlpha) {
      force = new Vec3((bx-ax)/distance,
                       (by-ay)/distance,
                       (bz-az)/distance);  // unit vector from alpha to beta.
      force.multiply(pullForceAmplitude);
    }
    else if (particleId==endParticleIdBeta) {
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
  
  
  void draw() {
    Vec3 alpha = particles.pos[endParticleIdAlpha];
    Vec3 beta  = particles.pos[endParticleIdBeta];
    
    float ax = mapx(alpha.x);
    float ay = mapy(alpha.y);
    float az = mapz(alpha.z);
    float bx = mapx( beta.x);
    float by = mapy( beta.y);
    float bz = mapz( beta.z);

    line(ax,ay,az,bx,by,bz);
  }
  
}


class Springs 
{  
  private final int N_SPRINGS = 3*N_TRIANGLES + 6*(N_TRIANGLES-1);
    // 3*N_TRIANGLES: In each "hirozntal" layer, a triangle has three edges.
    // 6*(N_TRIANGLES-1): Each vertex in the triangle has two springs 
    //                    connected to two vertices in the lower layer.
  
  private SpringElement[] element = new SpringElement[N_SPRINGS];
    //                           
    //               
    //                H                           
    //              x   x  "self triangle layer"
    //            x       x
    //          H  x  x  x  H
    //           \         /
    //            \       /
    //          L .\. . ./. L
    //            . \   / .
    //              .\ /.   "lower triangle layer"
    //                L

    //
    //             U     U     U
    //        .   / \   / \   /
    //         \ /   \ /   \ /
    //        . H x x H x x H .
    //         / \   / \   / \
    //        .   \ /   \ /   \
    //             L     L     L
    //
    // each edge in the layer 'H' has two
    // springs connecting to a vertex
    // in the lower layer 'L'.

  Springs(float characteristicPeriod)
  {
    float omega = PI*2 / characteristicPeriod;            
    float spc = PARTICLE_MASS * omega * omega;
              // spc = spring constant:  omega^2 = spc / mass

    int sCtr = 0; // spring counter
    for (int tl=0; tl<N_TRIANGLES; tl++) { // for "hirizontal" triangles.
      int pId0 = particles.id(tl,0); // 1st vertex in the triangle
      int pId1 = particles.id(tl,1); // 2nd 
      int pId2 = particles.id(tl,2); // 3rd

      register(spc, sCtr++, pId0, pId1);       
      register(spc, sCtr++, pId1, pId2);       
      register(spc, sCtr++, pId2, pId0);             
    }
    for (int tl=1; tl<N_TRIANGLES; tl++) { // skip the lowest layer.
      for (int me=0; me<3; me++) {
        //
        // when tl=even
        //
        // each vertex in the layer 'S' has two
        // springs connecting to two vertices
        // in the upper and lower layers 'U' and 'L'.
        //
        //    vertexId (0, 1, 2) in each layer.
        //
        //             0     1     2    upper layer
        //        .   / \   / \   /
        //         \ /   \ /   \ /
        //        . 0 x x 1 x x 2 .    tl (even layer)
        //         / \   / \   / \
        //        .   \ /   \ /   \
        //             0     1     2    lower layer
        //
        // Connection table. me and its counterparts.
        //
        //       same layer    lower layer
        //            /  \     /   \
        //   me  |  m1   m2   l1   l2
        //   ----+--------------------
        //    0  |   1    2    0    2
        //    1  |   2    0    1    0
        //    2  |   0    1    2    1
        //       +---------------------
        //       |  k1   k2   me   k2
        //
        //
        // when tl=odd
        //
        //       same layer    lower layer
        //            /  \     /   \
        //   me  |  m1   m2   l1   l2
        //   ----+--------------------
        //    0  |   1    2    0    1
        //    1  |   2    0    1    2
        //    2  |   0    1    2    0
        //       +--------------------
        //       |  k1   k2   me   k1
            
        int k1 = (me+1) % 3;
        int k2 = (me+2) % 3;

        int myPid = particles.id(tl,me);

        if ( tl%2==0 ) {
          register(spc, sCtr++, myPid, particles.id(tl-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(tl-1,k2)); // lower layer
        }
        else {
          register(spc, sCtr++, myPid, particles.id(tl-1,me)); // lower layer
          register(spc, sCtr++, myPid, particles.id(tl-1,k1)); // lower layer
        }
      }
    }
    println("debug: sCtr = ", sCtr);
    
  }
  
  

  
  private void register(float springConst, int springId, 
                        int alpha, int beta)
  {
    //
    // ids of particles on the both ends
    //           alpha         beta
    //             \           /
    //              O=========O 
    //

println(" springId = ", springId);
    element[springId] = new SpringElement(springConst,alpha,beta);

    particles.sixSpringsAppend(alpha, springId);
    particles.sixSpringsAppend(beta,  springId);
  }



  
  void draw()
  {
    stroke(150, 100, 70);

    for (int s=0; s<N_SPRINGS; s++) {
      element[s].draw();
    }
  }
  
  SpringElement getElement(int n)
  {
    return element[n];
  }
  

}

Springs springs = new Springs(SPRING_CHAR_PERIOD);



class ElasticString
{

  ElasticString()
  {
  }
  
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
  
    //step 2
    equationOfMotion(qwork, dq2, dt);
    rungeKuttaIncrement(NN, qwork, qprev, dq2, 0.5);
  
    //step 3
    equationOfMotion(qwork, dq3, dt);
    rungeKuttaIncrement(NN, qwork, qprev, dq3, 1.0);
  
    //step 4
    equationOfMotion(qwork, dq4, dt);

  
    //the result
    for (int tl=1; tl<N_TRIANGLES-1; tl++) { 
      // See boundaryCondition() for end points.
      for (int j=0; j<3; j++) { // three verteces in a triangle.
        int pid = particles.id(tl,j);        
        for (int i=0; i<6; i++) { // x,y,z,vx,vy,vz
          float newval = qprev[6*pid+i] + (
                                    ONE_SIXTH*dq1[6*pid+i]
                                  + ONE_THIRD*dq2[6*pid+i]
                                  + ONE_THIRD*dq3[6*pid+i]
                                  + ONE_SIXTH*dq4[6*pid+i]
                                  );
          if (i==0)
            particles.setPosX(pid,newval);
          else if (i==1)
            particles.setPosY(pid,newval);
          else if (i==2)
            particles.setPosZ(pid,newval);
          else if (i==3)
            particles.setVelX(pid,newval);
          else if (i==4)
            particles.setVelY(pid,newval);
          else if (i==5)
            particles.setVelZ(pid,newval);
        } 
      }
    }
  }
  


  void draw()
  {
    particles.draw();
    springs.draw();
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

ElasticString elasticString = new ElasticString();




float norma(float x) {
  float s= width / (x_coord_max-x_coord_min);
  return s*x;
}

float mapx(float x) {
//  x = min(x,x_coord_max);
//  x = max(x,x_coord_min);
  return norma(x);
}

float mapy(float y) {
//  y = min(y,y_coord_max);
//  y = max(y,y_coord_min);
  return -norma(y);
}

float mapz(float z) {
//  z = min(z,z_coord_max);
//  z =s max(z,z_coord_min);
  return -norma(z);
}



void draw_axes_xyz() {
        stroke(100, 100, 100);
        line(mapx(x_coord_min), 0, 0, mapx(x_coord_max), 0, 0);
        line(0, mapy(y_coord_min), 0, 0, mapy(y_coord_max), 0);
        line(0, 0, mapz(z_coord_min), 0, 0, mapz(z_coord_max));
}




void setup() {
  size(800,700,P3D);
  background(255);
  frameRate(60);
  
  //mouseCamera = new MouseCamera(10, 0, 0, (height/2.0)/tan(PI*30.0/180.0), 0, 0, 0, 0, -1, 0); // MouseCameraの生成
  //camera(x_coord_max, x_coord_max, x_coord_max, 0, 0, 0, 0, -1, 0);

  //
  //   +----                       y=0
  //   |  VERTICAL_MARGIN
  //   +----                       y=y1
  //   |
  //   |
  //   |
  //   |
  //   |
  //   +----                       y=y2
  //   |  VERTICAL_MARGIN
  //   +----                       y=height



 
}







//float total_energy()
//{
//    float x1 = part[0].pos.x;
//    float v1 = part[0].pos.vx;
//    float x2 = part[1].pos.x;
//    float v2 = part[1].pos.vx;

//    float v1sq = v1*v1;
//    float v2sq = v2*v2;

//    float y1 = parabola_func_upper(x1);
//    float y2 = parabola_func_lower(x2);
    
//    float y1dot = parabola_func_upper_derivative(x1)*v1;
//    float y2dot = parabola_func_lower_derivative(x2)*v2;
//    float y1dotsq = y1dot*y1dot;
//    float y2dotsq = y2dot*y2dot;
    
//    float s = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

//    float m1 = part[0].PARTICLE_MASS;
//    float m2 = part[1].PARTICLE_MASS;
//    float kinetic_e = 0.5*(m1*(v1sq+y1dotsq)
//                          +m2*(v2sq+y2dotsq));
//    float potential = 0.5*SPRING_K*(s-SPRING_L0)*(s-SPRING_L0);

//    return(kinetic_e + potential);
//}




void shoot() {
  
    for (int n=0; n<1; n++) { // to speed up the display
      elasticString.rungeKutta();
      //boundaryCondition();
      time += dt;
      step += 1;
      if ( step%10 == 0 ) {
        println("step=", step, " time=", time);
      }
    }
  
    background(255);
    pushMatrix();
      translate(width/2,height/2);
      rotateZ(rotor.rotz);
      rotateX(rotor.rotx);
      rotateY(rotor.roty);
      
      draw_axes_xyz();
      elasticString.draw();
    popMatrix();               
    
    //if ( step%1000 == 0 ) {
    //  println("step = ", step," time = ", time," energy = ",total_energy());
    //}
    

}


void draw() {

    //mouseCamera.update();
    
    rotor.update();    
    shoot(); 
}


void keyPressed() {
  switch (key) {
  case 'x':
    rotor.toggle('x');
    break;
  case 'y':
    rotor.toggle('y');
    break;
  case 'z':
    rotor.toggle('z');
    break;
  }
}


void keyReleased() {
  switch (key) {
  case 'x':
    rotor.toggle('x');
    break;
  case 'y':
    rotor.toggle('y');
    break;
  case 'z':
    rotor.toggle('z');
    break;
  }
}

//void mousePressed() {
//    mouseCamera.mousePressed();
//}
//void mouseDragged() {
//    mouseCamera.mouseDragged();
//}
//void mouseWheel(MouseEvent event) {
//    mouseCamera.mouseWheel(event);
//}
