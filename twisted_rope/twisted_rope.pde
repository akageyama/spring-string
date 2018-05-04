/* 

  twisted_rope.pde
 
 
 Developed 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.02 
 
 
 Usage:
 - Type u/d to speed up/down.
 - Type s to toggle start/stop.
 - Or mouse click to toggle start/stop. 
 
 */


//MouseCamera mouseCamera;

float time = 0.0;
int step = 0;

final int N_TRIANGLES = 6;
final float EDGE_LENGTH = 0.3;
final float MASS = 0.1;
final float SPRING_CHAR_PERIOD = 0.1; // second

// float dt = 0.0001;
float dt = SPRING_CHAR_PERIOD*0.1;

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
  
  void add(Vec3 v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
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


void rungeKuttaAdvance(int num, float[] p, float[] p1, float[] dp, float factor)
{
  for (int j=0; j<num; j++) {
    p[j] = p1[j] + factor*dp[j];
  }  
}


class Spring 
{
  final float SPRING_NATURAL_LENGTH = EDGE_LENGTH;
  final float SPRING_CHAR_OMEGA = PI*2 / SPRING_CHAR_PERIOD;
  final float SPRING_CHAR_OMEGA_SQ = SPRING_CHAR_OMEGA*SPRING_CHAR_OMEGA;
                // omega^2 = k/m
  final float SPRING_CONST = MASS * SPRING_CHAR_OMEGA_SQ;

  Spring() {
  }
  
  float getConst() {
    return SPRING_CONST;
  }
  
  Vec3 force(float xSelf, float ySelf, float zSelf,
             float xOther, float yOther, float zOther) {
    
    float distance = dist(xSelf,  ySelf,  zSelf,
                          xOther, yOther, zOther);

    float forceAmp = SPRING_CONST*(distance 
                                   - SPRING_NATURAL_LENGTH);

    float unitVectX = ( xOther - xSelf ) / distance;                                 
    float unitVectY = ( yOther - ySelf ) / distance;                                 
    float unitVectZ = ( zOther - zSelf ) / distance;

    float fx = forceAmp*unitVectX;
    float fy = forceAmp*unitVectY;
    float fz = forceAmp*unitVectZ;
          
    Vec3 force = new Vec3(fx,fy,fz);

    return force;                        
  }
}

Spring spring = new Spring();


class ElasticString {
  final int N_PARTICLES = N_TRIANGLES*3;
  Vec3[] pos = new Vec3[N_PARTICLES];
  Vec3[] vel = new Vec3[N_PARTICLES];

  //        
  //     i=5                       i=4         
  //       o x x x x x x x x x x x o
  //         x          i=2      x
  //           x         o     x
  //             x     .  .  x
  //               x .     x
  //               . x   x  .
  //             .     o     .
  //           .      i=3     .
  //         .                 .
  //       o . . . . . . . . . .o
  //      i=0                  i=1
  //
  
  ElasticString() {
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
        pos[3*n+0] = new Vec3(V0x,V0y,V0z); 
        pos[3*n+1] = new Vec3(V1x,V1y,V1z); 
        pos[3*n+2] = new Vec3(V2x,V2y,V2z); 
      }
      else {
        pos[3*n+0] = new Vec3(V3x,V3y,V3z); 
        pos[3*n+1] = new Vec3(V4x,V4y,V4z); 
        pos[3*n+2] = new Vec3(V5x,V5y,V5z); 
      }
      for (int i=0; i<3; i++) { // shift in z-direction.
        pos[3*n+i].z += 2*C3*(n/2);
      }
    }
    
    for (int i=0; i<N_PARTICLES; i++) {
      vel[i] = new Vec3(0.0, 0.0, 0.0);
    }
  }
  
  void drawBalls() {
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
  
  void drawSticksElement(Vec3 a, Vec3 b) {
      float ax = mapx(a.x);
      float ay = mapy(a.y);
      float az = mapz(a.z);
      float bx = mapx(b.x);
      float by = mapy(b.y);
      float bz = mapz(b.z);
      line(ax,ay,az,bx,by,bz);
  }


  void equationOfMotion(float q[], float dq[], float dt) 
  {
    for (int n=1; n<N_TRIANGLES-1; n++) { //  skip the ends.   
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
      
      for (int me=0; me<3; me++) {
        //
        // Connection table
        //
        // when n=even
        //
        //       same layer    upper     lower triangle
        //            /  \     /   \     /   \
        //   me  |  m1   m2   u1   u2   l1   l2
        //   ----+------------------------------
        //    0  |   1    2    0    2    0    2
        //    1  |   2    0    1    0    1    0
        //    2  |   0    1    2    1    2    1
        //       +------------------------------
        //       |  k1   k2   me   k2   me   k2
        //
        //
        // when n=odd
        //
        //       same layer    upper     lower triangle
        //            /  \     /   \     /   \
        //   me  |  m1   m2   u1   u2   l1   l2
        //   ----+------------------------------
        //    0  |   1    2    0    1    0    1
        //    1  |   2    0    1    2    1    2
        //    2  |   0    1    2    0    2    0
        //       +------------------------------
        //       |  k1   k2   me   k1   me   k1
            
        int k1 = (me+1) % 3;
        int k2 = (me+2) % 3;
        int myindex = 3*n + me;
        int[] connectIndex = new int[6];
        connectIndex[0] = 3*n + k1;  // same layer
        connectIndex[1] = 3*n + k2;
        if ( n%2==0 ) {
          connectIndex[2] = 3*(n+1) + me; // upper layer
          connectIndex[3] = 3*(n+1) + k2;
          connectIndex[4] = 3*(n-1) + me; // lower layer
          connectIndex[5] = 3*(n-1) + k2;
        }
        else {
          connectIndex[2] = 3*(n+1) + me; // upper layer
          connectIndex[3] = 3*(n+1) + k1;
          connectIndex[4] = 3*(n-1) + me; // lower layer
          connectIndex[5] = 3*(n-1) + k1;
        }

        Vec3 force = new Vec3(0.0, 0.0, 0.0); // spring force
        float xSelf = q[6*myindex+0];
        float ySelf = q[6*myindex+1];
        float zSelf = q[6*myindex+2];
        for (int j=0; j<6; j++) {
          int indexOther = connectIndex[j];
          float xOther = q[6*indexOther+0];
          float yOther = q[6*indexOther+1];
          float zOther = q[6*indexOther+2];
          Vec3 f = spring.force(xSelf,  ySelf,  zSelf,
                                xOther, yOther, zOther);
          force.add(f);
        }
        
        //float frictionCoeff = 0.0001;
        //float v_force_x = -frictionCoeff*q[4*i+2];
        //float v_force_y = -frictionCoeff*q[4*i+3];
        
        //float force_x = s_force12x - s_force01x + v_force_x;
        //float force_y = s_force12y - s_force01y + v_force_y + g_force_y;
    
        int i = myindex;
        dq[6*i+0] = ( q[6*i+3] ) * dt; // dx = vx * dt
        dq[6*i+1] = ( q[6*i+4] ) * dt; // dy = vy * dt
        dq[6*i+2] = ( q[6*i+5] ) * dt; // dz = vz * dt
        dq[6*i+3] = ( force.x ) / MASS * dt; // dvx = (fx/m)*dt 
        dq[6*i+4] = ( force.y ) / MASS * dt; // dvy = (fy/m)*dt 
        dq[6*i+5] = ( force.z ) / MASS * dt; // dvz = (fz/m)*dt 
      }
    }
  }
  

  void drawSticks() {
    stroke(0, 200, 50);

    for (int n=0; n<N_TRIANGLES; n++) {
      Vec3 v0 = pos[3*n+0];
      Vec3 v1 = pos[3*n+1];
      Vec3 v2 = pos[3*n+2];
      drawSticksElement(v0,v1);
      drawSticksElement(v1,v2);
      drawSticksElement(v2,v0);
    }
    
    for (int n=0; n<N_TRIANGLES-1; n++) {
      Vec3 v0 = pos[3*n+0];
      Vec3 v1 = pos[3*n+1];
      Vec3 v2 = pos[3*n+2];
      Vec3 v3 = pos[3*n+3]; // upper layer
      Vec3 v4 = pos[3*n+4];
      Vec3 v5 = pos[3*n+5];
      if ( n%2==0 ) {
        //        
        //     i=5                      i=4           
        //       o x x x x x x x x x x x o
        //         x          i=2      x
        //           x         o     x
        //             x     .  .  x
        //               x .     x
        //               . x   x  .
        //             .     o     .
        //           .      i=3     .
        //         .                 .
        //       o . . . . . . . . . .o
        //      i=0                  i=1
        drawSticksElement(v3,v0);
        drawSticksElement(v3,v1);
        drawSticksElement(v4,v1);
        drawSticksElement(v4,v2);      
        drawSticksElement(v5,v0);
        drawSticksElement(v5,v2);
      }
      else {
        //        
        //     i=2                      i=1           
        //       o x x x x x x x x x x x o
        //         x          i=5      x
        //           x         o     x
        //             x     .  .  x
        //               x .     x
        //               . x   x  .
        //             .     o     .
        //           .      i=0     .
        //         .                 .
        //       o . . . . . . . . . .o
        //      i=3                  i=4
        //   
        drawSticksElement(v0,v3);
        drawSticksElement(v0,v4);
        drawSticksElement(v1,v4);
        drawSticksElement(v1,v5);      
        drawSticksElement(v2,v3);
        drawSticksElement(v2,v5);
      }
    }
  }
  

  void rungeKutta()
  {
    final float ONE_SIXTH = 1.0/6.0;
    final float ONE_THIRD = 1.0/3.0;
    final int NN = 6*N_PARTICLES;  // pos.x,y,z and vel.x,y,z.
  
    float[] qprev = new float[NN];
    float[] qwork = new float[NN];
    float[] dq1 = new float[NN];
    float[] dq2 = new float[NN];
    float[] dq3 = new float[NN];
    float[] dq4 = new float[NN];
      
    for (int n=0; n<N_PARTICLES; n++) {
      qprev[6*n+0] = pos[n].x;
      qprev[6*n+1] = pos[n].y;
      qprev[6*n+2] = pos[n].z;
      qprev[6*n+3] = vel[n].x;
      qprev[6*n+4] = vel[n].y;
      qprev[6*n+5] = vel[n].z;
    }
  
    //step 1
    equationOfMotion(qprev, dq1, dt);
    rungeKuttaAdvance(NN, qwork, qprev, dq1, 0.5);
  
    //step 2
    equationOfMotion(qwork, dq2, dt);
    rungeKuttaAdvance(NN, qwork, qprev, dq2, 0.5);
  
    //step 3
    equationOfMotion(qwork, dq3, dt);
    rungeKuttaAdvance(NN, qwork, qprev, dq3, 1.0);
  
    //step 4
    equationOfMotion(qwork, dq4, dt);

  
    //the result
    for (int n=1; n<N_PARTICLES-1; n++) { 
      // See boundaryCondition() for end points.
      for (int i=0; i<6; i++) { // x,y,z,vx,vy,vz
        float newval = qprev[6*n+i] + (
                                  ONE_SIXTH*dq1[6*n+i]
                                + ONE_THIRD*dq2[6*n+i]
                                + ONE_THIRD*dq3[6*n+i]
                                + ONE_SIXTH*dq4[6*n+i]
                                );
        if (i==0)
          pos[n].x = newval;
        else if (i==1)
          pos[n].y = newval;
        else if (i==2)
          pos[n].z = newval;
        else if (i==3)
          vel[n].x = newval;
        else if (i==4)
          vel[n].y = newval;
        else if (i==5)
          vel[n].z = newval; 
      } 
    }
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
  //    sum_kinetic += 0.5*MASS*(velx*velx
  //                            +vely*vely
  //                            +velz*velz);
  //  }
  //  float kinetic = 0.5*MASS*(
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
//  z = max(z,z_coord_min);
  return -norma(z);
}



void draw_axes_xyz() {
        stroke(100, 100, 100);
        line(mapx(x_coord_min), 0, 0, mapx(x_coord_max), 0, 0);
        line(0, mapy(y_coord_min), 0, 0, mapy(y_coord_max), 0);
        line(0, 0, mapz(z_coord_min), 0, 0, mapz(z_coord_max));
}



void place_a_sphere(float x, float y, float z) {
    stroke(150, 100, 255);
    fill(0, 200,200);
    pushMatrix();
        translate(mapx(x), mapy(y), mapz(z));      
        sphere(2);
    popMatrix();
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

//    float m1 = part[0].mass;
//    float m2 = part[1].mass;
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
      elasticString.drawBalls();
      elasticString.drawSticks();
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
