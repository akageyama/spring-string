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
float dt = 0.001;

float x_coord_min = -3.0;
float x_coord_max =  3.0;
float y_coord_min = x_coord_min;
float y_coord_max = x_coord_max;
float z_coord_min = x_coord_min;
float z_coord_max = x_coord_max;

float vmax =  1.0;
float vmin = -vmax;

int speed = 10000;



class Rotor {
    float rotx = 0;
    float roty = 0;
    float rotz = 0;
    float delta = PI/500;
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


class GeneralCoords {
  float x;
  float vx;

  GeneralCoords(float x, float vx) {
    this.x = x;
    this.vx = vx;
  }

  GeneralCoords() {
    x = 0;
    vx = 0;
  }

  GeneralCoords(GeneralCoords copy) {
    x = copy.x;
    vx = copy.vx;
  }
}


class Vec3 {
  float x, y, z;
  Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
}

class ElasticString {
  final int N_TRIANGLES = 2;
  final float EDGE_LENGTH = 0.3;
  final float MASS = 0.01;
  final int N_PARTICLES = N_TRIANGLES*3;
  Vec3[] pos = new Vec3[N_PARTICLES];
  //        
  //     i=4                      i=5           
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
    final float V0x = -EDGE_LENGTH/2;
    final float V0y = -EDGE_LENGTH/(2*sqrt(3));
    final float V0z =  0;
    final float V1x =  EDGE_LENGTH/2;
    final float V1y = -EDGE_LENGTH/(2*sqrt(3));
    final float V1z =  0;
    final float V2x =  0;
    final float V2y =  EDGE_LENGTH/sqrt(3);
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
    final float V3x =  0;
    final float V3y = -EDGE_LENGTH/sqrt(3);
    final float V3z =  0.5;
    final float V4x =  EDGE_LENGTH/2;
    final float V4y =  EDGE_LENGTH/(2*sqrt(3));
    final float V4z =  0.5;
    final float V5x = -EDGE_LENGTH/2;
    final float V5y =  EDGE_LENGTH/(2*sqrt(3));
    final float V5z =  0.5;
    
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
    }
  }
  
  void placeBalls() {
    noStroke();
    fill(100,0,130);
    for (int i=0; i<N_PARTICLES; i++) {
      pushMatrix();
        float x = pos[i].x;
        float y = pos[i].y;
        float z = pos[i].z;
        translate(mapx(x), mapy(y), mapz(z));      
        sphere(10);
      popMatrix();
    }
  }
}

ElasticString elasticString = new ElasticString();

class Particle {
    float mass;
    GeneralCoords pos;
    
    Particle() {
        mass = STANDARD_PARTICLE_MASS;
        pos = new GeneralCoords();
    }
    
    Particle(Particle copy) {
        mass = copy.mass;
        pos = new GeneralCoords(copy.pos);
    }
    
    Particle(float mass) {
      this.mass = mass;
      pos =  new GeneralCoords();
    }
}

Particle[] part;
Particle[] part_prev;



final float SPRING_L0 = 0.0;
//final float SPRING_L0 = 1.0;  // x_coord_max * ;
final float SPRING_K  = 1.0;
final float STANDARD_PARTICLE_MASS = 1.0;





class DataSaver3 {
    float[] v0;
    float[] v1;
    float[] v2;
    int counter = 0;
    int size = 100;
    
    DataSaver3(int size_) {
        size = size_;
        counter = 0;
        v0 = new float[size];
        v1 = new float[size];
        v2 = new float[size];
    }
  
    void save(float v0_, float v1_, float v2_) {
        for (int i=size-1; i>0; i--) {
              v0[i] = v0[i-1];
              v1[i] = v1[i-1];
              v2[i] = v2[i-1];
        }
        v0[0] = v0_;
        v1[0] = v1_;
        v2[0] = v2_;
        if ( counter < size ) {
            counter++;
        } 
    }
    
    void reset() {
        counter = 0;
    }
}


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


DataSaver3 ds3 = new DataSaver3(100000);

void place_a_sphere(float x, float y, float z) {
    stroke(150, 100, 255);
    fill(0, 200,200);
    pushMatrix();
        translate(mapx(x), mapy(y), mapz(z));      
        sphere(2);
    popMatrix();
}


void place_balls() {
    float amp = 4;
    for (int i=0; i<ds3.counter; i++) {
        float red   = 255.0 * (amp*ds3.v0[i]-x_coord_min) / (x_coord_max-x_coord_min);
        float green = 255.0 * (amp*ds3.v1[i]-y_coord_min) / (y_coord_max-y_coord_min);
        float blue  = 255.0 * (amp*ds3.v2[i]-z_coord_min) / (z_coord_max-z_coord_min);
        float agefactor = 1-float(i)/ds3.counter;

        stroke(red*agefactor, green*agefactor, blue*agefactor);
            
        //pushMatrix();
        //    translate(mapx(ds3.v0[i]),mapy(ds3.v1[i]),mapz(ds3.v2[i]));
        //    sphere(5);
        //popMatrix();
        
        point(mapx(ds3.v0[i]),mapy(ds3.v1[i]),mapz(ds3.v2[i]));
    }
}


void place_a_sphere_with_tail() {
    for (int i=1; i<ds3.counter; i++) {
        float red   = 255.0 * (ds3.v0[i]-x_coord_min) / (x_coord_max-x_coord_min);
        float green = 255.0 * (ds3.v1[i]-y_coord_min) / (y_coord_max-y_coord_min);
        float blue  = 255.0 * (ds3.v2[i]-z_coord_min) / (z_coord_max-z_coord_min);
        float agefactor = 1-float(i)/ds3.counter;

        stroke(red*agefactor, green*agefactor, blue*agefactor);
            
        line(mapx(ds3.v0[i-1]),mapy(ds3.v1[i-1]),mapz(ds3.v2[i-1]),
             mapx(ds3.v0[i]  ),mapy(ds3.v1[i]  ),mapz(ds3.v2[i]  ));
    }
    stroke(150, 100, 255);
    fill(0, 200,200);
    pushMatrix();        
        translate(mapx(ds3.v0[0]),
                  mapy(ds3.v1[0]),    
                  mapy(ds3.v2[0]));     
        sphere(2);
    popMatrix();
}
    

//  void draw_orbit_in_phase_space(float x, float v) {
//      ps.save(x,v);
//      pushMatrix();
//          translate_origin();
//          for (int i=1; i<ps.counter; i++) {
//              float x0 = ps.v0[i-1];
//              float v0 = ps.v1[i-1];
//              float x1 = ps.v0[i];
//              float v1 = ps.v1[i];         
//              stroke((255.0*i)/ps.counter, 
//                     (255.0*i)/ps.counter,
//                     (255.0*i)/ps.counter);
//              line(mapx(x0), mapv(v0), mapx(x1), mapv(v1));
//          }
//          noStroke();
//          fill(255, 100, 100);
//          ellipse(mapx(x), mapv(v), 10, 10);
//      popMatrix();
//  }


  


void initialize_particles() {
    float x1 = x_coord_max*0.4; // chaos
    float x2 = x_coord_min*0.1; // chaos
//    float x1 = x_coord_max*0.4; // non-chaos
//    float x2 = x_coord_min*0.3; // non-chaos
    part[0].pos = new GeneralCoords(x1, 0);
    part[1].pos = new GeneralCoords(x2, 0);
    
    //for (int i=0; i<2; i++) {
    //    //float x  = random(x_coord_min*0.5, x_coord_max*0.5);
//  //      float vx = random(vmin*0.5, vmax*0.5);
    //    float vx = 0;
    //    part[i].pos = new GeneralCoords(x, vx);  
    //}    
    
    for (int i=0; i<2; i++) {
  //    part_prev[i] = part[i];   // This was a bug. Dangerous! 
        part_prev[i].mass   = part[i].mass;
        part_prev[i].pos.x  = part[i].pos.x;
        part_prev[i].pos.vx = part[i].pos.vx;
    }     
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


  part = new Particle[2];
  part_prev = new Particle[2];

  for (int i=0; i<2; i++) {  
    part[i] = new Particle();   
    part_prev[i] = new Particle();
    part[i].pos = new GeneralCoords();
    part_prev[i].pos = new GeneralCoords();
  }
  
  initialize_particles();


  for (int i=0; i<2; i++) {
//    part_prev[i] = part[i];   // This was a bug. Dangrous! 
      part_prev[i].mass   = part[i].mass;
      part_prev[i].pos.x  = part[i].pos.x;
      part_prev[i].pos.vx = part[i].pos.vx;
  }
}





float parabola_func_upper(float x) {
  // When you change this, revise its derivative
  // parabola_func_upper_derivative(), too.
  float y;
  y = x*x+1;
  return(y);
}

float parabola_func_upper_derivative(float x) {
  // When you change this, revise
  // parabola_func_upper(), too.
  float y;
  y = 2*x;
  return(y);
}


float parabola_func_lower(float x) {
  // When you change this, revise its derivative
  // parabola_func_lower_derivative(), too.
  float y;
  y = -x*x-1;
  return y;
}

float parabola_func_lower_derivative(float x) {
  // When you change this, revise
  // parabola_func_upper(), too.
  float y = -2*x;
  return y;
}



float total_energy()
{
    float x1 = part[0].pos.x;
    float v1 = part[0].pos.vx;
    float x2 = part[1].pos.x;
    float v2 = part[1].pos.vx;

    float v1sq = v1*v1;
    float v2sq = v2*v2;

    float y1 = parabola_func_upper(x1);
    float y2 = parabola_func_lower(x2);
    
    float y1dot = parabola_func_upper_derivative(x1)*v1;
    float y2dot = parabola_func_lower_derivative(x2)*v2;
    float y1dotsq = y1dot*y1dot;
    float y2dotsq = y2dot*y2dot;
    
    float s = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

    float m1 = part[0].mass;
    float m2 = part[1].mass;
    float kinetic_e = 0.5*(m1*(v1sq+y1dotsq)
                          +m2*(v2sq+y2dotsq));
    float potential = 0.5*SPRING_K*(s-SPRING_L0)*(s-SPRING_L0);

    return(kinetic_e + potential);
}

void rungekutta_advance(Particle[] p, Particle[] p1, Particle[] dp, float factor) {
  for (int i=0; i<2; i++) {
    p[i].pos.x  = p1[i].pos.x  + factor*dp[i].pos.x;
    p[i].pos.vx = p1[i].pos.vx + factor*dp[i].pos.vx;
  }
}


void equation_of_motion(Particle[] p, Particle[] dp, float dt) {

  //    Lagrangian
  //       L(x1,x2,x1',x2') = (m/2)*(x1'^2+4*x1^2*x1'^2)
  //                        + (m/2)*(x2'^2+4*x2^2*x2'^2)
  //                        - (k/2)*(s-L0)^2
  //    where
  //        s = sqrt(dx^2+dy^2), dx=x1-x2, dy=x1^2+x2^2+2
  // 
    float x1 = p[0].pos.x;
    float v1 = p[0].pos.vx;
    float x2 = p[1].pos.x;
    float v2 = p[1].pos.vx;

    float dx   = x1 - x2;
    float x1sq = x1*x1;
    float v1sq = v1*v1;
    float x2sq = x2*x2;
    float v2sq = v2*v2;
    float dy   = x1sq + x2sq + 2;
    float s    = sqrt(dx*dx + dy*dy);
    
    float m1 = p[0].mass;
    float m2 = p[1].mass;

    float f1 = (SPRING_K/m1)*(s-SPRING_L0)/s*( dx+2*x1*dy);
    float f2 = (SPRING_K/m2)*(s-SPRING_L0)/s*(-dx+2*x2*dy);

    dp[0].pos.x = ( v1 ) * dt;
    dp[0].pos.vx = ( -1.0/(1+4*x1sq)*(4*x1*v1sq + f1) ) * dt;
    dp[1].pos.x = ( v2 ) * dt;
    dp[1].pos.vx = ( -1.0/(1+4*x2sq)*(4*x2*v2sq + f2) ) * dt;
}




void rungeKutta4(Particle[] p, Particle[] prev, int nn)
{
  final float ONE_SIXTH = 1.0/6.0;
  final float ONE_THIRD = 1.0/3.0;
  
  Particle[] work = new Particle[nn];
  Particle[] dp01 = new Particle[nn];
  Particle[] dp02 = new Particle[nn];
  Particle[] dp03 = new Particle[nn];
  Particle[] dp04 = new Particle[nn];

  for (int i=0; i<nn; i++) {
      work[i] = new Particle(p[i].mass);
      dp01[i] = new Particle(p[i].mass);
      dp02[i] = new Particle(p[i].mass);
      dp03[i] = new Particle(p[i].mass);
      dp04[i] = new Particle(p[i].mass);
  
      prev[i].pos.x  = p[i].pos.x;
      prev[i].pos.vx = p[i].pos.vx;
  }

  //step 1
  equation_of_motion(prev, dp01, dt); 
  rungekutta_advance(work, prev, dp01, 0.5);

  //step 2
  equation_of_motion(work, dp02, dt);
  rungekutta_advance(work, prev, dp02, 0.5);

  //step 3
  equation_of_motion(work, dp03, dt);
  rungekutta_advance(work, prev, dp03, 1.0);

  //step 4
  equation_of_motion(work, dp04, dt);

  //the result
  for (int i=0; i<2; i++) {
    p[i].pos.x = prev[i].pos.x 
           + ( ONE_SIXTH*dp01[i].pos.x
             + ONE_THIRD*dp02[i].pos.x
             + ONE_THIRD*dp03[i].pos.x
             + ONE_SIXTH*dp04[i].pos.x 
             );
    p[i].pos.vx = prev[i].pos.vx 
            + ( ONE_SIXTH*dp01[i].pos.vx
              + ONE_THIRD*dp02[i].pos.vx
              + ONE_THIRD*dp03[i].pos.vx
              + ONE_SIXTH*dp04[i].pos.vx
            );
  }
  
}


class Poincare {
      //   |              .
      // va|____________.
      //   |          . |
      //   |        .   |
      // --+-xb---.-----xa------>x
      //   |  | .  \
      // vb|__.     x=xb-vb*(xa-xb)/(va-vb) 
      //   |.       (See below.)
      //
      // The equation of the linear function is  
      //       v(x) = (va-vb)/(xa-xb) * (x-xb) + vb.
      // Solving 
      //       v(x) = 0,
      // We get
      //       x = xb-vb*(xa-xb)/(va-vb)
      //         = xb+weight*(xa-xb)  [weight=-vb/(va-vb)]
      //         = weight*xa + (1-weight)*xb
    float wa, wb;
    boolean is_crossed;
    
    void check(float before, float after) {
        is_crossed = (before*after < 0);
        if (is_crossed) {
            wa = -before/(after-before);
            wb = 1 - wa;
        }
    }
    
    Poincare() {
        wa = 0.0;
        wb = 1.0;
        is_crossed = false;
    }
}


void reset() {
    time = 0.0;
    step = 0;
    initialize_particles();
}

  
Poincare poincare = new Poincare();

void shoot() {
  

    for (int icnt=0; icnt<speed; icnt++) {
        rungeKutta4(part, part_prev, 2);
        time += dt;
        step += 1;
        
        poincare.check(part_prev[1].pos.vx,
                            part[1].pos.vx);
                            
        if ( poincare.is_crossed  ) {
            float wb = poincare.wb;
            float wa = poincare.wa;
            
            float path_x = wa*part_prev[0].pos.x
                         + wb*     part[0].pos.x;
            float path_y = wa*part_prev[1].pos.x
                         + wb*     part[1].pos.x;
            float path_z = wa*part_prev[0].pos.vx
                         + wb*     part[0].pos.vx;
            ds3.save(path_x, path_y, path_z);
        }
    }

    //rotz += 0.01/60.0*log(speed);
    //rotx += 0.02/60.0*log(speed);
    //roty += 0.03/60.0*log(speed);
    //rotx %= 360;
    //roty %= 360;
    //rotz %= 360;

    background(255);
    pushMatrix();
      translate(width/2,height/2);
      rotateZ(rotor.rotz);
      rotateX(rotor.rotx);
      rotateY(rotor.roty);
      
      draw_axes_xyz();
      //place_a_sphere_with_tail();
      place_balls();
      elasticString.placeBalls();
    popMatrix();               
    
    //if ( step%1000 == 0 ) {
    //  println("step = ", step," time = ", time," energy = ",total_energy());
    //}
    

}


void draw() {

    //mouseCamera.update();
    
    rotor.update();
    
    shoot();
    
    println("speed=", nf(speed), " step=", nf(step), 
            "counter=",ds3.counter,
            "energy=",total_energy());

   
        
   // if ( step > 30000000 ) {
   //     reset();
   // }


    String str = "Speed = " + nf(speed) + "  (Type u/d to speed up/down)";
    str += "\nenergy = " + nf(total_energy(), 4, 3);
    //str += "  ang. mom. = " + nf(angular_momentum(), 4, 3);
    //str += "  mom. x = " + nf(momentum_x(), 4, 3);
    //str += "  mom. y = " + nf(momentum_y(), 4, 3);
    str += "\nt = " + nf(time, 6, 3);
    str += " (step = " + nf(step, 9) + ")";
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
  case 'u':
    speed *= 10;
    println("speed=", nf(speed), " step=", nf(step), "counter=",ds3.counter);
    break;
  case 'd':
    speed /= 10;
    if ( speed <= 0 ) speed = 1;
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
