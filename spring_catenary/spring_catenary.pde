/*
  spring_catenary.pde
 
 
 Developed started 
 - on 2018.05.01 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 
 Usage:  Start/Stop toggle by 
 - mouse click, or
 - keyboard type of "s"
 
 */

float time = 0.0;
int step = 0;
float dt = 0.005;

boolean running_state_toggle = true;

final float ROPE_MASS = 1.0;
final float ROPE_LENGTH = 1.0;
final int BALLS_NUM = 10;
//
//         when BALLS_NUM = 6
//         ......o............o........
//              0 \          / 5
//                 o        o
//                1 \      / 4
//                   o----o
//                  2      3 
//
final float BALLS_MASS = ROPE_MASS / (BALLS_NUM-2);

final float SPRING_NATURAL_LENGTH = ROPE_LENGTH / (BALLS_NUM-1);
//final float SPRING_NATURAL_LENGTH = 0.5;
//final float SPRING_CONST = BALLS_MASS; 
//// to keep the characteristic time scale O(1).
final float SPRING_CONST = 100.0; 

final float GRAVITY_ACCELERATION = 9.80665;  

float[] ballsCoord = new float[BALLS_NUM*4]; // (x,y,vx,vy)

float xmin = -3.0;
float xmax =  3.0;


float x = 0.0;

float[] particle_pos = new float[4];

float S0 = 0.0;
//float S0 = 2.9;
float K  = 8.0;
float M  = 1.0; 


/*
double double_sqrt_by_newton(double x) {
 int i = 0;
 double ans, guess = 1.0;
 float diff = 999.999; // any large num.
 ans = guess;
 while (i<20 && diff > 1.e-8) {
 ans = guess - (guess*guess-x)/(2*guess);
 diff = abs((float)(guess-ans));
 guess = ans;
 i++;
 }
 return ans; 
 }
 */


void initialize_particle()
{
  particle_pos[0] =   2.0;   // x1
  particle_pos[1] =   0.0;   // x1' = vx1 = (dx1/dt)
  particle_pos[2] =  -0.5;   // x2 
  particle_pos[3] =   0.0;   // x2' = vx2 = (dx2/dt)
}



void initialize()
{
  float footPointSeparation = ROPE_LENGTH * 0.75;
//
//          LeftX              RightX
//             |                  |
//             |<---separation--->|
//             |         |        |
//       ......o..................o........
//            0 \        |       / BALLS_NUM-1
//               o      x=0     o
//              1 \            / BALLS_NUM-2
//   
  float footPointRightX = footPointSeparation/2;
  float footPointLeftX = -footPointRightX;

//
//           angle
//           theta
//              \
//               \  |
//              \ \_|     separation/2
//               \/ |    /
//                \ |<------>|
//                 \|        |
//            ......o........+........o........
//                  |\               /    
//   natural length | o <--1        o <--NB-2
//              l0.....\           / 
//                  |   o         .            
//                  |    .       .               
//                  |     \     /
//            (NB-2)/2 --> o---o <--(NB-2)/2+1
//
//     (((NB-2)/2)*l0*sin(theta)+0.5*l0) \sim separation/2
//   or
//      sin(theta) \sim ((separation/2)-0.5*l0) / ((NB-2)/2)*l0
//   

  int nb = BALLS_NUM;
  float l0 = SPRING_NATURAL_LENGTH;
  float sinTheta = (footPointSeparation/2-0.5*l0) / ((nb-2)/2*l0);
  float cosTheta = sqrt(1-sinTheta*sinTheta);

  int i;

  i=0;
  ballsCoord[4*i+0] = footPointLeftX; // x coord
  ballsCoord[4*i+1] = 0.0; // y coord
  ballsCoord[4*i+2] = 0.0; // vx
  ballsCoord[4*i+3] = 0.0; // vy

  for (i=1; i<=(nb-2)/2; i++) {
    ballsCoord[4*i+0] = ballsCoord[4*(i-1)+0] + l0*sinTheta; // x
    ballsCoord[4*i+1] = ballsCoord[4*(i-1)+1] - l0*cosTheta; // y
    ballsCoord[4*i+2] = 0.0; // vx
    ballsCoord[4*i+3] = 0.0; // vy
  }

  i = nb-1;
  ballsCoord[4*i+0] = footPointRightX;
  ballsCoord[4*i+1] = 0.0;
  ballsCoord[4*i+2] = 0.0;
  ballsCoord[4*i+3] = 0.0;

  for (i=(nb-2); i>=(nb-2)/2+1; i--) {
    ballsCoord[4*i+0] = ballsCoord[4*(i+1)+0] - l0*sinTheta; // x
    ballsCoord[4*i+1] = ballsCoord[4*(i+1)+1] - l0*cosTheta; // y
    ballsCoord[4*i+2] = 0.0; // vx
    ballsCoord[4*i+3] = 0.0; // vy
  }
}


void setup() {
  size(600, 600);
  background(255);
  initialize_particle();
  initialize();
  frameRate(60);
}


float total_energy()
{
  float x1     = particle_pos[0];
  float x1_dot = particle_pos[1];
  float x2     = particle_pos[2];
  float x2_dot = particle_pos[3];

  float x1_sq     = x1*x1;
  float x1_dot_sq = x1_dot*x1_dot;
  float x2_sq     = x2*x2;
  float x2_dot_sq = x2_dot*x2_dot;

  float v1_sq = (4*x1_sq+1)*x1_dot_sq;
  float v2_sq = (4*x2_sq+1)*x2_dot_sq;

  float y1 =  x1_sq+1;
  float y2 = -x2_sq-1;

  float s = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

  float kinetic_e = 0.5*M*(v1_sq+v2_sq);
  float potential = 0.5*K*(s-S0)*(s-S0);

  return(kinetic_e + potential);
}


void runge_kutta4_advance(float[] p, float[] p1, float[] dp, float factor)
{
  for (int i=0; i<4; i++) {
    p[i] = p1[i] + factor*dp[i];
  }
}


void rungeKutta4Advance(float[] p, float[] p1, float[] dp, float factor)
{
  for (int j=0; j<BALLS_NUM*4; j++) {
    p[j] = p1[j] + factor*dp[j];
  }
}



void equationOfMotion(float q[], float dq[], float dt) 

{
  final float l0 = SPRING_NATURAL_LENGTH;
  final int NB = BALLS_NUM;

//    
//   when BALLS_NUM = 6
//   ......o............o........
//        0 \          / 5
//           o        o
//          1 \      / 4
//             o----o
//            2      3
//   

  for (int k=0; k<4; k++) {
    dq[4*0     +k] = 0.0;  // Ball i=1 is fixed.
    dq[4*(NB-1)+k] = 0.0;  // Ball i=NB-1 is also fixed.
  }

  for (int i=1; i<=NB-2; i++) {
    // 
    //     (x0,y0)          (x1,y1)
    //         i-1           i
    //          o------------o
    //         /            / \  ____dist12
    //        /<--dist01-->/   \/
    //                          \ 
    //                           o i+1
    //                         (x2,y2)
    //
    //    force_amp = k*(spring_length - l0)
    //
    float x0 = q[4*(i-1)+0];
    float x1 = q[4*(i  )+0];
    float x2 = q[4*(i+1)+0];
    float y0 = q[4*(i-1)+1];
    float y1 = q[4*(i  )+1];
    float y2 = q[4*(i+1)+1];
    float dist01 = dist(x0,y0,x1,y1);
    float dist12 = dist(x1,y1,x2,y2);


    float forceAmp01 = SPRING_CONST*(dist01-l0);
    float forceAmp12 = SPRING_CONST*(dist12-l0);
    
    float unitVec01x = (x1-x0)/dist01;
    float unitVec01y = (y1-y0)/dist01;
    float unitVec12x = (x2-x1)/dist12;
    float unitVec12y = (y2-y1)/dist12;
    float force01x = forceAmp01*unitVec01x;
    float force01y = forceAmp01*unitVec01y;
    float force12x = forceAmp12*unitVec12x;
    float force12y = forceAmp12*unitVec12y;

if (i==4) {
    println("i = ", i, "dist01 = ", dist01);
    println("    ", i, "dist12 = ", dist12);
    println("    ", i, "  l0 = ",  l0);
    println("   ", i, "forceAmp01 = ", forceAmp01);
    println("   ", i, "forceAmp12 = ", forceAmp12);
    println("   ", i, "unitVec01 = ", unitVec01x, unitVec01y);
    println("   ", i, "unitVec12 = ", unitVec12x, unitVec12y);
}

    dq[4*i+0] = ( q[4*i+2] ) * dt; // dx = vx * dt
    dq[4*i+1] = ( q[4*i+3] ) * dt; // dy = vy * dt
    dq[4*i+2] = ( force12x - force01x ) / BALLS_MASS * dt; // dvx = (fx/m)*dt 
    dq[4*i+3] = ( force12y - force01y ) / BALLS_MASS * dt
                              - GRAVITY_ACCELERATION * dt; // dvy = (fy/m- g)*dt
  }
}



void equation_of_motion(float pos[], float dpos[], float dt) 
{
  //    Lagrangian
  //       L(x1,x2,x1',x2') = (m/2)*(x1'^2+4*x1^2*x1'^2)
  //                        + (m/2)*(x2'^2+4*x2^2*x2'^2)
  //                        - (k/2)*(s-S0)^2
  //    where
  //        s = sqrt(dx^2+dy^2), dx=x1-x2, dy=x1^2+x2^2+2
  //
  float x1     = pos[0];
  float x1_dot = pos[1];
  float x2     = pos[2];
  float x2_dot = pos[3];

  float dx   = x1 - x2;
  float x1sq = x1*x1, x1_dot_sq = x1_dot*x1_dot;
  float x2sq = x2*x2, x2_dot_sq = x2_dot*x2_dot;
  float dy   = x1sq + x2sq + 2;
  float s    = sqrt(dx*dx + dy*dy);

  float f1 = (K/M)*(s-S0)/s*( dx+2*x1*dy);
  float f2 = (K/M)*(s-S0)/s*(-dx+2*x2*dy);

  dpos[0] = ( x1_dot ) * dt;
  dpos[1] = ( -1.0/(1+4*x1sq)*(4*x1*x1_dot_sq + f1) ) * dt;
  dpos[2] = ( x2_dot ) * dt;
  dpos[3] = ( -1.0/(1+4*x2sq)*(4*x2*x2_dot_sq + f2) ) * dt;
}


void runge_kutta4()
{
  final float ONE_SIXTH = 1.0/6.0;
  final float ONE_THIRD = 1.0/3.0;
  final int NB = BALLS_NUM;
  final int NB4 = NB*4;

  float[] pos_before = new float[4];
  float[] work_pos = new float[4];
  float[] dpos1 = new float[4];
  float[] dpos2 = new float[4];
  float[] dpos3 = new float[4];
  float[] dpos4 = new float[4];

  float[] qprev = new float[NB4];
  float[] qwork = new float[NB4];
  float[] dq1 = new float[NB4];
  float[] dq2 = new float[NB4];
  float[] dq3 = new float[NB4];
  float[] dq4 = new float[NB4];

  for (int j=0; j<NB4; j++) {
    qprev[j] = ballsCoord[j];
  }

  pos_before[0] = particle_pos[0];
  pos_before[1] = particle_pos[1];
  pos_before[2] = particle_pos[2];
  pos_before[3] = particle_pos[3];

  //step 1
  equation_of_motion(pos_before, dpos1, dt);
  runge_kutta4_advance(work_pos, pos_before, dpos1, 0.5);
  equationOfMotion(qprev, dq1, dt);
  rungeKutta4Advance(qwork, qprev, dq1, 0.5);
  //step 2
  equation_of_motion(work_pos, dpos2, dt);
  runge_kutta4_advance(work_pos, pos_before, dpos2, 0.5);
  equationOfMotion(qwork, dq2, dt);
  rungeKutta4Advance(qwork, qprev, dq2, 0.5);

  //step 3
  equation_of_motion(work_pos, dpos3, dt);
  runge_kutta4_advance(work_pos, pos_before, dpos3, 1.0);
  equationOfMotion(qwork, dq3, dt);
  rungeKutta4Advance(qwork, qprev, dq3, 1.0);


  //step 4
  equation_of_motion(work_pos, dpos4, dt);
  equationOfMotion(qwork, dq4, dt);

  //the result

  for (int j=0; j<NB*4; j++) {
    ballsCoord[j] = qprev[j] + (
      ONE_SIXTH*dq1[j]
      + ONE_THIRD*dq2[j]
      + ONE_THIRD*dq3[j]
      + ONE_SIXTH*dq4[j] 
      );
  }

  particle_pos[0] = pos_before[0] + (  
    ONE_SIXTH*dpos1[0]
    + ONE_THIRD*dpos2[0]
    + ONE_THIRD*dpos3[0]
    + ONE_SIXTH*dpos4[0] 
    );
  particle_pos[1] = pos_before[1] + (  
    ONE_SIXTH*dpos1[1]
    + ONE_THIRD*dpos2[1]
    + ONE_THIRD*dpos3[1]
    + ONE_SIXTH*dpos4[1] 
    );
  particle_pos[2] = pos_before[2] + (  
    ONE_SIXTH*dpos1[2]
    + ONE_THIRD*dpos2[2]
    + ONE_THIRD*dpos3[2]
    + ONE_SIXTH*dpos4[2] 
    );
  particle_pos[3] = pos_before[3] + (  
    ONE_SIXTH*dpos1[3]
    + ONE_THIRD*dpos2[3]
    + ONE_THIRD*dpos3[3]
    + ONE_SIXTH*dpos4[3] );
}


float parabola_func_upper(float x) {
  float y;
  y = x*x+1;
  return(y);
}


float parabola_func_lower(float x) {
  float y = -parabola_func_upper(x);
  return y;
}



float mapx(float x) {
  // (x,y) = physical unit coords. 
  // (map(x),map(y)) = pixel coords.
  float scale = width/(xmax-xmin);
  return map(x, xmin, xmax, scale*xmin, scale*xmax);
}


float mapy(float y) {
  // (x,y) = physical unit coords. 
  // (map(x),map(y)) = pixel coords.
  float ymax = parabola_func_upper(xmax);
  float ymin = -ymax;
  float scale = width/(xmax-xmin);
  return map(y, ymin, ymax, scale*ymin, scale*ymax);
}


void draw_parabolas() {
  stroke(100, 100, 100);

  int nx = 500;
  float dx = (xmax-xmin)/nx;
  float x, y;

  float x_prev = xmin;
  float y_prev = parabola_func_upper(x_prev);

  for (int i=1; i<=nx; i++) { // starts from i=1.
    x = xmin + dx*i;
    y = parabola_func_upper(x);
    line(mapx(x_prev), mapy(y_prev), mapx(x), mapy(y));
    x_prev = x;
    y_prev = y;
  }

  x_prev = xmin;
  y_prev = parabola_func_lower(x_prev);

  for (int i=1; i<=nx; i++) { // starts from i=1.
    x = xmin + dx*i;
    y = parabola_func_lower(x);
    line(mapx(x_prev), mapy(y_prev), mapx(x), mapy(y));
    x_prev = x;
    y_prev = y;
  }
}


void draw_text() {
  fill(0, 0, 0);
  scale(1, -1);
  text("Click to\nstart/stop.", -mapx(xmax*0.98), -mapx(xmax*0.1));
}



void draw_balls(float x1, float x2) {
  stroke(255, 100, 100);
  fill(255, 100, 100);

  float y1 = parabola_func_upper(x1);
  ellipse(mapx(x1), mapy(y1), 10, 10);

  stroke(50, 100, 255);
  fill(50, 100, 255);
  float y2 = parabola_func_lower(x2);
  ellipse(mapx(x2), mapy(y2), 10, 10);
}


void draw_spring() {
  stroke(0, 255, 0);
  float x1 = particle_pos[0];
  float y1 = parabola_func_upper(x1);
  float x2 = particle_pos[2];
  float y2 = parabola_func_lower(x2);
  line(mapx(x1), mapy(y1), mapx(x2), mapy(y2));
}

void drawRope() {
  stroke(50, 100, 200);

  int nb = BALLS_NUM;

  for (int i=0; i<nb-1; i++) {
    float x0 = ballsCoord[4*(i  )+0];
    float y0 = ballsCoord[4*(i  )+1];
    float x1 = ballsCoord[4*(i+1)+0];
    float y1 = ballsCoord[4*(i+1)+1];
    line(mapx(x0), mapy(y0), mapx(x1), mapy(y1));
  }

  for (int i=1; i<=nb-2; i++) {
    float x = ballsCoord[4*i+0];
    float y = ballsCoord[4*i+1];
    ellipse(mapx(x), mapy(y), 5, 5);
  }
}


void draw() {
  background(255);
  stroke(0, 0, 255);

  translate(width/2, height/2);
  scale(1, -1);

  draw_parabolas();
  drawRope();

  if ( keyPressed ) {
    if ( key == 's' ) {
      running_state_toggle = !running_state_toggle;
    }
  }
  if ( running_state_toggle ) {
    runge_kutta4();
    time += dt;
    step += 1;
    if ( step%10 == 0 ) {
      println("step = ", step, " time = ", time, " energy = ", total_energy());
    }
  }
  float x1 = particle_pos[0];
  float x2 = particle_pos[2];
  draw_spring();
  draw_balls(x1, x2);
  draw_text();
}

void mousePressed() {
  running_state_toggle = !running_state_toggle;
}
