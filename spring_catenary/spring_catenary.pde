/*
  spring_catenary.pde
 
 
 Developed started 
 - on 2018.05.01 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 
 Usage:  Start/Stop toggle by 
 - mouse click, or
 - keyboard type of "s"
 
 */

boolean RunningStateToggle = true;

final float ROPE_MASS = 0.01;
final float ROPE_LENGTH = 1.0;
final int BALLS_NUM = 50;
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
final float SPRING_CHAR_PERIOD = 0.06; // second
final float SPRING_CHAR_OMEGA = PI*2 / SPRING_CHAR_PERIOD;
final float SPRING_CHAR_OMEGA_SQ = SPRING_CHAR_OMEGA*SPRING_CHAR_OMEGA;
      // omega^2 = k/m
final float SPRING_CONST = BALLS_MASS * SPRING_CHAR_OMEGA_SQ; 
// to keep the characteristic time scale O(1).
// final float SPRING_CONST = 100.0; 

final float GRAVITY_ACCELERATION = 9.80665;  

float[] ballsCoord = new float[BALLS_NUM*4]; // (x,y,vx,vy)

float xmin = -1.0;
float xmax =  1.0;
float ymin = -1.8;
float ymax =  0.2;


float time = 0.0;
int step = 0;
float dt = SPRING_CHAR_PERIOD*0.05;



void initialize()
{
  float footPointSeparation = ROPE_LENGTH * 0.8;
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
  size(500, 500);
  background(255);
  initialize();
  frameRate(60);
}




float totalEnergy()
{
  
  float kineticEnergy = 0.0;
  float potentialEnergy = 0.0;
  
  for (int i=0; i<BALLS_NUM; i++) {
    float posx = ballsCoord[4*i+0];
    float posy = ballsCoord[4*i+1];
    float velx = ballsCoord[4*i+2];
    float vely = ballsCoord[4*i+3];
  
    kineticEnergy += 0.5*BALLS_MASS*(velx*velx+vely*vely);
    
    if ( i>0 ) {
      float posx0 = ballsCoord[4*(i-1)+0];
      float posy0 = ballsCoord[4*(i-1)+1];      
      float l = dist(posx,posy,posx0,posy0) - SPRING_NATURAL_LENGTH;
      float lsq = l*l;
      potentialEnergy += 0.5*SPRING_CONST*lsq; 
    }
    potentialEnergy += BALLS_MASS*GRAVITY_ACCELERATION*posy;
  }

  return(kineticEnergy + potentialEnergy);
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

    float s_forceAmp01 = SPRING_CONST*(dist01-l0);
    float s_forceAmp12 = SPRING_CONST*(dist12-l0);
    
    float unitVec01x = (x1-x0)/dist01;
    float unitVec01y = (y1-y0)/dist01;
    float unitVec12x = (x2-x1)/dist12;
    float unitVec12y = (y2-y1)/dist12;
    float s_force01x = s_forceAmp01*unitVec01x;
    float s_force01y = s_forceAmp01*unitVec01y;
    float s_force12x = s_forceAmp12*unitVec12x;
    float s_force12y = s_forceAmp12*unitVec12y;
    float  g_force_y = - BALLS_MASS*GRAVITY_ACCELERATION;
    
    float v_force_x = -0.001*q[4*i+2];
    float v_force_y = -0.001*q[4*i+3];
    
    float force_x = s_force12x - s_force01x + v_force_x;
    float force_y = s_force12y - s_force01y + v_force_y + g_force_y;

    dq[4*i+0] = ( q[4*i+2] ) * dt; // dx = vx * dt
    dq[4*i+1] = ( q[4*i+3] ) * dt; // dy = vy * dt
    dq[4*i+2] = ( force_x ) / BALLS_MASS * dt; // dvx = (fx/m)*dt 
    dq[4*i+3] = ( force_y ) / BALLS_MASS * dt;
                                 // dvy = (fy/m- g)*dt
  }
}





void rungeKutta4()
{
  final float ONE_SIXTH = 1.0/6.0;
  final float ONE_THIRD = 1.0/3.0;
  final int NB = BALLS_NUM;
  final int NB4 = NB*4;

  float[] qprev = new float[NB4];
  float[] qwork = new float[NB4];
  float[] dq1 = new float[NB4];
  float[] dq2 = new float[NB4];
  float[] dq3 = new float[NB4];
  float[] dq4 = new float[NB4];

  for (int j=0; j<NB4; j++) {
    qprev[j] = ballsCoord[j];
  }

  //step 1
  equationOfMotion(qprev, dq1, dt);
  rungeKutta4Advance(qwork, qprev, dq1, 0.5);

  //step 2
  equationOfMotion(qwork, dq2, dt);
  rungeKutta4Advance(qwork, qprev, dq2, 0.5);

  //step 3
  equationOfMotion(qwork, dq3, dt);
  rungeKutta4Advance(qwork, qprev, dq3, 1.0);

  //step 4
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

  ballsCoord[(NB-1)*4+0] += 0.001*sin(SPRING_CHAR_OMEGA*time);
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
  float scale = height/(ymax-ymin);
  return map(y, ymin, ymax, scale*ymin, scale*ymax);
}



void drawText() {
  fill(0, 0, 0);
  scale(1, -1);
  text("Click to\nstart/stop.", -mapx(xmax*0.98), -mapx(xmax*0.1));
}


void drawHorizontalLine() {
  stroke(200,0,200);
  line(mapx(xmin), mapy(0), mapx(xmax), mapy(0));
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

  for (int i=0; i<nb; i++) {
    float x = ballsCoord[4*i+0];
    float y = ballsCoord[4*i+1];
    ellipse(mapx(x), mapy(y), 5, 5);
  }
}


void draw() {
  background(255);
  stroke(0, 0, 255);

  translate(width/2, height*ymax/(ymax-ymin));
  scale(1, -1);

  drawRope();
  drawHorizontalLine();

  if ( keyPressed ) {
    if ( key == 's' ) {
      RunningStateToggle = !RunningStateToggle;
    }
  }
  if ( RunningStateToggle ) {
    rungeKutta4();
    time += dt;
    step += 1;
    if ( step%10 == 0 ) {
      println("step = ", step, " time = ", time, " energy = ", totalEnergy());
    }
  }
  drawText();
}

void mousePressed() {
  RunningStateToggle = !RunningStateToggle;
}
