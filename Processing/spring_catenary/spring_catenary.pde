/*
  spring_catenary.pde
 
 
 Developed started 
 - on 2018.05.01 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 
 Usage:  Start/Stop the end point shakes by the mouse click.
 
 */

boolean RunningStateToggle = true;

boolean shakeFlag = true;

final float ROPE_MASS = 0.005;
final float ROPE_LENGTH = 2.0;
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
final float SPRING_CHAR_PERIOD = 0.01; // second
final float SPRING_CHAR_OMEGA = PI*2 / SPRING_CHAR_PERIOD;
final float SPRING_CHAR_OMEGA_SQ = SPRING_CHAR_OMEGA*SPRING_CHAR_OMEGA;
      // omega^2 = k/m
final float SPRING_CONST = BALLS_MASS * SPRING_CHAR_OMEGA_SQ; 

final float GRAVITY_ACCELERATION = 9.80665;  

float[] ballPosX = new float[BALLS_NUM];
float[] ballPosY = new float[BALLS_NUM];
float[] ballVelX = new float[BALLS_NUM];
float[] ballVelY = new float[BALLS_NUM];


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
float footPointSeparation = ROPE_LENGTH * 0.5;
float footPointRightX = footPointSeparation/2;
float footPointLeftX = -footPointRightX;

float xmin = -1.0;
float xmax =  1.0;
float ymin = -1.5;
float ymax =  0.5;


float time = 0.0;
int step = 0;
float dt = SPRING_CHAR_PERIOD*0.05;



void initialize()
{


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
  ballPosX[i] = footPointLeftX; // x coord
  ballPosY[4*i+1] = 0.0; // y coord
  ballVelX[4*i+2] = 0.0; // vx
  ballVelY[4*i+3] = 0.0; // vy

  for (i=1; i<=(nb-2)/2; i++) {
    ballPosX[i] = ballPosX[i-1] + l0*sinTheta; // x
    ballPosY[i] = ballPosY[i-1] - l0*cosTheta; // y
    ballVelX[i] = 0.0; // vx
    ballVelY[i] = 0.0; // vy
  }

  i = nb-1;
  ballPosX[i] = footPointRightX;
  ballPosY[i] = 0.0;
  ballVelX[i] = 0.0;
  ballVelY[i] = 0.0;

  for (i=(nb-2); i>=(nb-2)/2+1; i--) {
    ballPosX[i] = ballPosX[i+1] - l0*sinTheta; // x
    ballPosY[i] = ballPosY[i+1] - l0*cosTheta; // y
    ballVelX[i] = 0.0; // vx
    ballVelY[i] = 0.0; // vy
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
    float posx = ballPosX[i];
    float posy = ballPosY[i];
    float velx = ballVelX[i];
    float vely = ballVelY[i];
  
    kineticEnergy += 0.5*BALLS_MASS*(velx*velx+vely*vely);
    
    if ( i>0 ) {
      float posx0 = ballPosX[i-1];
      float posy0 = ballPosY[i-1];      
      float l = dist(posx,posy,posx0,posy0) - SPRING_NATURAL_LENGTH;
      float lsq = l*l;
      potentialEnergy += 0.5*SPRING_CONST*lsq; 
    }
    potentialEnergy += BALLS_MASS*GRAVITY_ACCELERATION*posy;
  }

  return(kineticEnergy + potentialEnergy);
}


void rungeKutta4Advance(int num, 
                        float[] posx,
                        float[] posy,
                        float[] velx,
                        float[] vely,
                        float[] posx1,
                        float[] posy1,
                        float[] velx1,
                        float[] vely1,
                        float[] dposx,
                        float[] dposy,
                        float[] dvelx,
                        float[] dvely,
                        float factor)
{
  for (int j=0; j<num; j++) {
    posx[j] = posx1[j] + factor*dposx[j];
    posy[j] = posy1[j] + factor*dposy[j];
    velx[j] = velx1[j] + factor*dvelx[j];
    vely[j] = vely1[j] + factor*dvely[j];
  }
}



void equationOfMotion(float  posx[],
                      float  posy[],
                      float  velx[],
                      float  vely[],                      
                      float dposx[],
                      float dposy[],
                      float dvelx[],
                      float dvely[],
                      float dt) 

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

  for (int i=1; i<=NB-2; i++) {  // See boundaryCondition() for i=0 & NB-1.   
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
    float dtm = dt / BALLS_MASS;
    
    float x0 = posx[i-1];
    float x1 = posx[i  ];
    float x2 = posx[i+1];
    float y0 = posy[i-1];
    float y1 = posy[i  ];
    float y2 = posy[i+1];
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
    
    float frictionCoeff = 0.0001;
    float v_force_x = -frictionCoeff*velx[i];
    float v_force_y = -frictionCoeff*vely[i];
    
    float force_x = s_force12x - s_force01x + v_force_x;
    float force_y = s_force12y - s_force01y + v_force_y + g_force_y;

    dposx[i] = velx[i] * dt;  // dx = vx * dt
    dposy[i] = vely[i] * dt;  // dy = vy * dt
    dvelx[i] = force_x * dtm; // dvx = (fx/m)*dt 
    dvely[i] = force_y * dtm; // dvy = (fy/m- g)*dt
  }
}



void rungeKutta4()
{
  final float ONE_SIXTH = 1.0/6.0;
  final float ONE_THIRD = 1.0/3.0;
  final int NB = BALLS_NUM;

  float[] posxprev = new float[NB];
  float[] posxwork = new float[NB];
  float[]   dposx1 = new float[NB];
  float[]   dposx2 = new float[NB];
  float[]   dposx3 = new float[NB];
  float[]   dposx4 = new float[NB];
  float[] posyprev = new float[NB];
  float[] posywork = new float[NB];
  float[]   dposy1 = new float[NB];
  float[]   dposy2 = new float[NB];
  float[]   dposy3 = new float[NB];
  float[]   dposy4 = new float[NB];
  float[] velxprev = new float[NB];
  float[] velxwork = new float[NB];
  float[]   dvelx1 = new float[NB];
  float[]   dvelx2 = new float[NB];
  float[]   dvelx3 = new float[NB];
  float[]   dvelx4 = new float[NB];
  float[] velyprev = new float[NB];
  float[] velywork = new float[NB];
  float[]   dvely1 = new float[NB];
  float[]   dvely2 = new float[NB];
  float[]   dvely3 = new float[NB];
  float[]   dvely4 = new float[NB];

  for (int j=0; j<NB; j++) {
    posxprev[j] = ballPosX[j];
    posyprev[j] = ballPosY[j];
    velxprev[j] = ballVelX[j];
    velyprev[j] = ballVelY[j];
  }

  //step 1 
  equationOfMotion(posxprev,
                   posyprev,
                   velxprev,
                   velyprev,
                     dposx1,
                     dposy1,
                     dvelx1,
                     dvely1,
                         dt);
  rungeKutta4Advance(NB,
                     posxwork,
                     posywork,
                     velxwork,
                     velywork,
                     posxprev,
                     posyprev,
                     velxprev,
                     velyprev,
                       dposx1,
                       dposy1,
                       dvelx1,
                       dvely1,
                          0.5);                        
  boundaryCondition(time, posxwork, posywork);

  time += 0.5*dt;

  //step 2
  equationOfMotion(posxwork,
                   posywork,
                   velxwork,
                   velywork,
                     dposx2,
                     dposy2,
                     dvelx2,
                     dvely2,
                         dt);
  rungeKutta4Advance(NB,
                     posxwork,
                     posywork,
                     velxwork,
                     velywork,
                     posxprev,
                     posyprev,
                     velxprev,
                     velyprev,
                       dposx2,
                       dposy2,
                       dvelx2,
                       dvely2,
                          0.5);
  boundaryCondition(time, posxwork, posywork);
                          
  //step 3
  equationOfMotion(posxwork,
                   posywork,
                   velxwork,
                   velywork,
                     dposx3,
                     dposy3,
                     dvelx3,
                     dvely3,
                         dt);
  rungeKutta4Advance(NB,
                     posxwork,
                     posywork,
                     velxwork,
                     velywork,
                     posxprev,
                     posyprev,
                     velxprev,
                     velyprev,
                       dposx3,
                       dposy3,
                       dvelx3,
                       dvely3,
                          1.0);
  boundaryCondition(time, posxwork, posywork);

  time += 0.5*dt;

  //step 4
  equationOfMotion(posxwork,
                   posywork,
                   velxwork,
                   velywork,
                     dposx4,
                     dposy4,
                     dvelx4,
                     dvely4,
                         dt);
  
  //the result
  for (int j=1; j<NB-1; j++) { 
    posxwork[j] = posxprev[j] + (
                           ONE_SIXTH*dposx1[j]
                         + ONE_THIRD*dposx2[j]
                         + ONE_THIRD*dposx3[j]
                         + ONE_SIXTH*dposx4[j] 
                         );
    posywork[j] = posyprev[j] + (
                           ONE_SIXTH*dposy1[j]
                         + ONE_THIRD*dposy2[j]
                         + ONE_THIRD*dposy3[j]
                         + ONE_SIXTH*dposy4[j] 
                         );
    velxwork[j] = velxprev[j] + (
                           ONE_SIXTH*dvelx1[j]
                         + ONE_THIRD*dvelx2[j]
                         + ONE_THIRD*dvelx3[j]
                         + ONE_SIXTH*dvelx4[j] 
                         );
    velywork[j] = velyprev[j] + (
                           ONE_SIXTH*dvely1[j]
                         + ONE_THIRD*dvely2[j]
                         + ONE_THIRD*dvely3[j]
                         + ONE_SIXTH*dvely4[j] 
                         );
  }
  
  boundaryCondition(time, posxwork, posywork);
  
  for (int j=0; j<NB; j++) {
    ballPosX[j] = posxwork[j];
    ballPosY[j] = posywork[j];
    ballVelX[j] = velxwork[j];
    ballVelY[j] = velywork[j];
  }


}


void boundaryCondition(float t, float[] x, float y[]) 
{
  x[0] = footPointLeftX;   // x-coord if particle No.0.
  y[0] = 0.0;              // y-coord if particle No.0.
  x[BALLS_NUM-1] = footPointRightX;  // x-coord if the last particle.
  y[BALLS_NUM-1] = 0.0;              // y-coord if the last particle.
  if ( shakeFlag ) {
    float amp = SPRING_NATURAL_LENGTH*0.06;
    float omega = 0.80*SPRING_CHAR_OMEGA;
    x[0]           += amp*sin(omega*t);
    x[BALLS_NUM-1] += amp*sin(omega*t);
  }
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
  if (shakeFlag) {
    text("Shaing the end points. Click to\nstart/stop the shake.", 
          -mapx(xmax*0.98), -mapy(ymax*0.4));
  }
  else {
    text("No shaing. Click to\nstart/stop the shake.", 
          -mapx(xmax*0.98), -mapy(ymax*0.4));
  }
}


void drawHorizontalLine() {
  stroke(200,0,200);
  line(mapx(xmin), mapy(0), mapx(xmax), mapy(0));
}



void drawRope() {
  stroke(50, 100, 200);

  int nb = BALLS_NUM;

  for (int i=0; i<nb-1; i++) {
    float x0 = ballPosX[i];
    float y0 = ballPosY[i];
    float x1 = ballPosX[i+1];
    float y1 = ballPosY[i+1];
    line(mapx(x0), mapy(y0), mapx(x1), mapy(y1));
  }

  for (int i=0; i<nb; i++) {
    float x = ballPosX[i];
    float y = ballPosY[i];
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
      // RunningStateToggle = !RunningStateToggle;
      shakeFlag = !shakeFlag;
    }
  }
  if ( RunningStateToggle ) {
    for (int n=0; n<20; n++) { // to speed up the display
      rungeKutta4();
      step += 1;
      if ( step%10 == 0 ) {
        println("step=", step, " time=", time, " energy=", totalEnergy()," shake=",shakeFlag);
      }
    }
  }
  drawText();
}

void mousePressed() {
  // RunningStateToggle = !RunningStateToggle;
  shakeFlag = !shakeFlag;
}
