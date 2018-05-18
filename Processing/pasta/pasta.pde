/*

  pasta.pde


 Developed
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.16


 */


import peasy.PeasyCam;

PeasyCam cam;


final int N_TRIANGLES = 40;
final int N_PARTICLES = N_TRIANGLES*3;
final float STICK_MASS = 1.00;
final float PARTICLE_MASS = STICK_MASS / N_PARTICLES;
// final float SPRING_CHAR_PERIOD = 0.000000001; // second
final float SPRING_CHAR_PERIOD = 0.001; // second

final float STICK_LENGTH = 3.0;
final float TRIANGLE_NATURAL_SEPARATION = STICK_LENGTH / (N_TRIANGLES-1);
final float EDGE_LENGTH = TRIANGLE_NATURAL_SEPARATION * sqrt(3.0/2.0);

final float STICK_RADIUS = EDGE_LENGTH / sqrt(3.0);

final float EDGE_ELONGATION_CUT_LIMIT = EDGE_LENGTH*1.04;
final float EDGE_CONTRACTION_CUT_LIMIT = EDGE_LENGTH*0.96;

final float SOUND_SPEED = EDGE_LENGTH / SPRING_CHAR_PERIOD;
//final float STICK_END_POINT_MOVE_SPEED = SOUND_SPEED * 0.0001;
final float STICK_END_POINT_MOVE_SPEED = SOUND_SPEED * 0.001;

int drawTimeSkip = 128;

float time = 0.0;
int step = 0;
float dt = SPRING_CHAR_PERIOD*0.01;

boolean frictionFlag = true;
//final float FRICTION_COEFF = 0.1;
final float FRICTION_COEFF = 0.2;

float x_coord_min = -1.5;
float x_coord_max =  1.5;
float y_coord_min = x_coord_min;
float y_coord_max = x_coord_max;
float z_coord_min = x_coord_min;
float z_coord_max = x_coord_max;


Particles particles = new Particles();

Springs springs = new Springs(SPRING_CHAR_PERIOD);

Motion motion = new Motion();



float norma(float x) {
  float s = width / (x_coord_max-x_coord_min);
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
  return norma(y);
}

float mapz(float z) {
//  z = min(z,z_coord_max);
//  z =s max(z,z_coord_min);
  return norma(z); //<>// //<>//
}



void draw_axes_xyz() { //<>// //<>//
  stroke(100, 100, 100);
  line(mapx(x_coord_min), 0, 0, mapx(x_coord_max), 0, 0);
  line(0, mapy(y_coord_min), 0, 0, mapy(y_coord_max), 0);
  line(0, 0, mapz(z_coord_min), 0, 0, mapz(z_coord_max));
}




void setup() {
  size(800,800,P3D);
  background(255);
  frameRate(60);
  cam = new PeasyCam(this, 1000.0);
}



void integrate()
{
  motion.rungeKutta();
  step += 1;
}


void draw() {

    for (int i=0; i<drawTimeSkip; i++) {
      integrate();
    }

    if ( step%50 == 0 ) {
      println("step=", step, " time=", time,
              " friction=", frictionFlag,
              " timeSkip=", drawTimeSkip,
              " energy=", motion.totalEnergy());
    }

    background(255);
    pushMatrix();
      //translate(width/2,height/2);
      rotateX(PI/2);
      draw_axes_xyz();
      motion.display();
    popMatrix();

}


void keyPressed() {
  switch (key) {
  case 'a':
    drawTimeSkip *= 2;
    break;
  case 'd':
    drawTimeSkip /= 2;
    if ( drawTimeSkip<0 ) drawTimeSkip = 1;
    break;
  case 'f':
    frictionFlag = !frictionFlag;
    break;
  }
}


//void keyReleased() {
//  switch (key) {
//  }
//}
