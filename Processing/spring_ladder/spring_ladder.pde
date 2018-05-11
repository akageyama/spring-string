/*

  spring_ladder.pde


 Developed
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.07


 */


import peasy.PeasyCam;

PeasyCam cam;


final int N_PAIRS = 30;
final int N_PARTICLES = N_PAIRS*2;
final float LADDER_MASS = 1.0;
final float PARTICLE_MASS = LADDER_MASS / N_PARTICLES;
final float SPRING_CHAR_PERIOD = 0.01; // second

final float LADDER_LENGTH = 1.5;
final float SPRING_NATURAL_LENGTH = LADDER_LENGTH / (N_PAIRS-1);

final float SPRING_ELONGATION_CUT_LIMIT = SPRING_NATURAL_LENGTH*1.2;

final float SOUND_SPEED = SPRING_NATURAL_LENGTH / SPRING_CHAR_PERIOD;
final float SOUND_WAVE_TURN_OVER_TIME = LADDER_LENGTH / SOUND_SPEED;
final float LADDER_TWIST_TIME = SOUND_WAVE_TURN_OVER_TIME * 1;
final float LADDER_TWIST_RATE_OMEGA = TWO_PI / LADDER_TWIST_TIME;

float time = 0.0;
int step = 0;
float dt = SPRING_CHAR_PERIOD*0.01;

boolean frictionFlag = false;
final float FRICTION_COEFF = 0.001;

final float GRAVITY_ACCELERATION = 9.80665;

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
  return norma(z); //<>//
}



void draw_axes_xyz() { //<>//
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

    for (int i=0; i<100; i++) {
      integrate();
    }

    if ( step%100 == 0 ) {
      println("step=", step, " time=", time,
              " friction=", frictionFlag,
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
  case 'f':
    frictionFlag = !frictionFlag;
    break;
  }
}


//void keyReleased() {
//  switch (key) {
//  }
//}
