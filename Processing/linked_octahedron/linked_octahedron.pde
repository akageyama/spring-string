/*

  linked_octahedron.pde


 Developed
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.10


 */

import peasy.PeasyCam;

PeasyCam cam;

//final int N_TRIANGLES = 60;
//final int N_TRIANGLES = 50;
final int N_TRIANGLES = 30;
final int N_PARTICLES = N_TRIANGLES*3;
final float ROPE_MASS = 0.1;
final float PARTICLE_MASS = ROPE_MASS / N_PARTICLES;
//final float SPRING_CHAR_PERIOD = 0.00001; // second
final float SPRING_CHAR_PERIOD = 0.0001; // second
//final float SPRING_CHAR_PERIOD = 0.001; // second
final float SPRING_CHAR_OMEGA = TWO_PI / SPRING_CHAR_PERIOD;
final float SPRING_CHAR_OMEGA_SQ = SPRING_CHAR_OMEGA*SPRING_CHAR_OMEGA;
final float SPRING_CONST = PARTICLE_MASS * SPRING_CHAR_OMEGA_SQ;

final float ROPE_LENGTH = 5.0;
final float TRIANGLE_NATURAL_SEPARATION = ROPE_LENGTH / (N_TRIANGLES-1);
final float EDGE_LENGTH = TRIANGLE_NATURAL_SEPARATION * sqrt(3.0/2.0);
final float ROPE_RADIUS = EDGE_LENGTH / sqrt(3.0);

final float SPRING_CUT_LIMIT_LENGTH = 1.2*EDGE_LENGTH;
// final float SPRING_CUT_LIMIT_LENGTH = 1.3*EDGE_LENGTH;

float time = 0.0;
int step = 0;
final float DT_REF = SPRING_CHAR_PERIOD*0.01;
//final float DT_REF = SPRING_CHAR_PERIOD*0.1;
float dt = DT_REF;
int displayTimeSkip = 16;

boolean frictionFlag = true;
final float FRICTION_COEFF = 0.1;

final float SOUND_SPEED = EDGE_LENGTH / SPRING_CHAR_PERIOD;
final float SOUND_WAVE_TURN_OVER_TIME = 2*ROPE_LENGTH / SOUND_SPEED;
//final float EDGE_TWIST_TIME = SOUND_WAVE_TURN_OVER_TIME * 50;
 final float EDGE_TWIST_TIME = SOUND_WAVE_TURN_OVER_TIME * 20;
//
final float EDGE_TWIST_RATE_OMEGA = TWO_PI / EDGE_TWIST_TIME;
//final float EDGE_TWIST_RATE_OMEGA = 0.0;

//final float GRAVITY_ACCELERATION = 9.80665;
final float GRAVITY_ACCELERATION = 0.0;

float x_coord_min = -3.0;
float x_coord_max =  3.0;
float y_coord_min = x_coord_min;
float y_coord_max = x_coord_max;
float z_coord_min = x_coord_min;
float z_coord_max = x_coord_max;


Particles particles = new Particles();

Springs springs = new Springs();

Motion motion = new Motion();



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
} //<>//

float mapz(float z) {
//  z = min(z,z_coord_max);
//  z =s max(z,z_coord_min);
  return -norma(z); //<>// //<>//
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
  cam = new PeasyCam(this, 1000);
}


void integrate()
{
  motion.rungeKutta();

  step += 1;
}

void draw() {

    for (int i=0; i<displayTimeSkip; i++) {
//    for (int i=0; i<100; i++) {
//    for (int i=0; i<1; i++) {
      integrate();
      particles.resetDt();
    }

    if ( step%100 == 0 ) {
      println("step=", step, " time=", time,
              " friction=", frictionFlag,
              " energy=", motion.totalEnergy(),
              " skip=", displayTimeSkip);
    }

    background(255);
    pushMatrix();
//      translate(width/2,height/2);
//      rotateX(rotor.rotx);

      draw_axes_xyz();
      motion.display();
    popMatrix();

}


void keyPressed() {
  switch (key) {
  case 'f':
    frictionFlag = !frictionFlag;
    break;
  case 'a':
    displayTimeSkip *= 2;
    break;
  case 'd':
    if ( displayTimeSkip > 1 )
      displayTimeSkip /= 2;
    break;
  }
}

//void keyReleased() {
//  switch (key) {
//  }
//}
