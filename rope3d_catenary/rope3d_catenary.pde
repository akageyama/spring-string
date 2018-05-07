/*

  twisted_rope.pde


 Developed 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.06


 */


final int N_TRIANGLES = 60;
final int N_PARTICLES = N_TRIANGLES*3;
final float PARTICLE_MASS = 0.1;
final float SPRING_CHAR_PERIOD = 0.05; // second

final float ROPE_LENGTH = 5.0;
final float TRIANGLE_NATURAL_SEPARATION = ROPE_LENGTH / (N_TRIANGLES-1);
final float EDGE_LENGTH = TRIANGLE_NATURAL_SEPARATION * sqrt(3.0/2.0);

float time = 0.0;
int step = 0;
float dt = SPRING_CHAR_PERIOD*0.01;

boolean frictionFlag = true;
final float FRICTION_COEFF = 5.0;
boolean speedLimitFlag = false;

final float GRAVITY_ACCELERATION = 9.80665;

float x_coord_min = -3.0;
float x_coord_max =  3.0;
float y_coord_min = x_coord_min;
float y_coord_max = x_coord_max;
float z_coord_min = x_coord_min;
float z_coord_max = x_coord_max;


Rotor rotor = new Rotor(0,0,0); 

Particles particles = new Particles();

Springs springs = new Springs(SPRING_CHAR_PERIOD);

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
  return -norma(z); //<>//
}



void draw_axes_xyz() { //<>//
        stroke(100, 100, 100);
        line(mapx(x_coord_min), 0, 0, mapx(x_coord_max), 0, 0);
        line(0, mapy(y_coord_min), 0, 0, mapy(y_coord_max), 0);
        line(0, 0, mapz(z_coord_min), 0, 0, mapz(z_coord_max));
}




void setup() {
  size(800,700,P3D);
  background(255);
  frameRate(60);
}


void integrate()
{
  elasticString.rungeKutta();
  if ( speedLimitFlag ) particles.speedLimit();
  
  step += 1;
}


void draw() {

    rotor.update();

    for (int i=0; i<50; i++) {
      integrate();
    }

    if ( step%100 == 0 ) {
      println("step=", step, " time=", time,
              " friction=", frictionFlag,
              " speed limit=", speedLimitFlag,
              " energy=", elasticString.totalEnergy());
    }

    background(255);
    pushMatrix();
      translate(width/2,height/2);
      rotateZ(rotor.rotz);
      rotateX(rotor.rotx);
      rotateY(rotor.roty);

      draw_axes_xyz();
      elasticString.display();
    popMatrix();

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
  case 'f':
    frictionFlag = !frictionFlag;
    break;
  case 's':
    speedLimitFlag = !speedLimitFlag;
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
