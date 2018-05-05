/* 

  twisted_rope.pde
 
 
 Developed 
 - by Akira Kageyama (kage@port.kobe-u.ac.jp)
 - on 2018.05.05
 
 
 */


//MouseCamera mouseCamera;

float time = 0.0;
int step = 0;

final int N_TRIANGLES = 20;
final int N_PARTICLES = N_TRIANGLES*3;
final float EDGE_LENGTH = 0.3;
final float PARTICLE_MASS = 0.1;
final float SPRING_CHAR_PERIOD = 0.1; // second

float dt = SPRING_CHAR_PERIOD*0.001;

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
  return -norma(z);
}



void draw_axes_xyz() { //<>// //<>// //<>//
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





void shoot() {
  
    for (int n=0; n<20; n++) { // to speed up the display
      elasticString.rungeKutta();
      //boundaryCondition();
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
