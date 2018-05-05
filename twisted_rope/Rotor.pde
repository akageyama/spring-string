
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
