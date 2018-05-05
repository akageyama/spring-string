
class Vec3 {
  float x, y, z;
  
  Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  
  Vec3(Vec3 rhs) {
    this.x = rhs.x;
    this.y = rhs.y;
    this.z = rhs.z;
  }
  
  void add(Vec3 v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
  }
  
    
  void subtract(Vec3 v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
  }
  
  void multiply(float a) {
    this.x *= a;
    this.y *= a;
    this.z *= a;
  }
  
  void divide(float a) {
    this.x /= a;
    this.y /= a;
    this.z /= a;
  }
}
