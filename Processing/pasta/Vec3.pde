
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

  Vec3(float f) {
    this.x = f;
    this.y = f;
    this.z = f;
  }

  void add(Vec3 v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
  }

  void add(float x, float y, float z) {
    this.x += x;
    this.y += y;
    this.z += z;
  }

  void subtract(Vec3 v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
  }

  Vec3 ssubtract(Vec3 v) {
    Vec3 ans = new Vec3(0.0);
    ans.x = x - v.x;
    ans.y = y - v.y;
    ans.z = z - v.z;
    return ans;
  }

  void multiply(float a) {
    this.x *= a;
    this.y *= a;
    this.z *= a;
  }

  Vec3 mmultiply(float a) {
    Vec3 ans = new Vec3(0.0);
    ans.x = a*x;
    ans.y = a*y;
    ans.z = a*z;
    return ans;
  }

  void divide(float a) {
    this.x /= a;
    this.y /= a;
    this.z /= a;
  }
}
