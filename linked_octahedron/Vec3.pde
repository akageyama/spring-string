
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

  void add(float x, float y, float z) {
    this.x += x;
    this.y += y;
    this.z += z;
  }

  void subtract(Vec3 v) {
    this.x -= x;
    this.y -= y;
    this.z -= z;
  }

  void multiply(float a) {
    x *= a;
    y *= a;
    z *= a;
  }

  void divide(float a) {
    x /= a;
    y /= a;
    z /= a;
  }

  void normalize() {
    float amp = sqrt(x*x+y*y+z*z);
    x /= amp;
    y /= amp;
    z /= amp;
  }

  Vec3 crossProduct(Vec3 b) {
    Vec3 axb = new Vec3(0.0, 0.0, 0.0);

    axb.x = this.y*b.z - this.z*b.y;
    axb.y = this.z*b.x - this.x*b.z;
    axb.z = this.x*b.y - this.y*b.x;

    return axb;
  }
}
