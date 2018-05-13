
class Vec3 {
  float x, y, z;

  Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  Vec3() {
    this.x = 0.0;
    this.y = 0.0;
    this.z = 0.0;
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

  Vec3 aadd(Vec3 v) {
    Vec3 ans = new Vec3();
    ans.x = this.x + v.x;
    ans.y = this.y + v.y;
    ans.z = this.z + v.z;
    return ans;
  }

  void add(float x, float y, float z) {
    this.x += x;
    this.y += y;
    this.z += z;
  }

  Vec3 aadd(float x, float y, float z) {
    Vec3 ans = new Vec3();
    ans.x = this.x + x;
    ans.y = this.y + y;
    ans.z = this.z + z;
    return ans;
  }

  void subtract(Vec3 v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
  }

  Vec3 ssubtract(Vec3 v) {
    Vec3 ans = new Vec3();
    ans.x = this.x - v.x;
    ans.y = this.y - v.y;
    ans.z = this.z - v.z;
    return ans;
  }

  void multiply(float a) {
    x *= a;
    y *= a;
    z *= a;
  }

  Vec3 mmultiply(float a) {
    Vec3 ans = new Vec3();
    ans.x = this.x * a;
    ans.y = this.y * a;
    ans.z = this.z * a;
    return ans;
  }

  void divide(float a) {
    x /= a;
    y /= a;
    z /= a;
  }

  Vec3 ddivide(float a) {
    Vec3 ans = new Vec3();
    ans.x = this.x / a;
    ans.y = this.y / a;
    ans.z = this.z / a;
    return ans;
  }

  void normalize() {
    float amp = sqrt(x*x+y*y+z*z);
    x /= amp;
    y /= amp;
    z /= amp;
  }

  Vec3 crossProduct(Vec3 b) {
    Vec3 axb = new Vec3();

    axb.x = this.y*b.z - this.z*b.y;
    axb.y = this.z*b.x - this.x*b.z;
    axb.z = this.x*b.y - this.y*b.x;

    return axb;
  }

  float distance(Vec3 b)
  {
    return dist(this.x, this.y, this.z,
                   b.x,    b.y,    b.z);
  }
}
