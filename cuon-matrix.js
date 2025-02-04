///
// A more complete Matrix4 class in JavaScript, suitable for 3D transforms,
// camera operations, and projection matrices.
//
// Usage Example:
//   import { Matrix4 } from './cuon-matrix.js';
//   const M = new Matrix4();
//   M.setIdentity();
//   M.translate(1,2,3);
//   M.rotate(45, 0,1,0);
//   M.scale(2,2,2);
//   M.setPerspective(60, aspect, 0.1, 100);
//   M.lookAt(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ);
//   gl.uniformMatrix4fv(u_MvpMatrix, false, M.elements);
//

export class Matrix4 {
  constructor() {
    // We'll store the 4x4 matrix in a Float32Array, in column-major order
    this.elements = new Float32Array(16);
    this.setIdentity();
  }

  /**
   * Reset to identity matrix.
   */
  setIdentity() {
    const e = this.elements;
    e[0] = 1;  e[4] = 0;  e[8]  = 0;  e[12] = 0;
    e[1] = 0;  e[5] = 1;  e[9]  = 0;  e[13] = 0;
    e[2] = 0;  e[6] = 0;  e[10] = 1;  e[14] = 0;
    e[3] = 0;  e[7] = 0;  e[11] = 0;  e[15] = 1;
    return this;
  }

  /**
   * Copy from another Matrix4.
   */
  set(src) {
    const s = src.elements;
    const d = this.elements;
    for (let i = 0; i < 16; i++) {
      d[i] = s[i];
    }
    return this;
  }

  /**
   * Multiply this matrix by another matrix. (this = this * other)
   */
  multiply(other) {
    return this.concat(other);
  }

  /**
   * Concatenate (multiply) this matrix with another (on the right).
   * (this = this * other)
   */
  concat(other) {
    const a = this.elements;
    const b = other.elements;
    const e = new Float32Array(16);

    for (let i = 0; i < 4; i++) {
      const ai0 = a[i],   ai1 = a[i+4],  ai2 = a[i+8],   ai3 = a[i+12];
      e[i     ] = ai0*b[0]  + ai1*b[1]  + ai2*b[2]  + ai3*b[3];
      e[i + 4 ] = ai0*b[4]  + ai1*b[5]  + ai2*b[6]  + ai3*b[7];
      e[i + 8 ] = ai0*b[8]  + ai1*b[9]  + ai2*b[10] + ai3*b[11];
      e[i + 12] = ai0*b[12] + ai1*b[13] + ai2*b[14] + ai3*b[15];
    }

    this.elements = e;
    return this;
  }

  /**
   * Set this = matA * matB.
   */
  setConcat(matA, matB) {
    const a = matA.elements;
    const b = matB.elements;
    const e = this.elements;

    for (let i = 0; i < 4; i++) {
      const ai0 = a[i], ai1 = a[i+4], ai2 = a[i+8], ai3 = a[i+12];
      e[i     ] = ai0*b[0]  + ai1*b[1]  + ai2*b[2]  + ai3*b[3];
      e[i + 4 ] = ai0*b[4]  + ai1*b[5]  + ai2*b[6]  + ai3*b[7];
      e[i + 8 ] = ai0*b[8]  + ai1*b[9]  + ai2*b[10] + ai3*b[11];
      e[i + 12] = ai0*b[12] + ai1*b[13] + ai2*b[14] + ai3*b[15];
    }
    return this;
  }

  /**
   * Translate this matrix by (x, y, z).
   */
  translate(x, y, z) {
    const e = this.elements;
    e[12] += e[0]*x + e[4]*y + e[8]*z;
    e[13] += e[1]*x + e[5]*y + e[9]*z;
    e[14] += e[2]*x + e[6]*y + e[10]*z;
    e[15] += e[3]*x + e[7]*y + e[11]*z;
    return this;
  }

  /**
   * Rotate this matrix by `angle` (degrees) around (x, y, z).
   */
  rotate(angle, x, y, z) {
    const rad = Math.PI * angle / 180;
    const s = Math.sin(rad);
    const c = Math.cos(rad);

    let len = Math.sqrt(x*x + y*y + z*z);
    if (len === 0) return this;
    x /= len; y /= len; z /= len;

    const nc = 1 - c;
    const xy = x * y * nc;
    const yz = y * z * nc;
    const zx = z * x * nc;
    const xs = x * s;
    const ys = y * s;
    const zs = z * s;

    const rot = new Float32Array(16);
    rot[0]  = x*x*nc + c;   rot[4]  = xy - zs;      rot[8]  = zx + ys;      rot[12] = 0;
    rot[1]  = xy + zs;      rot[5]  = y*y*nc + c;   rot[9]  = yz - xs;      rot[13] = 0;
    rot[2]  = zx - ys;      rot[6]  = yz + xs;      rot[10] = z*z*nc + c;   rot[14] = 0;
    rot[3]  = 0;            rot[7]  = 0;            rot[11] = 0;            rot[15] = 1;

    return this.concat({ elements: rot });
  }

  /**
   * Scale this matrix by (sx, sy, sz).
   */
  scale(sx, sy, sz) {
    const e = this.elements;
    e[0] *= sx; e[4] *= sy; e[8]  *= sz;
    e[1] *= sx; e[5] *= sy; e[9]  *= sz;
    e[2] *= sx; e[6] *= sy; e[10] *= sz;
    e[3] *= sx; e[7] *= sy; e[11] *= sz;
    return this;
  }

  /**
   * Return a new Matrix4 that is a clone of this one.
   */
  clone() {
    const newMat = new Matrix4();
    newMat.set(this);
    return newMat;
  }

  /**
   * Transpose this matrix in-place.
   */
  transpose() {
    const e = this.elements;
    let t;
    t = e[1];  e[1]  = e[4];  e[4]  = t;
    t = e[2];  e[2]  = e[8];  e[8]  = t;
    t = e[3];  e[3]  = e[12]; e[12] = t;
    t = e[6];  e[6]  = e[9];  e[9]  = t;
    t = e[7];  e[7]  = e[13]; e[13] = t;
    t = e[11]; e[11] = e[14]; e[14] = t;
    return this;
  }

  /**
   * Invert this matrix in-place. Return this or null if not invertible.
   */
  invert() {
    const e = this.elements;
    const inv = new Float32Array(16);

    inv[0]  = e[5]*e[10]*e[15] - e[5]*e[11]*e[14] - e[9]*e[6]*e[15]
              + e[9]*e[7]*e[14] + e[13]*e[6]*e[11] - e[13]*e[7]*e[10];
    inv[4]  = -e[4]*e[10]*e[15] + e[4]*e[11]*e[14] + e[8]*e[6]*e[15]
              - e[8]*e[7]*e[14] - e[12]*e[6]*e[11] + e[12]*e[7]*e[10];
    inv[8]  = e[4]*e[9]*e[15] - e[4]*e[11]*e[13] - e[8]*e[5]*e[15]
              + e[8]*e[7]*e[13] + e[12]*e[5]*e[11] - e[12]*e[7]*e[9];
    inv[12] = -e[4]*e[9]*e[14] + e[4]*e[10]*e[13] + e[8]*e[5]*e[14]
              - e[8]*e[6]*e[13] - e[12]*e[5]*e[10] + e[12]*e[6]*e[9];

    inv[1]  = -e[1]*e[10]*e[15] + e[1]*e[11]*e[14] + e[9]*e[2]*e[15]
              - e[9]*e[3]*e[14] - e[13]*e[2]*e[11] + e[13]*e[3]*e[10];
    inv[5]  = e[0]*e[10]*e[15] - e[0]*e[11]*e[14] - e[8]*e[2]*e[15]
              + e[8]*e[3]*e[14] + e[12]*e[2]*e[11] - e[12]*e[3]*e[10];
    inv[9]  = -e[0]*e[9]*e[15] + e[0]*e[11]*e[13] + e[8]*e[1]*e[15]
              - e[8]*e[3]*e[13] - e[12]*e[1]*e[11] + e[12]*e[3]*e[9];
    inv[13] = e[0]*e[9]*e[14] - e[0]*e[10]*e[13] - e[8]*e[1]*e[14]
              + e[8]*e[2]*e[13] + e[12]*e[1]*e[10] - e[12]*e[2]*e[9];

    inv[2]  = e[1]*e[6]*e[15] - e[1]*e[7]*e[14] - e[5]*e[2]*e[15]
              + e[5]*e[3]*e[14] + e[13]*e[2]*e[7] - e[13]*e[3]*e[6];
    inv[6]  = -e[0]*e[6]*e[15] + e[0]*e[7]*e[14] + e[4]*e[2]*e[15]
              - e[4]*e[3]*e[14] - e[12]*e[2]*e[7] + e[12]*e[3]*e[6];
    inv[10] = e[0]*e[5]*e[15] - e[0]*e[7]*e[13] - e[4]*e[1]*e[15]
              + e[4]*e[3]*e[13] + e[12]*e[1]*e[7] - e[12]*e[3]*e[5];
    inv[14] = -e[0]*e[5]*e[14] + e[0]*e[6]*e[13] + e[4]*e[1]*e[14]
              - e[4]*e[2]*e[13] - e[12]*e[1]*e[6] + e[12]*e[2]*e[5];

    inv[3]  = -e[1]*e[6]*e[11] + e[1]*e[7]*e[10] + e[5]*e[2]*e[11]
              - e[5]*e[3]*e[10] - e[9]*e[2]*e[7]  + e[9]*e[3]*e[6];
    inv[7]  = e[0]*e[6]*e[11]  - e[0]*e[7]*e[10] - e[4]*e[2]*e[11]
              + e[4]*e[3]*e[10] + e[8]*e[2]*e[7]  - e[8]*e[3]*e[6];
    inv[11] = -e[0]*e[5]*e[11] + e[0]*e[7]*e[9]  + e[4]*e[1]*e[11]
              - e[4]*e[3]*e[9]  - e[8]*e[1]*e[7]  + e[8]*e[3]*e[5];
    inv[15] = e[0]*e[5]*e[10]  - e[0]*e[6]*e[9]  - e[4]*e[1]*e[10]
              + e[4]*e[2]*e[9]  + e[8]*e[1]*e[6]  - e[8]*e[2]*e[5];

    let det = e[0]*inv[0] + e[1]*inv[4] + e[2]*inv[8] + e[3]*inv[12];
    if (det === 0) {
      // not invertible
      return null;
    }
    det = 1 / det;
    for (let i = 0; i < 16; i++) {
      inv[i] *= det;
    }
    this.elements = inv;
    return this;
  }

  //=================================================
  //  CAMERA & PROJECTION
  //=================================================

  /**
   * Set this matrix to a perspective projection.
   * fovy: vertical field-of-view in degrees
   * aspect: width/height
   * near, far: distance to clipping planes
   */
  setPerspective(fovy, aspect, near, far) {
    this.setIdentity();
    const e = this.elements;

    const rad = Math.PI * fovy / 180.0;
    const tanHalfFovy = Math.tan(rad / 2);
    const nf = 1 / (near - far);

    e[0]  = 1 / (aspect * tanHalfFovy);
    e[5]  = 1 / tanHalfFovy;
    e[10] = (far + near) * nf;
    e[11] = -1;
    e[14] = (2 * far * near) * nf;
    e[15] = 0;
    return this;
  }

  /**
   * Multiply this matrix by a viewing matrix derived from eye, center, up.
   * This = this * lookAt(...)
   */
  lookAt(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ) {
    // f (forward)
    let fx = centerX - eyeX;
    let fy = centerY - eyeY;
    let fz = centerZ - eyeZ;
    // normalize f
    let rlf = 1/Math.sqrt(fx*fx + fy*fy + fz*fz);
    fx *= rlf;  fy *= rlf;  fz *= rlf;

    // s = f x up
    let sx = fy*upZ - fz*upY;
    let sy = fz*upX - fx*upZ;
    let sz = fx*upY - fy*upX;
    // normalize s
    let rls = 1/Math.sqrt(sx*sx + sy*sy + sz*sz);
    sx *= rls;  sy *= rls;  sz *= rls;

    // u = s x f
    let ux = sx*fy - sy*fx;
    let uy = sy*fz - sz*fy;
    let uz = sz*fx - sx*fz;

    const m = new Float32Array(16);
    m[0] = sx;   m[4] = sy;   m[8]  = sz;   m[12] = 0;
    m[1] = ux;   m[5] = uy;   m[9]  = uz;   m[13] = 0;
    m[2] = -fx;  m[6] = -fy;  m[10] = -fz;  m[14] = 0;
    m[3] = 0;    m[7] = 0;    m[11] = 0;    m[15] = 1;

    const viewMat = new Matrix4();
    viewMat.elements = m;
    // then translate(-eyeX, -eyeY, -eyeZ)
    viewMat.translate(-eyeX, -eyeY, -eyeZ);

    return this.concat(viewMat);
  }

  /**
   * Set this matrix to an orthographic projection.
   * left, right, bottom, top, near, far
   */
  ortho(left, right, bottom, top, near, far) {
    this.setIdentity();
    const e = this.elements;
    e[0] = 2 / (right - left);
    e[5] = 2 / (top - bottom);
    e[10] = -2 / (far - near);
    e[12] = - (right + left) / (right - left);
    e[13] = - (top + bottom) / (top - bottom);
    e[14] = - (far + near)   / (far - near);
    return this;
  }
}
