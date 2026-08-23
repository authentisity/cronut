#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

#define rad(x)                                                                 \
  (0.0174532925199432954743716805978692718781530857086181640625 * x)

#define CELL_RATIO 2 // terminal characters ~twice as tall as wide

#define W 80
#define H 24
#define FOV 72

#define delta 1e-4 // precision
#define MAX_STEPS 96

float f(float x, float y, float z) {
  // Define shape
  const float R1 = 3.0;
  const float R2 = 5.0;
  float inner = sqrt(x * x + y * y) - R1 * R1;
  return inner * inner + z * z - R2 * R2;
}

float grad(float *p) {
  float gx =
      (f(p[0] + delta, p[1], p[2]) - f(p[0] - delta, p[1], p[2])) / (2 * delta);
  float gy =
      (f(p[0], p[1] + delta, p[2]) - f(p[0], p[1] - delta, p[2])) / (2 * delta);
  float gz =
      (f(p[0], p[1], p[2] + delta) - f(p[0], p[1], p[2] - delta)) / (2 * delta);
  float g = sqrt(gx * gx + gy * gy + gz * gz);
  return g + 1e-9;
}

float sphere_trace(float *p, float *n) {
  const float thres = 1e-3;
  const float d_max = 100.0;

  float p1[3];
  memcpy(p1, p, 3 * sizeof(float));

  float d = 0.0;

  for (int i = 0; i < MAX_STEPS; i++) {
    float t = f(p1[0], p1[1], p1[2]);

    printf("iteration %i t: %f ", i + 1, t);

    if (t < thres) {
      break;
    }

    float c = t / grad(p1);

    p1[0] += n[0] * c;
    p1[1] += n[1] * c;
    p1[2] += n[2] * c;

    d += c;

    printf("c: %f dist: %f new point: (%f, %f, %f)\n", c, d, p1[0], p1[1],
           p1[2]);
  }

  return 1.0;
}

void ray(float i, float j, float *arr) {
  float scale = tan(rad(FOV / 2));
  float u = (2 * i + 1) / W - 1;
  float v = 1 - (2 * j + 1) / W;

  // normalize
  float m = sqrt(u * u + v * v + 1.0);

  arr[0] = u / m;
  arr[1] = v / m;
  arr[2] = 1 / m;
}

int main() {
  float ray1[3];
  ray(0, 0, ray1);
  for (int i = 0; i < 3; i++) {
    printf("%f ", ray1[i]);
  }
  printf("%f\n", grad((float[]){2.0, 0.0, 0.0}));
  printf("%f\n", grad(ray1));

  sphere_trace((float[]){0, 15.0, 5.0}, (float[]){0, -1, 0});
}
