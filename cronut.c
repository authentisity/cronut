#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

#define CELL_RATIO 2 // terminal characters ~twice as tall as wide
#define SCALE 30

#define W 79
#define H 23
#define FOV_H 2.0f * atanf((float)W / SCALE)
#define FOV_V 2.0f * atanf((float)(H * CELL_RATIO) / SCALE)

#define delta 1e-4f // precision
#define MAX_STEPS 96

#define CAMERA_Z 3

// convergence tolerance should be much smaller than character cell size
#define PIXEL_W (2.0f * CAMERA_Z / SCALE)
#define THRES (PIXEL_W / 256.0f)

void rot(float *p, float *p1, float x, float y) {
  float a = sinf(x);
  float b = cosf(x);
  float c = sinf(y);
  float d = cosf(y);

  p1[0] = p[0] * d + p[1] * a * c + p[2] * b * c;
  p1[1] = p[1] * b - p[2] * a;
  p1[2] = -p[0] * c + p[1] * a * d + p[2] * b * d;
}

float f(float x, float y, float z) {
  // lemniscate
  const float R = 0.35f;
  float s = x * x + y * y;
  float comp = (s * s) - 8 * (x * x - y * y);
  // adds thickness to thinnest part of curve
  float cx = 4 * x * (s - 4);
  float cy = 4 * y * (s + 4);
  float g = comp / (sqrtf(cx * cx + cy * cy) + 1e-6f);
  return g * g + z * z - R * R;

  // donut
  const float R1 = 3.0;
  const float R2 = 5.0;
  float inner = sqrt(x * x + y * y) - R1 * R1;
  return inner * inner + z * z - R2 * R2;
}

void grad(float *p, float *g) {
  g[0] =
      (f(p[0] + delta, p[1], p[2]) - f(p[0] - delta, p[1], p[2])) / (2 * delta);
  g[1] =
      (f(p[0], p[1] + delta, p[2]) - f(p[0], p[1] - delta, p[2])) / (2 * delta);
  g[2] =
      (f(p[0], p[1], p[2] + delta) - f(p[0], p[1], p[2] - delta)) / (2 * delta);
}

float mag(float *p) {
  float m = sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  return m + 1e-9f;
}

static float f0;

bool sphere_trace(float *p, float *p1, float *n) {
  const float thres = THRES; // world units, not units of f
  const float d_max = 100.0f;

  memcpy(p1, p, 3 * sizeof(float));

  float d = 0.0f;

  for (int i = 0; i < MAX_STEPS; i++) {
    if (d > d_max) {
      return false;
    }

    float t = f(p1[0], p1[1], p1[2]);

    float c = t / ({
                float g[3];
                grad(p1, g);
                mag(g);
              });

    if (fabsf(c) < thres && fabsf(t) <= f0) {
      return true;
    }

    p1[0] += n[0] * c;
    p1[1] += n[1] * c;
    p1[2] += n[2] * c;

    d += c;
  }

  return false;
}

float lambertian(float *p, float *nl) {
  float g[3];
  grad(p, g);
  float m = mag(g);
  float ng[3] = {g[0] / m, g[1] / m, g[2] / m};

  float Id = max(0.0f, ng[0] * nl[0] + ng[1] * nl[1] + ng[2] * nl[2]);

  return Id;
}

void ray(float i, float j, float *arr) {
  float u = ((2 * i + 1) / W - 1) * tanf(FOV_H / 2);
  float v = (1 - (2 * j + 1) / H) * tanf(FOV_V / 2);

  // normalize
  float m = sqrtf(u * u + v * v + 1.0f);

  arr[0] = u / m;
  arr[1] = v / m;
  arr[2] = 1 / m;
}

int main() {
  float rays[W * H][3];
  float rays_rot[W * H][3];

  float camera_rot[3];
  float light_rot[3];

  char b[(W + 1) * H];
  float p1[3];

  for (int i = 0; i < W; i++) {
    for (int j = 0; j < H; j++) {
      ray(i, j, rays[j * W + i]);
    }
  }

  for (int j = 0; j < H; j++) {
    b[(j + 1) * (W + 1) - 1] = '\n';
  }

  f0 = fabsf(f(0, 0, -CAMERA_Z));

  float alpha = 0.0f;
  float beta = 0.0f;

  while (1) {
    rot((float[]){0, 0, -CAMERA_Z}, camera_rot, alpha, beta);
    rot((float[]){0, 1, -1}, light_rot, alpha, beta);

    float lm = mag(light_rot);
    light_rot[0] /= lm;
    light_rot[1] /= lm;
    light_rot[2] /= lm;

    for (int i = 0; i < W; i++) {
      for (int j = 0; j < H; j++) {
        int idx = j * W + i;

        rot(rays[idx], rays_rot[idx], alpha, beta);

        if (sphere_trace(camera_rot, p1, rays_rot[idx])) {
          b[idx + j] = ".,-~:;=!*#$@"[(int)(lambertian(p1, light_rot) * 11)];
        } else {
          b[idx + j] = ' ';
        }
      }
    }

    fwrite(b, 1, (W + 1) * H, stdout);
    usleep(30000);

    alpha = fmodf(alpha + 0.07f, 6.28f);
    beta = fmodf(beta + 0.02f, 6.28f);

    fputs("\x1b[23A", stdout);
  }
}
