#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// terminal
#define W 79
#define H 23

// settings
#define SCALE 30
#define CAMERA_Z 3
#define A_SPEED 0.07f
#define B_SPEED 0.02f
#define DELAY 30000
#define LIGHT 0, 1, -1
#define LUM ".,-~:;=!*#$@"

// precision
#define dt 1e-4f
#define MAX_STEPS 96

// convergence tolerance should be much smaller than character cell size
#define THRES (CAMERA_Z / (SCALE * 128.0f))
#define D_MAX 100.0f

float A, B, f0;

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
  float inner = sqrt(x * x + y * y) - R1;
  return inner * inner + z * z - R2 * R2;
}

void rot(float *p, float *q) {
  float a = sinf(A), b = cosf(A), c = sinf(B), d = cosf(B);
  q[0] = p[0] * d + p[1] * a * c + p[2] * b * c;
  q[1] = p[1] * b - p[2] * a;
  q[2] = -p[0] * c + p[1] * a * d + p[2] * b * d;
}

void grad(float *p, float *g) {
  g[0] = (f(p[0] + dt, p[1], p[2]) - f(p[0] - dt, p[1], p[2])) / (2 * dt);
  g[1] = (f(p[0], p[1] + dt, p[2]) - f(p[0], p[1] - dt, p[2])) / (2 * dt);
  g[2] = (f(p[0], p[1], p[2] + dt) - f(p[0], p[1], p[2] - dt)) / (2 * dt);
}

float mag(float *p) {
  return sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) + 1e-9f;
}

char sphere_trace(float *p, float *p1, float *n) {
  memcpy(p1, p, 3 * sizeof(float));

  float d = 0.0f;

  for (int i = 0; i < MAX_STEPS; i++) {
    if (d > D_MAX) {
      return 0;
    }

    float t = f(p1[0], p1[1], p1[2]);

    float c = t / ({
                float g[3];
                grad(p1, g);
                mag(g);
              });

    if (fabsf(c) < THRES && fabsf(t) <= f0) {
      return 1;
    }

    p1[0] += n[0] * c;
    p1[1] += n[1] * c;
    p1[2] += n[2] * c;

    d += c;
  }

  return 1;
}

float lambertian(float *p, float *nl) {
  float g[3];
  grad(p, g);
  float m = mag(g);
  float ng[3] = {g[0] / m, g[1] / m, g[2] / m};

  float Id = fmaxf(0.0f, ng[0] * nl[0] + ng[1] * nl[1] + ng[2] * nl[2]);

  return Id;
}

void ray(float i, float j, float *arr) {
  float u = (2 * i + 1 - W) / SCALE;
  float v = 2 * (H - 2 * j - 1) / SCALE;

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

  for (;;) {
    rot((float[]){0, 0, -CAMERA_Z}, camera_rot);
    rot((float[]){0, 1, -1}, light_rot);

    float lm = mag(light_rot);
    light_rot[0] /= lm;
    light_rot[1] /= lm;
    light_rot[2] /= lm;

    for (int i = 0; i < W; i++) {
      for (int j = 0; j < H; j++) {
        int idx = j * W + i;

        rot(rays[idx], rays_rot[idx]);

        if (sphere_trace(camera_rot, p1, rays_rot[idx])) {
          b[idx + j] = LUM[(int)(lambertian(p1, light_rot) * 11)];
        } else {
          b[idx + j] = ' ';
        }
      }
    }

    fwrite(b, 1, (W + 1) * H, stdout);
    usleep(DELAY);

    A = fmodf(A + A_SPEED, 2 * M_PI);
    B = fmodf(B + B_SPEED, 2 * M_PI);

    fputs("\x1b[23A", stdout);
  }
}
