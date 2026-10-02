#include <math.h>
#include <stdio.h>
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

#define D(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define G(a, b, c)                                                             \
  (f(p[0] + a, p[1] + b, p[2] + c) - f(p[0] - a, p[1] - b, p[2] - c)) / (2 * dt)
#define T(p)                                                                   \
  {p[0] * y + p[1] * u * x + p[2] * v * x, p[1] * v - p[2] * u,                \
   -p[0] * x + p[1] * u * y + p[2] * v * y}

float A, B;
char b[(W + 1) * H];

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
}

int main() {
  float f0 = fabsf(f(0, 0, -CAMERA_Z));
  for (;;) {
    float u = sinf(A), v = cosf(A), x = sinf(B), y = cosf(B),
          l[] = T(((float[]){LIGHT}));
    char *o = b;
    for (int j = 0; j < H; j++, *o++ = '\n')
      for (int i = 0; i < W; i++) {
        float p[] = T(((float[]){0, 0, -CAMERA_Z})),
              n[] = T(((float[]){2 * i + 1 - W, 2 * (H - 2 * j - 1), SCALE})),
              m = sqrtf(D(n, n)), g[3], t, c, d = 0;
        for (int k = 0; k < MAX_STEPS && d <= D_MAX; k++) {
          t = f(p[0], p[1], p[2]);
          g[0] = G(dt, 0, 0), g[1] = G(0, dt, 0), g[2] = G(0, 0, dt);
          c = t / (sqrtf(D(g, g)) + 1e-9f);
          if (fabsf(c) < THRES && fabsf(t) <= f0)
            break;
          for (int q = 0; q < 3; q++)
            p[q] += n[q] / m * c;
          d += c;
        }
        *o++ = d > D_MAX
                   ? ' '
                   : LUM[(int)(fmaxf(0, D(g, l) / sqrtf(D(g, g) * D(l, l))) *
                               (sizeof LUM - 2))];
      }
    fwrite(b, 1, (W + 1) * H, stdout);
    usleep(DELAY);
    A = fmodf(A + A_SPEED, 2 * M_PI);
    B = fmodf(B + B_SPEED, 2 * M_PI);
    printf("\x1b[%dA", H);
  }
}
