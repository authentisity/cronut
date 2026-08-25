#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define max(a, b)           \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a > _b ? _a : _b;      \
  })

#define min(a, b)           \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a < _b ? _a : _b;      \
  })

#define rad(x) \
  (0.0174532925199432954743716805978692718781530857086181640625 * x)

#define CELL_RATIO 2 // terminal characters ~twice as tall as wide

#define W 80
#define H 24
#define FOV 72

#define delta 1e-4 // precision
#define MAX_STEPS 96

#define CAMERA_Z 21

float f(float x, float y, float z)
{
  // Define shape
  const float R1 = 3.0;
  const float R2 = 5.0;
  float inner = sqrt(x * x + y * y) - R1 * R1;
  return inner * inner + z * z - R2 * R2;
}

void grad(float *p, float *g)
{
  g[0] = (f(p[0] + delta, p[1], p[2]) - f(p[0] - delta, p[1], p[2])) / (2 * delta);
  g[1] = (f(p[0], p[1] + delta, p[2]) - f(p[0], p[1] - delta, p[2])) / (2 * delta);
  g[2] = (f(p[0], p[1], p[2] + delta) - f(p[0], p[1], p[2] - delta)) / (2 * delta);
}

float mag(float *p)
{
  float m = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  return m + 1e-9;
}

bool sphere_trace(float *p, float *p1, float *n)
{
  const float thres = 0.001;
  const float d_max = 100.0;

  memcpy(p1, p, 3 * sizeof(float));

  float d = 0.0;

  for (int i = 0; i < MAX_STEPS; i++)
  {
    if (d > d_max)
    {
      return false;
    }

    float t = f(p1[0], p1[1], p1[2]);

    // printf("iteration %i t: %f ", i + 1, t);

    if (t < thres)
    {
      // printf("converged: %f\n", d);
      return true;
    }

    float c = t / ({ float g[3]; grad(p, g); mag(g); });

    p1[0] += n[0] * c;
    p1[1] += n[1] * c;
    p1[2] += n[2] * c;

    d += c;

    // printf("c: %f dist: %f new point: (%f, %f, %f)\n", c, d, p1[0], p1[1],
    //  p1[2]);
  }

  return false;
}

float lambertian(float *p)
{
  const float nl[3] = {
      0.57735026919,
      0.57735026919,
      0.57735026919};

  float g[3];
  grad(p, g);
  float m = mag(g);
  float ng[3] = {g[0] / m, g[1] / m, g[2] / m};

  float Id = max(0, ng[0] * nl[0] + ng[1] * nl[1] + ng[1] * nl[1]);

  return min(1.0, Id + 0.1);
}

void ray(float i, float j, float *arr)
{
  float scale = tan(rad(FOV / 2));
  float u = (2 * i + 1) / W - 1;
  float v = 1 - (2 * j + 1) / H;

  // normalize
  float m = sqrt(u * u + v * v + 1.0);

  arr[0] = u / m;
  arr[1] = v / m;
  arr[2] = 1 / m;
}

int main()
{
  float rays[W * H][3];
  float buf[W * H];

  for (int i = 0; i < W; i++)
  {
    for (int j = 0; j < H; j++)
    {
      ray(i, j, rays[i * H + j]);
    }
  }

  float p1[3];

  for (int i = 0; i < W; i++)
  {
    for (int j = 0; j < H; j++)
    {
      if (sphere_trace((float[]){0, 0, -CAMERA_Z}, p1, rays[i * H + j]))
      {
        buf[i * H + j] = lambertian(p1);
      }
      else
      {
        buf[i * H + j] = 0.0;
      }
      printf("%f\n", buf[i * H + j]);
    }
  }

  printf("\x1b[H");
  for (int j = 0; j < H; j++)
  {
    for (int i = 0; i < W; i++)
    {
      putchar(" .,-~:;=!*#$@"[(int)(buf[i * H + j] * 12)]);
    }
    putchar('\n');
  }
}