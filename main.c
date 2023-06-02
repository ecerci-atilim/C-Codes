#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265359

float* rrcFilter(float beta, float T, float ts, int* lenRC);
float sinc(float x);

int main(){
    int lenRC = 10; // will be 9217 for this case
    float *rrc = rrcFilter(0.3,1.0,1.0/100.0,&lenRC);

    FILE *fid = fopen("rrc0001.txt", "w");

    for(int i = 0; i < lenRC; i++)
        fprintf(fid, "%f\n", *(rrc+i));

    fclose(fid);

    for(int i = 0; i < lenRC; i++)
      printf("Index %05d: %+1.10f\n", i+1, *(rrc+i));

    return 0;
}

float sinc(float x)
{
  return x == 0.0 ? 1.0 : sin(PI*x)/PI/x;
}

float* rrcFilter(float beta, float T, float ts, int* lenRC)
{
  float t;
  const int Nsymb = 12;
  const float sps = T/ts;
  unsigned long N = (unsigned long) Nsymb*sps + 1;
  (*lenRC) = N;
  float *rc = (float *) calloc(N,sizeof(N));

  float max = 0.0;
  float shift = -Nsymb*T/2.0;
  for(int i=0;i<N;i++)
  {
    t = shift+ts*((float) i);

    if(fabs(t) == T/2.0/beta)
    {
      rc[i] = PI*sinc(1.0/2.0/beta)/Nsymb/T/2.0;
      continue;
    }
    float tv = t/T;
    float tvb = beta*tv;
    float t1 = sinc(tv)/T;
    float t2 = cos(PI*tvb);
    float t3 = 1-pow(2.0*tvb,2);
    rc[i] = t1*t2/t3;
    if(rc[i]>max)
    {
      max = rc[i];
    }
  }
  for(int i=0;i<N;i++) rc[i] /= max;
  return rc;
}
