#include <stdio.h>
//#include <mkl.h>

extern void softposit_(float* rot, float* trans, int* foundPose, int* imagePts, float* worldPts, int* nbImagePts, int* nbWorldPts, float* beta0, float* noiseStd, float* initRot, float* initTrans, float* focalLength, int* center);

int main(void)
{
	//static float a[6] = { 1, 2, 3,4, 5,6 };
	//static float b[2][3] = { {1, 2, 4}, {3, 6, 9}};
	//static float c[3][3] ;
    //float a = 1.0;
	//float b = 2.0;
	//float c = 0.0;
    //float *pa = a;
	//float *pb = b;
	//float *pc = c;

	/* LARGE_INTEGER litmp; */
	/* LONGLONG QPart1,QPart2; */
	/* double dfFreq; */

	//static int imagePts[] = { 299, 641, 330, 522, 249, 254, 55, 436 };
	//static int imagePts1[] = { 407,593,344,481,235,408, 68, 83,172,182,280,411};
	//static int imagePts2[] = { 407,593,344,481,408, 235,68, 83,172,182,411,280};

	static int imagePts[] = {612, 486, 567, 476, 441, 117, 145, 206, 234, 329};
    //static float worldPts[] = { -3.75, 7.50, -3, 3, 0, 0, -5.0, 5.0, 0.5, 2.75, -2.0, -2.0};
	static float worldPts[] = { -3.75, 7.50, -3, 3, 0, 0, 0, 0, -5.0, 5.0, 2.25, -2.25, 0.5, 2.75, -2.0, -2.0, -0.75, -0.75};

	float beta0 = 2.0E-4f;
    float noiseStd = 10.0f;
    //static float initRot[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	static float initRot[] = {0, 1, 0, -1, 0, 0, 0, 0, 1};
    static float initTrans[] = {0.0, 0.0, 30.0};
    float focalLength = 982.1f;
	static int center[] = {376, 240};

	static float rot[3][3];
	static float trans[3];

	int counter = 0;
	int foundPose = 0;
	int nbImagePts = 5;
	int nbWorldPts = 6;

	/* QueryPerformanceFrequency(&litmp); */
	/* dfFreq = (double)litmp.QuadPart;// Get the clock frequency */

	/* QueryPerformanceCounter(&litmp); */
	/* QPart1 = litmp.QuadPart; */
    
	for(counter = 0; counter<1; counter++){
	  softposit_((float*) rot, trans, &foundPose, imagePts, worldPts, &nbImagePts, &nbWorldPts, &beta0, &noiseStd, initRot, initTrans, &focalLength, center);
    }

	/* QueryPerformanceCounter(&litmp); */
	/* QPart2 = litmp.QuadPart;//Get current value from counter */
	/* printf("\nTime difference = %lf\n", (double)(QPart2-QPart1) / dfFreq);// Calculate time difference */
	printf("Rotation:\n");
	printf("%f, %f, %f\n", rot[0][0], rot[1][0], rot[2][0]);
	printf("%f, %f, %f\n", rot[0][1], rot[1][1], rot[2][1]);
	printf("%f, %f, %f\n", rot[0][2], rot[1][2], rot[2][2]);

	printf("\n");
	printf("Translation:\n");
	printf("%f, %f, %f\n", trans[0], trans[1], trans[2]);



	return 0;

}
