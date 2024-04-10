/* 
 * Lab 5 - Triangle Rendering
 * Nitya Harikumar
 * nhariku
 * ECE 4680, Spring 2024
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define ROWS 256
#define COLS 256

typedef struct {
    double x;
    double y;
    double z;
} vector3;

void skip(int lines, FILE *fpt);
unsigned int asciiToInt(FILE *fpt);
vector3 rotateX(vector3 v, double theta);
vector3 rotateY(vector3 v, double theta);
vector3 rotateZ(vector3 v, double theta);

vector3 vectorCross(vector3 a, vector3 b);
double vectorDot (vector3 a, vector3 b);
vector3 vectorSub(vector3 a, vector3 b);
vector3 vectorAdd(vector3 a, vector3 b);
vector3 multiplyScalar(double i, vector3 v);
double magnitude(vector3 v);

int main (int argc, char *argv[]) {
	FILE *fpt, *output;
	unsigned char curr; 
	unsigned int vertices = 0;
	unsigned int faces = 0;

	if (argc != 5) {
		printf("Usage:lab5 file-name X_deg Y_deg Z_deg\n");
		exit(0);
	}

	char input_file[strlen(argv[1]) + 5]; //4 char for ".ply" + 1 null char
	sprintf(input_file,"%s.ply", argv[1]);
	char output_file[strlen(argv[1]) + 5]; //4 char for ".ply" + 1 null char
	sprintf(output_file,"%s.ppm", argv[1]);

	fpt = fopen(input_file, "r");
	if (fpt == NULL) {
		printf("Error in creating input file\n");
		exit(0);
	}
	output = fopen(output_file,"w");

	/* Step 1:i Parse PLY File header*/
	//Skip the first 2 lines
	skip(2, fpt);
	vertices = asciiToInt(fpt);
	
	//Skip the next 3 lines
	skip(3,fpt);
	faces = asciiToInt(fpt);
	skip(2,fpt);

	/* Step 2: Read PLY file*/
	/* Step 3 : Calculate Bounding Box */
	vector3 vertexArr[vertices];
	vector3 max;
	vector3 min;
	vector3 center;
	for (int i = 0; i < vertices; i++) {
		fscanf(fpt, "%lf %lf %lf", &vertexArr[i].x,&vertexArr[i].y,&vertexArr[i].z);
		
		// Set min and max arrays
		if (i == 0) {
			max.x = vertexArr[i].x; min.x = vertexArr[i].x;
			max.y = vertexArr[i].y; min.y = vertexArr[i].y;
			max.z = vertexArr[i].z; min.z = vertexArr[i].z; 		
		}
		if (vertexArr[i].x > max.x) {max.x = vertexArr[i].x;}
		if (vertexArr[i].x < min.x) {min.x = vertexArr[i].x;}

		if (vertexArr[i].y > max.y) {max.y = vertexArr[i].y;}
		if (vertexArr[i].y < min.y) {min.y = vertexArr[i].y;}

		if (vertexArr[i].z > max.z) {max.z = vertexArr[i].z;}
		if (vertexArr[i].z < min.z) {min.z = vertexArr[i].z;}

		// center.x += vertexArr[i].x;
		// center.y += vertexArr[i].y;
		// center.z += vertexArr[i].z;
	}

	// Center is the max-min
	center.x = (max.x + min.x)/ 2.0;
	center.y = (max.y + min.y)/ 2.0;
	center.z = (max.z + min.z)/ 2.0;


	double E = max.x - min.x; // maximum extent of bounding box
	if (E < (max.y - min.y)) {E = (max.y - min.y);} 
	if (E < (max.z - min.z)) {E = (max.z - min.z);} 


	int faceArr[faces][3]; 
	int dontcare; //First value is always a 3
	for (int i = 0; i < faces; i++) {
		fscanf(fpt, "%d %d %d %d",&dontcare,&faceArr[i][0],&faceArr[i][1],&faceArr[i][2]);
	}

	/* Step 4: Calculate the camera position and orientation using two vectors camera and up */
	vector3 up, camera;

	//Default
	camera.x = 1; camera.y = 0; camera.z = 0;
	up.x = 0; up.y = 0; up.z = 1;

	// //Calculations based on rotation matrices	
	double X_rad = atof(argv[2]) * M_PI / 180.0; //x rotation
    double Y_rad = atof(argv[3]) * M_PI / 180.0; //y rotation
    double Z_rad = atof(argv[4]) * M_PI / 180.0; //z rotation
	camera = rotateX(camera, X_rad);
    camera = rotateY(camera, Y_rad);
    camera = rotateZ(camera, Z_rad);
    up = rotateX(up, X_rad);
    up = rotateY(up, Y_rad);
    up = rotateZ(up, Z_rad);
	
	//Formula: <camera> = 1.5 * E * <camera> + <center>
	camera.x = (1.5 * E * camera.x) + center.x;
	camera.y = (1.5 * E * camera.y) + center.y;
	camera.z = (1.5 * E * camera.z) + center.z;
	//printf("center (%lf %lf %lf)\n",center.x,center.y,center.z);

	/* Step 5: Determine the 3D coordinates bounding the image */
	//<left> = <up> x <center - camera>
	vector3 left = vectorCross(up,vectorSub(center,camera));
	//a = ||<left>||
	double a = magnitude(left);
	//printf("%lf\n",a);
	//<left> = (E/(2a) <left>) + <center>
	left = vectorAdd(multiplyScalar(E/(2*a), left),center);	
	//<right> = <center - camera> x <up
	vector3 right = vectorCross(vectorSub(center,camera),up);
	//<right> = (E/(2a) <left>) + <right>
	right = vectorAdd(multiplyScalar(E/(2*a), right), center);
	//<top> = (E/2,<up>) + <center>
	vector3 top = vectorAdd(multiplyScalar(E/2, up), center);
	//<bottom> = (-E/2)<up> + <center>
	vector3 bottom = vectorAdd(multiplyScalar(-1*E/2, up), center);
	//<topleft> = (E/2)<up> + <left>
	vector3 topleft = vectorAdd(multiplyScalar(E/2, up), left);

	/* Step 6: Display Image */
	unsigned char PPM[ROWS][COLS] = {0}; //default image is black
	double z_buff = 999999; //default z-buffer is deep

	printf("Rendering...\n");

	vector3 image;
	for (int c = 0; c < COLS; c++) {
		for (int r = 0; r < ROWS; r++){
			image = multiplyScalar((double)c/(COLS-1),vectorSub(right,left));
			image = vectorAdd(topleft,image);
			image = vectorAdd(image, multiplyScalar((double)r/(ROWS -1), vectorSub(bottom,top)));
			//if (r == 0 && c ==0) printf("%lf %lf %lf",image.x,image.y,image.z);
			//if(r==103 && c ==103) printf("%lf %lf",c/(COLS-1),r/(ROWS -1) );
			//if (r==103 && c ==103) printf("image (%lf %lf %lf)\n",image.x,image.y,image.z);
			for (int i = 0; i < faces; i++) {
				//printf("%d ",r);
				vector3 v0, v1, v2;
				v0 = vertexArr[faceArr[i][0]];
				v1 = vertexArr[faceArr[i][1]];
				v2 = vertexArr[faceArr[i][2]];
				vector3 planeEq = vectorCross(vectorSub(v1,v0),vectorSub(v2,v0)); //plane array equations
				double D = -1*vectorDot(planeEq,v0);
				double n = (-1*vectorDot(planeEq,camera)) - D;
				double d = vectorDot(planeEq,vectorSub(image,camera));

				if (d < 0.01 && d > -0.01) continue; //if ray is parallel, skip the triangle
				vector3 intersect = vectorAdd(camera,multiplyScalar(n/d,vectorSub(image,camera)));
				//if (r==103 && c ==103 && i ==200 ) printf("image (%lf %lf %lf)\n",image.x,image.y,image.z);


				//Separate out dot product
				double dot1 = vectorDot(vectorCross(vectorSub(v2,v0),vectorSub(v1,v0)),vectorCross(vectorSub(intersect,v0),vectorSub(v1,v0)));
				double dot2 = vectorDot(vectorCross(vectorSub(v0,v1),vectorSub(v2,v1)),vectorCross(vectorSub(intersect,v1),vectorSub(v2,v1)));
				double dot3 = vectorDot(vectorCross(vectorSub(v1,v2),vectorSub(v0,v2)),vectorCross(vectorSub(intersect,v2),vectorSub(v0,v2)));

				if (dot1 < 0 || dot2 < 0 || dot3 < 0) continue; //intersection point lies outside of triangle and can be skipped
				//printf("%lf ", n/d);
				if( z_buff < (n/d)) z_buff = n/d; //distance of triangle n/d is > z_buff, then triangle lies behind a closer triangle, so skip
				PPM[r][c] = 155 + (i%100); 
				//printf("%d ",PPM[c][r]);
			}
			z_buff = 999999; //reset z buffer
		}
	}
	/* Step 7: Write PPM image */
	fprintf(output,"P5 %d %d 255\n",COLS,ROWS);
	fwrite(PPM,COLS*ROWS,sizeof(unsigned char),output);

	return 0;
}


/************************************/
/* Functions to read through header */
/************************************/

void skip(int lines, FILE *fpt) {
	int lineNum = 0;
	char curr;
	while (lineNum < lines) {
		fread(&curr,sizeof(char), 1,fpt);
		if (curr == '\n') {
			lineNum++;
		}
	}
}

unsigned int asciiToInt(FILE *fpt) {
	unsigned char curr;
	unsigned int digits[10],i,j,numChar,finalNum;
	finalNum = 0; i = 0;
	while (fread(&curr,sizeof(char), 1,fpt) && curr != '\n') {
		if (curr >= '0' && curr <= '9') {
			numChar = curr - '0';
			digits[i] = numChar;
			i++;
		}
	}
	for (j = 0 ;j < i; j++) {
		finalNum = finalNum * 10 + digits[j];
	}
	return finalNum;
}

/**********************/
/* Rotation Operations */
/*********************/

//Rotate a vector around the X-axis by angle theta
vector3 rotateX(vector3 v, double theta) {
    vector3 result;
    result.x = v.x * cos(theta) - v.y * sin(theta);
    result.y = v.x * sin(theta) + v.y * cos(theta);
    result.z = v.z;
    return result;
}

//Rotate a vector around the Y-axis by angle theta
vector3 rotateY(vector3 v, double theta) {
    vector3 result;
    result.x = v.x;
    result.y = v.y * cos(theta) - v.z * sin(theta);
    result.z = v.y * sin(theta) + v.z * cos(theta);
    return result;
}

//Rotate a vector around the Z-axis by angle theta
vector3 rotateZ(vector3 v, double theta) {
    vector3 result;
    result.x = v.x * cos(theta) + v.z * sin(theta);
    result.y = v.y;
    result.z = -v.x * sin(theta) + v.z * cos(theta);
    return result;
}

/*********************/
/* Vector Operations */
/*********************/
vector3 vectorCross(vector3 a, vector3 b) {
	vector3 result;
	result.x = (a.y * b.z) - (a.z * b.y);
	result.y = (a.z * b.x) - (a.x * b.z);
	result.z = (a.x * b.y) - (a.y * b.x);
	return result;
}

double vectorDot (vector3 a, vector3 b) {
	double result;
	result = (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
	return result;
}

vector3 vectorSub(vector3 a, vector3 b) {
	vector3 result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return result;
}
vector3 vectorAdd(vector3 a, vector3 b) {
	vector3 result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	return result;
}
vector3 multiplyScalar(double i, vector3 v){
	vector3 result; 
	result.x = v.x * i;
	result.y = v.y * i;
	result.z = v.z * i;
	return result;
}
double magnitude(vector3 v){
	double result; 
	result = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
	return result;
}
