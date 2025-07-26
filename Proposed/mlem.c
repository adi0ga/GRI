#include<omp.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double mu_a=0.0;
double mu_b=0.0;
double shit=0.0;
double   int_point[4];
double threshold=0.001;
//Returns the dot product for points in the arrays with indexes m,n,o,p
double length=1.0;
double length_z=1.0;
double dot(int m,int n,int o,int p,double * array[]){
	//Returns the dot product for points in the arrays with indexes m,n,o,p
	m=m-1;
	n=n-1;
	o=o-1;
	p=p-1;
	return ((array[m][0]-array[n][0])*(array[o][0]-array[p][0]))+((array[m][1]-array[n][1])*(array[o][1]-array[p][1]))+((array[m][2]-array[n][2])*(array[o][2]-array[p][2]));
	}


//Standard Linspace function
double * linspace(double x1, double x2,int n,double array[]){ //Remove array from here
	double * x=calloc(n,sizeof(double));
	double step=(x2-x1)/(double)(n-1);
	for (int i=0;i<n;i++){
		array[i]=x1+((double)i*step);
		x[i]=x1+((double)i*step);
	}
	return x;
}

//Calculates the distance of closest approach for lines passing through p_1 p_2 and p_3 p_4 and modifies temp to store those values.
//temp[0] is the x coordinate of the mid point of the shortest distance line segment
// temp[1] is the y coordinate of the mid point of the shortest distance line segment
//temp[2] is the z coordinate of the mid point of the shortest distance line segment
//temp[3] is the distance of closest approach
void points(double P_1[],double P_2[],double P_3[],double P_4[]){
	//dst is the distance of closest approach
	double dst=268*(((P_3[0]-P_1[0])*((P_2[1]-P_1[1])-(P_4[1]-P_3[1])))-((P_3[1]-P_1[1])*((P_2[0]-P_1[0])-(P_4[0]-P_3[0]))));
    //this check is to reduce the computations for a long loop later.We are interseted in only those lines that get closer than threshold
    if (dst>threshold){
		int_point[3]=999999; return;
		}
	
    double * array[4]={P_1,P_2,P_3,P_4};
    // Using the algorithm here The shortest line between two lines in 3D  https://paulbourke.net/geometry/pointlineplane/
    mu_a=((dot(1,3,4,3,array)*dot(4,3,2,1,array))-(dot(1,3,2,1,array)*dot(4,3,4,3,array)))/((dot(2,1,2,1,array)*dot(4,3,4,3,array))-(dot(4,3,2,1,array)*dot(4,3,2,1,array)));
    mu_b=(dot(1,3,4,3,array)+(mu_a*dot(4,3,2,1,array)))/dot(4,3,4,3,array);
    double P_a[3]={0,0,0};
    double P_b[3]={0,0,0};
    for (int i=0;i<3;i++){
        P_a[i]=P_1[i]+mu_a*(P_2[i]-P_1[i]);
        P_b[i]=P_3[i]+mu_b*(P_4[i]-P_3[i]);
    }
    //printf("%lf,%lf,%lf,\n",P_a[0],P_a[1],P_a[2]);
	dst=pow((pow((P_a[0]-P_b[0]),(double)(2))+pow((P_a[1]-P_b[1]),(double)(2))+pow(((P_a[2]-P_b[2])/10),(double)(2))),(0.5));
	//printf("%lf",dst);
    int_point[0]=(P_a[0]+P_b[0])/2;
	int_point[1]=(P_a[1]+P_b[1])/2;
	int_point[2]=(P_a[2]+P_b[2])/2;
	int_point[3]=dst;
	printf("%lf %lf %lf %lf\n",int_point[0],int_point[1],int_point[2],int_point[3]);
}
//Calculates the solid angle subtended by a rectangular plate(of size a*b) to a point at a distance d along the perpendicular  axis passing through the center of the plane. 

double omega(double a,double b,double d){
    return (a*b/d/d);
}






//Calculates the solid angle subtended by a rectangular plate(of size a*b) to a point at a perpendicular distance d from the plane and at a distance A from the nearest edge of length a and at a distance B from the nearest edge of length b.Formula from:

double omega_offaxis(double A,double B,double a, double b,double d){
	return (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
}






void main(){
//Initialize all arrays

double ** data=(double **)malloc(130746*sizeof( double *));
 for (int i=0;i<130746;i++){
	data[i]=(double *)malloc(6*sizeof(double));
}

//We have a space of size[-25,25]*[-25,25]*[-134,134] we divide it into small blocks of size  1*1*1
//Array Blocks represents the center of all blocks in our space.
double ** blocks=(double **)malloc(200000000*sizeof(double *));

for(int i=0;i<200000000;i++)
{
	blocks[i]=(double*)malloc(3*sizeof(double));
}

double ** intersect=(double **)malloc(6000000*sizeof(double *));
for(int i=0;i<6000000;i++)
{
	intersect[i]=(double*)malloc(4*sizeof(double));

}
//Detectors are small squares perpendicular to the x-axis. at z=-134 and z=134 and at  y=-134 and y=134
//The array detectors has the following structure
//detectors[i][0][0]represents lower bound for the x-coordinate of detector square i.
//detectors[i][1][0]represents upper bound for the x-coordinate of detector square i.
//detectors[i][0][1]represents lower bound for the y-coordinate of detector square i.
//detectors[i][1][1]represents upper bound for the y-coordinate of detector square i.
//detectors[i][0][2]represents lower bound for the z-coordinate of detector square i.
//detectors[i][1][2]represents upper bound for the z-coordinate of detector square i.
double *** detectors=(double ***)malloc(6000000*sizeof(double **));
for(int i=0;i<600000;i++)
{
	detectors[i]=(double**)malloc(4*sizeof(double*));
for(int j=0;j<3;j++){detectors[i][j]=(double*)malloc(3*sizeof(double));
}}


FILE * gilo;

gilo=fopen("output.txt","w");




//We want to divide our space between [-25,25]*[-25,25]*[-134,134] into blocks of size 1*1*1
//x_blocking,y_blocking,z_blocking represents the slicing of the space along x,y and z axes respectively.
double x_blocking[(int)(50/length)+1];
double y_blocking[(int)(50/length)+1];
double z_blocking[(int)(268/length_z)+1];
double z_blocking_1[(int)(50/length)+1];
double y_blocking_1[(int)(268/length_z/2)+1];
double y_blocking_2[(int)(268/length_z/2)+1];

linspace(-25.0,25.0,(int)(50/length)+1,x_blocking);
linspace(-25.0,25.0,(int)(50/length)+1,y_blocking);
linspace(-134.0,134.0,(int)(268/length_z)+1,z_blocking);
linspace(-25.0,25.0,(int)(50/length)+1,z_blocking_1);
linspace(-134.0,-25.0,(int)(109/length_z)+1,y_blocking_1);
linspace(25.0,134.0,(int)(109/length_z)+1,y_blocking_2);
//Reading input as x and y coordinates of points p1 and p2 on a line l.p1 is at z=0 and p2 is at z=268 for all lines.

//Here we initialize our blocks and detectors array by using x_blocking,y_blockiing,z_blocking.
//Variables for storing length of arrays,h1 for detectors and h for blocks
int h=0;
int h1=0;
for(int i=1;i<(int)(50/length)+1;i++){
    for (int j=1;j<(int)(50/length)+1;j++){
        for (int k=1;k<(int)(268/length_z)+1;k++){

            blocks[h][0]=(x_blocking[i-1]+x_blocking[i])/2.0;
            blocks[h][1]=(y_blocking[j-1]+y_blocking[j])/2.0;
            blocks[h][2]=(z_blocking[k-1]+z_blocking[k])/2.0;
            h++;
            if (k==1 || k==(int)(268/length_z)){
                if (k==(int)(268/length_z)) {k+=1;}
                for(int index=-1;index<=0;index++){
                detectors[h1][index+1][0]=x_blocking[i+index];
                detectors[h1][index+1][1]=y_blocking[j+index];
                detectors[h1][index+1][2]=z_blocking[k-1];
                }
                h1++;
                }
           }
}}
for(int i=1;i<(int)(50/length)+1;i++){
    for (int j=1;j<(int)(50/length)+1;j++){
        for (int k=1;k<(int)(109/length_z)+1;k++){

            blocks[h][0]=(x_blocking[i-1]+x_blocking[i])/2.0;
            blocks[h][2]=(z_blocking_1[j-1]+z_blocking_1[j])/2.0;
            blocks[h][1]=(y_blocking_1[k-1]+y_blocking_1[k])/2.0;
            h++;
            if (k==1){
                if (k==(int)(268/length_z)) {k+=1;}
                for(int index=-1;index<=0;index++){
                detectors[h1][index+1][0]=x_blocking[i+index];
                detectors[h1][index+1][2]=z_blocking_1[j+index];
                detectors[h1][index+1][1]=y_blocking_1[k-1];
                }
                h1++;
                }
           }
}}
for(int i=1;i<(int)(50/length)+1;i++){
    for (int j=1;j<(int)(50/length)+1;j++){
        for (int k=1;k<(int)(109/length_z)+1;k++){

            blocks[h][0]=(x_blocking[i-1]+x_blocking[i])/2.0;
            blocks[h][2]=(z_blocking_1[j-1]+z_blocking_1[j])/2.0;
            blocks[h][1]=(y_blocking_2[k-1]+y_blocking_2[k])/2.0;
            h++;
            if (k==(int)(109/length_z)){
                if (k==(int)(109/length_z)) {k+=1;}
                for(int index=-1;index<=0;index++){
                detectors[h1][index+1][0]=x_blocking[i+index];
                detectors[h1][index+1][2]=z_blocking_1[j+index];
                detectors[h1][index+1][1]=y_blocking_2[k-1];
                }
                h1++;
                }
           }
}}


double ** nstar=(double **)malloc(h1*sizeof(double *));

for(int i=0;i<h1;i++)
{
	nstar[i]=(double*)malloc(h1*sizeof(double));
}
for(int i=0;i<h1;i++){
	for(int j=0;j<h1;j++){
	nstar[i][j]=0.0;
}
}
double * lambda_0ld=(double *)malloc(h*sizeof(double ));
for (int j=0;j<h;j++){
	lambda_0ld[j]=0.0;
}
printf("Pass");
////Array to store new lambda values
double * Lambda=(double *)malloc(h*sizeof(double ));
double ** array_temp=(double **)malloc(h1*sizeof(double *));
for(int i=0;i<h1;i++){
	array_temp[i]=(double*)malloc(h1*sizeof(double));
}


printf("yay!!");


FILE * fila;
FILE * filb;

fila=fopen("inputa.txt","r");
filb=fopen("inputb.txt","r");
//int count=0;
int ind=0;

for(int i=0;i<81621;i++){
	//double y=200;
	fscanf(fila,"%lf %lf %lf %lf %lf \n",&data[ind][0],&data[ind][1],&data[ind][3],&data[ind][4],&shit);
	//if ((data[ind][0]<0.0)||(data[ind][1]<0.0)||(data[ind][3]<0.0)||(data[ind][4]<0.0)||(data[ind][0]>50.0)||(data[ind][1]>50.0)||(data[ind][3]>50.0)||(data[ind][4]>50.0)) {continue;}
	int n1=(int)data[ind][0];if (data[ind][0]<0) n1-=1;
	int n2 =(int)data[ind][1];if (data[ind][1]<0) n2-=1;
	int na=(100*n1)+(2*n2)+2550;
	int n3=(int)data[ind][3];if (data[ind][3]<0) n3-=1;
	int n4 =(int)data[ind][4];if (data[ind][4]<0) n4-=1;
	int nb=(100*n3)+(2*n4)+1+2550;
	data[ind][2]=-134;
	data[ind][5]=134;
	ind++;
	if ((na>=0)&&(nb>=0)&&(na<5000)&&(nb<5000)){
	nstar[na][nb]=nstar[na][nb]+1;
}
}
for(int i=0;i<49125;i++){
	//double y=200;
	fscanf(filb,"%lf %lf %lf %lf %lf\n",&data[ind][0],&data[ind][2],&data[ind][3],&data[ind][5],&shit);
	//if ((data[ind][0]<0.0)||(data[ind][1]<0.0)||(data[ind][3]<0.0)||(data[ind][4]<0.0)||(data[ind][0]>50.0)||(data[ind][1]>50.0)||(data[ind][3]>50.0)||(data[ind][4]>50.0)) {continue;}
	int n1=(int)data[ind][0];if (data[ind][0]<0) n1-=1;
	int n2 =(int)data[ind][2];if (data[ind][2]<0) n2-=1;
	int na=(50*n1)+(n2)+6275;
	int n3=(int)data[ind][3];if (data[ind][3]<0) n3-=1;
	int n4 =(int)data[ind][5];if (data[ind][5]<0) n4-=1;
	int nb=(50*n3)+(n4)+8775;
	data[ind][1]=-134;
	data[ind][4]=134;
	ind++;
	if ((na>=0)&&(nb>=0)&&(na<h1)&&(nb<h1)){
	nstar[na][nb]=nstar[na][nb]+1;
}
}

long kl =0;

//fprintf(gilo,"//////////////////////////////START\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
printf("Pass");



//The following for loop caculates all the intersection points betwen lines.
#pragma omp parallel for collapse(2)
for (int i=0;i<ind;i++){
for (int j=i+1;j<ind;j++){
            double P_1[3]={data[i][0],data[i][1],data[i][2]};
            double P_2[3]={data[i][3],data[i][4],data[i][5]};
            double P_3[3]={data[j][0],data[j][1],data[j][2]};
            double P_4[3]={data[j][3],data[j][4],data[j][5]};
         //   """   if project_xy(P_1,P_2,P_3,P_4) && project_xz(P_1,P_2,P_3,P_4) && ////////project_yz(P_1,P_2,P_3,P_4):"""
            //points(P_1,P_2,P_3,P_4);
			double dst;
			//dst=268*(((P_3[0]-P_1[0])*((P_2[1]-P_1[1])-(P_4[1]-P_3[1])))-((P_3[1]-P_1[1])*((P_2[0]-P_1[0])-(P_4[0]-P_3[0]))));
    //this check is to reduce the computations for a long loop later.We are interseted in only those lines that get closer than threshold
    //if (dst>=threshold) continue;

    double * array[4]={P_1,P_2,P_3,P_4};
    // Using the algorithm here The shortest line between two lines in 3D  https://paulbourke.net/geometry/pointlineplane/
    mu_a=((dot(1,3,4,3,array)*dot(4,3,2,1,array))-(dot(1,3,2,1,array)*dot(4,3,4,3,array)))/((dot(2,1,2,1,array)*dot(4,3,4,3,array))-(dot(4,3,2,1,array)*dot(4,3,2,1,array)));
    mu_b=(dot(1,3,4,3,array)+(mu_a*dot(4,3,2,1,array)))/dot(4,3,4,3,array);
    double P_a[3]={0,0,0};
    double P_b[3]={0,0,0};
    for (int q=0;q<3;q++){
        P_a[q]=P_1[q]+mu_a*(P_2[q]-P_1[q]);
        P_b[q]=P_3[q]+mu_b*(P_4[q]-P_3[q]);
    }
    //printf("%lf,%lf,%lf,\n",P_a[0],P_a[1],P_a[2]);
	dst=pow((pow((P_a[0]-P_b[0]),(double)(2))+pow((P_a[1]-P_b[1]),(double)(2))+pow(((P_a[2]-P_b[2])/10),(double)(2))),(0.5));
	//printf("%lf",dst);
    if (dst>=threshold) continue;
	intersect[kl][0]=(P_a[0]+P_b[0])/(double)2;
	intersect[kl][1]=(P_a[1]+P_b[1])/2;
	intersect[kl][2]=(P_a[2]+P_b[2])/2;
	intersect[kl][3]=dst;

                    kl++;
                    fprintf(gilo,"%d-%lf,%lf,%lf--%lf\n",kl,intersect[kl-1][0],intersect[kl-1][1],intersect[kl-1][2],intersect[kl-1][3]);
	}}

//Unused code for printing variables  to files
/*
fprintf(gilo,"||||||detectors||||||");
fprintf(gilo,"\n");
for (int i=0;i<h1;i++){
fprintf(gilo,"[%lf,%lf,%lf],[%lf,%lf,%lf]",detectors[i][0][0],detectors[i][0][1],detectors[i][0][2],detectors[i][1][0],detectors[i][1][1],detectors[i][1][2]);
fprintf(gilo,"\n");
}
fprintf(gilo,"||||||blocks||||||");
fprintf(gilo,"\n");
for (int i=0;i<h;i++){
fprintf(gilo,"[%lf,%lf,%lf]",blocks[i][0],blocks[i][1],blocks[i][2]);
fprintf(gilo,"\n");
}
*/

//Array to store old lambda values
printf("Pass");
printf("\n%d\n",h);
printf("\n%d\n",h1);









double calc_prob(int index,int i1,int i2, int indi){
if(indi==0){if (index>=670000){return 0.0;}
	double x_c_1=((detectors[i2][0][0]-detectors[i1][0][0])*(blocks[index][2]+134.0)/(268.0))+detectors[i1][0][0];
    double x_c_2=((detectors[i2][1][0]-detectors[i1][1][0])*(blocks[index][2]-134.0)/(268.0))+detectors[i1][1][0];
    double y_c_1=((detectors[i2][0][1]-detectors[i1][0][1])*(blocks[index][2]+134.0)/(268.0))+detectors[i1][0][1];
    double y_c_2=((detectors[i2][1][1]-detectors[i1][1][1])*(blocks[index][2]-134.0)/(268.0))+detectors[i1][1][1];
    //fprintf(gilo,"\n%lf,%lf,%lf\n",blocks[index][0],blocks[index][1],blocks[index][2]);
    //fprintf(gilo,"\n%lf,%lf,%lf,%lf\n",x_c_1,y_c_1,x_c_2,y_c_2);

    if (!((((blocks[index][0]-0.5)<=x_c_1)&&((blocks[index][0]+0.5)>=x_c_1) )&&(((blocks[index][1]-0.5)<=y_c_1)&&((blocks[index][1]+0.5)>=y_c_1) )|| (((blocks[index][0]-0.5)<=x_c_2)&&((blocks[index][0]+0.5)>=x_c_2) )&&(((blocks[index][1]-0.5)<=y_c_2)&&((blocks[index][1]+0.5)>=y_c_2) ))){
        return 0.0;
    }
    //fprintf(gilo,"\n True \n");
	int dex;
    double dist;
	if (blocks[index][2]<=0.0){
        	dist=268.0-blocks[index][2];
			dex=1;
	}
    else{
		dex=2;
        dist=blocks[index][2];
	}
/*
	double jeel1=blocks[index][0];
	double jeel2=blocks[index][1];
	if (jeel1>0.0) jeel1=25.0-jeel1;else jeel1=-jeel1;
	if (jeel2>0.0) jeel2=25.0-jeel2;else jeel2=-jeel2;

	double sangle1=((omega(2*(50.0-jeel1),2*(50.0-jeel2),dist))+(omega(2*(50.0),2*(50.0-jeel2),dist))+(omega(2*(50.0-jeel1),2*(50.0),dist))+(omega(2*(50.0),2*(50.0),dist)))/4.0;
*/
	//double sangle2=((omega(2*(50.0-blocks[index][0]),2*(50.0-blocks[index][1]),268.0-dist))+(omega(2*(50.0),2*(50.0-blocks[index][1]),268.0-dist))+(omega(2*(50.0-blocks[index][0]),2*(50.0),268.0-dist))+(omega(2*(50.0-blocks[index][0]),2*(50.0),268.0-dist)))/4.0;

    if ((detectors[i1][0][0]==detectors[i2][0][0])&&(detectors[i1][0][1]==detectors[i2][0][1])) return omega(length,length,dist)/(4*M_PI);
	double x_a,x_b,y_a,y_b;
	if (dex==1){
		double x_a=detectors[i1][0][0],y_a=detectors[i1][0][1];
		double x_b=detectors[i1][1][0],y_b=detectors[i1][1][1];
	}
	else{
		double x_a=detectors[i2][0][0],y_a=detectors[i2][0][1];
		double x_b=detectors[i2][1][0],y_b=detectors[i2][1][1];
	}
	double A1=x_a-blocks[index][0];double A2=x_b-blocks[index][0];if (A1<0) A1=-A1;if (A2<0) A2=-A2;
	double A;if (A1>A2) A=A1;else A=A2;
	double B1=y_a-blocks[index][1];double B2=y_b-blocks[index][1];if (B1<0) B1=-B1;if (B2<0) B2=-B2;
	double B;if (B1>B2) B=B1;else B=B2;
    return (omega_offaxis(A,B,length, length, dist)/(4*M_PI));
	}
	else{
	double x_c_1=((detectors[i2][0][0]-detectors[i1][0][0])*(blocks[index][2]+134.0)/(268.0))+detectors[i1][0][0];
    double x_c_2=((detectors[i2][1][0]-detectors[i1][1][0])*(blocks[index][2]-134.0)/(268.0))+detectors[i1][1][0];
    double z_c_1=((detectors[i2][0][1]-detectors[i1][0][1])*(blocks[index][2]+134.0)/(268.0))+detectors[i1][0][1];
    double z_c_2=((detectors[i2][1][1]-detectors[i1][1][1])*(blocks[index][2]-134.0)/(268.0))+detectors[i1][1][1];
    //fprintf(gilo,"\n%lf,%lf,%lf\n",blocks[index][0],blocks[index][1],blocks[index][2]);
    //fprintf(gilo,"\n%lf,%lf,%lf,%lf\n",x_c_1,y_c_1,x_c_2,y_c_2);

    if (!((((blocks[index][0]-0.5)<=x_c_1)&&((blocks[index][0]+0.5)>=x_c_1) )&&(((blocks[index][2]-0.5)<=z_c_1)&&((blocks[index][2]+0.5)>=z_c_1) )|| (((blocks[index][0]-0.5)<=x_c_2)&&((blocks[index][0]+0.5)>=x_c_2) )&&(((blocks[index][2]-0.5)<=z_c_2)&&((blocks[index][2]+0.5)>=z_c_2) ))){
        return 0.0;
    }
    //fprintf(gilo,"\n True \n");
	int dex;
    double dist;
	if (blocks[index][1]<=0.0){
        	dist=268.0-blocks[index][1];
			dex=1;
	}
    else{
		dex=2;
        dist=blocks[index][1];
	}
/*
	double jeel1=blocks[index][0];
	double jeel2=blocks[index][2];
	if (jeel1>0.0) jeel1=25.0-jeel1;else jeel1=-jeel1;
	if (jeel2>0.0) jeel2=25.0-jeel2;else jeel2=-jeel2;

	double sangle1=((omega(2*(50.0-jeel1),2*(50.0-jeel2),dist))+(omega(2*(50.0),2*(50.0-jeel2),dist))+(omega(2*(50.0-jeel1),2*(50.0),dist))+(omega(2*(50.0),2*(50.0),dist)))/4.0;
*/
	//double sangle2=((omega(2*(50.0-blocks[index][0]),2*(50.0-blocks[index][1]),268.0-dist))+(omega(2*(50.0),2*(50.0-blocks[index][1]),268.0-dist))+(omega(2*(50.0-blocks[index][0]),2*(50.0),268.0-dist))+(omega(2*(50.0-blocks[index][0]),2*(50.0),268.0-dist)))/4.0;

    if ((detectors[i1][0][0]==detectors[i2][0][0])&&(detectors[i1][0][2]==detectors[i2][0][2])) return omega(length,length,dist)/(4*M_PI);
	double x_a,x_b,z_a,z_b;
	if (dex==1){
		double x_a=detectors[i1][0][0],z_a=detectors[i1][0][2];
		double x_b=detectors[i1][1][0],z_b=detectors[i1][1][2];
	}
	else{
		double x_a=detectors[i2][0][0],z_a=detectors[i2][0][2];
		double x_b=detectors[i2][1][0],z_b=detectors[i2][1][2];
	}
	double A1=x_a-blocks[index][0];double A2=x_b-blocks[index][0];if (A1<0) A1=-A1;if (A2<0) A2=-A2;
	double A;if (A1<A2) A=A1;else A=A2;
	double B1=z_a-blocks[index][2];double B2=z_b-blocks[index][2];if (B1<0) B1=-B1;if (B2<0) B2=-B2;
	double B;if (B1<B2) B=B1;else B=B2;
    return (omega_offaxis(A,B,length, length, dist)/(4*M_PI));
	}
}



















//Assigns to each block the number of intersect points in it 
///void classify(){
long zero=0;
#pragma omp parallel for
//for (int j=0;j<h;j++){
	for (int i=0;i<kl;i++){
		int x1=(int)intersect[i][0];if (intersect[i][0]<0.0) x1-=1.0;
		int y1=(int)intersect[i][1];if (intersect[i][1]<0.0) y1-=1.0;
		int z1=(int)intersect[i][2];if (intersect[i][2]<0.0) z1-=1.0;
		if ((x1>=-25.0)&&(y1>=-25.0)&&(z1>=-134.0)&&(x1<25.0)&&(y1<25.0)&&(z1<134.0)){
		int nl=(x1*13400)+(y1*268)+z1+341834;

		if ((nl<h)&&(nl>=0)) 	lambda_0ld[nl]++;
	}
	if ((x1>=-25.0)&&(y1>=-134.0)&&(z1>=-25.0)&&(x1<25.0)&&(y1<-25.0)&&(z1<25.0)){
		int nl=(x1*5450)+(y1)+z1*109+809109;

		if ((nl<h)&&(nl>=0)) 	lambda_0ld[nl]++;
	}
	if ((x1>=-25.0)&&(y1>=25.0)&&(z1>=-25.0)&&(x1<25.0)&&(y1<134.0)&&(z1<25.0)){
		int nl=(x1*5450)+(y1)+z1*109+1081450;

		if ((nl<h)&&(nl>=0)) 	lambda_0ld[nl]++;
	}

	}
	//}


printf("%ld#%ld#",zero, kl);

//		 }
//classify();
printf("Pass");




















//Computation for revised lambda

double lambda(int b){
	if (lambda_0ld[b]==0){return 0.0;}
	double sum=0.0;
	double var;
	double numerator;
	#pragma omp parallel for collapse(2)
	for(int i1=0;i1<5000;i1+=2){
		for(int j1=1;j1<5000;j1+=2){
		double temp0=0.0;
		temp0=array_temp[i1][j1];
		numerator=(nstar[i1][j1]*calc_prob(b,i1,j1,0));
		if ((numerator==0.0)&&(temp0==0.0)){
			 continue;
		}
		
			sum=sum+(numerator/(double)temp0);
		}
	}
	for(int i1=5000;i1<7500;i1++){
		for(int j1=7500;j1<h1;j1++){
		double temp0=0.0;
		temp0=array_temp[i1][j1];
		numerator=(nstar[i1][j1]*calc_prob(b,i1,j1,1));
		if ((numerator==0.0)&&(temp0==0.0)){
			 continue;
		}

			sum=sum+(numerator/(double)temp0);
		}
	}
	return (double)(lambda_0ld[b]*(sum));
}


	int * array_index=(int *)malloc(h*sizeof(int));
	int new_var=0;

	for (int i=0;i<h;i++){

		if(lambda_0ld[i]!=0.00){
			array_index[new_var]=i;
			new_var++;
			}
	}
	//Calculates Lamda new as in the formula
	//We don't do the calculation for detectors with same z cooridinate.
	printf("**%d**",new_var);

	printf("Pass");
for(int y=0;y<5;y++){
	#pragma omp parallel for collapse(2)
	for(int i1=0;i1<h1;i1+=2){
		for(int j1=1;j1<h1;j1+=2){
			double temp=0.0;
			//#pragma omp parallel for
			for (int p=0;p<new_var;p++){
				int index=array_index[p];
				double var=0.0;
					var=(lambda_0ld[index]*calc_prob(index,i1,j1,0));
				temp=temp+var;
			}
		array_temp[i1][j1]=temp;
		}
	}
	for(int i1=5000;i1<7500;i1++){
		for(int j1=7500;j1<h1;j1++){
			double temp=0.0;
			//#pragma omp parallel for
			for (int p=0;p<new_var;p++){
				int index=array_index[p];
				double var=0.0;
					var=(lambda_0ld[index]*calc_prob(index,i1,j1,1));
				temp=temp+var;
			}
		array_temp[i1][j1]=temp;
		}
	}
	printf("Hello");



	fprintf(gilo,"//////////////////////////////MLEM\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");





	//#pragma omp parallel for

		fprintf(gilo,"//////////////////////////////ITERATION_%d\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n",y);
	for (int u=0;u<new_var;u++){
		Lambda[array_index[u]]=lambda(array_index[u]);
		fprintf(gilo,"%d&%lf,%lf,%lf&%lf\n",array_index[u],blocks[array_index[u]][0],blocks[array_index[u]][1],blocks[array_index[u]][2],Lambda[array_index[u]]);
	}
	for (int u=0;u<new_var;u++){
		lambda_0ld[array_index[u]]=Lambda[array_index[u]];
		Lambda[array_index[u]]=0.0;
		}

		}

fclose(fila);
fclose(filb);
fclose(gilo);


}
