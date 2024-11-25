#include<stdio.h>
#include<math.h>
#include<stdlib.h>
double mu_a=0.0;
double mu_b=0.0;
//Returns the dot product for points in the arrays with indexes m,n,o,p
//double data[49125][6];
float length=1;
double dot(int m,int n,int o,int p,double * array[]){
	m=m-1;
	n=n-1;
	o=o-1;
	p=p-1;
	return ((array[m][0]-array[n][0])*(array[o][0]-array[p][0]))+((array[m][1]-array[n][1])*(array[o][1]-array[p][1]))+((array[m][2]-array[n][2])*(array[o][2]-array[p][2]));
	}


//Standard Linspace function
double * linspace(double x1, double x2,int n,double array[]){
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
void points(double P_1[],double P_2[],double P_3[],double P_4[],double temp[]){
	//dst is the distance of closest approach
	double dst=280*(P_2[2]-P_1[2])*(((P_3[0]-P_1[0])*((P_2[1]-P_1[1])-(P_4[1]-P_3[1])))-((P_3[1]-P_1[1])*((P_2[0]-P_1[0])-(P_4[0]-P_3[0]))));
    //this check is to reduce the computations for a long loop later.We are interseted in only those lines that get closer than 0.001
    if (dst>0.001){
		temp[3]=999999; return;
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
    temp[0]=(P_a[0]+P_b[0])/2;
	temp[1]=(P_a[1]+P_b[1])/2;
	temp[2]=(P_a[2]+P_b[2])/2;
	temp[3]=dst;
}
//Calculates the solid angle subtended by a rectangular plate(of size a*b) to a point at a distance d along the perpendicular  axis passing through the center of the plane. 

double omega(double a,double b,double d){
	double alpha=a/(double)2/d;
    double beta=b/(double)2/d;
    return (a*b/pow(d,(double)2));
}






//Calculates the solid angle subtended by a rectangular plate(of size a*b) to a point at a perpendicular distance d from the plane and at a distance A from the nearest edge of length a and at a distance B from the nearest edge of length b.Formula from:  https://vixra.org/pdf/2001.0603v2.pdf

double omega_offaxis(double A,double B,double a, double b,double d){
	return (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/(double)4;
}

//Does a probability calculation for blocks for a given pair of detectors d1 and d2
double calc_prob(double b[],double ** d1,double ** d2){
	double d[2][3];
	double dist;
    	if (b[2]<=200.0){
        	dist=280.0-b[2];
		for(int i=0;i<2;i++){
			for (int j=0;j<3;j++){
				d[i][j]=d2[i][j];
			}
		}
	}
    else{

	for(int i=0;i<2;i++){
		for (int j=0;j<3;j++){
			d[i][j]=d1[i][j];
		}
	}	
        dist=b[2];
	}    
	double x1=d1[0][0],y1=d1[0][1],z1=d1[0][2];
	double x2=d1[1][0],y2=d1[1][1],z2=d1[1][2];
	double x3=d2[0][0],y3=d2[0][1],z3=d2[0][2];
	double x4=d2[1][0],y4=d2[1][1],z4=d2[1][2];

    double x_c_1=((x3-x1)*b[2]/(double)(280))+x1;
    double x_c_2=((x4-x2)*b[2]/(double)(280))+x2;
    double y_c_1=((y3-y1)*b[2]/(double)(280))+y1;
    double y_c_2=((y4-y2)*b[2]/(double)(280))+y2;
    //printf("\n%lf,%lf,%lf\n",b[0],b[1],b[2]);
    //printf("\n%lf,%lf,%lf,%lf\n",x1,y1,x4,y4);

    if ((b[0]<x_c_1) || (b[0]>x_c_2) || (b[1]<y_c_1) || (b[1]>y_c_2)){
        return 0;
    }
	double x_a=d[0][0],y_a=d[0][1],z_a=d[0][2];
	double x_b=d[1][0],y_b=d[1][1],z_b=d[1][2];
    double a1=length;
	double b1=length;
    return (omega_offaxis(x_a-b[0],y_a-b[1],a1, b1, dist)/(double)4*M_PI);
}



///////////////////////////////////////////////////
//Calculates n* as in the formula
int nstar(int i, int j,double ** data,double *** detectors,double h){
	int tempu=0;
	for(int k=0;k<49125;k++){
                        if( (data[k][0]>detectors[i][0][0]) &&( data[k][0]<detectors[i][1][0]) && (data[k][3]>detectors[i][0][1]) && (data[k][3]<detectors[i][1][1])){
                            if( (data[k][1]>detectors[j][0][0]) &&(data[k][1]<detectors[j][1][0]) && (data[k][4]>detectors[j][0][1] )&&( data[k][4]<detectors[j][1][1])){
                                tempu++;
}}
}

return tempu;
}










void main(){
printf("True");
double ** data=(double **)malloc(49125*sizeof( double *));
 for (int i=0;i<49125;i++){
	data[i]=(double *)malloc(6*sizeof(double));
}


//We have a space of size[0,50]*[0,50]*[0,400] we divide it into small blocks of size  1*1*1
//Array Blocks represents the center of all blocks in our space.
double ** blocks=(double **)malloc(100000000*sizeof(double *));
if (blocks==NULL){printf("Holy  !!!!!!");return;}
for(int i=0;i<100000000;i++)
{
	blocks[i]=(double*)malloc(3*sizeof(double));
	if (blocks[i]==NULL){printf("Holy  !!!!!!");return;}
}
double ** intersect=(double **)malloc(1000000*sizeof(double *));
if (intersect==NULL) {printf("Holy  !!!!!!");return;}
for(int i=0;i<1000000;i++)
{
	intersect[i]=(double*)malloc(3*sizeof(double));
	if (intersect[i]==NULL){printf("Holy  !!!!!!");return;}
}
//Detectors are small squares perpendicular to the z-axis. at z=0 and z=400
//The array detectors has the following structure
//detectors[i][0][0]represents lower bound for the x-coordinate of detector square i.
//detectors[i][1][0]represents upper bound for the x-coordinate of detector square i.
//detectors[i][0][1]represents lower bound for the y-coordinate of detector square i.
//detectors[i][1][1]represents upper bound for the y-coordinate of detector square i.
//detectors[i][.][2]represents z-coordinate of teh detector square i.
double *** detectors=(double ***)malloc(60000*sizeof(double **));
if (detectors==NULL){printf("Holy  !!!!!!");return;}
for(int i=0;i<6000;i++)
{
	detectors[i]=(double**)malloc(4*sizeof(double*));
	if (detectors[i]==NULL){printf("Holy  !!!!!!");return;}
for(int j=0;j<3;j++){detectors[i][j]=(double*)malloc(3*sizeof(double));
if (detectors[i][j]==NULL){printf("Holy  !!!!!!");return;}
}}

printf("%lf",M_PI);
printf("True");


//We want to divide our space between [0,50]*[0,50]*[0,400] into blocks of size 1*1*1
//x_blocking,y_blocking,z_blocking represents the slicing of the space along x,y and z axes respectively.
double x_blocking[(int)(50/length)+1];
double y_blocking[(int)(50/length)+1];
double z_blocking[(int)(280/length)+1];
linspace(0.0,50.0,(int)(50/length)+1,x_blocking);
linspace(0.0,50.0,(int)(50/length)+1,y_blocking);
linspace(0.0,280.0,(int)(280/length)+1,z_blocking);
//Reading input as x and y coordinates of points p1 and p2 on a line l.p1 is at z=0 and p2 is at z=400 for all lines.

FILE * filo;
filo=fopen("input.txt","r");
printf("True");
if (filo==NULL){
	printf("%s","noice");
}
double FrontSheild_2=0.0;
double FrontSheild_1=268.0;
double D_dist_1=2.5;
double D_dist_2=-2.5;
float rand1=(float)rand()/(float)RAND_MAX;
float rand2=(float)rand()/(float)RAND_MAX;

for(int i=0;i<49125;i++){
	//double y=200;
	double unused;
	fscanf(filo,"%lf %lf %lf %lf %lf \n",&data[i][0],&data[i][1],&data[i][3],&data[i][4],&unused);
	data[i][1]=data[i][0]+50;
	data[i][1]=data[i][1]+50;
	data[i][1]=data[i][4]+50;
	data[i][1]=data[i][3]+50;
	data[i][2]=((rand1-0.5)*6)+(1.0*FrontSheild_1+D_dist_1);
	data[i][5]=((rand2-0.5)*6)+(1.0*FrontSheild_2+D_dist_2);	
}

FILE * gilo;

gilo=fopen("output.txt","w");
int kl =0;

fprintf(gilo,"//////////////////////////////START\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");



//Here we initialize our blocks and detectors array by using x_blocking,y_blockiing,z_blocking.
//Variables for storing length of arrays,h1 for detectors and h for blocks
int h=0;
int h1=0;

for(int i=1;i<(int)(50/length)+1;i++){
    for (int j=1;j<(int)(50/length)+1;j++){
        for (int k=1;k<(int)(280/length)+1;k++){
        	
            blocks[h][0]=(x_blocking[i-1]+x_blocking[i])/2.0;
            blocks[h][1]=(y_blocking[j-1]+y_blocking[j])/2.0;
            blocks[h][2]=(z_blocking[k-1]+z_blocking[k])/2.0;
            h++;
            if (k==1 || k==(int)(280/length)){
                if (k==(int)(280/length)) {k+=1;}
                for(int index=-1;index<=0;index++){
                detectors[h1][index+1][0]=x_blocking[i+index];
                detectors[h1][index+1][1]=y_blocking[j+index];
                detectors[h1][index+1][2]=z_blocking[k-1];
                }
                h1++;
                }
           }
}}

//The following for loop caculates all the intersection points betwen lines.
for (int i=0;i<49125;i++){
for (int j=i;j<49125;j++){
		double temp[4];
            double P_1[3]={data[i][0],data[i][1],data[i][2]};
            double P_2[3]={data[i][3],data[i][4],data[i][5]};
            double P_3[3]={data[j][0],data[j][1],data[j][2]};
            double P_4[3]={data[j][3],data[j][4],data[j][5]};
         //   """   if project_xy(P_1,P_2,P_3,P_4) && project_xz(P_1,P_2,P_3,P_4) && ////////project_yz(P_1,P_2,P_3,P_4):"""
            points(P_1,P_2,P_3,P_4,temp);
            //printf("%lf\n",temp[3]);
          
           if (temp[3]<=0.0001){
               
                //printf("nice||||");
                    //#xs.append((temp[0][0]+temp[1][0])/2)
                    //#ys.append((temp[0][1]+temp[1][1])/2)
                    //#zs.append((temp[0][2]+temp[1][2])/2)
                    intersect[kl][0]=temp[0];
                    intersect[kl][1]=temp[1];
                    intersect[kl][2]=temp[2];
                    kl++;
                    //fprintf(gilo,"%lf,%lf,%lf",intersect[kl][0],intersect[kl][1],intersect[kl][2],intersect[kl][3]);
	}}}

//Unused code for printing variables  to files

//for (int i=0;i<h;i++){
//
//printf("i,%lf\n",calc_prob(blocks[i],detectors[1],detectors[49]));
//}
//
//fprintf(gilo,"||||||detectors||||||");
//fprintf(gilo,"\n");
//for (int i=0;i<h1;i++){
//fprintf(gilo,"[%lf,%lf,%lf],[%lf,%lf,%lf]",detectors[i][0][0],detectors[i][0][1],detectors[i][0][2],detectors[i][1][0],detectors[i][1][1],detectors[i][1][2]);
//fprintf(gilo,"\n");
//}
//fprintf(gilo,"||||||blocks||||||");
//fprintf(gilo,"\n");
//for (int i=0;i<h;i++){
//fprintf(gilo,"[%lf,%lf,%lf]",blocks[i][0],blocks[i][1],blocks[i][2]);
//fprintf(gilo,"\n");
//}

//Array to store old lambda values

printf("\n%d\n",h);
printf("\n%d\n",h1);

double * lambda_0ld=(double *)malloc(4*h*sizeof(double ));
if (lambda_0ld==NULL){printf("Holy  !!!!!!");return;}


////Array to store new lambda values
double * Lambda=(double *)malloc(4*h*sizeof(double ));
if (Lambda==NULL){printf("Holy  !!!!!!");return;}
//Assigns to each block the number of intersect points in it 
void classify(){
int zero=0;


	for (int j=0;j<h;j++){
	lambda_0ld[j]=0;
	for (int i=0;i<kl;i++){
		
	if((intersect[i][0]<(blocks[j][0]+(length/2)))&&(intersect[i][0]>(blocks[j][0]-(length/2)))&&(intersect[i][1]<(blocks[j][1]+(length/2)))&&(intersect[i][1]>(blocks[j][1]-(length/2)))&&(intersect[i][2]<(blocks[j][2]+(length/2)))&&(intersect[i][2]>(blocks[j][2]-(length/2)))){
		lambda_0ld[j]++;};
		zero++;
	}}

}

classify();
//for (int i=0;i<h;i++){
//fprintf(gilo,"%d,,,,%lf\n",i,lambda_0ld[i]);
//if(lambda_0ld[i]!=0) printf("%d,,,,%lf\n",i,lambda_0ld[i]);
//fflush(gilo);}
//Calculates Lamda new as in the formula
//We don't do the calculation for detectors with same z cooridinate.
printf("%d",h1);
double lambda(int b){
	if (lambda_0ld[b]==0){return 0.0;}
	//classify();	
	double sum=0;
	double var;
	double numerator;
	for(int i1=0;i1<h1;i1++){
		for(int j1=i1;j1<h1;j1++){
			
			if (detectors[i1][0][2]==detectors[j1][0][2]) continue;
			double temp=0.0;
			//if(i1%4==0) printf("----%d-%d-----\n",i1,j1);
			//printf("%d,%d----",i1,j);
			for (int p=0;p<h;p++){
				if (lambda_0ld[p]==0) continue;
				var=0.0;
				//if (detectors[i][2]!=detectors[j][2]){
					var=(lambda_0ld[p]*calc_prob(blocks[p],detectors[i1],detectors[j1]));
					//printf("-------%d-----%lf-----%lf---\n",p,lambda_0ld[p],calc_prob(blocks[p],detectors[i1],detectors[j1]));
						
				temp=temp+var;	
			}
		//if (temp==0.0) printf("\n%d,%d---%f\n",i1,j1,temp);
		numerator=(double)(nstar(i1,j1,data,detectors,h)*calc_prob(blocks[b],detectors[i1],detectors[j1]));
		if (numerator==0.0&&temp==0.0){
			 continue;  //printf("\n%d,%d---%f\n",i1,j1,temp);;
		}
		//if (temp==0.0) printf("\n%d,%d---%f\n",i1,j1,temp);
		sum=sum+(numerator/(double)temp); ////////printf("\nsum:%lf",sum);
		}
	}
	
	return (double)(lambda_0ld[b]*(sum));
}

fprintf(gilo,"//////////////////////////////MLEM\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
for (int u=0;u<h;u++){
	Lambda[u]=lambda(u);
	if (Lambda[u]!=0){
	printf("%d&%lf,%lf,%lf&%lf\n",u,blocks[u][0],blocks[u][1],blocks[u][2],Lambda[u]);
	fprintf(gilo,"%d&%lf,%lf,%lf&%lf\n",u,blocks[u][0],blocks[u][1],blocks[u][2],Lambda[u]);
	fflush(gilo);
}}
fclose(filo);


fclose(gilo);
}
