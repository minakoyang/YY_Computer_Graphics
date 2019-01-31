//
//		          Programming Assignment #3
//
//			        Yuzhe Yang
//
//
//
/***************************************************************************/

                                                   /* Include needed files */
//using namespace std;
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <GL/glut.h>        // The GL Utility Toolkit (Glut) Header
#include <limits>

#define WIDTH 500
#define HEIGHT 500

int x_last,y_last;

double dinf = std::numeric_limits<double>::infinity();

int pointCounter = 0;       //For the first 4 control points

int updatecounter = 0;      //For the updated control points

int update = 0;             //State register

int selectsize = 6;         //Size of the diamond (point selecting area)

double segmentation=0.06;    //For DDA line approximation, have tested the value for best result

double red=1.0, green=1.0, blue=1.0;

double zbuffer[WIDTH][HEIGHT];

double initz, zmax, zmin, zorigin;

const double PI = 3.14159265358979323846;

int modereg=0;

int drawreg=0;

static int seed;

int seedcounter=0;

int wireframemode=0;

int perspectivemode=0;

int phongSmooth = 0;

int bumpMapping = 0;

//static double t=dinf;

double kd = 1.0; // phong model diffuse weight
double ks = 0.8; // phong model specular weight
double ka = 0.5;
double shininess = 100;   // phong specular exponent

//double pvector[4]={};

class Point{
public: double x, y, z;
    void setxyz(double x1, double y1, double z1)
    { x=x1; y=y1; z=z1;}
};

Point points[4];

/*   OBJ Structures */

// Our vertex type
typedef struct{
	double x,y,z;
}vertexinfo;

// Our normal type
typedef struct{
	double i,j,k;
}norminfo;

// The polygon (triangle), 3 numbers that aim 3 vertices
typedef struct{
	int v[3],t[3],n[3];
}polygoninfo;

// The mapcoord type, 2 texture coordinates for each vertex
typedef struct{
	double u,v;
}textureinfo;

typedef struct{
	double red, green, blue;
}colorinfo;

typedef struct{
	std::vector<int> vv;
	std::vector<int> tt;
	std::vector<int> nn;
}faceinfo;


// The object type
class object{
	public:
	std::vector<vertexinfo> vertex;
	std::vector<textureinfo> texture;
	std::vector<norminfo> norm;
	std::vector<faceinfo> face;
	std::vector<polygoninfo> polygon;
	object(){}
	~object(){}
	void objloader(const char *filename);

};

typedef struct{
	Point lightposition;
	double red, green, blue, Ld, Ls, La;
}lightinfo;

class light{
	public:
	std::vector<lightinfo> lightpoint;
	light(){}
	~light(){}
};

light lightsource;

object obj;

object objbackup;

bool cw(int i, int j, int k)
{
    	double cwv0= (obj.vertex[j].x-obj.vertex[i].x)*(obj.vertex[j].x-obj.vertex[i].x)+(obj.vertex[j].y-obj.vertex[i].y)*(obj.vertex[j].y-obj.vertex[i].y)+(obj.vertex[j].z-obj.vertex[i].z)*(obj.vertex[j].z-obj.vertex[i].z);
	double cwv1= (obj.vertex[k].x-obj.vertex[i].x)*(obj.vertex[k].x-obj.vertex[i].x)+(obj.vertex[k].y-obj.vertex[i].y)*(obj.vertex[k].y-obj.vertex[i].y)+(obj.vertex[k].z-obj.vertex[i].z)*(obj.vertex[k].z-obj.vertex[i].z);
	double cwv2= (obj.vertex[j].x-obj.vertex[k].x)*(obj.vertex[j].x-obj.vertex[k].x)+(obj.vertex[j].y-obj.vertex[k].y)*(obj.vertex[j].y-obj.vertex[k].y)+(obj.vertex[j].z-obj.vertex[k].z)*(obj.vertex[j].z-obj.vertex[k].z);

	if (cwv0+cwv1<cwv2 || cwv0+cwv2<cwv1 || cwv2+cwv1<cwv0){
		return true;

		}
	else
		return false;
}

bool point_in_triangle(int verid, polygoninfo p)
{

	for (unsigned int i=0; i<3; i++)
		if (cw(p.v[i],p.v[(i+1)%3],verid)) return false;

	return true;
}

bool ear_Q(int i, int j, int k, faceinfo f)
{
	polygoninfo p;			/* coordinates for points i,j,k */
	p.v[0]=f.vv[i];
	p.v[1]=f.vv[j];
	p.v[2]=f.vv[k];



	if (cw(p.v[0],p.v[1],p.v[2])) return false;

	for (int m=0; m<(int)f.vv.size(); m++) {
		if ((m!=i) && (m!=j) && (m!=k))
			{

			if (point_in_triangle(f.vv[m],p))
			{

				return false;
			}
			}
	}

	return true;
}


Point crossProduct(Point vector1, Point vector2) {
	static Point cpresult;
        cpresult.setxyz(vector1.y * vector2.z - vector1.z * vector2.y, vector1.z * vector2.x - vector1.x * vector2.z, vector1.x * vector2.y - vector1.y * vector2.x);
	return cpresult;
}

void triangulate(faceinfo f)
{



	int l[100], r[100];	/* left/right neighbor indices */


	for (unsigned int i=0; i<f.vv.size(); i++) {	/* initialization */
		l[i] = ((i-1) + f.vv.size()) % f.vv.size();
		r[i] = ((i+1) + f.vv.size()) % f.vv.size();
	}


	int tcounter=0;
	int i = f.vv.size()-1;
	while (tcounter < (int)(f.vv.size()-2)) {
		i = r[i];

		if (cw(l[i],i,r[i])==false) {


			polygoninfo p;
			p.v[0]=f.vv[l[i]];
			p.v[1]=f.vv[i];
			p.v[2]=f.vv[r[i]];
			p.n[0]=f.nn[l[i]];
			p.n[1]=f.nn[i];
			p.n[2]=f.nn[r[i]];
			p.t[0]=f.tt[l[i]];
			p.t[1]=f.tt[i];
			p.t[2]=f.tt[r[i]];

			obj.polygon.push_back(p);
			tcounter=tcounter+1;

			l[ r[i] ] = l[i];
			r[ l[i] ] = r[i];
		}
	}
}

void triangulation(faceinfo f)
{
	  	polygoninfo p;

		p.v[0]=f.vv[0];
		p.v[1]=f.vv[1];
		p.v[2]=f.vv[2];
		p.n[0]=f.nn[0];
		p.n[1]=f.nn[1];
		p.n[2]=f.nn[2];
		p.t[0]=f.tt[0];
		p.t[1]=f.tt[1];
		p.t[2]=f.tt[2];

	        obj.polygon.push_back(p);

		std::vector<polygoninfo> tempp;
		tempp.push_back(p);
	for(unsigned int i=3; i<f.vv.size(); i++)
	{
		polygoninfo p;
		for(unsigned int j=0; j<tempp.size(); j++)
		{

		    for(int m=0; m<3; m++)
		    {
			Point vgiven, vside1;
			vgiven.setxyz(obj.vertex[tempp[j].v[m]].x-obj.vertex[tempp[j].v[(m+1)%3]].x, obj.vertex[tempp[j].v[m]].y -obj.vertex[tempp[j].v[(m+1)%3]].y, obj.vertex[tempp[j].v[m]].z-obj.vertex[tempp[j].v[(m+1)%3]].z);
			vside1.setxyz(obj.vertex[tempp[j].v[m]].x-obj.vertex[f.vv[i]].x, obj.vertex[tempp[j].v[m]].y -obj.vertex[f.vv[i]].y, obj.vertex[tempp[j].v[m]].z-obj.vertex[f.vv[i]].z);
			int oneside=0;
			for(unsigned int k=0; k<tempp.size(); k++)
			{

				for (int p=0; p<3; p++)
				{
					Point vside2, vP1, vP2;
					vside2.setxyz(obj.vertex[tempp[j].v[m]].x-obj.vertex[tempp[k].v[p]].x, obj.vertex[tempp[j].v[m]].y -obj.vertex[tempp[k].v[p]].y, obj.vertex[tempp[j].v[m]].z-obj.vertex[tempp[k].v[p]].z);
					vP1=crossProduct(vgiven, vside1);
					vP2=crossProduct(vgiven, vside2);

					if(abs(vP1.x+vP2.x)>abs(vP1.x-vP2.x) || abs(vP1.y+vP2.y)>abs(vP1.y-vP2.y) || abs(vP1.z+vP2.z)>abs(vP1.z-vP2.z))
					{
						oneside=oneside+1;



					}
				}
			}

				if(oneside==0)
				{

		printf("here");
						p.v[0]=tempp[j].v[m];
						p.v[2]=tempp[j].v[(m+1)%3];
						p.v[1]=f.vv[i];
						p.n[0]=tempp[j].n[m];
						p.n[2]=tempp[j].n[(m+1)%3];
						p.n[1]=f.nn[i];
						p.t[0]=tempp[j].t[m];
						p.t[2]=tempp[j].t[(m+1)%3];
						p.t[1]=f.tt[i];


				}
		    }
		}
		obj.polygon.push_back(p);
		tempp.push_back(p);

	}

}



void object::objloader(const char *filename)
{
    std::ifstream in(filename, std::ios::in);
    if (!in)
    {
        std::cerr << "Cannot open " << filename << std::endl; exit(1);
    }
	srand(time(NULL));
	seed=rand();

    std::string line;
    while (getline(in, line))
    {
        if (line.substr(0,2) == "v ")
        {
            std::istringstream s(line.substr(2));
            vertexinfo v;
 	    s>>v.x; s>>v.y; s>>v.z;
            vertex.push_back(v);
        }
        else if (line.substr(0,2) == "vt")
        {
            std::istringstream s(line.substr(2));
            textureinfo t;
 	    s>>t.u; s>>t.v;
            texture.push_back(t);
        }
	else if (line.substr(0,2) == "vn")
        {
            std::istringstream s(line.substr(2));
            norminfo n;
 	    s>>n.i; s>>n.j; s>>n.k;
            norm.push_back(n);

        }
        else if (line.substr(0,2) == "f ")
        {
            std::istringstream s(line.substr(2));
            faceinfo f;






		do
		{
		std::string strsface;
		s>>strsface;
 		std::stringstream s1(strsface);

		std::string token;
		int counter=0;
			while (std::getline(s1, token,'/')){

			if (counter==0) f.vv.push_back(atoi(token.c_str())-1);
			if (counter==1) f.tt.push_back(atoi(token.c_str())-1);
			if (counter==2) f.nn.push_back(atoi(token.c_str())-1);
			counter++;

			}


 		} while (s);

	     	face.push_back(f);
		triangulation(f);

        }
        else
        {
            /* ignoring this line */
        }
    }


}


/***************************************************************************/

void init_window()
                 /* Clear the image area, and set up the coordinate system */
{

        					       /* Clear the window */
        glClearColor(0.0,0.0,0.0,0.0);
	    glShadeModel(GL_SMOOTH);
        glOrtho(0,WIDTH,0,HEIGHT,-1.0,1.0);
}

/***************************************************************************/

int Round(double number)               // Round function
{
        return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

int min(int a, int b)
{
        int x = a < b ? a:b;
        return (x);
}



void write_pixel(int x, int y, double intensity)//, double red, double green, double blue)
                                         /* Turn on the pixel found at x,y */
{

    glColor3f (red, green, blue);
    glBegin(GL_POINTS);
    glVertex3i( x, y, 0);
    glEnd();

}


void write_pixel_color(int x, int y, colorinfo c)//, double red, double green, double blue)
                                         /* Turn on the pixel found at x,y */
{

    glColor3f (c.red, c.green, c.blue);
    glBegin(GL_POINTS);
    glVertex3i( x, y, 0);
    glEnd();

}


/***************************************************************************/


/***************************************************************************/
void DDA(Point A, Point B)                //DDA algorithm following the updating concept from the lecture
{

    double xi;
    double yi;
    double zi;
    double x1=A.x;
    double x2=B.x;
    double y1=A.y;
    double y2=B.y;
    double z1=A.z;
    double z2=B.z;
    double deltaX=x2-x1;
    double deltaY=y2-y1;
    double deltaZ=z2-z1;
    double m, mz;
    double swapx, swapy, swapz;
    int xr, yr;




if(abs(deltaY)<=abs(deltaX))
{
	if (int(x1)!=int(x2))
	{
      if(x1>x2)
      {
          swapx = x1;
          x1 = x2;
          x2= swapx;

          swapy = y1;
          y1 = y2;
          y2 = swapy;

	  swapz = z1;
          z1 = z2;
          z2 = swapz;
      }

        xi=x1;
        yi=y1;
	zi=z1;
        deltaX = x2-x1;
        deltaY = y2-y1;
	deltaZ = z2-z1;
        m = deltaY / deltaX;
	mz = deltaZ / deltaX;

	xr=Round(xi);
	yr=Round(yi);
       while (xr<Round(x2))
       {



           if (zi<=zbuffer[xr][yr])
	   {
           	write_pixel(xr, yr, 1.0);
		zbuffer[xr][yr]=zi;

	   }
	           yi = yi + m;

           yr = Round(yi);
	   zi = zi + mz;

           xr++;
        }
    }




    else
    {

      if(y1>y2)
      {
          swapx = x1;
          x1 = x2;
          x2= swapx;

          swapy = y1;
          y1 = y2;
          y2 = swapy;

	  swapz = z1;
          z1 = z2;
          z2 = swapz;
      }

        xi=x1;
        yi=y1;
	zi=z1;
        deltaX = x2-x1;
        deltaY = y2-y1;
	deltaZ = z2-z1;

        m = deltaX / deltaY;
	mz = deltaZ / deltaY;
	yr=Round(yi);
	xr=Round(xi);
       while (yr<Round(y2))
      {


          if (zi<=zbuffer[xr][yr])
	  {
           	write_pixel(xr, yr, 1.0);
		zbuffer[xr][yr]=zi;
	   }
          xi = xi + m;
          xr = Round(xi);

          zi = zi + mz;
          yr++;
        }
    }

}

else
{

if(abs(deltaY)>abs(deltaX))
    {

      if(y1>y2)
      {
          swapx = x1;
          x1 = x2;
          x2= swapx;

          swapy = y1;
          y1 = y2;
          y2 = swapy;

	        swapz = z1;
          z1 = z2;
          z2 = swapz;
      }

        xi=x1;
        yi=y1;
	      zi=z1;
        deltaX = x2-x1;
        deltaY = y2-y1;
	      deltaZ = z2-z1;

        m = deltaX / deltaY;
	      mz = deltaZ / deltaY;
	      yr=Round(yi);
	      xr=Round(xi);
        while (yr<Round(y2))
      {

          if (zi<=zbuffer[xr][yr])
	  {
           	write_pixel(xr, yr, 1.0);
		        zbuffer[xr][yr]=zi;
	   }
	          xi = xi + m;
          xr = Round(xi);

          zi = zi + mz;
          yr++;
        }
    }
}


}





/* normalize */

object normalize(object obj)
{
	double xmax=obj.vertex[0].x;
	double xmin=obj.vertex[0].x;
	double ymax=obj.vertex[0].y;
	double ymin=obj.vertex[0].y;
	zmax=obj.vertex[0].z;
	zmin=obj.vertex[0].z;
	double xsum=0;
	double ysum=0;
	double zsum=0;

        for (unsigned int i=0;i<obj.vertex.size();i++)
	{
		xmax=xmax>obj.vertex[i].x?xmax:obj.vertex[i].x;
		ymax=ymax>obj.vertex[i].y?ymax:obj.vertex[i].y;
		zmax=zmax>obj.vertex[i].z?zmax:obj.vertex[i].z;
		xmin=xmin<obj.vertex[i].x?xmin:obj.vertex[i].x;
		ymin=ymin<obj.vertex[i].y?ymin:obj.vertex[i].y;
		zmin=zmin<obj.vertex[i].z?zmin:obj.vertex[i].z;
		xsum=xsum+obj.vertex[i].x;
		ysum=ysum+obj.vertex[i].y;
		zsum=zsum+obj.vertex[i].z;
	}

	double xorigin=xsum/obj.vertex.size();
	double yorigin=ysum/obj.vertex.size();
	zorigin=zsum/obj.vertex.size();

	double scale=300/((xmax-xmin)>(ymax-ymin)?(xmax-xmin):(ymax-ymin));
	initz=(abs(zmax)+abs(zmin))*100;
	for (unsigned int i=0;i<obj.vertex.size();i++)
	{
	obj.vertex[i].x=(obj.vertex[i].x-xorigin)*scale;
	obj.vertex[i].y=(obj.vertex[i].y-yorigin)*scale;
	obj.vertex[i].z=(obj.vertex[i].z-zorigin);

	}

	return obj;
}

double *matrixmul(double matrix[4][4], double pvector[4]) {
	static double newpvector[4]={};

    for (int i = 0; i < 4; i++)
	{
	newpvector[i]=0;
        for (int j = 0; j < 4; j++)
            newpvector[i]+= matrix[i][j] * pvector[j];

	}
    	return newpvector;

}

void drawTriangle(Point A, Point B, Point C)
{
	srand(seed+seedcounter);
	red=rand()/((double)RAND_MAX+1);
	green=rand()/((double)RAND_MAX+1);
	blue=rand()/((double)RAND_MAX+1);


	Point tmpP1, tmpP2, tmpP3, tmpA, tmpB;





if (wireframemode==1)
{
	      DDA(A, B);
        DDA(A, C);
	      DDA(B, C);
	printf("A.x=%f, A.y=%f, B.x=%f, B.y=%f, C.x=%f, C.y=%f\n", A.x, A.y, B.x, B.y, C.x, C.y);
}
else
{
	if (A.y<=B.y && A.y<=C.y)
	{
		tmpP1.setxyz(A.x, A.y, A.z);
		if (B.y<C.y)
		{
			tmpP2.setxyz(B.x, B.y, B.z);
			tmpP3.setxyz(C.x, C.y, C.z);
		}
		else
		{
			tmpP3.setxyz(B.x, B.y, B.z);
			tmpP2.setxyz(C.x, C.y, C.z);
		}
	}
	else if (B.y<=A.y && B.y<=C.y)
	{
		tmpP1.setxyz(B.x, B.y, B.z);
		if (A.y<C.y)
		{
			tmpP2.setxyz(A.x, A.y, A.z);
			tmpP3.setxyz(C.x, C.y, C.z);
		}
		else
		{
			tmpP3.setxyz(A.x, A.y, A.z);
			tmpP2.setxyz(C.x, C.y, C.z);
		}
	}
	else if (C.y<=B.y && C.y<=A.y)
	{
		tmpP1.setxyz(C.x, C.y, C.z);
		if (B.y<A.y)
		{
			tmpP2.setxyz(B.x, B.y, B.z);
			tmpP3.setxyz(A.x, A.y, A.z);
		}
		else
		{
			tmpP3.setxyz(B.x, B.y, B.z);
			tmpP2.setxyz(A.x, A.y, A.z);
		}
	}



	if (Round(tmpP1.y)!=Round(tmpP3.y))
	{
		if (Round(tmpP1.y)!=Round(tmpP2.y))
		{
			int scany=0;
			while (scany<Round(tmpP2.y)-Round(tmpP1.y))
			{
				tmpA.setxyz(tmpP1.x+scany*(tmpP2.x-tmpP1.x)/(tmpP2.y-tmpP1.y),Round(tmpP1.y)+scany+0.0,tmpP1.z+scany*(tmpP2.z-tmpP1.z)/(tmpP2.y-tmpP1.y));
				tmpB.setxyz(tmpP1.x+scany*(tmpP3.x-tmpP1.x)/(tmpP3.y-tmpP1.y),Round(tmpP1.y)+scany+0.0,tmpP1.z+scany*(tmpP3.z-tmpP1.z)/(tmpP3.y-tmpP1.y));
				DDA(tmpA, tmpB);
				scany=scany+1;

			}

			int scany1=0;
			while (scany<Round(tmpP3.y)-Round(tmpP1.y))
			{
				tmpA.setxyz(tmpP2.x+scany1*(tmpP3.x-tmpP2.x)/(tmpP3.y-tmpP2.y),Round(tmpP1.y)+scany+0.0,tmpP2.z+scany1*(tmpP3.z-tmpP2.z)/(tmpP3.y-tmpP2.y));
				tmpB.setxyz(tmpP1.x+scany*(tmpP3.x-tmpP1.x)/(tmpP3.y-tmpP1.y),Round(tmpP1.y)+scany+0.0,tmpP1.z+scany*(tmpP3.z-tmpP1.z)/(tmpP3.y-tmpP1.y));
				DDA(tmpA, tmpB);
				scany=scany+1;
				scany1=scany1+1;

			}


		}
		else
		{
			int scany2=0;
			while (scany2<int(tmpP3.y)-int(tmpP1.y))
			{
				tmpA.setxyz(tmpP2.x+scany2*(tmpP3.x-tmpP2.x)/(tmpP3.y-tmpP2.y),Round(tmpP2.y)+scany2+0.0,tmpP2.z+scany2*(tmpP3.z-tmpP2.z)/(tmpP3.y-tmpP2.y));
				tmpB.setxyz(tmpP1.x+scany2*(tmpP3.x-tmpP1.x)/(tmpP3.y-tmpP1.y),Round(tmpP1.y)+scany2+0.0,tmpP1.z+scany2*(tmpP3.z-tmpP1.z)/(tmpP3.y-tmpP1.y));
				DDA(tmpA, tmpB);

				scany2=scany2+1;
			}

		}
	}
	else
	{
	      DDA(A, B);
        DDA(A, C);
	      DDA(B, C);
	}
}
}

void drawPolygon(object obj)
{
	seedcounter=0;
	Point A, B, C;

	if (perspectivemode==1)
	{
	object tmpobj=obj;
	double perspectivematrix[4][4];
    	perspectivematrix[0][0] = 1.0;
    	perspectivematrix[1][1] = 1.0;
    	perspectivematrix[2][2] = 1.0;
    	perspectivematrix[3][2] = 1.0/(abs(zmax)+abs(zmin));


		for (unsigned int i=0;i<obj.vertex.size();i++)
		{

 		double pvector[4];
 		pvector[0]=obj.vertex[i].x;
		pvector[1]=obj.vertex[i].y;
		pvector[2]=obj.vertex[i].z-zmin+abs(zmax)+abs(zmin)+zorigin;
		pvector[3]=1.0;

		double *pvectorafterperspective=matrixmul(perspectivematrix, pvector);

		tmpobj.vertex[i].x=pvectorafterperspective[0]/pvectorafterperspective[3];
		tmpobj.vertex[i].y=pvectorafterperspective[1]/pvectorafterperspective[3];


		}
		for (unsigned int j=0; j<tmpobj.polygon.size();j++)
		{


		A.setxyz(tmpobj.vertex[tmpobj.polygon[j].v[0]].x+WIDTH/2.0, tmpobj.vertex[tmpobj.polygon[j].v[0]].y+HEIGHT/2.0, tmpobj.vertex[tmpobj.polygon[j].v[0]].z);
		B.setxyz(tmpobj.vertex[tmpobj.polygon[j].v[1]].x+WIDTH/2.0, tmpobj.vertex[tmpobj.polygon[j].v[1]].y+HEIGHT/2.0, tmpobj.vertex[tmpobj.polygon[j].v[1]].z);
		C.setxyz(tmpobj.vertex[tmpobj.polygon[j].v[2]].x+WIDTH/2.0, tmpobj.vertex[tmpobj.polygon[j].v[2]].y+HEIGHT/2.0, tmpobj.vertex[tmpobj.polygon[j].v[2]].z);
		drawTriangle(A, B, C);

		seedcounter++;
		}
	}
	else
	{
		for (unsigned int j=0; j<obj.polygon.size();j++)
		{


		A.setxyz(obj.vertex[obj.polygon[j].v[0]].x+WIDTH/2.0, obj.vertex[obj.polygon[j].v[0]].y+HEIGHT/2.0, obj.vertex[obj.polygon[j].v[0]].z);
		B.setxyz(obj.vertex[obj.polygon[j].v[1]].x+WIDTH/2.0, obj.vertex[obj.polygon[j].v[1]].y+HEIGHT/2.0, obj.vertex[obj.polygon[j].v[1]].z);
		C.setxyz(obj.vertex[obj.polygon[j].v[2]].x+WIDTH/2.0, obj.vertex[obj.polygon[j].v[2]].y+HEIGHT/2.0, obj.vertex[obj.polygon[j].v[2]].z);
		drawTriangle(A, B, C);
		seedcounter++;
		}
	}
}



object scale(object obj, double scalar)
{
	object newobj=obj;

	double scalematrix[4][4];
	scalematrix[0][0] = scalar;
    	scalematrix[1][1] = scalar;
    	scalematrix[2][2] = scalar;

    	scalematrix[3][3] = 1.0;
	int overflow=0;
        for (unsigned int i=0;i<obj.vertex.size();i++)
	{
 		double pvector[4];
 		pvector[0]=obj.vertex[i].x;
		pvector[1]=obj.vertex[i].y;
		pvector[2]=obj.vertex[i].z;
		pvector[3]=1.0;


		double *pvectorafterscale=matrixmul(scalematrix, pvector);

		newobj.vertex[i].x=pvectorafterscale[0];
		newobj.vertex[i].y=pvectorafterscale[1];
		newobj.vertex[i].z=pvectorafterscale[2];
		if (abs(pvectorafterscale[0])>WIDTH/2.0 || abs(pvectorafterscale[1])>HEIGHT/2.0)
			overflow=1;

	}

	if (overflow)
	{
		printf("Overflow!\n");
		return obj;
	}
	else
        return newobj;
}

object translation(object obj, double horizontalmove, double verticalmove)
{
	object newobj=obj;
	double translatematrix[4][4];
    	translatematrix[0][0] = 1.0;
    	translatematrix[1][1] = 1.0;
    	translatematrix[2][2] = 1.0;
    	translatematrix[3][3] = 1.0;
    	translatematrix[0][3] = horizontalmove;
    	translatematrix[1][3] = verticalmove;
	int overflow=0;
	for (unsigned int i=0;i<obj.vertex.size();i++)
	{
 		double pvector[4];
 		pvector[0]=obj.vertex[i].x;
		pvector[1]=obj.vertex[i].y;
		pvector[2]=obj.vertex[i].z;
		pvector[3]=1.0;

		double *pvectoraftertranslate=matrixmul(translatematrix, pvector);

		newobj.vertex[i].x=pvectoraftertranslate[0];
		newobj.vertex[i].y=pvectoraftertranslate[1];
		newobj.vertex[i].z=pvectoraftertranslate[2];
		if (abs(pvectoraftertranslate[0])>WIDTH/2.0 || abs(pvectoraftertranslate[1])>HEIGHT/2.0)
			overflow=1;

	}

	if (overflow)
	{
		printf("Overflow!\n");
		return obj;
	}
	else
        return newobj;
}

object rotate(object obj, double rotatex, double rotatey)
{
	object newobj=obj;
	double cosrotatex = cos(rotatex);
    	double sinrotatex = sin(rotatex);
	double cosrotatey = cos(rotatey);
    	double sinrotatey = sin(rotatey);
	double rotatematrix[4][4];
    	rotatematrix[0][0] = cosrotatey;
    	rotatematrix[1][1] = cosrotatex;
    	rotatematrix[2][2] = cosrotatey+cosrotatex-1.0;
    	rotatematrix[3][3] = 1.0;
    	rotatematrix[0][2] = sinrotatey;
    	rotatematrix[2][0] = -1.0*sinrotatey;
	rotatematrix[1][2] = -1.0*sinrotatex;
	rotatematrix[2][1] = sinrotatex;
	int overflow=0;
	for (unsigned int i=0;i<obj.vertex.size();i++)
	{
 		double pvector[4];
 		pvector[0]=obj.vertex[i].x;
		pvector[1]=obj.vertex[i].y;
		pvector[2]=obj.vertex[i].z;
		pvector[3]=1.0;

		double *pvectorafterrotate=matrixmul(rotatematrix, pvector);

		newobj.vertex[i].x=pvectorafterrotate[0];
		newobj.vertex[i].y=pvectorafterrotate[1];
		newobj.vertex[i].z=pvectorafterrotate[2];
		if (abs(pvectorafterrotate[0])>WIDTH/2.0 || abs(pvectorafterrotate[1])>HEIGHT/2.0)
			overflow=1;

	}

        if (overflow)
	{
		printf("Overflow!\n");
		return obj;
	}
	else
        return newobj;
}

/*************************************************************************/
//dot product function
double dotProduct(Point vector1, Point vector2) {
	double dpresult=vector1.x*vector2.x+vector1.y*vector2.y+vector1.z*vector2.z;
    	return dpresult;
}

double findmag(Point vector){
        double magnitude=sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
	return magnitude;
}

//vector normalization function

Point vectornormalize(Point vector){
	static Point nvector;
	double magnitude=findmag(vector);
	if (magnitude>0)
		nvector.setxyz(vector.x/magnitude, vector.y/magnitude, vector.z/magnitude);
	return nvector;
}

//intersection
Point internormal(int polygonid, Point interpoint) {
    Point v0v1;
    Point v0v2;
    Point v0P;

    v0v1.setxyz(obj.vertex[obj.polygon[polygonid].v[1]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[1]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[1]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    v0v2.setxyz(obj.vertex[obj.polygon[polygonid].v[2]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[2]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[2]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    v0P.setxyz(interpoint.x - obj.vertex[obj.polygon[polygonid].v[0]].x, interpoint.y - obj.vertex[obj.polygon[polygonid].v[0]].y, interpoint.z - obj.vertex[obj.polygon[polygonid].v[0]].z);

    Point v2CrossP = crossProduct(v0v2, v0P);
    Point v1CrossP = crossProduct(v0v1, v0P);
    Point v1CrossV2 = crossProduct(v0v1, v0v2);

    double u = findmag(v2CrossP) / findmag(v1CrossV2);
    double v = findmag(v1CrossP) / findmag(v1CrossV2);
    Point n;
    n.setxyz(u * obj.vertex[obj.polygon[polygonid].v[1]].x + v * obj.vertex[obj.polygon[polygonid].v[2]].x + (1 - u - v) * obj.vertex[obj.polygon[polygonid].v[0]].x, u * obj.vertex[obj.polygon[polygonid].v[1]].y + v * obj.vertex[obj.polygon[polygonid].v[2]].y + (1 - u - v) * obj.vertex[obj.polygon[polygonid].v[0]].y, u * obj.vertex[obj.polygon[polygonid].v[1]].z + v * obj.vertex[obj.polygon[polygonid].v[2]].z + (1 - u - v) * obj.vertex[obj.polygon[polygonid].v[0]].z);
    n=vectornormalize(n);
    return n;
}

//intersection
bool intersect(int polygonid, Point rayorig, Point raydir) {
    Point v0v1;
    Point v0v2;
    v0v1.setxyz(obj.vertex[obj.polygon[polygonid].v[1]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[1]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[1]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    v0v2.setxyz(obj.vertex[obj.polygon[polygonid].v[2]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[2]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[2]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    Point P = crossProduct(raydir, v0v2);
    double det = dotProduct(v0v1, P);

    Point T;
    if (det>0)
    {
	T.setxyz(rayorig.x - obj.vertex[obj.polygon[polygonid].v[0]].x, rayorig.y - obj.vertex[obj.polygon[polygonid].v[0]].y, rayorig.z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    }
    else
    {
	T.setxyz(-1.0*rayorig.x + obj.vertex[obj.polygon[polygonid].v[0]].x, -1.0*rayorig.y + obj.vertex[obj.polygon[polygonid].v[0]].y, -1.0*rayorig.z + obj.vertex[obj.polygon[polygonid].v[0]].z);
	det=-1.0*det;
    }
    if (det < 0.0001) return false;

    double u = dotProduct(T, P);
    if (u < 0.0 || u > det) return false;

    Point Q = crossProduct(T, v0v1);
    double v = dotProduct(raydir, Q);
    if (v < 0.0 || u + v > det) return false;


    return true;
}

double tafterintersect(int polygonid, Point rayorig, Point raydir) {
    Point v0v1;
    Point v0v2;
    v0v1.setxyz(obj.vertex[obj.polygon[polygonid].v[1]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[1]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[1]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    v0v2.setxyz(obj.vertex[obj.polygon[polygonid].v[2]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[2]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[2]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    Point P = crossProduct(raydir, v0v2);
    double det = dotProduct(v0v1, P);
    //printf("v0v2=%f, %f, %f\n", v0v2.x, v0v2.y, v0v2.z);
    Point T;
    if (det>0)
    {
	T.setxyz(rayorig.x - obj.vertex[obj.polygon[polygonid].v[0]].x, rayorig.y - obj.vertex[obj.polygon[polygonid].v[0]].y, rayorig.z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    }
    else
    {
	T.setxyz(-1.0*rayorig.x + obj.vertex[obj.polygon[polygonid].v[0]].x, -1.0*rayorig.y + obj.vertex[obj.polygon[polygonid].v[0]].y, -1.0*rayorig.z + obj.vertex[obj.polygon[polygonid].v[0]].z);
	det=-1.0*det;
    }
    Point Q = crossProduct(T, v0v1);



    double invDet = 1.0 / det;

    double t = dotProduct(v0v2, Q)*invDet;

    return t;
}

//bump mapping function
Point bumpfunc(Point beforebump, int polygonid, Point UA, Point UB){
	double Ha=0.2; //magnitude of the dimple
	double Hb=15.0;	//size of the dimple
	double realx, realy;
	Point perturb;
	perturb.setxyz(0.0,0.0,0.0);
	double centerx = (obj.vertex[obj.polygon[polygonid].v[0]].x+obj.vertex[obj.polygon[polygonid].v[1]].x+obj.vertex[obj.polygon[polygonid].v[2]].x)/3.0;
	double centery = (obj.vertex[obj.polygon[polygonid].v[0]].y+obj.vertex[obj.polygon[polygonid].v[1]].y+obj.vertex[obj.polygon[polygonid].v[2]].y)/3.0;
	//perturb.setxyz(sin(Hb*beforebump.x)+1.1, sin(Hb*beforebump.y)+1.1, 0.0);
	if (abs(UB.x)>0.0001)
		realx=(beforebump.x-centerx)*sqrt(UB.x*UB.x+UB.z*UB.z)/abs(UB.x);
	else
		realx=(beforebump.x-centerx);
	if (abs(UA.y)>0.0001)
		realy=(beforebump.y-centery)*sqrt(UA.y*UA.y+UA.z*UA.z)/abs(UA.y);
	else
		realy=(beforebump.y-centery);
	double Rdist=sqrt(realx*realx+realy*realy);
	if (Rdist<Hb)
		perturb.setxyz(-Ha*sin(PI*Rdist*0.5/Hb)*realx/Rdist, -Ha*sin(PI*Rdist*0.5/Hb)*realy/Rdist, 0.0);

	return perturb;

}

colorinfo phong(int polygonid, Point interpoint, Point rayorig, lightinfo curlight) {
    Point v0v1;
    Point v0v2;
    v0v1.setxyz(obj.vertex[obj.polygon[polygonid].v[1]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[1]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[1]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    v0v2.setxyz(obj.vertex[obj.polygon[polygonid].v[2]].x - obj.vertex[obj.polygon[polygonid].v[0]].x, obj.vertex[obj.polygon[polygonid].v[2]].y - obj.vertex[obj.polygon[polygonid].v[0]].y, obj.vertex[obj.polygon[polygonid].v[2]].z - obj.vertex[obj.polygon[polygonid].v[0]].z);
    Point norm;
    if (phongSmooth) {
        norm = internormal(polygonid, interpoint);

    }
    else {
        norm = crossProduct(v0v1, v0v2);
    }

    norm = vectornormalize(norm);
    if (bumpMapping)
    {


        Point UA, UB;
	Point xplane, yplane;
	xplane.setxyz(1.0, 0.0, 0.0);
	yplane.setxyz(0.0, 1.0, 0.0);
	UA=crossProduct(norm, xplane);
        UB=crossProduct(norm, yplane);

        Point perturb=bumpfunc(interpoint, polygonid, UA, UB);
	norm.setxyz(norm.x+perturb.x, norm.y+perturb.y, norm.z);

    }

    norm = vectornormalize(norm);
    Point l, r, v;
    l.setxyz(curlight.lightposition.x - interpoint.x, curlight.lightposition.y - interpoint.y, curlight.lightposition.z - interpoint.z);
    l=vectornormalize(l);
    r.setxyz(2.0*dotProduct(l, norm)*norm.x - l.x, 2.0*dotProduct(l, norm)*norm.y - l.y, 2.0*dotProduct(l, norm)*norm.z - l.z);
    r=vectornormalize(r);
    v.setxyz(rayorig.x-interpoint.x, rayorig.y-interpoint.y, rayorig.z-interpoint.z);
    v=vectornormalize(v);

    colorinfo color;
    color.red = ka * curlight.La * curlight.red;
    color.green = ka * curlight.La * curlight.green;
    color.blue = ka * curlight.La * curlight.blue;

    double lDotN = dotProduct(l, norm);
    //printf("lDotN=%f",lDotN);
    if (lDotN < 0) {
        lDotN = 0;
    }
    color.red += kd * curlight.Ld * curlight.red * lDotN;
    color.green += kd * curlight.Ld * curlight.green * lDotN;
    color.blue += kd * curlight.Ld * curlight.blue * lDotN;

    double rDotV = dotProduct(r, v);
    if (rDotV < 0) {
        rDotV = 0;
    }
    color.red += ks * curlight.Ls * curlight.red * pow(rDotV, shininess);
    color.green += ks * curlight.Ls * curlight.green * pow(rDotV, shininess);
    color.blue += ks * curlight.Ls * curlight.blue * pow(rDotV, shininess);

    if(color.red<0.0) color.red=0.0;
    if(color.red>1.0) color.red=1.0;
    if(color.green<0.0) color.green=0.0;
    if(color.green>1.0) color.green=1.0;
    if(color.blue<0.0) color.blue=0.0;
    if(color.blue>1.0) color.blue=1.0;

    return color;

}

colorinfo raytrace(Point rayorig, Point raydir) {
    double distance = dinf; // closest intersection, set to INFINITY to start with
    unsigned int interpolygonid = obj.polygon.size()+1;

    Point interpoint;
    double t=dinf;
    for (unsigned int j=0; j<obj.polygon.size();j++) {
        bool isintersect = intersect(j, rayorig, raydir);

        if (isintersect) {
	t = tafterintersect(j, rayorig, raydir);

            if (t < distance) {

                distance = t;
                interpolygonid = j;

                interpoint.setxyz(rayorig.x + t * raydir.x, rayorig.y + t * raydir.y, rayorig.z + t * raydir.z);

            }
        }
    }

    colorinfo color, phongcolor;
    color.red = 0.0;
    color.green = 0.0;
    color.blue = 0.0;

    if(t<dinf)
    {
        for (unsigned int i = 0; i < lightsource.lightpoint.size(); i++) {
            bool notShadow = true;
	    Point shadowraydir;
            shadowraydir.setxyz(interpoint.x - lightsource.lightpoint[i].lightposition.x, interpoint.y - lightsource.lightpoint[i].lightposition.y, interpoint.z - lightsource.lightpoint[i].lightposition.z);

              double shadowdistance = dinf; // closest intersection, set to INFINITY to start with
    		unsigned int shadowpolygonid = obj.polygon.size()+1;


    		double shadowt=dinf;

            for (unsigned int j=0; j<obj.polygon.size();j++) {
        	bool shadowintersect = intersect(j, lightsource.lightpoint[i].lightposition, shadowraydir);

        	if (shadowintersect) {
			shadowt = tafterintersect(j, lightsource.lightpoint[i].lightposition, shadowraydir);

           	 	if (shadowt < shadowdistance) {

                		shadowdistance = shadowt;
                		shadowpolygonid = j;

                	}
            	}
            }

	    if (shadowpolygonid!=interpolygonid){
		notShadow = false;
	    }

            if (notShadow) {

                phongcolor = phong(interpolygonid, interpoint, rayorig, lightsource.lightpoint[i]);
                color.red += phongcolor.red;
                color.green += phongcolor.green;
                color.blue += phongcolor.blue;
            }
            else {
                color.red += ka * lightsource.lightpoint[i].La * lightsource.lightpoint[i].red;
                color.green += ka * lightsource.lightpoint[i].La * lightsource.lightpoint[i].green;
                color.blue += ka * lightsource.lightpoint[i].La * lightsource.lightpoint[i].blue;
            }
        }

    }

    return color;
}

void render() {
    Point rayorig;
    rayorig.setxyz(0.0, 0.0, 1000.0);

    for (int y = -HEIGHT/2; y <= HEIGHT/2; y++) {
        for (int x = -WIDTH/2; x <= WIDTH/2; x++) {
            Point raydir;
	    raydir.setxyz(x-rayorig.x, y-rayorig.y, 0.0-rayorig.z);
            colorinfo color;
            color = raytrace(rayorig, raydir);

            write_pixel_color(x+WIDTH/2, y+HEIGHT/2, color);

        }
    }
}

/*************************************************************************/

void display ( void )   // Create The Display Function
{
         //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if(drawreg==1)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		render();
		drawreg=0;
   	}
  //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen
        //if (drawreg==1)
	//{
	//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//	for (int i=0; i<WIDTH; i++)
	//	for (int j=0; j<HEIGHT; j++)
	//		zbuffer[i][j]=initz;
	//printf("%f", zbuffer[499][499]);

        //	drawPolygon(obj);
	//	drawreg=0;
	//}
//printf("display zbuffer = %f, x=%f", zbuffer[300][300], obj.vertex[0].z);
	glutSwapBuffers();                                        // Draw Frame Buffer
}

/***************************************************************************/
void mouse(int button, int state, int x, int y)
{
/* This function I finessed a bit, the value of the printed x,y should
   match the screen, also it remembers where the old value was to avoid multiple
   readings from the same mouse click.  This can cause problems when trying to
   start a line or curve where the last one ended */


}

/***************************************************************************/
void keyboard ( unsigned char key, int x, int y )  // Create Keyboard Function
{

	switch ( key ) {
		case 27:              // When Escape Is Pressed...
			exit ( 0 );       // Exit The Program
			break;
	        case '1':         // stub for new screen
		        printf("New screen\n");
            		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	      // Clear Screen

            		//obj.objloader("model2.obj");
			//obj=normalize(obj);
			obj=objbackup;
			srand(time(NULL));
			seed=rand();
			modereg=0;
			wireframemode=0;
			perspectivemode=0;
			drawreg=1;
			break;
		case 'x':
			if (phongSmooth==0)
			{
				printf("Phong Smoothing Render\n");
				phongSmooth=1;
			}
			else
			{
				printf("Flat Shaded Render\n");
				phongSmooth=0;
			}

			//glutSwapBuffers();
			drawreg=1;
			break;
		case 'b':
			if (bumpMapping==0)
			{
				printf("Bump Mapping with a Dimple in the Center of Each Triangle\n");
				bumpMapping=1;
			}
			else
			{
				printf("Normal Mapping\n");
				bumpMapping=0;
			}

			//glutSwapBuffers();
			drawreg=1;
			break;
		case 'z':
			if (wireframemode==0)
			{
				printf("Wireframe mode\n");
				wireframemode=1;
			}
			else
			{
				printf("Z-buffer mode\n");
				wireframemode=0;
			}

			//glutSwapBuffers();
			drawreg=1;
			break;
		case 'v':
			if (perspectivemode==0)
			{
				printf("Perspective mode\n");
				perspectivemode=1;
			}
			else
			{
				printf("Orthogonal mode\n");
				perspectivemode=0;
			}

			//glutSwapBuffers();
			drawreg=1;
			break;
		case 'e':
			printf("Scale mode\n");
			modereg=1;
			break;
		case 't':
			printf("Translation mode\n");
			modereg=2;
			break;
		case 'r':
			printf("Rotation mode\n");
			modereg=3;
			break;
		case 'a':
			switch (modereg)
				{
				case 1:
				printf("Scale Up\n");
				obj=scale(obj,1.25);
            			drawreg=1;
				break;
				case 2:
				printf("Move Left\n");
				obj=translation(obj,-10.0,0.0);
            			drawreg=1;
				break;
				case 3:
				printf("Rotate about X, Clockwise\n");
				obj=rotate(obj,-0.1,0.0);
            			drawreg=1;
				break;
				}
			break;
		case 'd':
			switch (modereg)
				{
				case 1:
				printf("Scale Up\n");
				obj=scale(obj,1.25);
            			drawreg=1;
				break;
				case 2:
				printf("Move Right\n");
				obj=translation(obj,10.0,0.0);
            			drawreg=1;
				break;
				case 3:
				printf("Rotate about X, Counterclockwise\n");
				obj=rotate(obj,0.1,0.0);
            			drawreg=1;
				break;
				}
			break;
		case 's':
			switch (modereg)
				{
				case 1:
				printf("Scale Down\n");
				obj=scale(obj,0.8);
            			drawreg=1;
				break;
				case 2:
				printf("Move Down\n");
				obj=translation(obj,0.0,-10.0);
            			drawreg=1;
				break;
				case 3:
				printf("Rotate about Y, Clockwise\n");
				obj=rotate(obj,0.0,-0.1);
            			drawreg=1;
				break;
				}
			break;
		case 'w':
			switch (modereg)
				{
				case 1:
				printf("Scale Down\n");
				obj=scale(obj,0.8);
            			drawreg=1;
				break;
				case 2:
				printf("Move Up\n");
				obj=translation(obj,0.0,10.0);
            			drawreg=1;
				break;
				case 3:
				printf("Rotate about Y, Counterclockwise\n");
				obj=rotate(obj,0.0,0.1);
            			drawreg=1;
				break;
				}
			break;

		default:
            break;
	}
}
/***************************************************************************/


int main (int argc, char *argv[])
{
/* This main function sets up the main loop of the program and continues the
   loop until the end of the data is reached.  Then the window can be closed
   using the escape key.

						  */



	glutInit            ( &argc, argv );
    	glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize  ( WIDTH,HEIGHT );
	glutCreateWindow    ( "Computer Graphics" );


    	if (argc < 2)
        	obj.objloader("model1.obj");
    	else
        	obj.objloader(argv[1]);

        init_window();				             //create_window

	lightinfo l1, l2;
	l1.lightposition.setxyz(1000.0, 0.0, 2000.0);
	l1.red=0.5;
	l1.green=0.0;
	l1.blue=0.8;
  l1.La=0.5;
	l1.Ld=0.5;
	l1.Ls=0.5;
	l2.red=1.0;
	l2.green=0.6;
	l2.blue=0.0;
  l2.La=0.5;
	l2.Ld=0.5;
	l2.Ls=0.5;
	l2.lightposition.setxyz(-1000.0, -0.0, 2000.0);
	lightsource.lightpoint.push_back(l1);
	lightsource.lightpoint.push_back(l2);
	obj=normalize(obj);
	objbackup=obj;

        drawreg=1;

	glutDisplayFunc     ( display );
	glutIdleFunc	    ( display );
	glutMouseFunc       ( mouse );
	glutKeyboardFunc    ( keyboard );


    	glutMainLoop        ( );                // Initialize The Main Loop


}
