//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

struct Light
{
	string name;
	vec4 color;
	vec4 position;
};

struct Sphere
{
	string name;
	vec4 position;
	vec3 scale;
	vec4 color;
	float ka; //ambient reflection
	float kd; //diffuse reflection
	float ks; //specular reflection
	float kr; //reflection
	float n; //specularity
	mat4 inverse_transformed;
	mat4 inverse_squared;
};

const int max_spheres = 10;
int sphere_number = 0;
Sphere spheres[max_spheres];

const int max_lights = 10;
int light_number = 0;
Light lights[max_lights];

vector<vec4> g_colors;

int g_width;
int g_height;

float g_near;
float g_left;
float g_right;
float g_bottom;
float g_top;

vec4 g_bg_color;
vec4 g_ambient_color;

string g_output_file;

// TODO: add structs for spheres, lights and anything else you may need.

// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

vec3 toVec3(const string& s1, const string& s2, const string& s3)
{
	stringstream ss(s1 + " " + s2 + " " + s3);
	vec3 result;
	ss >> result.x >> result.y >> result.z;
	//result.w = 1.0f;
	return result;
}

vec3 vec4to3(vec4 v4){
	vec3 v3;
	v3.x = v4.x;
	v3.y = v4.y;
	v3.z = v4.z;
	return v3;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
	if (vs[0] == "NEAR")
	{
		g_near = toFloat(vs[1]);
	}

	else if (vs[0] == "LEFT")
	{
		g_left = toFloat(vs[1]);
	}

	else if (vs[0] == "RIGHT")
	{
		g_right = toFloat(vs[1]);
	}

	else if (vs[0] == "BOTTOM")
	{
		g_bottom = toFloat(vs[1]);
	}

	else if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }

	else if (vs[0] == "SPHERE")
	{
		spheres[sphere_number].name = vs[1];
		spheres[sphere_number].position = toVec4(vs[2], vs[3], vs[4]);
		spheres[sphere_number].scale = toVec3(vs[5], vs[6], vs[7]);
		spheres[sphere_number].color = toVec4(vs[8], vs[9], vs[10]);
		spheres[sphere_number].ka = toFloat(vs[11]);
		spheres[sphere_number].kd = toFloat(vs[12]);
		spheres[sphere_number].ks = toFloat(vs[13]);
		spheres[sphere_number].kr = toFloat(vs[14]);
		spheres[sphere_number].n = toFloat(vs[15]);
		sphere_number++;
	}

	else if (vs[0] == "LIGHT")
	{
		lights[light_number].name = vs[1];
		lights[light_number].position = toVec4(vs[2], vs[3], vs[4]);
		lights[light_number].color = toVec4(vs[5], vs[6], vs[7]);
		light_number++;
	}

	else if (vs[0] == "BACK")
	{
		g_bg_color = toVec4(vs[1], vs[2], vs[3]);
	}

	else if (vs[0] == "AMBIENT")
	{
		g_ambient_color = toVec4(vs[1], vs[2], vs[3]);
	}

	else if (vs[0] == "OUTPUT")
	{
		g_output_file = vs[1];
	}

}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

bool intersect(Sphere& hit_sphere, Ray& ray, vec4& intersection_point, float& closest_t, float min_t, float max_t)
{
	//for each sphere, see if the ray intersects the sphere
	//return the closest intersection position and sphere
	bool hit = false;
	float temp_t = 0;

	for (int i = 0; i < sphere_number; i++){

		vec3 S = vec4to3(spheres[i].inverse_transformed * ray.origin);
		vec3 C = vec4to3(normalize(spheres[i].inverse_transformed * ray.dir));

		float a = dot(C, C);
		float b = 2 * dot(S, C);
		float c = dot(S, S) - 1;

		float discrimant = b * b - 4 * a * c;
		
		if (discrimant == 0){
  			temp_t = -b / (2 * a);

			if ((temp_t < closest_t) && (temp_t > min_t) && (temp_t < max_t)){
				hit = true;
				closest_t = temp_t;
				hit_sphere = spheres[i];
				intersection_point = ray.origin + closest_t * ray.dir;
			}
		} else if (discrimant > 0){

			float root1 = (-b - sqrtf(discrimant)) / (2 * a);
			float root2 = (-b + sqrtf(discrimant)) / (2 * a);

			if (root1 > root2){
				temp_t = root2;
				if ((temp_t < closest_t) && (temp_t > min_t) && (temp_t < max_t)){
					hit = true;
					closest_t = temp_t;
					hit_sphere = spheres[i];
					intersection_point = ray.origin + closest_t * ray.dir;
				}
			} else {
				temp_t = root1;
				if ((temp_t < closest_t) && (temp_t > min_t) && (temp_t < max_t)){
					hit = true;
					closest_t = temp_t;
					hit_sphere = spheres[i];
					intersection_point = ray.origin + closest_t * ray.dir;
				}
			}
		}
	}

	return hit;
}

vec4 trace(Ray& ray, int step)
{
	vec4 pixel_color;
	bool intersection;
	vec4 intersection_point;
	Sphere hit_sphere;
	float closest_t = 1000.0f;
	float min_t = 1;
	float max_t = 1000;
	vec4 N;
	vec4 L;
	vec4 V;
	vec4 R;

	intersection = intersect(hit_sphere, ray, intersection_point, closest_t, min_t, max_t);

	if (!intersection && step == 0) {
		return g_bg_color;
	} 

	else if (!intersection && step > 0){
		return (0, 0, 0, 0);
	}
	
	else {
		pixel_color += hit_sphere.ka * hit_sphere.color * g_ambient_color; // ambient contribution
		
		for (int i = 0; i < light_number; i++){ //check if lights visible from intersection point
			Ray s_ray;
			Sphere s_sphere;
			vec4 s_intersection_point;
			float s_closest_t = 1000;
			float s_min_t = .00001;
			float s_max_t = 1000;

			s_ray.origin = intersection_point;
			s_ray.dir = normalize(lights[i].position - intersection_point);
	
			bool shadow = intersect(s_sphere, s_ray, s_intersection_point, s_closest_t, s_min_t, s_max_t);

			L = normalize(lights[i].position - intersection_point);
			N = normalize(hit_sphere.inverse_squared * (intersection_point - hit_sphere.position));
			V = normalize(ray.origin - intersection_point);
			R = 2 * dot(N, L) * N - L;

			float N_L = dot(N, L);
			if (N_L < 0){
				N_L = 0.0f;
			}

			float R_V = dot(R, V);
			if (R_V < 0){
				R_V = 0.0f;
			}

			if ( !shadow ){ //if shadow ray doesn't hit anything, light visible
				pixel_color += hit_sphere.kd * lights[i].color * N_L * hit_sphere.color; //diffuse contribution
				pixel_color += hit_sphere.ks * lights[i].color * pow(R_V, hit_sphere.n); //specular contribution
			}
		}

		if ( (hit_sphere.kr > 0) && (step < 3) ){
			Ray r_ray;
			r_ray.origin = intersection_point;
			r_ray.dir = normalize(ray.dir - (2 * dot(ray.dir, N) * N));

			vec4 reflection = trace(r_ray, step + 1) * hit_sphere.kr;
		//	cout << "reflection: " << reflection << endl;
			pixel_color += reflection;
		}
	}

    return pixel_color;
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
	
	float x = (g_right - g_left) * ((float)ix / g_width - .5f);
	float y = (g_top - g_bottom) * ((float)iy / g_height - .5f) * 2;

    vec4 dir = vec4(x, y, -g_near, 0.0f);
	dir = normalize(dir);
    
	return dir;
}

void renderPixel(int ix, int iy) //sets a color for a pixel
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f); //eye
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray, 0);
    setColor(ix, iy, color);
}

void setMatrices()
{
	for (int i = 0; i < sphere_number; i++){
		mat4 translate = Translate(spheres[i].position);
		mat4 scale = Scale(spheres[i].scale);
		mat4 transformed_sphere = translate * scale;
		mat4 inverse_squared = scale * scale;
		InvertMatrix(transformed_sphere, spheres[i].inverse_transformed);
		InvertMatrix(inverse_squared, spheres[i].inverse_squared);
	}
}

void render() // iterates through all pixels in the view area
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
			for (int i = 0; i < 3; i++){
				float value = ((float*)g_colors[y*g_width + x])[i];
				if (value > 1.0)
					value = 1.0;
				buf[y*g_width * 3 + x * 3 + i] = (unsigned char)(value * 255.9f);
			}

    // TODO: change file name based on input file name.
	char *file_name = new char[g_output_file.length() + 1];
	strcpy(file_name, g_output_file.c_str());

	savePPM(g_width, g_height, file_name, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
	setMatrices();
    render();
	saveFile();
	return 0;
}

