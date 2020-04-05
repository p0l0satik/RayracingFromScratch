#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// #include "geometry.h"
#define EPS    0.001
#define C_EPS  0.00001
#define REQ    4
#define INF    1000000000

struct Primitive {
    virtual bool ray_intersect(glm::vec3 &orig, glm::vec3 &dir, glm::vec3 &cross_point) = 0;
    virtual glm::vec3 get_normal(glm::vec3& point) = 0;

};

struct Sphere {
    glm::vec3 center;
    float radius;
    glm::vec3 color;
    float blink_density;
    float mirror ;
    Sphere(const glm::vec3 &c, const float &r, const glm::vec3 &col, float bd, float m) 
    : center(c), radius(r), color(col), blink_density(bd), mirror(m) {}

    Sphere(){
        center = glm::vec3(0, 0, 0);
        radius = 0;
        blink_density = 0;
        mirror = 0;
        color = glm::vec3(0, 0, 0);
    };

    bool ray_intersect(glm::vec3 &orig, glm::vec3 &dir, glm::vec3 &cross_point) {
        float cross_dist;
        glm::vec3 orig_center = center - orig;
        if (glm::dot(orig_center, dir) <= 0) return false;

        glm::vec3 center_ray_proj = orig + dir * (glm::dot(orig_center, dir) / glm::length(dir));

        float center_dist = glm::distance(center_ray_proj, center);
        if (center_dist > radius) return false;
        if (center_dist < radius){
            float incide_part = sqrt(radius * radius - center_dist * center_dist);
            cross_dist = glm::distance(center_ray_proj, orig) - incide_part;
            cross_point = orig + dir * cross_dist;
        } else {
            cross_point = center_ray_proj;

        }
        return true;
    }

    glm::vec3 get_normal(glm::vec3& point){
        return glm::normalize(point - center);
    }
};

struct Triangle {
    glm::vec3 v0, v1, v2;
    glm::vec3 color;
    glm::vec3 N;
    Triangle(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 c) : v0(p0), v1(p1), v2(p2), color(c){
        N = glm::normalize(glm::cross(p1 - p0, p2 - p0));
    }
    Triangle() {}

    bool ray_intersect(glm::vec3 &orig, glm::vec3 &dir, glm::vec3 &cross_point) {
        glm::vec3 E1 = v1 - v0;
        glm::vec3 E2 = v2 - v0;
        glm::vec3 T  = orig - v0;
        glm::vec3 P  = glm::cross(dir, E2);
        glm::vec3 Q  = glm::cross(T, E1);
        float betha  = glm::dot(P, E1);
        
        if (abs(betha) < C_EPS) return false;

        float alpha = 1 / betha;

        float u = alpha * glm::dot(P, T);
        if (u < 0.0 || u > 1.0) return false;


        float v = alpha * glm::dot(Q, dir);
        if (v < 0.0 || u + v > 1.0) return false;

        float t = alpha * glm::dot(Q, E2);

        if (t > C_EPS){
            cross_point = orig + dir * t;
            return true;
        }
        return false;
    }

    glm::vec3 get_normal(glm::vec3& point){
        return N;
    }
};

struct Rect {
    glm::vec3 v0, v1, v2, v3; //vertexes clocwise
    Triangle t1, t2;
    glm::vec3 N;
    glm::vec3 color;
    Rect() {}
    Rect(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 c) {
        t1 = Triangle(p0, p1, p2, c);
        t2 = Triangle(p0, p3, p2, c);
        N = glm::normalize(glm::cross(p1 - p0, p2 - p0));
        color = c;
    }

    glm::vec3 get_normal(glm::vec3& point){
        return N;
    }

};

struct Light {
    Light(const glm::vec3 &p, const float &i) : position(p), intensity(i) {}
    glm::vec3 position;
    float intensity;
};

struct Objects {
    std::vector<Sphere> spheres; 
    std::vector<Triangle> triangles;
    std::vector<Rect> rects;
};


glm::vec3 cast_ray(glm::vec3 &orig,  
                   glm::vec3 &dir, 
                   Objects& objects,
                   std::vector<Light>& lights, 
                   int req = REQ) {
    // float sphere_dist = std::numeric_limits<float>::max();
    if (!req) return glm::vec3(0, 0, 0);

    Sphere nearest_one;
    Triangle n_tr;
    float min_dist = INF;
    int intersect = 0;
    for (auto sphere : objects.spheres){
        glm::vec3 point;
        if (sphere.ray_intersect(orig, dir, point)){
            intersect = 1;
            if (glm::distance(orig, point) < min_dist) {
                min_dist = glm::distance(orig, point);
                nearest_one = sphere;
            }
        } 
    }
    for (auto triangle : objects.triangles) {
        glm::vec3 point;
        if (triangle.ray_intersect(orig, dir, point)){
            if (glm::distance(orig, point) < min_dist) {
                min_dist = glm::distance(orig, point);
                n_tr = triangle;
                intersect = 2;
            }
        }
    }
    if (intersect == 1){

        glm::vec3 point;
        nearest_one.ray_intersect(orig, dir, point);
        glm::vec3 N = nearest_one.get_normal(point);
        float diff_light_intens = 0;
        float blinks = 0;
        for(auto light : lights){
            glm::vec3 light_dir = glm::normalize(light.position - point);
            bool use_light = true;
            for (auto cross_sp : objects.spheres) {
                glm::vec3 rev_point;
                glm::vec3 moved_point = point * float(1 + EPS) ;
                glm::vec3 rev_light_dir = light_dir ;
                if (cross_sp.ray_intersect(moved_point, rev_light_dir, rev_point)) {
                    if (glm::distance(rev_point, point) < glm::distance(light.position, point)){
                        use_light = false;
                        break;
                    }
                }
            }
            if (!use_light) continue;
            diff_light_intens += light.intensity * std::max(0.f, glm::dot(light_dir, N));
            glm::vec3 reflect_dir = glm::reflect(light_dir * float(-1), N);
            blinks += light.intensity * std::max(0.0, pow(glm::dot(glm::normalize(reflect_dir), (dir * float(-1))),
                                                    nearest_one.blink_density));
            
        }
        glm::vec3 reflect_dir = glm::normalize(N * float(2.0) + dir);
        glm::vec3 mirror_color = cast_ray(point, reflect_dir, objects, lights, req - 1);
        return nearest_one.color * (diff_light_intens + blinks) * (1 - nearest_one.mirror) + 
        mirror_color * nearest_one.mirror;
    }
    if (intersect == 2){
        return n_tr.color;      
    }
    return glm::vec3(0.2, 0.7, 0.8);
}

template<typename T>
void get_color(T& obj){
    std::cout << obj.color.x << std::endl;
}



glm::vec3 red(1, 0, 0), green(1, 1, 1), blue(0, 0, 1);
glm::vec3 yellow(1, 1, 0);

void render() {
    const int width    = 1024;
    const int height   = 768;
    std::vector<glm::vec3> framebuffer(width*height);
    
    float fov = M_PI / 2;
    Objects objects;
    objects.spheres= { Sphere(glm::vec3(0, 4, -45 ), 9, red, 60, 0), 
                                   Sphere(glm::vec3(-20, 4, -50), 9, blue, 700, 0),
                                   Sphere(glm::vec3(-10, -10, -30), 9, green, 1000, 0.9),
                                   Sphere(glm::vec3(20, 20, -35), 12, green, 1000, 0.8)};
    std::vector<Light> lights = {Light(glm::vec3(40, 20, 0), 0.7), Light(glm::vec3(0, 0, 0), 0.0),
     Light(glm::vec3(-10, -40, -20), 0)};
    objects.triangles = {Triangle(glm::vec3(-40, 0, -40), 
                                                glm::vec3(-30, 0, -40),
                                                glm::vec3(-30, -10, -40), red),
                                    Triangle(glm::vec3(-40, 0, -40), 
                                             glm::vec3(-40, -10, -40),
                                             glm::vec3(-30, -10, -40), red)};

    // get_color(triangles[0]);
    // get_color(spheres[0]);

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)* width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            glm::vec3 dir = glm::normalize(glm::vec3(x, y, -1));
            glm::vec3 orig = glm::vec3(0, 0, 0);
            framebuffer[i+j*width] = cast_ray(orig, dir, objects, lights);
        }
    }


    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    render();
    return 0;
}