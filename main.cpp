#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// #include "geometry.h"
#define EPS    0.0001
#define C_EPS  0.00001
#define REQ    4
#define INF    100
using Vec = glm::vec3;

struct Primitive {
    virtual bool ray_intersect(Vec &orig, Vec &dir, Vec &cross_point) = 0;
    virtual Vec get_normal(Vec& point) = 0;

};

struct Material{
    Vec color;
    float blink_density;
    float mirror;
    Material(Vec c, float bl, float mirr) : color(c), blink_density(bl), mirror(mirr){}
    Material(){}
};


struct Sphere {
    Vec center;
    float radius;
    Material props;
    Sphere(Vec c, float r, Material m) 
    : center(c), radius(r), props(m) {}

    Sphere(){}

    bool ray_intersect(Vec &orig, Vec &dir, Vec &cross_point) {
        float cross_dist;
        Vec orig_center = center - orig;
        if (glm::dot(orig_center, dir) <= 0) return false;

        Vec center_ray_proj = orig + dir * (glm::dot(orig_center, dir) / glm::length(dir));

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

    Vec get_normal(Vec point, Vec dir, Vec orig){
        return glm::normalize(point - center);
    }
};

struct Triangle {
    Vec v0, v1, v2;
    Vec N;
    Material props;

    Triangle(Vec p0, Vec p1, Vec p2, Material m) : v0(p0), v1(p1), v2(p2), props(m){
        N = glm::normalize(glm::cross(p1 - p0, p2 - p0));
    }
    Triangle() {}

    bool ray_intersect(Vec &orig, Vec &dir, Vec &cross_point) {
        Vec E1 = v1 - v0;
        Vec E2 = v2 - v0;
        Vec T  = orig - v0;
        Vec P  = glm::cross(dir, E2);
        Vec Q  = glm::cross(T, E1);
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

    Vec get_normal(Vec point, Vec dir, Vec orig){
        if (glm::dot(dir, N) > 0) return N * (-1.0F);  
        return  N;
    }
};

struct Rect {
    Vec v0, v1, v2, v3; //vertexes clocwise
    Triangle t1, t2;
    Vec N;
    Material props;
    
    Rect() {}
    Rect(Vec p0, Vec p1, Vec p2, Vec p3, Material m) {
        t1 = Triangle(p0, p1, p2, m);
        t2 = Triangle(p0, p3, p2, m);
        N = glm::normalize(glm::cross(p1 - p0, p2 - p0));
        props = m;
    }

    bool ray_intersect(Vec &orig, Vec &dir, Vec &cross_point) {
        if(t1.ray_intersect(orig, dir, cross_point)) return true;
        if(t2.ray_intersect(orig, dir, cross_point)) return true;
        return false;
    }

    Vec get_normal(Vec point, Vec dir, Vec orig){
        if (glm::dot(dir, N) > 0) return N * (-1.0F);  
        return  N;
    }

};

struct Piramid {
    Vec v0, v1, v2, v3;
    std::vector<Triangle> T;
    Material props;

    Piramid(Vec p0, Vec p1, Vec p2, Vec p3, Material m){
        T.push_back(Triangle(p0, p1, p2, m));
        T.push_back(Triangle(p0, p3, p2, m));
        T.push_back(Triangle(p0, p1, p3, m));
        T.push_back(Triangle(p1, p3, p2, m));
        props = m;
    }

    Piramid(){}
    bool ray_intersect(Vec &orig, Vec &dir, Vec &cross_point) {
        Triangle res;
        bool intersect = false;
        float min_dist = INF;
        for (auto triangle : T) {
            if (triangle.ray_intersect(orig, dir, cross_point) && min_dist > glm::distance(orig, cross_point)){
                min_dist = glm::distance(orig, cross_point);
                intersect = true;
                res = triangle;
            }
        }
        res.ray_intersect(orig, dir, cross_point);
        return intersect;
    }
    Vec get_normal(Vec point, Vec dir, Vec orig){
        Triangle res;
        Vec cross_point;
        float min_dist = INF;
        for (auto triangle : T) {
            if (triangle.ray_intersect(orig, dir, cross_point) && min_dist > glm::distance(orig, cross_point)){
                min_dist = glm::distance(orig, cross_point);
                res = triangle;
            }
        }
        return res.get_normal(point, dir, orig);
    }
};

struct Light {
    Light(const Vec &p, const float &i) : position(p), intensity(i) {}
    Vec position;
    float intensity;
};

struct Objects {
    std::vector<Sphere> spheres; 
    std::vector<Triangle> triangles;
    std::vector<Rect> rects;
    std::vector<Piramid> piramids;
};

template<typename T>
bool get_shadow(std::vector<T>& objects, Vec point, Light light, Vec orig){
    
    Vec rev_light_dir = glm::normalize(light.position - point);
    
    for (auto obj : objects){
        Vec rev_point;
        Vec moved_point = point + obj.get_normal(point, rev_light_dir * (-1.0F), orig) * float(EPS);
        if (obj.ray_intersect(moved_point, rev_light_dir, rev_point)) {
            if (glm::distance(rev_point, point) < glm::distance(light.position, point)){
                return true;
            }
        }
    }
    
    return false;
}
Vec cast_ray(Vec &orig,  
             Vec &dir, 
             Objects& objects,
             std::vector<Light>& lights, 
             int req = REQ);

template<typename T>
Vec get_color(Vec &orig, Vec &dir, T& obj, Objects& objects, std::vector<Light>& lights, int req){
    Vec point;
    obj.ray_intersect(orig, dir, point);
    Vec N = obj.get_normal(point, dir, orig);
    float diff_light_intens = 0;
    float blinks = 0;
    for(auto light : lights){
        Vec light_dir = glm::normalize(light.position - point);
        if (get_shadow(objects.spheres, point, light, orig)) continue;
        if (get_shadow(objects.triangles, point, light, orig)) continue;
        // if (get_shadow(objects., point, light, orig)) continue;

        diff_light_intens += light.intensity * std::max(0.f, glm::dot(light_dir, N));
        Vec reflect_dir = glm::reflect(light_dir * float(-1), N);
        float cos_val = glm::dot(glm::normalize(reflect_dir), (dir * float(-1)));
        blinks += light.intensity * std::max(0.0, pow(cos_val, obj.props.blink_density));
        
    }
    Vec reflect_dir = glm::normalize(N * float(2.0) + dir);
    Vec mirror_color = cast_ray(point, reflect_dir, objects, lights, req - 1);
    return obj.props.color * (diff_light_intens + blinks) * (1 - obj.props.mirror) + mirror_color * obj.props.mirror;
}

template <typename T> 
bool nearest_dist(std::vector<T> objects, T& nearest, Vec &orig, Vec &dir, float& min_dist) {
    bool intersect = false;
    for (auto obj : objects){
        Vec point;
        if (obj.ray_intersect(orig, dir, point)){
            if (glm::distance(orig, point) < min_dist) {

            intersect = true;
                min_dist = glm::distance(orig, point);
                nearest = obj;
            }
        } 
    }
    return intersect;
}

Vec cast_ray(Vec &orig,  
             Vec &dir, 
             Objects& objects,
             std::vector<Light>& lights, 
             int req) {

    if (!req) return Vec(0, 0, 0);
    Sphere n_sph;
    Triangle n_tr;
    Rect n_rect;
    Piramid n_pir;
    float min_dist = INF;
    int intersect = 0;
    if (nearest_dist(objects.spheres, n_sph, orig, dir, min_dist)) intersect = 1;
    if (nearest_dist(objects.triangles, n_tr, orig, dir, min_dist)) intersect = 2;
    if (nearest_dist(objects.rects, n_rect, orig, dir, min_dist)) intersect = 3;
    if (nearest_dist(objects.piramids, n_pir, orig, dir, min_dist)) intersect = 4;
    if (intersect == 1){
        return get_color(orig, dir, n_sph, objects, lights, req);
    }
    if (intersect == 2){
        return get_color(orig, dir, n_tr, objects, lights, req);;      
    }
    if (intersect == 3){
        return get_color(orig, dir, n_rect, objects, lights, req);;      
    }
    if (intersect == 4){
        return get_color(orig, dir, n_pir, objects, lights, req);;      
    }
    return Vec(0, 0, 0);
}

Vec red(1, 0, 0), green(1, 1, 1), blue(0, 0, 1);
Vec yellow(1, 1, 0), white(1, 1, 1), purple(1, 0, 1);

Material glass_gr(green, 1000, 0.9), resin_red(red, 30, 0), metal_blue(blue, 700, 0.75), metal_yellow(yellow, 700, 0.75);
Material purple_wood_polished(purple, 1000, 0.8);

void render() {
    const int width    = 1024;
    const int height   = 768;
    std::vector<Vec> framebuffer(width*height);
    
    float fov = M_PI / 2;
    Objects objects;
    objects.spheres= { Sphere(Vec(0, 4, -45 ), 9.0, resin_red), 
                                   Sphere(Vec(-60, -5, -50), 9.0, metal_blue),
                                   Sphere(Vec(-10, -10, -30), 9.0, glass_gr),
                                   Sphere(Vec(20, 20, -35), 12, glass_gr)};

    std::vector<Light> lights = {Light(Vec(0, 10, 0), 0.9), Light(Vec(35, -17, -30), 0.6),
     Light(Vec(-10, -40, -20), 0.7), };

    objects.triangles = {Triangle(Vec(-40, 10, -40), 
                                  Vec(-30, 10, -40),
                                  Vec(-30, 0, -40), resin_red),

                         Triangle(Vec(-40, 10, -40), 
                                  Vec(-40, 0, -45),
                                  Vec(-30, 20, -40), resin_red),

                         Triangle(Vec(-40, 30, -40), 
                                  Vec(-30, 30, -40),
                                  Vec(-40, 20, -45), resin_red),

                         Triangle(Vec(-30, 20, -40), 
                                  Vec(-30, 30, -40),
                                  Vec(-40, 20, -45), resin_red)};
    
    objects.rects = {Rect(Vec(-INF, -20, INF), Vec(INF, -20, INF), Vec(INF, -20, -INF), Vec(-INF, -20, -INF), metal_yellow)};

    objects.piramids = {Piramid(Vec(40, 0, -40), Vec(50, -15, -30), Vec(30, -15, -40), Vec(20, -15, -20), purple_wood_polished)};

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)* width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec dir = glm::normalize(Vec(x, y, -1));
            Vec orig = Vec(0, 0, 0);
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