#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "geometry.h"

struct Sphere {
    Vec3f center;
    float radius;
    Vec3f color;
    Sphere(const Vec3f &c, const float &r, const Vec3f &col) : center(c), radius(r), color(col) {}

    bool ray_intersect(const Vec3f &orig, Vec3f &dir, Vec3f &cross_point, Vec3f& N) const {
        float cross_dist;
        Vec3f orig_center = center - orig;
        if (orig_center * dir <= 0 && (orig - center).norm() > radius) return false;
        Vec3f center_ray_proj = orig + dir * ((orig_center * dir) / (dir.norm()));
        float center_dist = (center_ray_proj - center).norm();
        if (center_dist > radius) return false;
        if (center_dist < radius){
            float incide_part = radius * radius - center_dist * center_dist;
            if ((orig - center).norm() > radius){
                cross_dist = (center_ray_proj - orig).norm() - incide_part;
            } else {
                if (orig_center * dir <= 0)
                    cross_dist = incide_part - (center_ray_proj - orig).norm();
                else 
                    cross_dist = (center_ray_proj - orig).norm() + incide_part;
            }
            cross_point = orig + dir * cross_dist;
        } else {
            cross_point = center_ray_proj;

        }
        
        N = (cross_point - center).normalize();
        return true;
        
    }

    bool ray_intersect2(const Vec3f &orig, const Vec3f &dir, float &t0) const {
        Vec3f L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

struct Light {
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position;
    float intensity;
};

Vec3f cast_ray(const Vec3f &orig,  Vec3f &dir, std::vector<Sphere>& spheres, std::vector<Light>& lights) {
    // float sphere_dist = std::numeric_limits<float>::max();
    for (auto sphere : spheres){
        Vec3f point, N;
        if (sphere.ray_intersect(orig, dir, point, N)){
            float diff_light_intens = 0;
            for(auto light : lights){
                Vec3f light_dir = (point - light.position).normalize();
                diff_light_intens += light.intensity * std::max(0.f, light_dir * N);
            }
            return sphere.color * diff_light_intens;
        } 
    }
    
    return Vec3f(0.2, 0.7, 0.8);
}

Vec3f red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1);

void render() {
    const int width    = 1024;
    const int height   = 768;
    std::vector<Vec3f> framebuffer(width*height);
    
    float fov = M_PI / 2;

    std::vector<Sphere> spheres= { Sphere(Vec3f( 1.5, -0.5, -18), 2, red), 
    Sphere(Vec3f( 7,    5,   -18), 5, green), Sphere(Vec3f(10, 0, -20), 5, blue)};
    std::vector<Light> lights = {Light(Vec3f(-10, 20, 20), 1.5)};

    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)* width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, spheres, lights);
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