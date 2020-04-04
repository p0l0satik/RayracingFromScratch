#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "geometry.h"
#define EPS  0.001
#define REQ  2
struct Sphere {
    Vec3f center;
    float radius;
    Vec3f color;
    float blink_density;
    float mirror ;
    Sphere(const Vec3f &c, const float &r, const Vec3f &col, float bd, float m) 
    : center(c), radius(r), color(col), blink_density(bd), mirror(m) {}

    bool ray_intersect(const Vec3f &orig, Vec3f &dir, Vec3f &cross_point, Vec3f& N) const {
        float cross_dist;
        Vec3f orig_center = center - orig;
        if (orig_center * dir <= 0) return false;

        Vec3f center_ray_proj = orig + dir * ((orig_center * dir) / (dir.norm()));

        float center_dist = (center_ray_proj - center).norm();
        if (center_dist > radius) return false;
        if (center_dist < radius){
            float incide_part = sqrt(radius * radius - center_dist * center_dist);
            cross_dist = (center_ray_proj - orig).norm() - incide_part;
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

Vec3f cast_ray(const Vec3f &orig,  Vec3f &dir, std::vector<Sphere>& spheres, std::vector<Light>& lights, int req = REQ) {
    // float sphere_dist = std::numeric_limits<float>::max();
    if (!req) return Vec3f(0, 0, 0);

    std::vector < std::pair <float, Vec3f>> col_dist;
    for (auto sphere : spheres){
        Vec3f point, N;
        if (sphere.ray_intersect(orig, dir, point, N)){
            float diff_light_intens = 0;
            float blinks = 0;
            for(auto light : lights){
                Vec3f light_dir = (light.position - point).normalize();
                bool use_light = true;
                for (auto cross_sp : spheres) {
                    Vec3f rev_point, n;
                    Vec3f moved_point = point * (1 + EPS ) ;
                    Vec3f rev_light_dir = light_dir ;
                    if (cross_sp.ray_intersect(moved_point, rev_light_dir, rev_point, n)) {
                        // if ((rev_point - point).norm() > (orig - point).norm()){
                            
                        // }
                        use_light = false;
                        break;
                    }
                }
                if (!use_light) continue;
                diff_light_intens += light.intensity * std::max(0.f, light_dir * N);
                Vec3f reflect_dir = N * 2.0 + light_dir * (-1);
                blinks += light.intensity * std::max(0.0, pow((reflect_dir).normalize() * (dir * (-1)),
                                                     sphere.blink_density));
                
            }
            Vec3f reflect_dir = N * 2.0 - dir;
            Vec3f mirror_color = cast_ray(point, reflect_dir.normalize(), spheres, lights, req - 1);
            Vec3f col = sphere.color * (diff_light_intens + blinks) * (1 - sphere.mirror) + mirror_color * sphere.mirror;

            col_dist.push_back(std::make_pair((orig - point).norm(), col));
        } 
    }
    float min_dist = 1000000000;
    Vec3f col;
    if (col_dist.size()){
        for(auto dist : col_dist){
            if (dist.first < min_dist){
                col = dist.second;
                min_dist = dist.first;
            }

        }
        return col;
    }
    // if (req ==  REQ){
    //     return ;
    // }
    // return Vec3f(0, 0, 0);
    return Vec3f(0.2, 0.7, 0.8);
}

Vec3f red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1);
Vec3f yellow(1, 1, 0);

void render() {
    const int width    = 1024;
    const int height   = 768;
    std::vector<Vec3f> framebuffer(width*height);
    
    float fov = M_PI / 2;

    std::vector<Sphere> spheres= { Sphere(Vec3f(10, 4, -30 ), 8, red, 600, 0.5), 
                                   Sphere(Vec3f(-10, 4, -30), 8, blue, 700, 0.5),
                                   Sphere(Vec3f(0, -10, -25), 8, green, 1000, 0.9),
                                   Sphere(Vec3f(0, -10000000, -30), 9000, yellow, 1000, 0)};
    std::vector<Light> lights = {Light(Vec3f(20, 10, 0), 1), Light(Vec3f(0, 0, 0), 0.0),
     Light(Vec3f(-10, -40, -20), 0)};


    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)* width/(float)height;
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            Vec3f dir = Vec3f(x, y, -1).normalize();
            framebuffer[i+j*width] = cast_ray(Vec3f(0, 0, 0), dir.normalize(), spheres, lights);
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