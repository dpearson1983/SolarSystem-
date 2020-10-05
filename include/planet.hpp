#ifndef _PLANET_HPP_
#define _PLANET_HPP_

#include <vector>
#include <string>
#include <cmath>

template<typename T>
std::vector<T>& operator+=(std::vector<T> &lhs, const std::vector<T> &rhs) {
    if (lhs.size() != rhs.size())
        throw std::length_error("vectors must be same size to add");

    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

template<typename T>
std::vector<T> operator+ (std::vector<T> lhs, const std::vector<T> &rhs) {
    if (lhs.size() != rhs.size())
        throw std::length_error("vectors must be same size to add");
    return lhs += rhs;
}

template<typename T>
std::vector<T> operator/=(std::vector<T> &lhs, double rhs) {
    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] /= rhs;
    return lhs;
}

template<typename T>
std::vector<T> operator/(std::vector<T> lhs, double rhs) {
    return lhs /= rhs;
}

class planet{
    std::vector<double> pos;
    std::vector<double> vel;
    std::string name;
    double a, e, mass, M_Sun, G;
    
    std::vector<double> f(double h, double t, std::vector<double> &y, std::vector<planet> &planets) {
        double r_mag = std::sqrt(y[0]*y[0] + y[2]*y[2] + y[4]*y[4]);
        std::vector<double> r_hat = {-y[0]/r_mag, -y[2]/r_mag, -y[4]/r_mag};
        double F_g = this->G*this->M_Sun/(r_mag*r_mag);
        std::vector<double> acc = {F_g*r_hat[0], F_g*r_hat[1], F_g*r_hat[2]};
        for (size_t i = 0; i < planets.size(); ++i) {
            if (planets[i].getName() != this->name) {
                //std::cout << this->name << " " << planets[i].getName() << "\n";
                std::vector<double> pos = planets[i].getPos();
                std::vector<double> r_pp = {pos[0] - y[0], pos[1] - y[2], pos[2] - y[4]};
                r_mag = std::sqrt(r_pp[0]*r_pp[0] + r_pp[1]*r_pp[1] + r_pp[2]*r_pp[2]);
                F_g = this->G*planets[i].getMass()/(r_mag*r_mag);
                r_hat[0] = r_pp[0]/r_mag;
                r_hat[1] = r_pp[1]/r_mag;
                r_hat[2] = r_pp[2]/r_mag;
                acc[0] += F_g*r_hat[0];
                acc[1] += F_g*r_hat[1];
                acc[2] += F_g*r_hat[2];
            }
        }
        std::vector<double> F = {h*y[1], h*acc[0], h*y[3], h*acc[1], h*y[5], h*acc[2]};
        return F;
    }
    
    void rk4(double t, std::vector<double> &y, double h, std::vector<planet> &planets) {
        std::vector<double> k_1 = f(h, t, y, planets);
        std::vector<double> yt = y + k_1/2.0;
        std::vector<double> k_2 = f(h, t + h/2.0, yt, planets);
        yt = y + k_2/2.0;
        std::vector<double> k_3 = f(h, t + h/2.0, yt, planets);
        yt = y + k_3;
        std::vector<double> k_4 = f(h, t + h, yt, planets);
        for (size_t i = 0; i < y.size(); ++i) {
            y[i] += (1.0/6.0)*(k_1[i] + 2.0*k_2[i] + 2.0*k_3[i] + k_4[i]);
        }
    }
    
public:
    planet(std::string name, std::vector<double> pos_0, std::vector<double> vel_0, double mass, double M_Sun, double G) {
        this->pos = pos_0;
        this->vel = vel_0;
        this->mass = mass;
        this->name = name;
        this->M_Sun = M_Sun;
        this->G = G;
    }
    
    void update(double t, double h, std::vector<planet> &planets) {
        std::vector<double> y;
        for (size_t i = 0; i < this->pos.size(); ++i) {
            y.push_back(this->pos[i]);
            y.push_back(this->vel[i]);
        }
        this->rk4(t, y, h, planets);
        for (size_t i = 0; i < this->pos.size(); ++i) {
            this->pos[i] = y[2*i];
            this->vel[i] = y[2*i + 1];
        }
    }
    
    std::vector<double> getPos() {
        return this->pos;
    }
    
    std::vector<double> getVel() {
        return this->vel;
    }
    
    std::string getName() {
        return this->name;
    }
    
    double getMass() {
        return this->mass;
    }
    
};

#endif
