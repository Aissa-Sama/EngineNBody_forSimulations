// initial_conditions.cpp
#include "initial_conditions.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>

// Constante matemática portable (funciona en cualquier compilador)
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 6.28318530717958647692;

// ============================================================================
// Sistema solar simplificado (Sol + Tierra)
// ============================================================================
NBodySystem InitialConditions::solar_system() {
    NBodySystem system;
    system.G = 1.0;
    
    // Sol en el centro
    system.bodies.push_back({{0, 0, 0}, {0, 0, 0}, 1.0});
    
    // Tierra en órbita circular (masa 3e-6 masas solares)
    double a = 1.0;  // 1 UA en nuestras unidades
    double v = std::sqrt(system.G * 1.0 / a);  // v = sqrt(G*M/a)
    system.bodies.push_back({{a, 0, 0}, {0, v, 0}, 3e-6});
    
    return system;
}

// ============================================================================
// Binaria + estrella de campo (para probar el integrador híbrido)
// ============================================================================
NBodySystem InitialConditions::binary_with_field() {
    NBodySystem system;
    system.G = 1.0;
    
    // Binaria de igual masa en órbita circular
    double a = 1.0;  // separación
    double v = 0.5 * std::sqrt(system.G * (1.0 + 1.0) / a);  // v = 0.5 * sqrt(G*M/a)
    
    system.bodies.push_back({{-a/2, 0, 0}, {0, -v, 0}, 1.0});
    system.bodies.push_back({{ a/2, 0, 0}, {0,  v, 0}, 1.0});
    
    // Estrella de campo (lejana)
    system.bodies.push_back({{10, 0, 0}, {0, -0.1, 0}, 1.0});
    
    return system;
}

// ============================================================================
// Órbita figura-8 de 3 cuerpos (Chenciner & Montgomery, 2000)
// ============================================================================
NBodySystem InitialConditions::figure_eight() {
    NBodySystem system;
    system.G = 1.0;

    // Solución periódica figura-8 de Chenciner & Montgomery (2000)
    // Referencia: Ann. of Math. 152, 881-901 (2000)
    // Valores numéricos de Simó (2002), precisión 16 dígitos.
    // CM en reposo exacto: sum(m_i * r_i) = 0, sum(m_i * v_i) = 0
    // Período T = 6.3259 (unidades G=1, m=1)
    const double m = 1.0;

    system.bodies.push_back({
        {  0.9700436926041021, -0.2430865994534989, 0.0 },
        {  0.4662036850003824,  0.4323657300939067, 0.0 },
        m
    });
    system.bodies.push_back({
        { -0.9700436926041021, -0.2430865994534989, 0.0 },
        {  0.4662036850003824, -0.4323657300939067, 0.0 },
        m
    });
    system.bodies.push_back({
        {  0.0,                 0.4861731989069978, 0.0 },
        { -0.9324073700007648,  0.0,                0.0 },
        m
    });

    return system;
}
// ============================================================================
// Sistema aleatorio para testing
// ============================================================================
NBodySystem InitialConditions::random_system(int n_bodies, double radius) {
    NBodySystem system;
    system.G = 1.0;
    
    // Generador de números aleatorios (semilla fija para reproducibilidad)
    std::mt19937 gen(42);  // Semilla fija para resultados consistentes
    std::uniform_real_distribution<> pos_dist(-radius, radius);
    std::uniform_real_distribution<> vel_dist(-0.5, 0.5);
    std::uniform_real_distribution<> mass_dist(0.5, 1.5);
    
    for (int i = 0; i < n_bodies; ++i) {
        Body b;
        b.position = {pos_dist(gen), pos_dist(gen), pos_dist(gen)};
        b.velocity = {vel_dist(gen), vel_dist(gen), vel_dist(gen)};
        b.mass = mass_dist(gen);
        system.bodies.push_back(b);
    }
    
    return system;
}

// ============================================================================
// Cúmulo de Plummer (distribución realista de estrellas)
// ============================================================================
NBodySystem InitialConditions::plummer_cluster(int n_bodies) {
    NBodySystem system;
    system.G = 1.0;
    
    // Generador de números aleatorios
    std::mt19937 gen(42);
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    
    double total_mass = 0.0;
    std::vector<Body> bodies;
    
    for (int i = 0; i < n_bodies; ++i) {
        Body b;
        
        // 1. Generar radio según distribución de Plummer: f(r) ∝ r^2 / (1 + r^2)^(5/2)
        // Usamos el método de transformada inversa
        double u = uniform(gen);
        double r = 1.0 / std::sqrt(std::pow(u, -2.0/3.0) - 1.0);
        
        // 2. Dirección aleatoria en esfera unitaria
        double theta = TWO_PI * uniform(gen);
        double phi = std::acos(2.0 * uniform(gen) - 1.0);
        
        b.position.x = r * std::sin(phi) * std::cos(theta);
        b.position.y = r * std::sin(phi) * std::sin(theta);
        b.position.z = r * std::cos(phi);
        
        // 3. Velocidad según distribución de Plummer (simplificada)
        // En una implementación real, aquí iría la función de distribución de velocidades
        b.velocity = {0, 0, 0};  // Por ahora, partículas quietas
        
        // 4. Masa: todas iguales para simplificar
        b.mass = 1.0;
        total_mass += b.mass;
        
        bodies.push_back(b);
    }
    
    // Normalizar masas para que la masa total sea 1
    for (auto& b : bodies) {
        b.mass /= total_mass;
        system.bodies.push_back(b);
    }
    
    return system;
}

// ============================================================================
// Binaria Kepleriana pura (para tests de KS)
// ============================================================================
NBodySystem InitialConditions::kepler_binary(double a, double e) {
    NBodySystem system;
    system.G = 1.0;
    
    // Masas iguales para simplificar
    double m1 = 1.0;
    double m2 = 1.0;
    double M = m1 + m2;
    
    // Velocidad para órbita circular (e=0) o elíptica simplificada
    double v_circular = std::sqrt(system.G * M / a);
    
    // Para excentricidad e, en periastro: r_peri = a*(1-e), v_peri = v_circular * sqrt((1+e)/(1-e))
    double r_peri = a * (1.0 - e);
    double v_peri = v_circular * std::sqrt((1.0 + e) / (1.0 - e));
    
    // Posiciones en periastro (a lo largo del eje X)
    system.bodies.push_back({{-r_peri * m2/M, 0, 0}, {0, -v_peri * m2/M, 0}, m1});
    system.bodies.push_back({{ r_peri * m1/M, 0, 0}, {0,  v_peri * m1/M, 0}, m2});
    
    return system;
}