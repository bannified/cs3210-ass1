#ifdef _WIN64
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include "thrust/sort.h"
#include "thrust/device_ptr.h"
#include <thrust/execution_policy.h>
#include "device_vector2.h"
#include "math.h"
#include <vector>

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

struct particle_t {
    double i;
    vector2 position;
    vector2 velocity;
    int p_collisions;
    int w_collisions;
};

__constant__ int l, r, s;

__managed__ particle_t* particles;
__constant__ int n;
int host_n;

struct Collision
{
    __device__ Collision()
    {
    }

    __device__ Collision(int index1, int index2, double stepValue)
        : index1(index1), index2(index2), stepValue(stepValue)
    {
    }

    int index1;
    int index2;
    double stepValue;

    __device__ bool operator<(Collision& rhs)
    {
        return stepValue < rhs.stepValue;
    }
};

__host__ __device__ bool operator<(const Collision& lhs, const Collision& rhs)
{
    return lhs.stepValue < rhs.stepValue;
}

__managed__ int* numCollisions; // numCollisions for each 
__managed__ Collision* collisionSteps;

__device__ bool isStepValid(double step)
{
    return 0 <= step && step < 1;
}

// called per particle
__device__ double detectParticleCollision_cuda(particle_t a, particle_t b)
{
    double distance = dist(b.position, a.position);
    double sumRadii = r + r;
    distance -= sumRadii;

    vector2 resultVector = a.velocity - b.velocity;
    double resultMag = magnitude(resultVector);

    // Early escape: Can't reach just based on maximum travel possible
    if (resultMag < distance) {
        return -1;
    }

    vector2 unitResultVector = resultVector;
    unitResultVector.normalize();

    vector2 c = b.position - a.position;
    double d = unitResultVector * c;

    // Early escape: result vector does not cause the particles to move closer together.
    // since one is not moving in the direction of the other.
    if (d <= 0) {
        return -1;
    }

    double lengthC = magnitude(c);
    double fSquared = lengthC * lengthC - d * d;

    double sumRadiiSquared = sumRadii * sumRadii;

    // Escape: closest that a will get to b.
    if (fSquared >= sumRadiiSquared) {
        return -1;
    }

    double tSquared = sumRadiiSquared - fSquared;

    // negative tSquared. Probably don't have to do this check because the one preceding 
    // this one already ensures that tSquared isn't negative.
    if (tSquared < 0) {
        return -1;
    }

    double distanceToCollide = d - std::sqrt(tSquared);

    // Ensure that distance A has to move to touch B
    // is not greater than the magnitude of the movement vector
    if (resultMag < distanceToCollide) {
        return -1;
    }

    // the final displacement that the particle would have just before the collision.
    // can also choose to return this in a result vector;
    vector2 finalVector = unitResultVector * distanceToCollide;

    return magnitude(finalVector) / resultMag;
}

__device__ void checkParticleCollisions(int particleIndex, const int numParticles)
{
    numCollisions[particleIndex] = 0;
    const particle_t& current = particles[particleIndex];
    for (int j = particleIndex + 1; j < n; j++) {
        const particle_t& target = particles[j];
        double step = detectParticleCollision_cuda(current, target);
        if (isStepValid(step)) {
            collisionSteps[particleIndex * numParticles + numCollisions[particleIndex]++] = Collision(particleIndex, j, step);
        }
    }
}

__device__ Collision detectWallCollision_cuda(const particle_t p)
{
    vector2 end_pos = p.position + p.velocity;
    Collision result(0, 0, 2.0); // stepValue > 1 means no collision
    // TODO: reduce branching
    if (end_pos.x - r <= 0) { // left, -1
        Collision temp = Collision(p.i, -1, (r - p.position.x) / p.velocity.x);
        if (temp < result) {
            result = temp;
        }
    }
    if (end_pos.x + r >= l) { // right, -2
        Collision temp = Collision(p.i, -2, (l - r - p.position.x) / p.velocity.x);
        if (temp < result) {
            result = temp;
        }
    }
    if (end_pos.y - r <= 0) { // bottom, -3
        Collision temp = Collision(p.i, -3, (r - p.position.y) / p.velocity.y);
        if (temp < result) {
            result = temp;
        }
    }
    if (end_pos.y + r >= l) { // top, -4
        Collision temp = Collision(p.i, -4, (l - r - p.position.y) / p.velocity.y);
        if (temp < result) {
            result = temp;
        }
    }

    return result;
}

__device__ void checkWallCollisions(int particleIndex, const int numParticles)
{
    const particle_t& current = particles[particleIndex];
    Collision result = detectWallCollision_cuda(current);
    if (isStepValid(result.stepValue)) {
        collisionSteps[particleIndex * numParticles + numCollisions[particleIndex]++] = result;
    }
}

__host__ double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

__global__ void runCollisionChecks(int numParticles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    checkParticleCollisions(i, numParticles);

    checkWallCollisions(i, numParticles);
    /* Dummy code that does not check for collision or walls */
    //particles[i].x += particles[i].vx;
    //particles[i].y += particles[i].vy;
}

void sortCollisions(Collision* unsortedColls, const int numColls)
{
    thrust::device_ptr<Collision> ptr(unsortedColls);
    // TODO: Change to use thrust::sort_by_key (radix sort)
    thrust::sort(thrust::device, ptr, ptr + numColls);
}

__device__ double clamp(double d, double min, double max)
{
    const double t = d < min ? min : d;
    return t > max ? max : t;
}

// keep a particle within bounds
__device__ void clampParticleBounds(particle_t& p)
{
    double x = p.position.x;
    double y = p.position.y;
    p.position.x = clamp(x, r, s - r);
    p.position.y = clamp(y, r, s - r);
}

__device__ void resolveWallCollision(particle_t& p, int wall, double stepProportion)
{
    if (wall == -1) {
        p.position += p.velocity * stepProportion;
        p.velocity.x *= -1;
    }
    else if (wall == -2) {
        p.position += p.velocity * stepProportion;
        p.velocity.x *= -1;
    }
    else if (wall == -3) {
        p.position += p.velocity * stepProportion;
        p.velocity.y *= -1;
    }
    else {
        p.position += p.velocity * stepProportion;
        p.velocity.y *= -1;
    }
    p.position += p.velocity * (1 - stepProportion);
}

__device__ void resolveParticleCollision(particle_t& a, particle_t& b, double stepProportion)
{
    vector2 aImpact = a.position + a.velocity * stepProportion;
    vector2 bImpact = b.position + b.velocity * stepProportion;

    double d = dist(aImpact, bImpact);

    vector2 n = vector2((bImpact.x - aImpact.x) / d, (bImpact.y - aImpact.y) / d);
    double p = 2 * (a.velocity * n - b.velocity * n) / 2;

    a.velocity = a.velocity - n * p * 1;
    b.velocity = b.velocity + n * p * 1;

    a.position = aImpact + a.velocity * (1.0f - stepProportion);
    b.position = bImpact + b.velocity * (1.0f - stepProportion);
}

__global__ void resolveCollisions(Collision* validCollisions)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    Collision res = validCollisions[i];
    if (res.index2 < 0) {
        resolveWallCollision(particles[res.index1], res.index2, res.stepValue);
        clampParticleBounds(particles[res.index1]);
        particles[res.index1].w_collisions++;
    }
    else {
        resolveParticleCollision(particles[res.index1], particles[res.index2], res.stepValue);
        clampParticleBounds(particles[res.index1]);
        clampParticleBounds(particles[res.index2]);
        particles[res.index1].p_collisions++;
        particles[res.index2].p_collisions++;
    }
}

__global__ void moveUnresolvedParticles(int* resolvedArr)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i > n) {
        return;
    }

    if (!resolvedArr[i]) {
        particles[i].position += particles[i].velocity;
    }
}

__host__ void print_particles(int step)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", step, i, particles[i].position.x, particles[i].position.y,
            particles[i].velocity.x, particles[i].velocity.y);
    }
}

__host__ void print_statistics(int num_step)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f %d %d\n", num_step, i, particles[i].position.x,
            particles[i].position.y, particles[i].velocity.x, particles[i].velocity.y,
            particles[i].p_collisions, particles[i].w_collisions);
    }
}

int main(int argc, char** argv)
{
    int i;
    double x, y, vx, vy;
    int num_blocks, num_threads;
    int step;
    int host_l, host_r, host_s;
    simulation_mode_t mode;
    char mode_buf[6];

    if (argc != 3) {
        printf("Usage:\n%s num_blocks numParticles\n", argv[0]);
        return 1;
    }

    num_blocks = atoi(argv[1]);
    num_threads = atoi(argv[2]);

    scanf("%d", &host_n);
    scanf("%d", &host_l);
    scanf("%d", &host_r);
    scanf("%d", &host_s);
    scanf("%5s", mode_buf);

    cudaMallocManaged(&numCollisions, sizeof(int) * host_n);
    cudaMallocManaged(&particles, sizeof(particle_t) * host_n);
    cudaMallocManaged(&collisionSteps, sizeof(double) * host_n * host_n); // n x n matrix of results

    for (i = 0; i < host_n; i++) {
        particles[i].i = -1;
        particles[i].p_collisions = 0;
        particles[i].w_collisions = 0;
    }

    while (scanf("%d %lf %lf %lf %lf", &i, &x, &y, &vx, &vy) != EOF) {
        particles[i].i = i;
        particles[i].position = vector2(x, y);
        particles[i].velocity = vector2(vx, vy);
    }

    if (particles[0].i == -1) {
        srand(time(NULL));
        double minVelocity = 1 / 4;
        double maxVelocity = 1 / (8 * host_r);
        for (int i = 0; i < host_n; i++) {
            int sign1 = (rand() % 2) ? 1 : -1;
            int sign2 = (rand() % 2) ? 1 : -1;
            particles[i].position = vector2(fRand(host_r, 1 - host_r), fRand(host_r, 1 - host_r));
            particles[i].velocity = vector2(sign1 * fRand(minVelocity, maxVelocity), sign2 * fRand(minVelocity, maxVelocity));
        }
    }

    mode = strcmp(mode_buf, "print") == 0 ? MODE_PRINT : MODE_PERF;

    /* Copy to GPU constant memory */
    cudaMemcpyToSymbol(n, &host_n, sizeof(n));
    cudaMemcpyToSymbol(l, &host_l, sizeof(l));
    cudaMemcpyToSymbol(r, &host_r, sizeof(r));
    cudaMemcpyToSymbol(s, &host_s, sizeof(s));

    for (step = 0; step < host_s; step++) {
        if (mode == MODE_PRINT || step == 0) {
            print_particles(step);
        }

        /* Check collisions */
        runCollisionChecks<<<num_blocks, num_threads>>>(host_n);

        /* Barrier */
        cudaDeviceSynchronize();

        /* Sort with thrust */
        sortCollisions(collisionSteps, *numCollisions);

        // Collision Validation
        std::vector<int> resolved(host_n, 0);

        std::vector<Collision> validCollisions; // Stores all valid collision results to be resolved
        validCollisions.reserve(host_n / 2);

        for (int i = 0; i < *numCollisions; i++) {
            Collision res = collisionSteps[i];
            if (resolved[res.index1]) continue;
            if (res.index2 < 0) {
                validCollisions.push_back(res);
                resolved[res.index1] = true;
            }
            else {
                if (resolved[res.index2]) continue;
                validCollisions.push_back(res);
                resolved[res.index1] = true;
                resolved[res.index2] = true;
            }
        }

        resolveCollisions<<<num_blocks, num_threads >>>(&validCollisions[0]);

        moveUnresolvedParticles<<<num_blocks, num_threads >>>(&resolved[0]);
    }

    print_statistics(host_s);

    return 0;
}
