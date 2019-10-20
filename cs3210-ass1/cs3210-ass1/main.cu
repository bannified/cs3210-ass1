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
#include "thrust/host_vector.h"
#include <thrust/execution_policy.h>
#include "device_vector2.h"
#include "math.h"
#include <vector>
#include <cuda_profiler_api.h>

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

struct Collision {
    __host__ __device__ Collision() {}
    __device__ Collision(int i, int j, double stepValue) : i(i), j(j), stepValue(stepValue) {}

    int i;
    int j;
    double stepValue;
};

__host__ __device__ bool operator<(const Collision& lhs, const Collision& rhs) {
    return lhs.stepValue < rhs.stepValue;
}

__managed__ int* numCollisions;
Collision* collisions, *collisions_map;
__managed__ bool* resolved;
__managed__ Collision* validCollisions;

__device__ bool isStepValid(double step) {
    return 0 <= step && step < 1;
}

// 2 is returned as "infinity"
__device__ double detectParticleCollision_cuda(particle_t a, particle_t b) {
    double invalidationAdd = 0;

    double distance = dist(b.position, a.position);
    double sumRadii = r + r;
    distance -= sumRadii;

    vector2 resultVector = a.velocity - b.velocity;
    double resultMag = magnitude(resultVector);

    if (resultMag < distance) return 2;

    vector2 unitResultVector = resultVector;
    unitResultVector.normalize();

    vector2 c = b.position - a.position;
    double d = unitResultVector * c;

    if (d <= 0) return 2;

    double lengthC = magnitude(c);
    double fSquared = lengthC * lengthC - d * d;
    double sumRadiiSquared = sumRadii * sumRadii;

    // Escape: closest that a will get to b.
    if (fSquared >= sumRadiiSquared) return 2;

    double tSquared = sumRadiiSquared - fSquared;

    // negative tSquared. Probably don't have to do this check because the one preceding 
    // this one already ensures that tSquared isn't negative.
    if (tSquared < 0) return 2;

    double distanceToCollide = d - std::sqrt(tSquared);

    // Ensure that distance A has to move to touch B
    // is not greater than the magnitude of the movement vector
    if (resultMag < distanceToCollide) return 2;

    // the final displacement that the particle would have just before the collision.
    // can also choose to return this in a result vector;
    vector2 finalVector = unitResultVector * distanceToCollide;

    return magnitude(finalVector) / resultMag;
}

__device__ void checkParticleCollisions(Collision* coll_map, int i, const int max_collisions) {
    // i: particle index
    const particle_t& current = particles[i];
    for (int j = 0; j < i + 1; j++) {
        coll_map[i * max_collisions + j] = Collision(i, j, 2.0);
    }
    for (int j = i + 1; j < max_collisions - 1; j++) {
        const particle_t& target = particles[j];
        double step = detectParticleCollision_cuda(current, target);
        coll_map[i * max_collisions + j] = Collision(i, j, step);
    }
}

__device__ Collision detectWallCollision_cuda(const particle_t& p) {
    vector2 end_pos = p.position + p.velocity;
    Collision result(0, 0, 2.0); // stepValue > 1 means no collision
    // TODO: reduce branching
    if (end_pos.x - r <= 0) { // left, -1
        Collision temp = Collision(p.i, -1, (r - p.position.x) / p.velocity.x);
        if (temp < result) result = temp;
    }
    if (end_pos.x + r >= l) { // right, -2
        Collision temp = Collision(p.i, -2, (l - r - p.position.x) / p.velocity.x);
        if (temp < result) result = temp;
    }
    if (end_pos.y - r <= 0) { // bottom, -3
        Collision temp = Collision(p.i, -3, (r - p.position.y) / p.velocity.y);
        if (temp < result) result = temp;
    }
    if (end_pos.y + r >= l) { // top, -4
        Collision temp = Collision(p.i, -4, (l - r - p.position.y) / p.velocity.y);
        if (temp < result) result = temp;
    }

    return result;
}

__device__ void checkWallCollisions(Collision* coll_map, int i, const int max_collisions) {
    coll_map[i * max_collisions + max_collisions-1] = detectWallCollision_cuda(particles[i]);
}

void gatherCollisions(thrust::host_vector<Collision>& resultVector, const int numParticles) {
    for (int i = 0; i < numParticles; i++) {
        int numColl = numCollisions[i];
        for (int j = 0; j < numColl; j++) {
            resultVector.push_back(collisions_map[i * numParticles + j]);
        }
    }
}

__host__ double fRand(double fMin, double fMax) {
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

__global__ void runCollisionChecks(Collision* coll_map, int numParticles, int threadsTotal, int chunkNo)
{
    int i = chunkNo * threadsTotal + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numParticles) return;

    checkParticleCollisions(coll_map, i, numParticles + 1);

    checkWallCollisions(coll_map, i, numParticles + 1);
}

__device__ double clamp(double d, double min, double max) {
    const double t = d < min ? min : d;
    return t > max ? max : t;
}

// keep a particle within bounds
__device__ void clampParticleBounds(particle_t& p) {
    double x = p.position.x;
    double y = p.position.y;
    p.position.x = clamp(x, r, s - r);
    p.position.y = clamp(y, r, s - r);
}

__device__ void resolveWallCollision(particle_t& p, int wall, double stepProportion) {
    switch (wall) {
        case -1:
            p.position += p.velocity * stepProportion;
            p.velocity.x *= -1;
            break;
        case -2:
            p.position += p.velocity * stepProportion;
            p.velocity.x *= -1;
            break;
        case -3:
            p.position += p.velocity * stepProportion;
            p.velocity.y *= -1;
            break;
        case -4:
            p.position += p.velocity * stepProportion;
            p.velocity.y *= -1;
    }
    p.position += p.velocity * (1 - stepProportion);
}

__device__ void resolveParticleCollision(particle_t& a, particle_t& b, double stepProportion) {
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

__global__ void resolveCollisions(int size, int threadsTotal, int chunkNo) {
    int i = chunkNo * threadsTotal + blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= size) return;

    Collision res = validCollisions[i];
    if (res.j < 0) {
        resolveWallCollision(particles[res.i], res.j, res.stepValue);
        clampParticleBounds(particles[res.i]);
        particles[res.i].w_collisions++;
    } else {
        resolveParticleCollision(particles[res.i], particles[res.j], res.stepValue);
        clampParticleBounds(particles[res.i]);
        clampParticleBounds(particles[res.j]);
        particles[res.i].p_collisions++;
        particles[res.j].p_collisions++;
    }
}

__global__ void moveUnresolvedParticles(bool* resolved, int numParticles, int threadsTotal, int chunkNo) {
    int i = chunkNo * threadsTotal + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numParticles) return;

    if (!resolved[i]) {
        particles[i].position += particles[i].velocity;
    }
}

__host__ void print_particles(int step) {
    for (int i = 0; i < host_n; i++) {
        printf("%d %d %10.8f %10.8f %10.8f %10.8f\n", step, i, particles[i].position.x, particles[i].position.y,
                particles[i].velocity.x, particles[i].velocity.y);
    }
}

__host__ void print_statistics(int num_step) {
    for (int i = 0; i < host_n; i++) {
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
    cudaMallocManaged(&collisions, sizeof(Collision) * host_n * (host_n+1)); // [particle][collided object (host_n=wall)]
    
    cudaSetDeviceFlags(cudaDeviceMapHost);
    cudaHostAlloc(&collisions, sizeof(Collision) * host_n * (host_n+1), cudaHostAllocMapped);
    cudaHostGetDevicePointer(&collisions_map, collisions, 0);
    
    cudaMallocManaged(&resolved, sizeof(bool) * host_n);
    cudaMallocManaged(&validCollisions, sizeof(Collision) * host_n);

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

    int threadsTotal = num_blocks * num_threads;

    cudaProfilerStart();
    for (step = 0; step < host_s; step++) {
        if (mode == MODE_PRINT || step == 0) {
            print_particles(step);
        }

        int numChunks = ceil((double)host_n / (double)threadsTotal);
        for (int chunkNo = 0; chunkNo < numChunks; chunkNo++) {
            /* Check collisions */
            runCollisionChecks<<<num_blocks, num_threads>>>(collisions_map, host_n, threadsTotal, chunkNo);
        }

        /* Barrier */
        cudaDeviceSynchronize();

        thrust::sort(collisions, collisions + host_n * (host_n + 1));

        // Collision Validation
        cudaMemset(resolved, 0, host_n);

        cudaMemset(validCollisions, 0, host_n * sizeof(Collision));

        int next_idx = 0;
        for (int i = 0; collisions[i].stepValue < 1.5; i++) {
            Collision res = collisions[i];
            if (resolved[res.i]) continue;
            if (res.j < 0) {
                validCollisions[next_idx] = res;
                ++next_idx;
                resolved[res.i] = true;
            }
            else {
                if (resolved[res.j]) continue;
                validCollisions[next_idx] = res;
                ++next_idx;
                resolved[res.i] = true;
                resolved[res.j] = true;
            }
        }

        int numValidCollisionChunks = ceil(next_idx / (double)threadsTotal);

        for (int chunkNo = 0; chunkNo < numValidCollisionChunks; chunkNo++) {
            resolveCollisions<<<num_blocks, num_threads>>>(next_idx, threadsTotal, chunkNo);
        }
        cudaDeviceSynchronize();

        for (int chunkNo = 0; chunkNo < numChunks; chunkNo++) {
            moveUnresolvedParticles<<<num_blocks, num_threads>>>(resolved, host_n, threadsTotal, chunkNo);
        }
        cudaDeviceSynchronize();
    }
    cudaProfilerStop();

    print_statistics(host_s);

    return 0;
}
