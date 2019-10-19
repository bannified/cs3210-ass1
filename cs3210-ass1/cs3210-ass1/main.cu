#ifdef _WIN64
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "thrust/sort.h"

#include <cmath>

#pragma region vector2

struct vector2;

double magnitude(const vector2& p);

struct vector2
{
	double x, y;

	vector2()
	{
	}

	vector2(double x, double y)
		: x(x), y(y)
	{
	};

	__device__ vector2& operator+=(const vector2& rhs)
	{
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	__device__ vector2& operator-=(const vector2& rhs)
	{
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}
	__device__ vector2& operator*=(const double& rhs)
	{
		x *= rhs;
		y *= rhs;
		return *this;
	}
	__device__ vector2& operator/=(const double& rhs)
	{
		x /= rhs;
		y /= rhs;
		return *this;
	}
	__device__ vector2 operator-() const
	{
		return { -x, -y };
	}

	__device__ void normalize()
	{
		double mag = magnitude(*this);
		x /= mag;
		y /= mag;
	}
};

__device__ inline vector2 operator+(vector2 lhs, const vector2& rhs)
{
	lhs += rhs;
	return lhs;
}
__device__ inline vector2 operator-(vector2 lhs, const vector2& rhs)
{
	lhs -= rhs;
	return lhs;
}
__device__ inline vector2 operator*(vector2 lhs, const double& rhs)
{
	lhs *= rhs;
	return lhs;
}
__device__ inline vector2 operator/(vector2 lhs, const double& rhs)
{
	lhs /= rhs;
	return lhs;
}
__device__ inline double operator*(const vector2& lhs, const vector2& rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y;
}

__device__ double distSq(const vector2& p, const vector2& q)
{
	return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
}
__device__ double dist(const vector2& p, const vector2& q)
{
	return std::sqrt(distSq(p, q));
}

__device__ double dist(const vector2& p)
{
	return std::sqrt((p.x) * (p.x) + (p.y) * (p.y));
}

__device__ double magnitudeSq(const vector2& p)
{
	return p * p;
}

__device__ double magnitude(const vector2& p)
{
	return std::sqrt(p * p);
}

__device__ double dotProduct(const vector2& a, const vector2& b)
{
	return a * b;
}

#pragma endregion

typedef enum {
    MODE_PRINT,
    MODE_PERF
} simulation_mode_t;

typedef struct {
    double i;
	vector2 position;
    /*double x;
    double y;
    double vx;
    double vy;*/
	vector2 velocity;
    int p_collisions;
    int w_collisions;
} particle_t;

__constant__ int l, r, s;

__managed__ particle_t* particles;
__constant__ int n;
int host_n;

struct Collision
{
	__device__ Collision(uint32_t index1, uint32_t index2, double stepValue)
		: index1(index1), index2(index2), stepValue(stepValue)
	{
	}

	int index1;
	int index2;
	double stepValue;

	bool operator<(const Collision& rhs) const
	{
		return stepValue < rhs.stepValue;
	}
};

__managed__ int* numCollisions; // numCollisions for each 
__managed__ Collision* collisionSteps; // 2D grid of collision's step sizes

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

__device__ void checkParticleCollisions(int particleIndex)
{
	numCollisions[particleIndex] = 0;
	const particle_t& current = particles[particleIndex];
	for (int j = particleIndex + 1; j < n; j++) {
		const particle_t& target = particles[j];
		double step = detectParticleCollision_cuda(current, target);
		if (isStepValid(step)) {
			collisionSteps[numCollisions[particleIndex]++] = Collision(particleIndex, j, step);
		}
	}
}


__host__ double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

__global__ void simulate_step(int num_threads)
{
    int i = blockIdx.x * num_threads + threadIdx.x;

	checkParticleCollisions(i);

    /* Dummy code that does not check for collision or walls */
    //particles[i].x += particles[i].vx;
    //particles[i].y += particles[i].vy;
}

__host__ void randomise_particles()
{
    srand(time(NULL));
	double minVelocity = 1 / 4;
	double maxVelocity = 1 / (8 * r);
    /* TODO Implement randomisation */
    for (int i = 0; i < host_n; i++) {
		int sign1 = (rand() % 2) ? 1 : -1;
		int sign2 = (rand() % 2) ? 1 : -1;
		particles[i].position = vector2(fRand(r, 1 - r), fRand(r, 1 - r));
		particles[i].velocity = vector2(sign1 * fRand(minVelocity, maxVelocity), sign2 * fRand(minVelocity, maxVelocity));
		//particles[i].x = fRand(r, l - r);
		//particles[i].y = fRand(r, 1 - r);
		//particles[i].vx = sign1 * fRand(minVelocity, maxVelocity);
		//particles[i].vy = sign2 * fRand(minVelocity, maxVelocity);
    }
}

__host__ void print_particles(int step)
{
    int i;
    for (i = 0; i < host_n; i++) {
        printf("%d %d %d %d %d %d\n", step, i, particles[i].position.x, particles[i].position.y,
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
    int i, x, y, vx, vy;
    int num_blocks, num_threads;
    int step;
    int host_l, host_r, host_s;
    simulation_mode_t mode;
    char mode_buf[6];

    if (argc != 3) {
        printf("Usage:\n%s num_blocks num_threads\n", argv[0]);
        return 1;
    }

    num_blocks = atoi(argv[1]);
    num_threads = atoi(argv[2]);

    scanf("%d", &host_n);
    scanf("%d", &host_l);
    scanf("%d", &host_r);
    scanf("%d", &host_s);
    scanf("%5s", mode_buf);

    cudaMallocManaged(&particles, sizeof(particle_t) * host_n);
	cudaMallocManaged(&collisionSteps, sizeof(double) * host_n * host_n); // n x n matrix of results

    for (i = 0; i < host_n; i++) {
        particles[i].i = -1;
        particles[i].p_collisions = 0;
        particles[i].w_collisions = 0;
    }

    while (scanf("%d %d %d %d %d", &i, &x, &y, &vx, &vy) != EOF) {
        particles[i].i = i;
		particles[i].position = vector2(x, y);
        //particles[i].x = x;
        //particles[i].y = y;
		particles[i].velocity = vector2(vx, vy);
        //particles[i].vx = vx;
        //particles[i].vy = vy;
    }

    if (particles[0].i == -1) {
        randomise_particles();
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

        /* Call the kernel */
        simulate_step<<<num_blocks, num_threads>>>(num_threads);

        /* Barrier */
        cudaDeviceSynchronize();
    }

    print_statistics(host_s);

    return 0;
}
