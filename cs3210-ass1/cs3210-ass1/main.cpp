/**
	CS3210 Assignment 1 Part 1 (OpenMP)
	Authors:	Kuan Wei Heng (A.......)
				Chong Jun Hong, Dominic (A0121237U)
	std=C++11
*/

#include <vector>
#include <string>
#include <iostream>

Vector2 gStageSize;
float gStepSize; 
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;

int gRadius;

std::vector<Particle> particles;

inline void PrintParticle(const Particle particle);

inline void PrintParticle(const Particle particle)
{
	printf("%i %i %10.8f %10.8f %10.8f\n %10.8f\n\n",
		   particle.index,
		   gStepNumber,
		   particle.position.x,
		   particle.position.y,
		   particle.velocity.x,
		   particle.velocity.y
	);
}

struct Vector2
{
	Vector2(double x, double y) 
		: x(x), y(y) {};

	double x;
	double y;
};

struct Particle
{
	Particle(Vector2 pos, Vector2 vel, uint32_t index)
		: position(pos), velocity(vel), index(index) {};

	Vector2 position;
	Vector2 velocity;

	uint32_t index;

	uint32_t numWallCollisions;
	uint32_t numParticleCollisions;
};

int main(int argc, char *argv[])
{
	// Num of particles, Size of square, Radius of particle, and number of steps
	int N, L;

	scanf("%i", &N);
	scanf("%i", &L);
	gStageSize = Vector2(L, L);
	scanf("%i", &gRadius);
	scanf("%i", &gNumSteps);

	std::string inputBuffer;
	scanf("%s", inputBuffer);
	if (inputBuffer.compare("print")) {
		// print for every timestep
		gPrintAll = true;
	}

	int runningIndex = 0;
	particles.reserve(N);

	// Generate particles based on data
	while (getline(std::cin, inputBuffer)) {
		// TODO: Create particle
	}

	// Generate random particles

	for (int i = runningIndex; i < N; i++) {
		// TODO: Create random particles
	}

	return 0;
}