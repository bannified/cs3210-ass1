/**
    CS3210 Assignment 1 Part 1 (OpenMP)
    Authors: Kuan Wei Heng (A0121712X)
             Chong Jun Hong, Dominic (A0121237U)
    std=C++11
*/
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "vector2.h"
#include "collision.h"
#include "Particle.h"
#include <mpi.h>

int procs;
vector2 gStageSize;
double gStepSize;
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;
// std::vector<int> gCollisionCounts;
// Num of particles, Size of square, Radius of particle, and number of steps
int N;
int L;
double r;

int nprocs;
int myid;
int workers;
int size;
#define MASTER_ID 0
int numParticlesPerWorker;

MPI_Datatype collisionDataType;
MPI_Datatype vector2DataType;
MPI_Datatype particleDataType;
int vector2MPISize;

inline void PrintParticle(const Particle p);

inline void PrintParticleFull(const Particle p);

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

void master()
{
    /* ------ Setup ------ */
    std::cin >> N >> L >> r >> gNumSteps;

    int buffer[3];
    double rBuffer;
    buffer[0] = N;
    buffer[1] = L;
    buffer[2] = gNumSteps;
    rBuffer = r;

    fprintf(stderr, "Rank: %d // N: %d//L: %d//r: %lf//gNumSteps: %d\n", myid, N, L, r, gNumSteps);

    for (int i = 1; i < nprocs; i++) {
        MPI_Send(&buffer[0], 3, MPI_INT, i, i + 100, MPI_COMM_WORLD);
        MPI_Send(&rBuffer, 1, MPI_DOUBLE, i, i + 101, MPI_COMM_WORLD);
    }

    std::vector<Particle> particles;
    particles.reserve(N);

    gStageSize = vector2(L, L);

    std::string mode;
    std::cin >> mode;
    if (mode == "print") {
        gPrintAll = true;
    }
    
    {
        int particleIndex;
        while (std::cin >> particleIndex) {
            vector2 initialPosition, initialVelocity;
            std::cin >> initialPosition.x >> initialPosition.y >> initialVelocity.x >> initialVelocity.y;
            particles.emplace_back(Particle(initialPosition, initialVelocity, r, particleIndex));
        }
    }

    // Generate up to N random particles
    srand(time(NULL));
    double minVelocity = L / 4;
    double maxVelocity = L / (8 * r);
    while ((int)particles.size() < N) {
        int sign = (rand() % 2) ? 1 : -1;
        vector2 initialPosition(fRand(r, L - r), fRand(r, L - r));
        vector2 initialVelocity(sign * fRand(minVelocity, maxVelocity), sign * fRand(minVelocity, maxVelocity));

        particles.push_back(Particle(initialPosition, initialVelocity, r, particles.size()));
    }

    /* ------ End Particles Setup ------ */

    numParticlesPerWorker = N / workers;

    // Start simulation for gNumSteps
    for (; gStepNumber < gNumSteps; gStepNumber++) {
        std::vector<Collision> collisions;
        collisions.reserve(N);

        // Print all particles
        if (gPrintAll) {
            for (const Particle p : particles) PrintParticle(p);
        }

        int workerId;
        // send to all worker processes the particles
        for (workerId = 0; workerId < workers; workerId++) {
            MPI_Send(&(particles[0]), N, particleDataType, workerId, 1, MPI_COMM_WORLD);
        }

        // get all collision results
        for (int i = 0; i < workers; i++) {
            std::vector<Collision> buffer;
            int numCollisions = 0;
            MPI_Status stat;
            MPI_Recv(&numCollisions, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &stat);
            buffer.resize(numCollisions);
            MPI_Recv(&(buffer[0]), numCollisions, collisionDataType, i, 1, MPI_COMM_WORLD, &stat);
            
            // transfer from buffer to master's collisions
            for (int j = 0; j < numCollisions; j++) {
                collisions.push_back(buffer[j]);
            }
        }

        // sort all collisions by increasing time
        std::sort(collisions.begin(), collisions.end());

        std::vector<bool> resolved(N);

        std::vector<Collision> validCollisions; // Stores all valid collision results to be resolved
        validCollisions.reserve(N / 2);

        // Pick out the collisions that are valid starting from the smallest step value
        // and not allowing for repeated collisions for any one particle.
        for (int i = 0; i < collisions.size(); i++) {
            Collision res = collisions[i];
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

        // resolve valid collisions
        //#pragma omp parallel for default(none) shared(particles, validCollisions, gStageSize)
        for (int i = 0; i < validCollisions.size(); i++) {
            Collision res = validCollisions[i];
            if (res.index2 < 0) {
                resolveWallCollision(particles[res.index1], res.index2, res.stepValue, gStageSize);
                clamp(particles[res.index1], gStageSize);
                particles[res.index1].numWallCollisions++;
            }
            else {
                resolveParticleCollision(particles[res.index1], particles[res.index2], res.stepValue);
                clamp(particles[res.index1], gStageSize);
                clamp(particles[res.index2], gStageSize);
                particles[res.index1].numParticleCollisions++;
                particles[res.index2].numParticleCollisions++;
            }
        }
        // move remaining particles that were not involved in collisions
        //#pragma omp parallel for default(none) shared(particles, resolved)
        for (int i = 0; i < particles.size(); i++) {
            if (!resolved[i]) {
                particles[i].position += particles[i].velocity;
            }
        }
    }

    // Print all particles
    if (gPrintAll) {
        for (const Particle p : particles) PrintParticle(p);
    }

    for (const Particle p : particles) PrintParticleFull(p);
    // for (auto i : gCollisionCounts) std::cout << i << '\n';
}

/************** Worker Code **************/

void workerReceiveParticles(std::vector<Particle> &particles)
{
    MPI_Status stat;

    MPI_Recv(&(particles[0]), N, particleDataType, MASTER_ID, 1, MPI_COMM_WORLD, &stat);

    fprintf(stderr, " --- WORKER %d: Received particles\n", myid);
}

void workerComputeCollisions(std::vector<Particle> &particles, std::vector<Collision> &collisions)
{
    fprintf(stderr, "Worker %d start computing collisions\n", myid);
    int startIdx = (myid - 1) * numParticlesPerWorker;
    int endIdx = startIdx + numParticlesPerWorker;
    endIdx = std::min(endIdx, N);

    // particle-to-particle collision detection
    // Each process works on N / nProcs number of particles
    for (int i = startIdx; i < endIdx; i++) {
        const Particle& particle = particles[i];
        for (int j = i + 1; j < particles.size(); j++) {
            const Particle& target = particles[j];
            double step = detectParticleCollision(particle, target);
            if (isStepValid(step)) {
                //#pragma omp critical
                collisions.push_back({ particle.index, target.index, step });
            }
        }
    }

    // particle-to-particle collision detection
        // Each process works on N / nProcs number of particles
    for (int i = startIdx; i < endIdx; i++) {
        const Particle& particle = particles[i];
        for (int j = i + 1; j < particles.size(); j++) {
            const Particle& target = particles[j];
            double step = detectParticleCollision(particle, target);
            if (isStepValid(step)) {
                //#pragma omp critical
                collisions.push_back({ particle.index, target.index, step });
            }
        }
    }
}

void worker()
{
    MPI_Status Stat;
    int buffer[3] = {0, 0, 0};
    double rBuffer = 0;

    MPI_Recv(buffer, 3, MPI_INT, MASTER_ID, myid + 100, MPI_COMM_WORLD, &Stat);
    MPI_Recv(&rBuffer, 1, MPI_DOUBLE, MASTER_ID, myid + 101, MPI_COMM_WORLD, &Stat);

    N = buffer[0];
    L = buffer[1];
    gNumSteps = buffer[2];
    r = rBuffer;

    for (int step = 0; step < gNumSteps; step++) {
        std::vector<Collision> resultCollisions;
        std::vector<Particle> particles;
        resultCollisions.reserve(N);
        particles.resize(N);

        workerReceiveParticles(particles);
        workerComputeCollisions(particles, resultCollisions);

        // send collision results
        int numCollisions = resultCollisions.size();
        MPI_Send(&numCollisions, 1, MPI_INT, MASTER_ID, 1, MPI_COMM_WORLD);
        MPI_Send(&(resultCollisions[0]), numCollisions, collisionDataType, MASTER_ID, 1, MPI_COMM_WORLD);

    }
}

/****************************************/

int main(int argc, char *argv[])
{
    //std::cout << "Usage: " << argv[0] << " <procs>\n";

    if (argc >= 2)
        procs = atoi(argv[1]);
    else
        procs = -1;

    // MPI Setup
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    int collisionBlockLengths[3] = { 1, 1, 1 };
    MPI_Aint collisionDisplacements[3] = {
                                            offsetof(Collision, index1),
                                            offsetof(Collision, index2),
                                            offsetof(Collision, stepValue)
    };
    MPI_Datatype collisionDataTypes[3] = { MPI_INT, MPI_INT, MPI_DOUBLE };

    // init collision DataType
    MPI_Type_create_struct(3, collisionBlockLengths, collisionDisplacements, collisionDataTypes, &collisionDataType);
    MPI_Type_commit(&collisionDataType);

    // init vector2 DataType
    int vector2BlockLengths[2] = { 1, 1 };
    MPI_Aint vector2Displacements[2] = {
                                        offsetof(vector2, x),
                                        offsetof(vector2, y)
    };
    MPI_Datatype vector2DataTypes[2] = { MPI_DOUBLE, MPI_DOUBLE };
    MPI_Type_create_struct(2, vector2BlockLengths, vector2Displacements, vector2DataTypes, &vector2DataType);
    MPI_Type_commit(&vector2DataType);

    MPI_Type_size(vector2DataType, &vector2MPISize);

    // init particle DataType
    int particleBlockLengths[6] = { 1, 1, 1, 1, 1 };
    MPI_Aint particleDisplacements[6] = {
                                            offsetof(Particle, position),
                                            offsetof(Particle, velocity),
                                            offsetof(Particle, radius),
                                            offsetof(Particle, index),
                                            offsetof(Particle, numWallCollisions),
                                            offsetof(Particle, numParticleCollisions)
    };
    MPI_Datatype particleDataTypes[6] = { vector2DataType, vector2DataType, MPI_DOUBLE, MPI_UINT32_T, MPI_UINT32_T, MPI_UINT32_T };
    MPI_Type_create_struct(6, particleBlockLengths, particleDisplacements, particleDataTypes, &particleDataType);
    MPI_Type_commit(&particleDataType);

    workers = nprocs - 1;

    if (myid == MASTER_ID) {
        master();
    }
    else {
        worker();
    }

    MPI_Finalize();
    return 0;

}

inline void PrintParticle(const Particle p)
{
    std::cout << std::fixed << std::setprecision(8);
    std::cout << gStepNumber << ' ';
    std::cout << p.index << ' ';
    std::cout << p.position.x << ' ';
    std::cout << p.position.y << ' ';
    std::cout << p.velocity.x << ' ';
    std::cout << p.velocity.y << '\n';
}

inline void PrintParticleFull(const Particle p)
{
    std::cout << std::fixed << std::setprecision(8);
    std::cout << gStepNumber << ' ';
    std::cout << p.index << ' ';
    std::cout << p.position.x << ' ';
    std::cout << p.position.y << ' ';
    std::cout << p.velocity.x << ' ';
    std::cout << p.velocity.y << ' ';
    std::cout << p.numParticleCollisions << ' ';
    std::cout << p.numWallCollisions << '\n';
}
