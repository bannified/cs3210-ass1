import subprocess

BINARY = './a.out'
IN_FILE_FORMAT = 'input-%d'

THREAD_COUNTS = [1, 2, 4, 8, 16, 32]
PARTICLE_COUNTS = [100, 200, 400, 800, 1600, 3200, 6400, 12800]

def run(command, particles, samples):
    perf = '(perf stat -e cache-references -e cache-misses -e cycles -e instructions -- %s < %s > /dev/null) 2>&1' \
        % (command, IN_FILE_FORMAT % particles)
    output = subprocess.check_output(perf, shell=True)
    _,_,_,a,b,c,d,_,e,_,_ = output.decode('utf-8').split('\n')
    results = list(zip(*[[eval(s.strip().split()[0].replace(',','')) for s in (a,b,c,d,e)] for _ in range(samples)]))
    cache_references = max(results[0])
    cache_misses     = min(results[1])
    cycles           = min(results[2])
    instructions     = max(results[3])
    time             = min(results[4]) 
    return cache_references, cache_misses, cycles, instructions, time

for threads in THREAD_COUNTS:
    for particles in PARTICLE_COUNTS:
        command = '%s %d' % (BINARY, threads)
        cache_references, cache_misses, cycles, instructions, time = run(command, particles, 5)
        print('Threads:', threads, 'Particles:', particles)
        print(cache_references, cache_misses, cycles, instructions, time)
