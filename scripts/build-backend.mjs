import { spawnSync } from 'node:child_process';
import { mkdirSync, copyFileSync } from 'node:fs';
import { basename, join } from 'node:path';

mkdirSync('.simulators', { recursive: true });
const root = 'backend/model_potts_homo';
const result = spawnSync('g++', [
  '-O3', '-std=gnu++17', '-Wl,-rpath,$ORIGIN/lib', `-I${root}/include`,
  `${root}/main.cpp`, `${root}/functions.cpp`, `${root}/rand_gen.cpp`,
  '-o', '.simulators/homo'
], { stdio: 'inherit' });
if (result.error) throw result.error;
if (result.status !== 0) process.exit(result.status ?? 1);

// Vercel's build image has the compiler's shared C++ libraries, but its Node
// runtime need not. Ship those exact libraries beside the executable.
mkdirSync('.simulators/lib', { recursive: true });
for (const library of ['libstdc++.so.6', 'libgcc_s.so.1']) {
  const lookup = spawnSync('g++', [`-print-file-name=${library}`], { encoding: 'utf8' });
  if (lookup.status !== 0) throw new Error(`Could not locate ${library}`);
  const path = lookup.stdout.trim();
  if (path === library) throw new Error(`Missing compiler runtime: ${library}`);
  copyFileSync(path, join('.simulators/lib', basename(path)));
}
