import { spawnSync } from 'node:child_process';
import { mkdirSync } from 'node:fs';

mkdirSync('.simulators', { recursive: true });
const root = 'backend/model_potts_homo';
const result = spawnSync('g++', [
  '-O3', '-std=gnu++17', '-static-libstdc++', '-static-libgcc', `-I${root}/include`,
  `${root}/main.cpp`, `${root}/functions.cpp`, `${root}/rand_gen.cpp`,
  '-o', '.simulators/homo'
], { stdio: 'inherit' });
if (result.error) throw result.error;
if (result.status !== 0) process.exit(result.status ?? 1);
