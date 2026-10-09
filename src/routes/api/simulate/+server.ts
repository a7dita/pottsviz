import { randomInt } from 'node:crypto';
import { resolve } from 'node:path';
import manifest from '$lib/topologies/manifest.json';
import { json } from '@sveltejs/kit';
import type { RequestHandler } from './$types';
import { streamSimulation } from '$lib/server/simulation';

export const config = { runtime: 'nodejs22.x', maxDuration: 300 };

export const POST: RequestHandler = async ({ request }) => {
  let body;
  try { body = await request.json(); } catch { return json({ message: 'Invalid JSON.' }, { status: 400 }); }
  const { model, S, w, tau2, steps = 5000, posterior, frontal, lambda, topology = 'one-to-one' } = body ?? {};
  const number = (value: unknown, min: number, max: number) =>
    typeof value === 'number' && Number.isFinite(value) && value >= min && value <= max;
  const region = (value: { S?: number; w?: number; tau2?: number } | undefined, maxTau: number) =>
    value && Number.isInteger(value.S) && number(value.S, 3, 11) &&
    number(value.w, 0, 2) && number(value.tau2, 100, maxTau);
  if (!Number.isInteger(steps) || !number(steps, 1, 5000)) {
    return json({ message: 'Please choose between 1 and 5000 sweeps.' }, { status: 400 });
  }
  let args: (number | string)[];
  let metadata = {};
  if (model === 'homo' && Number.isInteger(S) && number(S, 3, 11) && number(w, 0, 2) && number(tau2, 100, 800)) {
    args = [S, w, tau2, steps];
  } else if (model === 'hybrid' && region(posterior, 800) && region(frontal, 1600) &&
    number(lambda, 0, 1) && topology === 'one-to-one') {
    args = [posterior.w, frontal.w, posterior.S, frontal.S, posterior.tau2, frontal.tau2, lambda, steps];
    metadata = { memories: 50 };
  } else if (model === 'topology' && region(posterior, 800) && region(frontal, 1600)) {
    const selected = manifest.topologies.find(item => item.id === topology);
    const global = body.global;
    const cue = (value: { cue?: number; randomCue?: boolean }) =>
      Number.isInteger(value.cue) && number(value.cue, 0, 48) && typeof value.randomCue === 'boolean';
    if (!selected || !number(posterior.lambda, 0, 1) || !number(frontal.lambda, 0, 1) ||
      !global || !number(global.U, 0, 0.6) || !number(global.beta, 1, 21) ||
      !number(global.a, 0.1, 0.4) || !number(global.tau1, 5, 35) || !number(global.density, 0.05, 0.25) ||
      !cue(posterior) || !cue(frontal)) {
      return json({ message: 'Choose a saved topology, valid parameters and regional cues.' }, { status: 400 });
    }
    // Fresh independent seeds are generated on the server for every Start.
    // This also refreshes dilution and update tables; the topology stays frozen.
    const patternSeed = randomInt(1, 2147482001), runtimeSeed = randomInt(1, 2147482001);
    args = [resolve('.simulators/topologies', selected.file), posterior.w, frontal.w,
      posterior.S, frontal.S, posterior.tau2, frontal.tau2, posterior.lambda, frontal.lambda,
      global.U, global.beta, global.a, global.tau1, global.density,
      patternSeed, runtimeSeed, posterior.cue, frontal.cue,
      Number(posterior.randomCue), Number(frontal.randomCue), 256, steps];
    metadata = { memories: 49, topology: selected.id, patternSeed, runtimeSeed,
      posteriorCue: posterior.randomCue ? 'random' : posterior.cue,
      frontalCue: frontal.randomCue ? 'random' : frontal.cue };
  } else {
    return json({ message: 'Please choose valid model parameters.' }, { status: 400 });
  }
  return new Response(streamSimulation(model, args, request.signal, metadata), {
    headers: { 'Content-Type': 'application/x-ndjson', 'Cache-Control': 'no-store', 'X-Accel-Buffering': 'no' }
  });
};
